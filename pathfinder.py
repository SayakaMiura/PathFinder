import sys
import os
import time
import shutil
import glob
import subprocess
import math
import copy
import pydot
import networkx
import random
import pylab
import datetime
import argparse
import random
import math
import types
from Bio import Phylo
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from io import StringIO
from tree_permute import permute_unique_trees

# read aln
# generate initial MP tree from aln and root on normal
# run ML ancestral reconstruction
# read EP data
# map EP data to nodes in initial tree
# assign each internal node to a tumor
# make edge list
# read migration file
# calculate accuracy
# generate outputs

#megacc_app = "megacc_200422.exe"
#megacc_app = "megacc.exe"
megacc_app = "megacc_11210415.exe"
mp_tree_infer_mao = "infer_NJ_amino_acid.mao"
ancestral_seqs_mao = "ancestral_seqs_ML_protein.mao"
outgroup_file = "outgroup.txt"
# resolve_polytomy_to_most_diverse_node = True

print_megacc_cmd = False
tryto_fix_anc_seq_inference = True
fix_0length_issues = True
# relax_node_threshold = True  # If node has no AA with probability >threshold, fall back to selecting single AA with the highest probability if it exists

aa_label_list = ['C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
tumor_label_dict = {'A': 'Normal'}
rev_tumor_label_dict = {'Normal': 'A'}

parser = argparse.ArgumentParser(description="PathFinder tumor migration path solver.")
parser.add_argument("aln", help="Clone sequence alignment file.", type=str)
parser.add_argument("clone_locations", help="Per tumor clone frequency/presence matrix file.", type=str)
parser.add_argument("--primary", help="Specify name of primary tumor.", type=str, default="Primary")
parser.add_argument("--anc_tumor_threshold", help="Lower limit on ancestral node tumor threshold.", type=float, default=0.15)
parser.add_argument("--mig_event_threshold", help="Migration edge detection threshold.", type=float, default=0.5)
parser.add_argument("--max_permutations", help="Maximum number of tree polytomy permutations to use.", type=int, default=30)
parser.add_argument("--max_graphs_per_tree", help="Maximum number of node-tumor membership states to analyze per initial tree polytomy permutation.", type=int, default=-1)
parser.add_argument("--verbosity", help="Status message verbosity from 0-4.", type=int, default=1)
parser.add_argument("--use_all_weighted_outputs", help="Make final weighted edge list from all outputs instead of count-based selection.", action='store_true', default=False)
# parser.add_argument("--machina_inputs", help="Generate mutation tree and label files suitable for MACHINA input.", action='store_true', default=False)
parser.add_argument("--infer_primary_seq", help="Infer primary clone sequence based on Normal and all clone sequences.", action='store_true', default=False)
parser.add_argument("--default_normal_char", help="Default character to use for Normal sequence, if Normal sequence is not specified in input alignment.", type=str, default=None)
parser.add_argument("--draw_all_outputs", help="Generate migration graph and phylogeny images for all sampled permutations/configurations.", action='store_true', default=False)
parser.add_argument("--phy_scale_x", help="Horizontal scaling factor for tree image outputs.", type=float, default=1.0)
parser.add_argument("--phy_scale_y", help="Vertical scaling factor for tree image outputs.", type=float, default=1.0)
parser.add_argument("--acc_by_edge_type", help="Break accuracy counts into P->M, M->M, M->P.", action='store_true', default=False)
parser.add_argument("-t", "--true_paths", help="List of true migration paths.", type=str, default=None)
parser.add_argument("-o", "--output", help="Output directory to put results in.", type=str, default=".")
parser.add_argument("--log_ancestral_inference", help="Keep inputs/outputs of MEGA ancestral sequence inference calculations for debugging purposes.", action='store_true', default=False)
parser.add_argument("--relax_threshold", help="If a node has no site with probability>threshold, fall back to selecting the most likely site(s).", action='store_true', default=False)
parser.add_argument("--use_select_weighted_outputs", help="Make final ancestral seqs probability weighted edge list from outputs with minimized count-based selection.", action='store_true', default=False)
parser.add_argument("--keep_ambiguous_results", help="Treat nodes with fully ambiguous inferences as belonging equally to all sites, instead of discarding the ancestral inference result.", action='store_true', default=False)

args = parser.parse_args()

v_level = args.verbosity

mega_io_logging = args.log_ancestral_inference

if args.output != "." and not os.path.exists(args.output):
	os.mkdir(args.output)

if mega_io_logging and not os.path.exists(os.path.join(args.output, "ancestral_inference_logging")):
	os.mkdir(os.path.join(args.output, "ancestral_inference_logging"))

dump_edges = True
edge_dump_file = None
if edge_dump_file is None: dump_edges = False

tumor_membership_cutoff = args.anc_tumor_threshold
avg_edge_weight_cutoff = args.mig_event_threshold

count_keys = ["TP", "FP", "FN"]

if args.acc_by_edge_type:
	count_keys = ['TP(P->M)', 'FP(P->M)', 'FN(P->M)', 'TP(M->M)', 'FP(M->M)', 'FN(M->M)', 'TP(M->P)', 'FP(M->P)', 'FN(M->P)']

if args.primary == "Primary":
	args.infer_primary_seq = True

def parse_clone_freqs(file_name):
	# Read tumor/clone frequencies file into dict of dicts
	clones = {}
	if args.primary == "Primary":
		clones.update({"Primary": {}})
		tumor = args.primary
		tumor_label_dict[aa_label_list[0]] = tumor
		rev_tumor_label_dict[tumor] = aa_label_list[0]
		aa_label_list.remove(aa_label_list[0])
		clones["Primary"].update({tumor: 1.0})
	with open(file_name) as file:
		line_num = 0
		for line in file.readlines():
			line_num += 1
			data = line.rstrip().split('\t')
			if line_num == 1:
				if data[0] == "Tumor":
					clone_list = data[1:]
					for clone in clone_list:
						clones.update({clone: {}})
				else:
					raise Exception("Unexpected header format")
			else:
				tumor = data[0]
				if tumor in ['PSec1', 'PSec2']:
					tumor = args.primary
				if len(aa_label_list) < 1:
					raise Exception("More than the maximum number of tumor sites (19) specified")
				if tumor not in tumor_label_dict.values():
					tumor_label_dict[aa_label_list[0]] = tumor
					rev_tumor_label_dict[tumor] = aa_label_list[0]
					aa_label_list.remove(aa_label_list[0])
				freqs = data[1:]
				for i in range(len(freqs)):
					if float(freqs[i]) == 0.0:
						pass
					else:
						# Note this currently erases all but the last non-zero frequency for a tumor/clone
						clones[clone_list[i]].update({tumor: freqs[i]})
	return clones


def make_machina_labels(clones, basename):
	labels_filename = basename + ".labeling"
	with open(os.path.join(args.output, labels_filename), 'w') as labels_file:
		for clone in clones.keys():
			for tumor in clones[clone].keys():
				if clones[clone][tumor] > 0.0:
					labels_file.write("{}_{} {}\n".format(clone, tumor, tumor))


def make_machina_trees(in_tree, basename):
	tree = copy.deepcopy(in_tree)
	labels_filename = basename + ".labeling"
	btree_edges = []
	mtree_edges = []
	clades = lookup_by_names(tree)
	for key in clades.keys():
		for child in clades[key].root.clades:
			btree_edges.append((clades[key].name, child.name.replace(":","_")))
	with open(os.path.join(args.output, labels_filename), 'r') as labels_file:
		for line in labels_file:
			values = line.strip().split(' ')[0].split('_')
			btree_edges.append((values[0], "{}_{}".format(values[0], values[1])))
	updated = True
	while updated:
		updated = False
		nodes = lookup_by_names(tree)
		for key in nodes.keys():
			for child in nodes[key].root.clades:
				if child.branch_length == 0.0 and len(child.root.clades) > 0 and "root_0" not in child.name:
					for gchild in nodes[child.name].root.clades:
						nodes[key].root.clades.append(gchild)
					nodes[key].root.clades.remove(child)
					updated = True
	clades = lookup_by_names(tree)
	for key in clades.keys():
		for child in clades[key].root.clades:
			mtree_edges.append((clades[key].name, child.name.replace(":", "_")))
	with open(os.path.join(args.output, labels_filename), 'r') as labels_file:
		for line in labels_file:
			values = line.strip().split(' ')[0].split('_')
			mtree_edges.append((values[0], "{}_{}".format(values[0], values[1])))
	btree_edges = [(x[0].replace("_root_0", ""), x[1].replace("_root_0", "")) for x in btree_edges]
	mtree_edges = [(x[0].replace("_root_0", ""), x[1].replace("_root_0", "")) for x in mtree_edges]
	with open(os.path.join(args.output, basename + ".btree"), 'w') as tree_file:
		for btree_edge in sorted(list(set(btree_edges))):
			if btree_edge[1] in ["Normal", "Primary"]: continue
			tree_file.write("{}\n".format(" ".join(btree_edge)))
	with open(os.path.join(args.output, basename + ".mtree"), 'w') as tree_file:
		for mtree_edge in sorted(list(set(mtree_edges))):
			if mtree_edge[1] in ["Normal", "Primary"]: continue
			tree_file.write("{}\n".format(" ".join(mtree_edge)))


def parse_input_aln(aln_file_in, aln_file_out):
	seqs = {}
	with open(aln_file_out, 'w') as out_file:
		with open(aln_file_in, 'r') as file:
			seqname = ""
			out_file.write("#MEGA\n")
			out_file.write("!Title SNVs;\n")
			out_file.write("!Format datatype=Protein;\n")
			for line in file:
				line = line.strip()
				if len(line) == 0: continue
				if line == "#MEGA" or line[0] == "!":
					continue
				if line[0] == "#" or line[0] == ">":
					seqname = line[1:]
				else:
					seqs[seqname] = line
		if 'Normal' not in seqs.keys():
			if args.default_normal_char == None:
				raise Exception("No Normal sequence found, Normal must be specified, or a default Normal sequence character must be specified with --default_normal_char option.")
			print("No Normal sequence found, using user specified default char {} * {}".format(args.default_normal_char,len(next(iter(seqs.values())))))
			seqs['Normal'] = 'A' * len(next(iter(seqs.values())))
		seq_len = len(seqs[next(iter(seqs.keys()))])
		target_seq_len = 3
		repeat_count = int(math.ceil(float(target_seq_len) / float(seq_len)))
		primary = ""
		for i in range(seq_len):
			normal_allele = seqs['Normal'][i]
			mut_alleles = list(set([seqs[seqid][i] for seqid in seqs.keys() if seqid != 'Normal']))
			if len(mut_alleles) == 1 and mut_alleles[0] != normal_allele:
				allele = mut_alleles[0]
			else:
				allele = normal_allele
			primary = primary + allele
		for key in seqs.keys(): seqs[key] = seqs[key] * repeat_count
		primary = primary * repeat_count
		if args.infer_primary_seq and primary not in seqs.values():
			seqs['Primary'] = seqs.get('Primary', primary)
		for key in seqs.keys():
			out_file.write("#" + key + '\n')
			out_file.write(seqs[key] + '\n')
	return seqs


def read_paths(filename):
	paths = []
	with open(filename, 'r') as file:
		for line in file:
			data = line.strip().split('->')
			paths.append((data[0], data[1]))
	return paths


def lookup_by_names(tree):
	names = {}
	for clade in tree.find_clades():
		if clade.name:
			if clade.name in names:
				raise ValueError("Duplicate key: %s" % clade.name)
			names[clade.name] = clade
	return names


def get_parent(tree, child_clade):
	node_path = tree.get_path(child_clade)
	parent = None
	try:
		parent = node_path[-2]
	except:  # Deal with nodes that have a parent not found by get_path, for some reason
		for node in lookup_by_names(tree):
			children = list(tree.find_any(node).__iter__())
			for child in children:
				if child.name == child_clade:
					new_clade = Phylo.Newick.Clade()
					new_clade.name = node
					new_clade.branch_length = tree.distance(node, child)
					parent = new_clade
	return parent

def get_label(leaf):
	return leaf.name

def pre_cleanup(mega_aln_filename):
	base_filename = os.path.splitext(os.path.basename(mega_aln_filename))[0]
	try:
		os.remove(base_filename + '.csv')
	except:
		pass
	try:
		os.remove(base_filename + '_nodeMap.txt')
	except:
		pass
	try:
		os.remove(base_filename + '.nwk')
	except:
		pass
	try:
		os.remove(base_filename + '_summary.txt')
	except:
		pass


def infer_mp_tree(mega_aln_filename):
	base_filename = os.path.splitext(os.path.basename(mega_aln_filename))[0]
	tree_filename = os.path.join(scratch_dir, base_filename + ".nwk")
	megacc_cmd = "{} -a {} -d {} -o {}".format(megacc_app, mp_tree_infer_mao, mega_aln_filename, tree_filename)
	if print_megacc_cmd: print(megacc_cmd)
	FNULL = open(os.devnull, 'w')
	return_code = subprocess.call(megacc_cmd, stdout=FNULL, stderr=subprocess.STDOUT)
	if return_code != 0:
		raise ValueError('MEGACC returned error code', return_code)
	trees = Phylo.parse(tree_filename, 'newick')
	tree = Phylo.BaseTree.Tree.from_clade(trees.__next__().clade)
	tree.root_with_outgroup({'name': 'Normal'})
	anc_id = 0
	clade_count = 0
	for clade in tree.find_clades():
		clade_count += 1
		if str(clade.name) == "None":
			clade.name = 'anc_node_' + str(anc_id)
			anc_id += 1
	anc_states_files = glob.glob(base_filename + '_ancestral_states_*.txt')
	for file in anc_states_files:
		try:
			os.remove(file)
		except:
			print("Error while deleting ancestral states file : ", file)
	changes_list_files = glob.glob(base_filename + '_changes_list_*.txt')
	for file in changes_list_files:
		try:
			os.remove(file)
		except:
			print("Error while deleting ancestral states file : ", file)
	for clade in tree.find_clades():
		if clade.branch_length is None:
			clade.branch_length = 0.0
	return tree, clade_count

def parse_nodemap_file(nodemap_filename):
	nodemap_table = []
	new_nodemap_table = []
	node_ids = {}
	leaf_names = {}
	with open(nodemap_filename, 'r') as file:
		for line in file.readlines():
			line = line.replace('anc node ', 'anc_node_')
			values = line.rstrip().split()
			if values[0] == 'NodeLabel':
				continue
			if values[0][0:5] != "Node_":
				leaf_names[values[0]] = "Node_" + values[1]
			nodemap_table.append((values[0], values[1], values[2], values[3]))
			node_ids[values[1]] = values[0]
	for node in nodemap_table:
		if node[2] in node_ids.keys() and node[3] in node_ids.keys():
			new_nodemap_table.append((node[0], node_ids[node[2]], node_ids[node[3]]))
		else:
			new_nodemap_table.append((node[0], '-', '-'))
	return new_nodemap_table, leaf_names

def map_mega_node_names(nodemap_filename, tree):
	node_label_pairs = {}
	nodemap_table, leaf_names = parse_nodemap_file(nodemap_filename)
	mega_node_names = [val[0] for val in nodemap_table]
	incomplete = True
	while incomplete:
		incomplete = False
		for node in lookup_by_names(tree):
			if node in node_label_pairs.keys():
				parent = get_parent(tree, node)
				if parent == None:
					continue
				parent_name = parent.name
				for mega_node in nodemap_table:
					if mega_node[1] == node_label_pairs[node] or mega_node[2] == node_label_pairs[node]:
						node_label_pairs[parent_name] = mega_node[0]
			else:
				incomplete = True
				#print(node)
				if node in mega_node_names:
					node_label_pairs[node] = node
	for key in node_label_pairs.keys():
		if key in leaf_names.keys():
			node_label_pairs[key] = leaf_names[key]
	return node_label_pairs

def parse_ep_outputs(basename, tree, tumor_site_map):
	temp_map = map_mega_node_names("{}eps_nodeMap.txt".format(basename), tree)
	node_label_pairs = {}
	for key in temp_map.keys():
		node_label_pairs[temp_map[key]] = key
	ep_filename = "{}eps.csv".format(basename)
	anc_snp_probs = {}
	with open(ep_filename, 'r') as file:
		for line in file.readlines()[0:len(tumor_label_dict)+1]:
			values = line.rstrip().split(',')
			if values[0] == 'site':
				header = copy.deepcopy(values)
				for node_label in values[2:]:
					anc_snp_probs.update({node_label_pairs[node_label]: {}})
				continue
			for node_label, prob in zip(header[2:], values[2:]):
				anc_snp_probs[node_label_pairs[node_label]][tumor_site_map[(values[1])]] = float(prob)
	return anc_snp_probs

def clear_ep_files(basename):
	os.remove("{}eps.nwk".format(basename))
	os.remove("{}eps.csv".format(basename))
	try:
		os.remove("{}eps_avg_blens.csv".format(basename))
		os.remove("{}eps_data_coverage.csv".format(basename))
		os.remove("{}eps_ML_data.csv".format(basename))
	except:
		pass
	os.remove("{}eps_nodeMap.txt".format(basename))
	os.remove("{}eps_summary.txt".format(basename))

def get_eps(tree, aln_filename, site_labels):
	temp_id = random.randint(100000, 999999)
	temp_basename = os.path.join(scratch_dir, "temp_{}_".format(temp_id))
	tree_filename = os.path.join(scratch_dir, "{}_anc_seqs_in.nwk".format(os.path.splitext(os.path.basename(aln_filename))[0]))
	temp_tree = copy.deepcopy(tree)
	Phylo.write(temp_tree, tree_filename, 'newick')
	anc_seqs_filename = "{}eps.csv".format(temp_basename)
	node_map_filename = "{}eps_nodeMap.txt".format(temp_basename)
	megacc_cmd = "{} -a {} -d {} -t {} --keep-tree-blens -o {} -g {}".format(megacc_app, ancestral_seqs_mao, aln_filename, tree_filename, anc_seqs_filename, outgroup_file)
	if print_megacc_cmd: print(megacc_cmd)
	if mega_io_logging:
		shutil.copy(aln_filename, os.path.join(args.output, "ancestral_inference_logging", "temp_{}_".format(temp_id) + os.path.basename(aln_filename)))
		shutil.copy(tree_filename, os.path.join(args.output, "ancestral_inference_logging", "temp_{}_".format(temp_id) + os.path.basename(tree_filename)))
	FNULL = open(os.devnull, 'w')
	return_code = subprocess.call(megacc_cmd, stdout=FNULL, stderr=subprocess.STDOUT)
	if return_code != 0:
		raise ValueError('MEGACC returned error code', return_code)
	if mega_io_logging:
		try:
			shutil.copy(anc_seqs_filename, os.path.join(args.output, "ancestral_inference_logging"))
			shutil.copy(node_map_filename, os.path.join(args.output, "ancestral_inference_logging"))
		except:
			print("Could not copy ancestral sequence inference result file to logging directory.")
	eps = parse_ep_outputs(temp_basename, tree, site_labels)
	clear_ep_files(temp_basename)
	return eps

def get_node_tumors(eps, threshold):
	new_tumor_membership = {}
	for key in eps.keys():
		new_tumor_membership[key] = {}
		for tumor in eps[key].keys():
			if eps[key][tumor] >= threshold:
				new_tumor_membership[key].update({tumor: eps[key][tumor]})
	return new_tumor_membership

def tree_to_digraph(tree, start_node = 'Normal'):
	new_nwx_graph = Phylo.to_networkx(tree)
	new_pydot_graph = networkx.drawing.nx_pydot.to_pydot(new_nwx_graph)
	for edge in new_pydot_graph.get_edge_list():
		edge.obj_dict["attributes"]["weight"] = tree.distance(edge.get_source(), edge.get_destination())
	# Remove cyclical edges
	for node in new_pydot_graph.get_nodes():
		if new_pydot_graph.del_edge(node, node):
			pass
			print("Deleted circular edge on node: {}".format(node))

	# Convert non-directed graph into directed graph by starting with normal and radiating outward
	node_list = [start_node]
	new_edge_list = []
	while len(node_list) > 0:
		new_node_list = []
		for edge in new_pydot_graph.get_edge_list():
			if edge.get_destination() in node_list:
				new_edge_list.append((edge.get_destination(), edge.get_source(), edge.obj_dict['attributes']['weight']))
				new_node_list.append(edge.get_source())
				new_pydot_graph.del_edge(edge.get_source(), edge.get_destination())
			if edge.get_source() in node_list:
				new_edge_list.append((edge.get_source(), edge.get_destination(), edge.obj_dict['attributes']['weight']))
				new_node_list.append(edge.get_destination())
				new_pydot_graph.del_edge(edge.get_source(), edge.get_destination())
		node_list = set(new_node_list)
	new_edge_list = set(new_edge_list)
	# Generate directed graph
	new_pydot_digraph = pydot.Dot(graph_type='digraph')
	for edge in new_edge_list:
		new_pydot_digraph.add_edge(pydot.Edge(edge[0], edge[1], length=edge[2]))
		#new_pydot_digraph.add_edge(pydot.Edge(edge[0], edge[1]))
	return new_pydot_digraph

def draw_tumor_map(static_tree, map_tuple, filename, timed_tree):
	tumor_map = map_tuple[0]
	tree = copy.deepcopy(static_tree)
	nodes = lookup_by_names(tree)
	timed_nodes = lookup_by_names(timed_tree)
	node_times = {}
	for key in timed_nodes.keys():
		node_times[key] = timed_nodes[key].branch_length
	for node_name in nodes.keys():
		if len(tumor_map[node_name].keys()) == 1:
			for key in tumor_map[node_name].keys():
				if tumor_map[node_name][key] >= 0.995:
					name_str = key
				else:
					name_str = "{}_{}".format(key, round(tumor_map[node_name][key], 2))
				if node_name[0:8] != "anc_node":
					name_str = name_str + "({})".format(node_name)
				if nodes[node_name].branch_length is None:
					if node_times[node_name] is not None:
						name_str = "{}\n{}".format(name_str, round(node_times[node_name],4))
				nodes[node_name].name = name_str
				#nodes[node_name].name = "{}_{}({})".format(key, tumor_map[node_name][key], node_name)
		else:
			raise ValueError("Node {} has {} tumors assigned.".format(key, len(tumor_map[node_name].keys())))
	pylab.rcParams['figure.figsize'] = (len(nodes)*args.phy_scale_x)/4.0, (len(nodes)*args.phy_scale_y)/4.0
	Phylo.draw(tree, label_func=get_label, do_show=False)
	pylab.suptitle("Probability: {}".format(map_tuple[1]))
	pylab.savefig(filename)
	pylab.close('all')

def get_gv_terminals(graph):
	edges = graph.get_edge_list()
	terminals = set()
	for edge in edges:
		terminals.add(edge.get_source())
		terminals.add(edge.get_destination())
	for edge in edges:
		if edge.get_source() in terminals:
			terminals.remove(edge.get_source())
	return terminals

def get_gv_path(graph, terminal, root):
	child = terminal
	path = []
	parents = { child:parent for (child, parent) in [ (edge.get_destination(), edge.get_source()) for edge in graph.get_edge_list() ] }
	while child != root:
		path.append((parents[child], child))
		child = parents[child]
	path.reverse()
	return path


def reverse_edge(edge_in):
	edge_tup = edge_in.split("->")
	return "{}->{}".format(edge_tup[1],edge_tup[0])


def generate_edge_list(static_tree, tumor_map):
	tree = copy.deepcopy(static_tree)
	temp_digraph = tree_to_digraph(tree, "Normal")
	temp_edge_list = temp_digraph.get_edge_list()
	temp_edge_list.sort(key=lambda tup: str(tup))
	temp_migration_edges = []
	temp_migration_edge_probabilities = []
	temp_migration_edge_lens = []

	# Group nodes into polytomies
	ungrouped_nodes = list(tumor_map.keys())
	grouped_nodes = {}
	node_groups = {}
	while len(ungrouped_nodes) > 0:
		key_node = ungrouped_nodes[0]
		node_group = []
		for node in ungrouped_nodes:
			if initial_tree.distance(key_node, node) == 0.0:
				node_group.append(node)
		for node in node_group:
			ungrouped_nodes.remove(node)
		grouped_nodes[key_node] = node_group
	for key_node in grouped_nodes.keys():
		for node in grouped_nodes[key_node]:
			node_groups[node] = key_node
	temp_node_edges = []

	for edge in temp_edge_list:
		src_node = tumor_map[edge.get_source()]
		dst_node = tumor_map[edge.get_destination()]
		if 'Normal' in [next(iter(x.keys())) for x in [src_node, dst_node]]: continue
		#migration_edge = "{}->{}".format(src_node.keys()[0], dst_node.keys()[0])
		if next(iter(src_node.keys())) != next(iter(dst_node.keys())):
			temp_migration_edges.append("{}->{}".format(next(iter(src_node.keys())), next(iter(dst_node.keys()))))
			temp_migration_edge_probabilities.append(next(iter(src_node.values())) * next(iter(dst_node.values())))
			temp_migration_edge_lens.append(edge.obj_dict["attributes"]["length"])
			temp_node_edges.append((node_groups[edge.get_source()], node_groups[edge.get_destination()], static_tree.distance(edge.get_source(), edge.get_destination())))
	temp_edge_list = list(zip(temp_migration_edges, temp_migration_edge_probabilities, temp_migration_edge_lens))

	# If a 0-length edge contained in a polytomy duplicates or reverses another edge attached to that polytomy, the 0-length edge should be dropped
	i = 0
	nonzero_edges = [record for record in zip(temp_edge_list, temp_node_edges) if record[1][2] > 0.0]
	bad_edges = []
	# for record in zip(temp_edge_list, temp_node_edges):
	# 	if record[1][2] == 0.0:
	# 		for edge in nonzero_edges:
	# 			if edge[0][0] in [record[0][0], reverse_edge(record[0][0])]:
	# 				if record[1][0] in [edge[1][0], edge[1][1]]:
	# 					bad_edges.append((temp_migration_edges[i], temp_migration_edge_probabilities[i]))
	# 					break
	# 	i += 1
	mig_count1 = len(temp_edge_list)
	source_count1 = len(set([temp_edge[0].split('->')[0] for temp_edge in temp_edge_list]))
	for bad_edge in bad_edges:
		temp_edge_list.remove(bad_edge)
	mig_count2 = len(temp_edge_list)
	source_count2 = len(set([temp_edge[0].split('->')[0] for temp_edge in temp_edge_list]))

	temp_edge_list.sort(key=lambda tup: tup[1], reverse=True)
	temp_migration_edges = []
	temp_migration_edge_probabilities = []
	temp_migration_edge_lens = []
	for edge in temp_edge_list:
		j = 1
		while "{}[{}]".format(edge[0], j) in temp_migration_edges: j += 1
		temp_migration_edges.append("{}[{}]".format(edge[0], j))
		temp_migration_edge_probabilities.append(edge[1])
		temp_migration_edge_lens.append(float(edge[2]))

	# Count comigrations
	terminals = get_gv_terminals(temp_digraph)
	mig_edge_counts = {}
	for terminal in terminals:
		path_counts = {}
		path = get_gv_path(temp_digraph, terminal, "Normal")
		for edge in path:
			mig_edge = (next(iter(tumor_map[edge[0]].keys())), next(iter(tumor_map[edge[1]].keys())))
			if mig_edge[0] != mig_edge[1]:
				path_counts[mig_edge] = path_counts.get(mig_edge, 0) + 1
		for key in path_counts.keys():
			mig_edge_counts[key] = max(mig_edge_counts.get(key, 0), path_counts[key])
	comigration_count = sum(mig_edge_counts.values())

	return temp_migration_edges, temp_migration_edge_probabilities, temp_migration_edge_lens, [comigration_count, mig_count1, mig_count2, source_count1, source_count2]

def analyze_edge_list(true_paths, edge_list):
	temp_true_paths = copy.deepcopy(true_paths)
	counts = {x: 0 for x in count_keys}
	for edge in edge_list:
		edge_type = ""
		if args.acc_by_edge_type:
			edge_type = "({})".format('->'.join(['P' if x == args.primary else 'M' for x in edge.split('[')[0].split('->')]))
			#edge_type = "({})".format('->'.join([x[0] for x in edge.split('->')]))
		if edge in temp_true_paths:
			counts["TP{}".format(edge_type)] += 1
			del temp_true_paths[temp_true_paths.index(edge)]
		else:
			counts["FP{}".format(edge_type)] += 1
	for edge in temp_true_paths:
		edge_type = ""
		if args.acc_by_edge_type:
			edge_type = "({})".format('->'.join(['P' if x == args.primary else 'M' for x in edge.split('[')[0].split('->')]))
			#edge_type = "({})".format('->'.join([x[0] for x in edge.split('->')]))
		counts["FN{}".format(edge_type)] += 1
	return counts

def make_pydot_seeding_graph(edges, probabilities, lens):
	tumors = set()
	for edge in edges:
		tumors.add(edge.split('[')[0].split('->')[0])
		tumors.add(edge.split('[')[0].split('->')[1])
	tumors = list(tumors)
	new_graph = pydot.Dot(graph_type='digraph')
	new_graph.add_node(pydot.Node(name=args.primary, rank="min"))
	for tumor in tumors:
		if tumor == args.primary: continue
		new_graph.add_node(pydot.Node(name=tumor))
	for record in zip(edges, probabilities, lens):
		edge = record[0]
		edge_color = 'grey'
		edge_label = "{}".format(round(record[2], 4))
		if record[1] == '-':
			edge_color = 'black'
		else:
			edge_label = "{}\n({})".format(round(record[1], 2), round(record[2], 4))
			if record[1] > args.mig_event_threshold:
				edge_color = 'black'
		new_graph.add_edge(pydot.Edge(edge.split('[')[0].split('->')[0], edge.split('[')[0].split('->')[1], label=edge_label, color=edge_color))
	return new_graph

def make_proto_graph(edges):
	tumors = set()
	for edge in edges:
		tumors.add(edge.split('[')[0].split('->')[0])
		tumors.add(edge.split('[')[0].split('->')[1])
	tumors = list(tumors)
	new_graph = pydot.Dot(graph_type='digraph')
	new_graph.add_node(pydot.Node(name=args.primary, rank="min"))
	for tumor in tumors:
		if tumor == args.primary: continue
		new_graph.add_node(pydot.Node(name=tumor))
	for edge in edges:
		new_graph.add_edge(pydot.Edge(edge.split('[')[0].split('->')[0], edge.split('[')[0].split('->')[1], id=edge, style="invis"))
	return new_graph

def parse_anc_seqs_times(initial_tree, anc_seqs_tree):
	anc_seqs_tree.root_with_outgroup({'name': 'Normal'})
	nodes = lookup_by_names(initial_tree)
	skipped = True
	while skipped:
		clades = anc_seqs_tree.find_clades()
		skipped = False
		for clade in clades:
			if clade.name != None:
				continue
			if clade[0].name == None or clade[1].name == None:
				skipped = True
				continue
			for key in nodes.keys():
				if len(nodes[key]) == 2 and nodes[key][0].name in [clade[1].name, clade[0].name] and nodes[key][1].name in [clade[1].name, clade[0].name]:
					clade.name = nodes[key].name
	return anc_seqs_tree

def get_avg_mig_probs(mig_file_list):
	edge_lists = []
	weight_list = []
	avg_edge_weights = {}
	for mig_filename in mig_file_list:
		with open(os.path.join(args.output, mig_filename), 'r') as mig_file:
			edge_list = []
			for line in mig_file:
				line = line.strip()
				if line[0:30] == "Tumor Membership Probability: ":
					weight_list.append(float(line[30:]))
					continue
				if line == "":
					continue
				edge_list.append(line.split("]")[0] + "]")
			edge_lists.append(edge_list)
	weight_factor = 1.0 / sum(weight_list)
	adj_weight_list = [weight_factor * weight for weight in weight_list]
	for (weight, edge_list) in zip(adj_weight_list, edge_lists):
		for edge in edge_list:
			avg_edge_weights[edge] = avg_edge_weights.get(edge, 0.0) + weight
	return avg_edge_weights

def split_mt_leaves(tree, clones, basename, mut_seqs=None): # Split multi-tumor leaf nodes in tree and generate tumor-only alignment file
	out_filename = basename + "_tumors.meg"
	tumor_seqs = {}
	with open(out_filename, 'w') as file:
		file.write("#MEGA\n")
		file.write("!Title SNVs;\n")
		file.write("!Format datatype=Protein;\n")
		for clone in clones.keys():
			if len(clones[clone].keys()) == 1:
				tumor_seqs[clone] = rev_tumor_label_dict[next(iter(clones[clone].keys()))]
			elif len(clones[clone].keys()) >= 1:
				nodes = lookup_by_names(tree)
				tumor_list = list(clones[clone].keys())
				nodes[clone].name = "{}:{}".format(clone, tumor_list[0])
				tumor_seqs["{}:{}".format(clone, tumor_list[0])] = rev_tumor_label_dict[tumor_list[0]]
				if mut_seqs is not None:
					mut_seqs["{}:{}".format(clone, tumor_list[0])] = mut_seqs[clone]
				node_to_split = "{}:{}".format(clone, tumor_list[0])
				i = 0
				for tumor in tumor_list[1:]:
					nodes = lookup_by_names(tree)
					tumor_seqs["{}:{}".format(clone, tumor)] = rev_tumor_label_dict[tumor]
					if mut_seqs is not None:
						mut_seqs["{}:{}".format(clone, tumor)] = mut_seqs[clone]
					new_root_name = "{}_root_{}".format(clone, i)
					new_clade1 = Phylo.Newick.Clade(0, node_to_split)
					new_clade2 = Phylo.Newick.Clade(0, "{}:{}".format(clone, tumor))
					nodes[node_to_split].root.clades = [new_clade1, new_clade2]
					nodes[node_to_split].name = new_root_name
					node_to_split = new_clade2.name
					i += 1
			else:
				raise ValueError("Something went wrong...")
		if args.infer_primary_seq and "Primary" in mut_seqs.keys():
			tumor_seqs["Primary"] = rev_tumor_label_dict[args.primary]
		tumor_seqs["Normal"] = rev_tumor_label_dict[args.primary]
		if mut_seqs is not None:
			for key in tumor_seqs:
				tumor_seqs[key] = tumor_seqs[key] + mut_seqs[key]
		for key in tumor_seqs:
			file.write("#{}\n".format(key))
			file.write("{}\n".format(tumor_seqs[key]))
	return out_filename


def group_polytomies(input_tree):
	tree_nodes = lookup_by_names(input_tree)
	node_groups = []
	ungrouped_nodes = list(lookup_by_names(input_tree).keys())
	polytomies = {}
	# Make polytomy node groups
	while len(ungrouped_nodes) > 0:
		key_node = ungrouped_nodes[0]
		node_group = []
		for node in ungrouped_nodes:
			if input_tree.distance(key_node, node) == 0.0:
				node_group.append(node)
		for node in node_group:
			ungrouped_nodes.remove(node)
		node_groups.append(node_group)
	ptmy_nodes = [group for group in node_groups if len(group) > 2 and 'Normal' not in group]
	for ptmy in ptmy_nodes:
		top_node = None
		for node in ptmy:
			if get_parent(input_tree, node).name not in ptmy:
				top_node = node
				break
		children = []
		for node in ptmy:
			for child in tree_nodes[node].root.clades:
				if child.name not in ptmy or child.is_terminal():
					children.append(child.name)
		polytomies[top_node] = (children, ptmy)
	return polytomies


def permute_polytomies(input_tree):
	pmt_downsample_thresh = 10000
	trees = [copy.deepcopy(input_tree)]
	input_nodes = lookup_by_names(input_tree)
	polytomies = group_polytomies(input_tree)
	depths = input_tree.depths(True)
	depths = {key.name:depths[key] for key in depths.keys()}
	for top_node in polytomies.keys():
		new_trees = []
		internal_names = [node for node in polytomies[top_node][1] if node not in polytomies[top_node][0]]
		internal_name_depths = [depths[node] for node in internal_names]
		sorted_internal_names = [x[0] for x in sorted(zip(internal_names, internal_name_depths), key=lambda row: row[1])]
		ptmy_tree_strings = permute_unique_trees(polytomies[top_node][0], 1000)
		if len(trees)*len(ptmy_tree_strings) > pmt_downsample_thresh:
			trees_sample_size = int(round(math.sqrt(pmt_downsample_thresh)))
			ptmy_strings_sample_size = int(round(math.sqrt(pmt_downsample_thresh)))
			if len(trees) < math.sqrt(pmt_downsample_thresh):
				ptmy_strings_sample_size = int(round(pmt_downsample_thresh/len(trees)))
			if len(ptmy_tree_strings) < math.sqrt(pmt_downsample_thresh):
				trees_sample_size = int(round(pmt_downsample_thresh/len(ptmy_tree_strings)))
			random.shuffle(ptmy_tree_strings)
			ptmy_tree_strings = ptmy_tree_strings[0:ptmy_strings_sample_size]
			random.shuffle(trees)
			trees = trees[0:trees_sample_size]
		for ptmy_tree_string in ptmy_tree_strings:
			handle = StringIO(str(ptmy_tree_string))
			ptmy_tree = Phylo.read(handle, "newick")
			ptmy_nodes = lookup_by_names(ptmy_tree)
			ptmy_depths = ptmy_tree.depths(True)
			ptmy_depth_sorted_internal = [x[0] for x in sorted([(key, ptmy_depths[key]) for key in ptmy_depths.keys()],
															key=lambda row: row[1]) if x[0].name is None]
			for node, name in zip(ptmy_depth_sorted_internal, sorted_internal_names):
				node.name = name
				node.branch_length = 0.0
			for edge_node in polytomies[top_node][0]:
				ptmy_nodes[edge_node].root.clades = copy.deepcopy(input_nodes[edge_node].root.clades)
				ptmy_nodes[edge_node].branch_length = input_nodes[edge_node].branch_length
			for tree in trees:
				new_tree = copy.deepcopy(tree)
				tree_nodes = lookup_by_names(new_tree)
				tree_nodes[ptmy_tree.root.name].clades = ptmy_tree.root.clades
				new_trees.append(new_tree)
		trees = copy.deepcopy(new_trees)
	return trees


def cleanup(basename):
	suffixes = ["_summary.txt","_tumors.meg","_tumors_anc_seqs_in.nwk","_tumors_anc_seqs_out.nwk",".nwk",".meg", "_data_coverage.csv"]
	for filename in [basename + x for x in suffixes]:
		try:
			if os.path.exists(os.path.join(args.output, filename)):
				os.remove(os.path.join(args.output, filename))
			shutil.move(filename, args.output)
		except:
			pass
	try:
		os.remove(basename + "_tumors_anc_seqs_out.nwk")
	except:
		pass


def make_scratch_dir(basename):
	idx = random.randint(1, 1000)
	while os.path.exists("{}_{}".format(basename, idx)):
		idx = random.randint(1, 1000)
	os.mkdir("{}_{}".format(basename, idx))
	return "{}_{}".format(basename, idx)


def draw_tumor_labeled_phylogeny(input_tree, node_tumors, filename):
	new_tree = copy.deepcopy(input_tree)
	nodes = lookup_by_names(new_tree)
	for node_name in nodes.keys():
		name_str = node_name
		if len(node_tumors[node_name].keys()) == 1:
			for key in node_tumors[node_name].keys():
				if node_tumors[node_name][key] >= 0.995:
					name_str = key
				else:
					name_str = "{}_{}".format(key, round(node_tumors[node_name][key], 2))
				if node_name[0:8] != "anc_node":
					name_str = name_str + "({})".format(node_name)
		# name_str = "{}_{}({})".format(key, node_tumors[node_name][key], node_name)
		else:
			for key in node_tumors[node_name].keys():
				name_str = name_str + "\n{}: {}".format(key, round(node_tumors[node_name][key], 2))
		nodes[node_name].name = name_str
	pylab.rcParams['figure.figsize'] = (len(nodes)*args.phy_scale_x) / 2.0, (len(nodes)*args.phy_scale_y) / 2.0
	Phylo.draw(new_tree, label_func=get_label, do_show=False)
	pylab.savefig(filename)
	pylab.close()


def fix_anc_seq_inference(tree, tumor_map, node):
	parent = get_parent(tree, node)
	try:
		node_name = node.name
	except:
		node_name = node
	if parent is not None:
		try:
			distance = tree.distance(parent.name, node)
		except:
			distance = tree.distance(parent, node)
	if parent is None or distance > 0.0:
		if "Normal" in tumor_map[node_name].keys():
			return None
		else:
			return tumor_map[node_name]
	else:
		return fix_anc_seq_inference(tree, tumor_map, parent)


def get_top_ptmy_node(tree, node):
	if not isinstance(node, str):
		node = node.name
	parent = get_parent(tree, node)
	if parent is None:
		return node
	else:
		parent = parent.name
	if tree.distance(parent, node) == 0.0:
		return get_top_ptmy_node(tree, parent)
	else:
		return node


def correct_zero_branch_eps(tree, node_membership):
	terminals = tree.get_terminals()
	nonterminals = tree.get_nonterminals()
	zero_threshold = (tree.total_branch_length() / (len(terminals) + len(nonterminals))) / 50.0
	for node in nonterminals:
		membership = {}
		colocation_count = 0.0
		for terminal in terminals:
			if tree.distance(node, terminal) < zero_threshold:
				colocation_count += 1.0
				for key in node_membership[terminal.name]:
					membership[key] = membership.get(key, 0) + node_membership[terminal.name][key]
		if colocation_count > 0.0:
			for key in membership.keys():
				membership[key] = membership[key]/colocation_count
			node_membership[node.name] = membership


class PermutedMembership:
	def __init__(self, membership, groups):
		self.membership = membership
		self.groups = groups
		self.group_anchors = list(self.groups.keys())
		self.anchor_membership_sizes = {x: len(self.membership[x]) for x in self.group_anchors}
		self.anchor_membership = {x: list(self.membership[x].keys()) for x in self.group_anchors}
		# self.group_sizes = [len(self.groups[x]) for x in self.group_anchors]
		# self.idx_matrix = [list(self.membership[x].keys()) for x in self.membership.keys()]
		# self.idx_matrix_sizes = [len(x) for x in idx_matrix]

	def perm_count(self):
		pcount = 1
		for val in self.anchor_membership_sizes.values():
			pcount = pcount * val
		return pcount

	def get_perm_by_idx(self, idx):
		new_tumor_map = (copy.deepcopy(self.membership), 1.0)
		idx_val = idx
		for anchor_node in self.group_anchors:
			tumor_idx = idx_val % self.anchor_membership_sizes[anchor_node]
			idx_val = int((idx_val - tumor_idx)/self.anchor_membership_sizes[anchor_node])
			tumor = self.anchor_membership[anchor_node][tumor_idx]
			new_tumor_map = (new_tumor_map[0], new_tumor_map[1] * self.membership[anchor_node][tumor])
			for ptmy_node in self.groups[anchor_node]:
				new_tumor_map[0][ptmy_node] = {tumor: self.membership[anchor_node][tumor]}
		return new_tumor_map


def derive_mut_scale(tree, seqs):
	keys = list(seqs.keys())
	keys.remove('Normal')
	min_distance = None
	while len(keys) > 0:
		distance = 0.0
		key = keys.pop()
		for position in zip(seqs[key], seqs['Normal']):
			if position[0].upper() != position[1].upper():
				distance += 1
		if min_distance is None or min_distance > distance:
			if distance > 0:
				min_distance = distance
				mut_scale = distance / tree.distance(key, 'Normal')
	return mut_scale


scratch_dir = make_scratch_dir(os.path.join(args.output, "scratch"))

read_true_paths = False

if args.true_paths is not None:
	read_true_paths = True
	temp_paths = read_paths(args.true_paths)
	true_paths = []
	for path in temp_paths:
		path = "{}->{}".format(path[0], path[1])
		i = 1
		while "{}[{}]".format(path, i) in true_paths: i += 1
		true_paths.append("{}[{}]".format(path, i))

start_time = datetime.datetime.now()

aln_file_out = os.path.join(scratch_dir, os.path.splitext(os.path.basename(args.aln))[0] + ".meg")

mut_seqs = parse_input_aln(args.aln, aln_file_out)
clones = parse_clone_freqs(args.clone_locations)

if "Primary" in clones.keys() and "Primary" not in mut_seqs.keys():
	del clones["Primary"]

mega_aln_filename = aln_file_out

initial_tree, node_count = infer_mp_tree(mega_aln_filename)

mut_scale = derive_mut_scale(initial_tree, mut_seqs)

for key in mut_seqs.keys():
	mut_seqs[key] = ''
mega_aln_filename = split_mt_leaves(initial_tree, clones, os.path.join(scratch_dir, os.path.splitext(os.path.basename(args.aln))[0]), mut_seqs)


# detect polytomies in initial tree
# for each polytomy, permute each unique bifurcating resolution
# make list of trees as cartesian product of all possibilities for each polytomy
# for every tree, run ancestral seqs, parse eps, generate all possible path sets
# at this point, you've basically got a hierarchical set of results, structured like:
#	the base tree (i.e. what was input to the ancestral seqs analysis)
#	the tumor probabilities for each node from anc_seqs analysis
#		the combined probability of one possible set of tumor assignments
#		the set of clone tumor assignments itself
#		the set of migration edges produced by this set of tumor assignments

permuted_trees = permute_polytomies(initial_tree)

if len(permuted_trees) > args.max_permutations:
	print("Found {} possible permutations of initial tree, using a sample of size {}.".format(len(permuted_trees), args.max_permutations))
	random.shuffle(permuted_trees)
	permuted_trees = permuted_trees[0:args.max_permutations]
else:
	print("Found {} possible permutations of initial tree.".format(len(permuted_trees)))

data_by_tree = {} # Structured as {tree_permutation: (tumor_probabilities, [(joint_probability, node_tumor_membership, migration_edges), ...])}

update_interval = 10
update_user = False
last_update_time = datetime.datetime.now()
trees_processed = 0.0

for pmt_tree in permuted_trees:
	eps = get_eps(pmt_tree, mega_aln_filename, tumor_label_dict)
#	for clade in pmt_tree.find_clades():
#		clade.branch_length = round(clade.branch_length * mut_scale)
	digraph = tree_to_digraph(pmt_tree)
	tumor_membership = get_node_tumors(eps, tumor_membership_cutoff)
	bad_result = False
	nodes_below_threshold = []
	# min_max_ep = maximum ep of node with lowest maximum ep below threshold
	min_max_ep = 1.0
	for key in tumor_membership.keys():
		if len(tumor_membership[key]) == 0:
			# print("Warning: tree permutation yields node({}) with maximum ancestral probability({}) below threshold value({}), skipping...".format(key, max(eps[key].values()),tumor_membership_cutoff))
			nodes_below_threshold.append(key)
			min_max_ep = min(min_max_ep, max(eps[key].values()))
			bad_result = True
			# break
			if args.relax_threshold:
				multi_select = True
				max_select = not multi_select
				if max_select:
					max_ep = max(eps[key].values())
					if len([val for val in eps[key].values() if val >= max_ep]) == 1:
						selected = [val for val in eps[key].keys() if eps[key][val] >= max_ep][0]
						if len(selected) == 0:
							print("Warning: tree permutation yields invalid ancestral sequence inference results, skipping...")
							bad_result = True
							break
						tumor_membership[key][selected] = max_ep
						print("No sites exceed probability threshold {} for node {}, found exactly one site {} with maximum probability {}, selecting it.".format(tumor_membership_cutoff, key, selected, max_ep))
						bad_result = False
				if multi_select:
					min_ep = min(eps[key].values())
					if len([val for val in eps[key].values() if val > min_ep]) <= 19 - len(aa_label_list):
						selected = [val for val in eps[key].keys() if eps[key][val] > min_ep]
						if len(selected) == 0 and args.keep_ambiguous_results:
							selected = [val for val in tumor_label_dict.values() if val != "Normal"]
						if len(selected) == 0:
							print("No sites exceed probability threshold {} for node {}, found {} sites with probability greater than minimum value({}).".format(
									tumor_membership_cutoff, key, len(selected), min_ep))
							print("Warning: tree permutation yields invalid ancestral sequence inference results, skipping...")
							bad_result = True
							break
						total_prob = sum([eps[key][val] for val in selected])
						for val in selected:
							tumor_membership[key][val] = eps[key][val] / total_prob
						print("No sites exceed probability threshold {} for node {}, found {} sites with probability greater than minimum value({}), selecting them and normalizing probability.".format(tumor_membership_cutoff, key, len(selected), min_ep))
						bad_result = False
		if "Normal" in tumor_membership[key].keys():
			if tryto_fix_anc_seq_inference:
				tumor_membership[key] = fix_anc_seq_inference(pmt_tree, tumor_membership, key)
				if tumor_membership[key] is None:
					print("Warning: tree permutation yields invalid ancestral sequence inference results, skipping...")
					bad_result = True
					break
			else:
				print("Warning: tree permutation yields invalid ancestral sequence inference results, skipping...")
				bad_result = True
				break
	if bad_result:
		if min_max_ep < tumor_membership_cutoff:
			print("Warning: tree permutation yields nodes({}) with maximum ancestral probability({}) below threshold value({}), skipping...".format(",".join(nodes_below_threshold), min_max_ep, tumor_membership_cutoff))
		continue
	mt_nodes = []
	for key in tumor_membership.keys():
		if len(tumor_membership[key].keys()) != 1:
			mt_nodes.append(key)
	if fix_0length_issues:
		correct_zero_branch_eps(pmt_tree, tumor_membership)
	tumor_map_list = [(copy.deepcopy(tumor_membership), 1.0)]
	# Group nodes in polytomies
	ungrouped_nodes = [x for x in lookup_by_names(pmt_tree).keys() if x not in [y.name for y in pmt_tree.get_terminals()]]
	grouped_nodes = {}
	while len(ungrouped_nodes) > 0:
		key_node = ungrouped_nodes[0]
		node_group = []
		for node in ungrouped_nodes:
			if pmt_tree.distance(key_node, node) == 0.0:
				node_group.append(node)
		for node in node_group:
			ungrouped_nodes.remove(node)
		top_node = get_top_ptmy_node(pmt_tree, key_node)
		if not isinstance(top_node, str):
			top_node = top_node.name
		grouped_nodes[top_node] = node_group
	node_membership_permutation_count = 1
	for key_node in grouped_nodes.keys():
		node_membership_permutation_count = node_membership_permutation_count * len(tumor_membership[key_node])
	if v_level > 2:
		print("grouped_nodes len: {}".format(len(grouped_nodes)))
		print("grouped_nodes key membership counts: {}".format([len(tumor_membership[key_node]) for key_node in grouped_nodes.keys()]))
		print("Total grouped node membership combinations: {}".format(node_membership_permutation_count))
	perm_obj = PermutedMembership(tumor_membership, grouped_nodes)
	# Generate sample set from all possible combinations of multi tumor node assignments with probabilities
	tumor_map_list = []
	if args.max_graphs_per_tree >= 1:
		sample_size = min(perm_obj.perm_count(), args.max_graphs_per_tree)
	else:
		sample_size = perm_obj.perm_count()
	if v_level > 1:
		print("Generating sample of {} out of a possible {} tumor resolution graphs...".format(sample_size, perm_obj.perm_count()))
	sample_idxs = random.sample(list(range(0, perm_obj.perm_count())), sample_size)
	for i in range(0, len(sample_idxs)):
		tumor_map_list.append(perm_obj.get_perm_by_idx(sample_idxs[i]))
	# # Generate list of all possible combinations of multi tumor node assignments with probabilities
	# for key_node in grouped_nodes.keys():
	# 	temp_tumor_map_list = []
	# 	if v_level > 3:
	# 		print("tumor_map_list len: {}".format(len(tumor_map_list)))
	# 	# if len(tumor_map_list) > 1000:
	# 	# 	continue
	# 	for tumor_map in tumor_map_list:
	# 		for key in tumor_membership[key_node].keys():
	# 			new_tumor_map = (copy.deepcopy(tumor_map[0]), tumor_map[1] * tumor_membership[key_node][key])
	# 			for polytomy_node in grouped_nodes[key_node]:
	# 				new_tumor_map[0][polytomy_node] = {key: tumor_membership[key_node][key]}
	# 			temp_tumor_map_list.append(new_tumor_map)
	# 	tumor_map_list = temp_tumor_map_list
	if v_level > 1:
		print("Processing sample of {} tumor resolution graphs...".format(sample_size))
	edge_lists = []
	edge_lists_probabilities = []
	edge_lists_lens = []
	graphic_file_list = []
	migration_counts = []
	graphs_processed = 0.0
	for record in tumor_map_list:
		graphs_processed += 1
		if (datetime.datetime.now() - last_update_time).seconds > update_interval:
			est_completion = (trees_processed + (graphs_processed/len(tumor_map_list))) / float(len(permuted_trees))
			print("Finished processing {}% of possible graphs...".format(round(est_completion * 100.0,1)))
			last_update_time = datetime.datetime.now()
		edge_list, edge_list_probabilities, edge_list_lens, counts = generate_edge_list(pmt_tree, record[0])
		edge_lists.append(edge_list)
		edge_lists_probabilities.append(edge_list_probabilities)
		edge_lists_lens.append([round(x * mut_scale) for x in edge_list_lens])
		migration_counts.append(counts)
	trees_processed += 1.0
	temp_sorted = sorted(list(zip(tumor_map_list, edge_lists, edge_lists_probabilities, edge_lists_lens, migration_counts)),
						 key=lambda row: row[2], reverse=True)
	tumor_map_list, edge_lists, edge_lists_probabilities, edge_lists_lens, migration_counts = list(zip(*temp_sorted))
	data_by_tree[pmt_tree] = (tumor_membership, list(zip([x[1] for x in tumor_map_list], [x[0] for x in tumor_map_list], edge_lists, migration_counts, edge_lists_lens)))

print("Finished processing graphs, generating outputs...")

if len(data_by_tree) == 0:
	tree_idx = 0
	for tree in permuted_trees:
		tree_idx += 1
		Phylo.write(tree, os.path.join(scratch_dir, os.path.splitext(os.path.basename(args.aln))[0] + "_tree_{}.nwk".format(tree_idx)), 'newick')
	raise Exception("None of the tree permutations sampled produced usable output, try lowering anc_tumor_threshold({}), or enabling --relax_threshold option.".format(tumor_membership_cutoff))

key_trees = list(data_by_tree.keys())

tree_idx = 0
drawn_graphs = set()

for tree in key_trees:
	tree_idx += 1
	if args.draw_all_outputs:
		draw_tumor_labeled_phylogeny(tree, data_by_tree[tree][0], os.path.join(scratch_dir, os.path.splitext(os.path.basename(args.aln))[0] + "_tree_{}.png".format(tree_idx)))
		anc_states_idx = 0
		for record in data_by_tree[tree][1]:
			graph_string = ";".join(sorted(record[2]))
			if graph_string in drawn_graphs:
				continue
			drawn_graphs.add(graph_string)
			anc_states_idx += 1
			seeding_graph = make_pydot_seeding_graph(record[2], ['-' for x in record[2]], record[4])
			seeding_graph.set('label', "Probability: {}".format(record[0]))
			seeding_graph.set('labelloc', 't')
			try:
				seeding_graph.write_png(os.path.join(scratch_dir, os.path.splitext(os.path.basename(args.aln))[0] + "_{}_{}_migration.png".format(tree_idx, anc_states_idx)))
			except OSError:
				stuff = sys.exc_info()
				if stuff[1].strerror in ['"dot.exe" not found in path.', '"dot" not found in path.']:
					seeding_graph.write(os.path.join(scratch_dir, os.path.splitext(os.path.basename(args.aln))[0] + "_{}_{}_migration.png".format(tree_idx, anc_states_idx)))
					print("Warning: 'dot.exe' not found in path, migration paths graphic not generated, see README file...")
				else:
					raise (stuff[0], stuff[1], stuff[2])
	Phylo.write(tree, os.path.join(scratch_dir, os.path.splitext(os.path.basename(args.aln))[0] + "_tree_{}.nwk".format(tree_idx)), 'newick')

composite_weighted_edges = {}
composite_weighted_edge_lens = {}

with open(os.path.join(args.output, os.path.splitext(os.path.basename(args.aln))[0] + "_all_output_counts.txt"), 'w') as output_accs_file:
	fields = ["dataset", "probability", "migration_count", "comigration_count", "source_count", "tree_idx"]
	if read_true_paths:
		fields = fields + count_keys
	fields = fields + ["paths"]
	output_accs_file.write("{}\n".format("\t".join(fields)))
	if args.use_all_weighted_outputs:
		tree_idx = 0
		for tree in key_trees:
			tree_idx += 1
			pmt_weighted_edges = {}
			pmt_weighted_edge_lens = {}
			probability_sum = sum([record[0] for record in data_by_tree[tree][1]])
			for record in data_by_tree[tree][1]:
				# format of record is (joint_probability, node_tumor_membership, migration_edges, migration_counts)
				accs_file_record = [os.path.splitext(os.path.basename(args.aln))[0], record[0], record[3][2], record[3][0], record[3][4], tree_idx]
				for edge, edge_len in zip(record[2], record[4]):
					pmt_weighted_edges[edge] = pmt_weighted_edges.get(edge, 0) + (record[0] / probability_sum)
					pmt_weighted_edge_lens[edge] = pmt_weighted_edge_lens.get(edge, 0) + ((record[0] * edge_len) / probability_sum)
				if read_true_paths:
					acc_counts = analyze_edge_list(true_paths, record[2])
					accs_file_record = accs_file_record + [acc_counts[x] for x in count_keys]
				accs_file_record = accs_file_record + [";".join(sorted(record[2]))]
				output_accs_file.write("{}\n".format("\t".join([str(x) for x in accs_file_record])))
			for edge in pmt_weighted_edges.keys():
				composite_weighted_edges[edge] = composite_weighted_edges.get(edge, 0) + pmt_weighted_edges[edge]
				composite_weighted_edge_lens[edge] = composite_weighted_edge_lens.get(edge, 0) + pmt_weighted_edge_lens[edge]
		for edge in composite_weighted_edges.keys():
			composite_weighted_edges[edge] = composite_weighted_edges[edge] / len(data_by_tree)
			composite_weighted_edge_lens[edge] = composite_weighted_edge_lens[edge] / len(data_by_tree)

	else:
		flat_results = [] # list of tuples with contents: (tree, tumor_map, probability, migration_count, comigration_count, source_count, edge_list)
		tree_idx = 0
		for tree in key_trees:
			tree_idx += 1
			for record in data_by_tree[tree][1]:
				flat_results.append((tree, record[1], record[0], record[3][2], record[3][0], record[3][4], record[2], record[4]))
				accs_file_record = [os.path.splitext(os.path.basename(args.aln))[0], record[0], record[3][2], record[3][0], record[3][4], tree_idx]
				if read_true_paths:
					acc_counts = analyze_edge_list(true_paths, record[2])
					accs_file_record = accs_file_record + [acc_counts[x] for x in count_keys]
				accs_file_record = accs_file_record + [";".join(sorted(record[2]))]
				output_accs_file.write("{}\n".format("\t".join([str(x) for x in accs_file_record])))

		mig_min = min([x[3] for x in flat_results])
		comig_min = min([x[4] for x in flat_results if x[3] == mig_min])
		src_min = min([x[5] for x in flat_results if x[3] == mig_min and x[4] == comig_min])
		prob_max = max([x[2] for x in flat_results if x[3] == mig_min and x[4] == comig_min and x[5] == src_min])
		print("Minimum migration count:{}   Results with minimum migration count: {}".format(mig_min, len([x[4] for x in flat_results if x[3] == mig_min])))
		print("Minimum comigration count:{}   Results with minimum comigration count: {}".format(comig_min, len([x[5] for x in flat_results if x[3] == mig_min and x[4] == comig_min])))
		print("Minimum source count:{}   Results with minimum source count: {}".format(src_min, len([x[2] for x in flat_results if x[3] == mig_min and x[4] == comig_min and x[5] == src_min])))
		print("Average probability of filtered result sets:{}".format(float(sum([x[2] for x in flat_results if x[3] == mig_min and x[4] == comig_min and x[5] == src_min]))/len([x[2] for x in flat_results if x[3] == mig_min and x[4] == comig_min and x[5] == src_min])))
		if args.use_select_weighted_outputs:
			picked_results = [x for x in flat_results if x[3] == mig_min and x[4] == comig_min and x[5] == src_min]
			temp_total_prob = 0
			for result in picked_results:
				temp_total_prob += result[2]
				for edge, edge_len in zip(result[6], result[7]):
					composite_weighted_edges[edge] = composite_weighted_edges.get(edge, 0.0) + result[2]
					composite_weighted_edge_lens[edge] = composite_weighted_edge_lens.get(edge, 0.0) + (edge_len * result[2])
			for edge in composite_weighted_edges.keys():
				composite_weighted_edges[edge] = composite_weighted_edges[edge] / temp_total_prob
				composite_weighted_edge_lens[edge] = composite_weighted_edge_lens[edge] / temp_total_prob
		else:
			picked_results = [x for x in flat_results if x[3] == mig_min and x[4] == comig_min and x[5] == src_min and x[2] == prob_max]
			print("Maximum probability of filtered result sets:{}   Results with maximum probability:{}".format(prob_max, len(picked_results)))
			for result in picked_results:
				for edge, edge_len in zip(result[6], result[7]):
					composite_weighted_edges[edge] = composite_weighted_edges.get(edge, 0.0) + 1.0
					composite_weighted_edge_lens[edge] = composite_weighted_edge_lens.get(edge, 0.0) + edge_len
			for edge in composite_weighted_edges.keys():
				composite_weighted_edges[edge] = composite_weighted_edges[edge] / len(picked_results)
				composite_weighted_edge_lens[edge] = composite_weighted_edge_lens[edge] / len(picked_results)

edge_list = sorted(composite_weighted_edges.keys())

with open(os.path.join(args.output, os.path.splitext(os.path.basename(args.aln))[0] + "_Mig.txt"), 'w') as file:
	file.write("Edge\tProbability\tMutation_Count\n")
	for edge in edge_list:
		file.write("{}\t{}\t{}\n".format(edge, round(composite_weighted_edges[edge], 4), round(composite_weighted_edge_lens[edge], 4)))

seeding_graph = make_pydot_seeding_graph(edge_list, [composite_weighted_edges[x] for x in edge_list], [composite_weighted_edge_lens[x] for x in edge_list])

try:
	seeding_graph.write_png(
		os.path.join(args.output, os.path.splitext(os.path.basename(args.aln))[0] + "_migration.png"))
except OSError:
	stuff = sys.exc_info()
	if stuff[1].strerror in ['"dot.exe" not found in path.', '"dot" not found in path.']:
		seeding_graph.write(os.path.join(args.output, os.path.splitext(os.path.basename(args.aln))[0] + "_migration.dot"))
		print("Warning: 'dot.exe' not found in path, migration paths graphic not generated, see README file...")
	else:
		raise (stuff[0], stuff[1], stuff[2])

if read_true_paths:
	accuracy_list = []
	all_edges = set()
	for edge in true_paths:
		all_edges.add(edge)
	for edge in composite_weighted_edges.keys():
		all_edges.add(edge)
	new_graph = make_proto_graph(all_edges)
	edge_probabilities = composite_weighted_edges
	temp_edge_list = list(composite_weighted_edges.keys())
	with open(os.path.join(args.output, os.path.splitext(os.path.basename(args.aln))[0] + "_Mig.txt"), 'w') as avg_mig_file:
		if dump_edges: edges_file = open(edge_dump_file,'a+')
		dataname = os.path.splitext(os.path.basename(args.aln))[0]
		avg_mig_file.write("Edge\tProbability\tClassification\n")
		for graph_edge in new_graph.get_edge_list():
			edge = graph_edge.obj_dict['attributes']['id']
			if edge in edge_probabilities.keys():
				graph_edge.obj_dict['attributes']['label'] = round(edge_probabilities[edge], 2)
			else:
				graph_edge.obj_dict['attributes']['label'] = ""
			if edge in true_paths and edge in temp_edge_list:	# TP edge
				if dump_edges: edges_file.write("{}\t{}\t{}\t{}\n".format(dataname, edge, edge_probabilities.get(edge, 0.0), "True"))
				graph_edge.obj_dict['attributes']['style'] = ""
				if edge_probabilities[edge] > 0.9:
					graph_edge.obj_dict['attributes']['color'] = "forestgreen"
					avg_mig_file.write("{}\t{}\t{}\n".format(edge, edge_probabilities[edge], "True Positive"))
				elif edge_probabilities[edge] > 0.5:
					graph_edge.obj_dict['attributes']['color'] = "chartreuse"
					avg_mig_file.write("{}\t{}\t{}\n".format(edge, edge_probabilities[edge], "True Positive"))
				elif edge_probabilities[edge] >= avg_edge_weight_cutoff:
					graph_edge.obj_dict['attributes']['color'] = "yellow"
					avg_mig_file.write("{}\t{}\t{}\n".format(edge, edge_probabilities[edge], "True Positive"))
				else:
					graph_edge.obj_dict['attributes']['color'] = "blue"
					avg_mig_file.write("{}\t{}\t{}\n".format(edge, edge_probabilities[edge], "False Negative"))
			elif edge in true_paths:							# FN edge
				if dump_edges: edges_file.write("{}\t{}\t{}\t{}\n".format(dataname, edge, edge_probabilities.get(edge, 0.0), "True"))
				graph_edge.obj_dict['attributes']['style'] = ""
				graph_edge.obj_dict['attributes']['color'] = "blue"
				avg_mig_file.write("{}\t{}\t{}\n".format(edge, edge_probabilities.get(edge, 0.0), "False Negative"))
			elif edge in temp_edge_list:						# FP edge
				if dump_edges: edges_file.write("{}\t{}\t{}\t{}\n".format(dataname, edge, edge_probabilities.get(edge, 0.0), "False"))
				if edge_probabilities[edge] >= avg_edge_weight_cutoff:
					graph_edge.obj_dict['attributes']['style'] = ""
					graph_edge.obj_dict['attributes']['color'] = "red"
					avg_mig_file.write("{}\t{}\t{}\n".format(edge, edge_probabilities[edge], "False Positive"))
				else:
					graph_edge.obj_dict['attributes']['style'] = ""
					graph_edge.obj_dict['attributes']['color'] = "grey"
					avg_mig_file.write("{}\t{}\t{}\n".format(edge, edge_probabilities[edge], "True Negative"))
		if dump_edges: edges_file.close()
	new_graph.set('label', "Composite Migration Graph")
	new_graph.set('labelloc', 't')
	try:
		new_graph.write_png(os.path.join(args.output, os.path.splitext(os.path.basename(args.aln))[0] + "_migration.png"))
	except OSError:
		stuff = sys.exc_info()
		if stuff[1].strerror in ['"dot.exe" not found in path.', '"dot" not found in path.']:
			new_graph.write(
				os.path.join(args.output, os.path.splitext(os.path.basename(args.aln))[0] + "_migration.dot"))
			print("Warning: 'dot.exe' not found in path, migration paths graphic not generated, see README file...")
		else:
			raise (stuff[0], stuff[1], stuff[2])

	avg_acc_counts = analyze_edge_list(true_paths, [x for x in composite_weighted_edges.keys() if composite_weighted_edges[x] >= avg_edge_weight_cutoff])
	print("{}\t{}".format(os.path.splitext(os.path.basename(args.aln))[0], '\t'.join([str(avg_acc_counts[x]) for x in count_keys])))

if args.output != ".":
	if not os.path.exists(os.path.join(args.output, "scratch")):
		shutil.move(scratch_dir, os.path.join(args.output, "scratch"))
	else:
		for file in os.listdir(scratch_dir):
			if os.path.exists(os.path.join(args.output, "scratch", os.path.basename(file))):
				os.remove(os.path.join(args.output, "scratch", os.path.basename(file)))
			shutil.move(os.path.join(scratch_dir, file), os.path.join(args.output, "scratch"))
		os.rmdir(scratch_dir)

print("runtime: {}".format(datetime.datetime.now() - start_time))

