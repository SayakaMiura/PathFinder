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
from Bio import Phylo
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

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

mp_tree_infer_mao = "infer_NJ_amino_acid.mao"
ancestral_seqs_mao = "ancestral_seqs_ML_protein.mao"
#megacc_app = "megacc_200327.exe"
megacc_app = "megacc_200324-2.exe"
#megacc_app = "megacc_bak.exe"
outgroup_file = "outgroup.txt"
specify_outgroup = False
anc_seqs_input_topology_only = False

if specify_outgroup:
	megacc_app = "megacc_200420.exe"

print_megacc_cmd = True

aa_label_list = ['C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
tumor_label_dict = {'A': 'Normal'}
rev_tumor_label_dict = {'Normal': 'A'}

parser = argparse.ArgumentParser(description="PathFinder tumor migration path solver.")
parser.add_argument("aln", help="Clone sequence alignment file.", type=str)
parser.add_argument("clone_locations", help="Per tumor clone frequency/presence matrix file.", type=str)
parser.add_argument("-t", "--true_paths", help="List of true migration paths.", type=str, default=None)
parser.add_argument("--anc_tumor_threshold", help="Lower limit on ancestral node tumor threshold.", type=float, default=0.15)
parser.add_argument("--mig_event_threshold", help="Migration edge detection threshold.", type=float, default=0.5)
parser.add_argument("-mp", help="Use {argument}th maximum parsimony tree instead of neighbor-joining tree.", type=int, default=None)
parser.add_argument("--machina_inputs", help="Generate mutation tree and label files suitable for MACHINA input.", action='store_true', default=False)
parser.add_argument("--edges_file", help="File to dump/append edges to.", type=str, default=None)
parser.add_argument("--infer_primary_seq", help="Infer primary clone sequence based on Normal and all clone sequences.", action='store_true', default=False)
parser.add_argument("--default_normal_char", help="Infer primary clone sequence based on Normal and all clone sequences.", type=str, default=None)
parser.add_argument("-o", "--output", help="Output directory to put results in.", type=str, default=None)


args = parser.parse_args()

if args.output is not None:
	try:
		os.mkdir(args.output)
	except:
		pass
else:
	args.output = "."

dump_edges = True
if args.edges_file is None: dump_edges = False

tumor_membership_cutoff = args.anc_tumor_threshold
avg_edge_weight_cutoff = args.mig_event_threshold
if args.mp is not None:
	mp_tree_infer_mao = "infer_MP_protein.mao"


def reverse_edge(edge_in):
	edge_tup = edge_in.split("->")
	return "{}->{}".format(edge_tup[1],edge_tup[0])


def parse_clone_freqs(file_name):
	# Read tumor/clone frequencies file into dict of dicts
	clones = {}
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
					tumor = 'P'
				if len(aa_label_list) < 1:
					raise Exception("More than the maximum number of tumor sites (19) specified")
				tumor_label_dict[aa_label_list[0]] = tumor
				rev_tumor_label_dict[tumor] = aa_label_list[0]
				aa_label_list.remove(aa_label_list[0])
				freqs = data[1:]
				for i in range(len(freqs)):
					if float(freqs[i]) == 0.0:
						pass
					else:
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
#	rev_tumor_label_dict = {'Normal': 'A', 'P': 'C', 'M1': 'D', 'M2': 'E', 'M3': 'F', 'M4': 'G', 'M5': 'H', 'M6': 'I',
#							'M7': 'K', 'M8': 'L', 'M9': 'M', 'M10': 'N', 'M11': 'P', 'M12': 'Q', 'M13': 'R', 'M14': 'S',
#							'M15': 'T', 'M16': 'V', 'M17': 'W', 'M18': 'Y'}
	# clones = parse_clone_freqs(clone_sites_file)
#	for key in clones.keys():
#		clones['#' + key] = clones[key]
#		del (clones[key])
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
			print("No Normal sequence found, using user specified default char {} * {}".format(args.default_normal_char,len(seqs.values()[0])))
			seqs['Normal'] = 'A' * len(seqs.values()[0])
		seq_len = len(seqs[seqs.keys()[0]])
		target_seq_len = 1000
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
		if args.infer_primary_seq:
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

def strip_tumor_chars(mega_aln_filename):
	stripped_filename = mega_aln_filename.replace(".meg", "_no_tumors.meg")
	with open(stripped_filename, 'w') as out_file:
		with open(mega_aln_filename, 'r') as in_file:
			for line in in_file:
				line = line.strip()
				if line[0] in ['!', '#']:
					out_file.write(line + "\n")
				else:
					out_file.write(line[1:] + "\n")
	return stripped_filename

def extract_tumor_chars(mega_aln_filename):
	stripped_filename = mega_aln_filename.replace(".meg", "_tumors_only.meg")
	with open(stripped_filename, 'w') as out_file:
		with open(mega_aln_filename, 'r') as in_file:
			for line in in_file:
				line = line.strip()
				if line[0] in ['!', '#']:
					out_file.write(line + "\n")
				else:
					out_file.write(line[0:1] + "\n")
	return stripped_filename

def infer_mp_tree(mega_aln_filename):
	base_filename = os.path.splitext(os.path.basename(mega_aln_filename))[0]
	tree_filename = base_filename + ".nwk"
	# stripped_filename = strip_tumor_chars(mega_aln_filename)
	megacc_cmd = "{} -a {} -d {} -o {}".format(megacc_app, mp_tree_infer_mao, mega_aln_filename, tree_filename)
	if print_megacc_cmd: print(megacc_cmd)
	FNULL = open(os.devnull, 'w')
	return_code = subprocess.call(megacc_cmd, stdout=FNULL, stderr=subprocess.STDOUT)
	if return_code != 0:
		raise ValueError('MEGACC returned error code', return_code)
	#tree = Phylo.BaseTree.Tree.from_clade(Phylo.parse(tree_filename, 'newick').next().clade)
	trees = Phylo.parse(tree_filename, 'newick')
	if args.mp is not None:
		for i in range(1, args.mp):
			trees.next()
		print("Using maximum parsimony tree #{}...".format(args.mp))
	tree = Phylo.BaseTree.Tree.from_clade(trees.next().clade)
	#print("Number of phylogenies found: {}".format(len(testobj)))
	#Phylo.write(tree, "pre_root.nwk", 'newick')
	tree.root_with_outgroup({'name': 'Normal'})
	#Phylo.write(tree, "post_root.nwk", 'newick')
	#print(tree)
	# Phylo.write(tree, tree_filename, 'newick')
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
#	for clade in tree.find_clades():
#		clade.branch_length = None
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
		# print(tree)
		# print(nodemap_table)
		# print(node_label_pairs)
		# print(mega_node_names)
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
	tumor_count = len(tumor_site_map)
	anc_snp_probs = {}
	#print(tumor_site_map)
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
	#os.remove("{}eps.nwk".format(basename))
	os.remove("{}eps.csv".format(basename))
	os.remove("{}eps_nodeMap.txt".format(basename))
	os.remove("{}eps_summary.txt".format(basename))

def get_eps(tree, aln_filename, site_labels):
	temp_basename = "temp_{}_".format(random.randint(100000, 999999))
	#tree_filename = "{}tree.nwk".format(temp_basename)
	tree_filename = "{}_anc_seqs_in.nwk".format(os.path.splitext(os.path.basename(aln_filename))[0])
	temp_tree = copy.deepcopy(tree)
	nodes = lookup_by_names(temp_tree)
	if anc_seqs_input_topology_only:
		for key in nodes.keys():
			nodes[key].branch_length = None
	Phylo.write(temp_tree, tree_filename, 'newick')
	anc_seqs_filename = "{}eps.csv".format(temp_basename)
	megacc_cmd = "{} -a {} -d {} -t {} -o {}".format(megacc_app, ancestral_seqs_mao, aln_filename, tree_filename, anc_seqs_filename)
	if specify_outgroup:
		megacc_cmd = megacc_cmd + " -g {}".format(outgroup_file)
	if print_megacc_cmd: print(megacc_cmd)
	FNULL = open(os.devnull, 'w')
	return_code = subprocess.call(megacc_cmd, stdout=FNULL, stderr=subprocess.STDOUT)
	if return_code != 0:
		raise ValueError('MEGACC returned error code', return_code)
	eps = parse_ep_outputs(temp_basename, tree, site_labels)
	#os.remove(tree_filename)
	clear_ep_files(temp_basename)
	shutil.move("{}eps.nwk".format(temp_basename), "{}_anc_seqs_out.nwk".format(os.path.splitext(os.path.basename(aln_filename))[0]))
	return eps

def get_node_tumors_max(eps, threshold):
	new_tumor_membership = {}
	for key in eps.keys():
		max_ep = max(eps[key].values())
		if max_ep < threshold:
			continue
		max_ep_tumors = []
		for tumor in eps[key].keys():
			if eps[key][tumor] == max_ep:
				max_ep_tumors.append(tumor)
		new_tumor_membership[key] = (max_ep_tumors, max_ep)
	return new_tumor_membership

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
				new_edge_list.append((edge.get_destination(), edge.get_source()))
				new_node_list.append(edge.get_source())
				new_pydot_graph.del_edge(edge.get_source(), edge.get_destination())
			if edge.get_source() in node_list:
				new_edge_list.append((edge.get_source(), edge.get_destination()))
				new_node_list.append(edge.get_destination())
				new_pydot_graph.del_edge(edge.get_source(), edge.get_destination())
		node_list = set(new_node_list)
	new_edge_list = set(new_edge_list)

	# Generate directed graph
	new_pydot_digraph = pydot.Dot(graph_type='digraph')
	for edge in new_edge_list:
		new_pydot_digraph.add_edge(pydot.Edge(edge[0], edge[1]))

	return new_pydot_digraph

def draw_tumor_map(static_tree, map_tuple, filename, timed_tree):
	tumor_map = map_tuple[0]
	tree = copy.deepcopy(static_tree)
	nodes = lookup_by_names(tree)
	timed_nodes = lookup_by_names(timed_tree)
	node_times = {}
	#tumor_map['anc_node_0'] = tumor_map.get('anc_node_0', tumor_map['anc_node_1'])
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
				#print(name_str)
				nodes[node_name].name = name_str
				#nodes[node_name].name = "{}_{}({})".format(key, tumor_map[node_name][key], node_name)
		else:
			raise ValueError("Node {} has {} tumors assigned.".format(key, len(tumor_map[node_name].keys())))
	pylab.rcParams['figure.figsize'] = len(nodes)/4.0, len(nodes)/4.0
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

def generate_edge_list(static_tree, tumor_map):
	tree = copy.deepcopy(static_tree)
	temp_digraph = tree_to_digraph(tree, "Normal")
	temp_edge_list = sorted(temp_digraph.get_edge_list())
	temp_migration_edges = []
	temp_migration_edge_probabilities = []

	# Group nodes into polytomies
	ungrouped_nodes = tumor_map.keys()
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

	#tumor_map['anc_node_0'] = tumor_map.get('anc_node_0', tumor_map['anc_node_1'])
	#print(tumor_map)
	for edge in temp_edge_list:
		src_node = tumor_map[edge.get_source()]
		dst_node = tumor_map[edge.get_destination()]
		if 'Normal' in [x.keys()[0] for x in [src_node, dst_node]]: continue
		#migration_edge = "{}->{}".format(src_node.keys()[0], dst_node.keys()[0])
		if src_node.keys()[0] != dst_node.keys()[0]:
			temp_migration_edges.append("{}->{}".format(src_node.keys()[0], dst_node.keys()[0]))
			temp_migration_edge_probabilities.append(src_node.values()[0] * dst_node.values()[0])
			temp_node_edges.append((node_groups[edge.get_source()], node_groups[edge.get_destination()], static_tree.distance(edge.get_source(), edge.get_destination())))
	temp_edge_list = zip(temp_migration_edges, temp_migration_edge_probabilities)

	# If a 0-length edge contained in a polytomy duplicates or reverses another edge attached to that polytomy, the 0-length edge should be dropped
	i = 0
	nonzero_edges = [record for record in zip(temp_edge_list, temp_node_edges) if record[1][2] > 0.0]
	bad_edges = []
	for record in zip(temp_edge_list, temp_node_edges):
		if record[1][2] == 0.0:
			for edge in nonzero_edges:
				if edge[0][0] in [record[0][0], reverse_edge(record[0][0])]:
					if record[1][0] in [edge[1][0], edge[1][1]]:
						bad_edges.append((temp_migration_edges[i], temp_migration_edge_probabilities[i]))
						break
		i += 1
	mig_count1 = len(temp_edge_list)
	source_count1 = len(set([temp_edge[0].split('->')[0] for temp_edge in temp_edge_list]))
	for bad_edge in bad_edges:
		temp_edge_list.remove(bad_edge)
	mig_count2 = len(temp_edge_list)
	source_count2 = len(set([temp_edge[0].split('->')[0] for temp_edge in temp_edge_list]))

	temp_edge_list.sort(key=lambda tup: tup[1], reverse=True)
	temp_migration_edges = []
	temp_migration_edge_probabilities = []
	for edge in temp_edge_list:
		j = 1
		while "{}[{}]".format(edge[0], j) in temp_migration_edges: j += 1
		temp_migration_edges.append("{}[{}]".format(edge[0], j))
		temp_migration_edge_probabilities.append(edge[1])

	# Count comigrations
	terminals = get_gv_terminals(temp_digraph)
	mig_edge_counts = {}
	for terminal in terminals:
		path_counts = {}
		path = get_gv_path(temp_digraph, terminal, "Normal")
		for edge in path:
			mig_edge = (tumor_map[edge[0]].keys()[0], tumor_map[edge[1]].keys()[0])
			if mig_edge[0] != mig_edge[1]:
				path_counts[mig_edge] = path_counts.get(mig_edge, 0) + 1
		for key in path_counts.keys():
			mig_edge_counts[key] = max(mig_edge_counts.get(key, 0), path_counts[key])
	comigration_count = sum(mig_edge_counts.values())

	return temp_migration_edges, temp_migration_edge_probabilities, [comigration_count, mig_count1, mig_count2, source_count1, source_count2]

def analyze_edge_list(true_paths, edge_list):
	temp_true_paths = copy.deepcopy(true_paths)
	counts = {'TP(P->M)': 0, 'FP(P->M)': 0, 'FN(P->M)': 0, 'TP(M->M)': 0, 'FP(M->M)': 0, 'FN(M->M)': 0, 'TP(M->P)': 0,
			  'FP(M->P)': 0, 'FN(M->P)': 0}
	for edge in edge_list:
		edge_type = '->'.join([x[0] for x in edge.split('->')])
		if edge in temp_true_paths:
			counts["TP({})".format(edge_type)] += 1
			del temp_true_paths[temp_true_paths.index(edge)]
		else:
			counts["FP({})".format(edge_type)] += 1
	for edge in temp_true_paths:
		edge_type = '->'.join([x[0] for x in edge.split('->')])
		counts["FN({})".format(edge_type)] += 1
	return counts

def make_pydot_seeding_graph(edges, probabilities):
	tumors = set()
	for edge in edges:
		tumors.add(edge.split('[')[0].split('->')[0])
		tumors.add(edge.split('[')[0].split('->')[1])
	tumors = list(tumors)
	new_graph = pydot.Dot(graph_type='digraph')
	new_graph.add_node(pydot.Node(name='P', rank="min"))
	for tumor in tumors:
		if tumor == 'P': continue
		new_graph.add_node(pydot.Node(name=tumor))
	for record in zip(edges, probabilities):
		edge = record[0]
		edge_color = 'grey'
		edge_label = ''
		if record[1] == '-':
			edge_color = 'black'
		else:
			edge_label = round(record[1], 2)
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
	new_graph.add_node(pydot.Node(name='P', rank="min"))
	for tumor in tumors:
		if tumor == 'P': continue
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
#	rev_tumor_label_dict = {'Normal': 'A', 'P': 'C', 'M1': 'D', 'M2': 'E', 'M3': 'F', 'M4': 'G', 'M5': 'H', 'M6': 'I',
#							'M7': 'K', 'M8': 'L', 'M9': 'M', 'M10': 'N', 'M11': 'P', 'M12': 'Q', 'M13': 'R', 'M14': 'S',
#							'M15': 'T', 'M16': 'V', 'M17': 'W', 'M18': 'Y'}
	out_filename = basename + "_tumors.meg"
	tumor_seqs = {}
	with open(out_filename, 'w') as file:
		file.write("#MEGA\n")
		file.write("!Title SNVs;\n")
		file.write("!Format datatype=Protein;\n")
		for clone in clones.keys():
			if len(clones[clone].keys()) == 1:
				tumor_seqs[clone] = rev_tumor_label_dict[clones[clone].keys()[0]]
			elif len(clones[clone].keys()) >= 1:
				nodes = lookup_by_names(tree)
				tumor_list = clones[clone].keys()
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
		if args.infer_primary_seq:
			tumor_seqs["Primary"] = rev_tumor_label_dict['P']
		tumor_seqs["Normal"] = rev_tumor_label_dict['P']
		if mut_seqs is not None:
			for key in tumor_seqs:
				tumor_seqs[key] = tumor_seqs[key] + mut_seqs[key]
		for key in tumor_seqs:
			file.write("#{}\n".format(key))
			file.write("{}\n".format(tumor_seqs[key]))
	return out_filename


def cleanup(basename):
	suffixes = ["_summary.txt","_tumors.meg","_tumors_anc_seqs_in.nwk","_tumors_anc_seqs_out.nwk",".nwk",".meg"]
	for filename in [basename + x for x in suffixes]:
		try:
			shutil.move(filename, args.output)
		except:
			pass
	try:
		os.remove(basename + "_tumors_anc_seqs_out.nwk")
	except:
		pass


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

#aln_file_out = os.path.splitext(os.path.basename(str(sys.argv[1])))[0] + "_with_tumors.meg"
aln_file_out = os.path.splitext(os.path.basename(args.aln))[0] + ".meg"

#parse_input_aln(sys.argv[1], sys.argv[2], aln_file_out)
mut_seqs = parse_input_aln(args.aln, aln_file_out)
clones = parse_clone_freqs(args.clone_locations)

mega_aln_filename = aln_file_out
pre_cleanup(mega_aln_filename)

initial_tree, node_count = infer_mp_tree(mega_aln_filename)
# mega_aln_filename = split_mt_leaves(initial_tree, clones, os.path.splitext(os.path.basename(args.aln))[0])
mega_aln_filename = split_mt_leaves(initial_tree, clones, os.path.splitext(os.path.basename(args.aln))[0], mut_seqs)
temp_nodes = lookup_by_names(initial_tree)
for key in temp_nodes.keys():
	if not temp_nodes[key].branch_length is None:
		temp_nodes[key].branch_length = temp_nodes[key].branch_length * 10.0
tree = copy.deepcopy(initial_tree)

if args.machina_inputs:
	make_machina_labels(clones, os.path.splitext(os.path.basename(args.aln))[0])
	make_machina_trees(initial_tree, os.path.splitext(os.path.basename(args.aln))[0])

eps = get_eps(tree, mega_aln_filename, tumor_label_dict)

timed_anc_seqs_trees = Phylo.parse("{}_anc_seqs_out.nwk".format(os.path.splitext(os.path.basename(mega_aln_filename))[0]), "newick")
timed_anc_seqs_tree = parse_anc_seqs_times(initial_tree, Phylo.BaseTree.Tree.from_clade(timed_anc_seqs_trees.next().clade))

digraph = tree_to_digraph(initial_tree, "Primary")

mt_nodes = []

tumor_membership = get_node_tumors(eps, tumor_membership_cutoff)
#print(initial_tree)
#print(tumor_membership)

for key in tumor_membership.keys():
	if len(tumor_membership[key].keys()) != 1:
		mt_nodes.append(key)
		#raise ValueError("Node {} has multiple tumors assigned.".format(key))


#tumor_map_list is a list of (dict, int) tuples
tumor_map_list = [(copy.deepcopy(tumor_membership), 1.0)]
for node in mt_nodes:
	del tumor_map_list[0][0][node]

#Group multi-tumor nodes in polytomies
ungrouped_nodes = copy.deepcopy(mt_nodes)
grouped_nodes = {}
while len(ungrouped_nodes) > 0:
	key_node = ungrouped_nodes[0]
	node_group = []
	for node in ungrouped_nodes:
		if initial_tree.distance(key_node, node) == 0.0:
			node_group.append(node)
	for node in node_group:
		ungrouped_nodes.remove(node)
	grouped_nodes[key_node] = node_group

#Generate list of all possible combinations of multi tumor node assignments with probabilities
for key_node in grouped_nodes.keys():
	temp_tumor_map_list = []
	for tumor_map in tumor_map_list:
		for key in tumor_membership[key_node].keys():
			new_tumor_map = (copy.deepcopy(tumor_map[0]), tumor_map[1] * tumor_membership[key_node][key])
			for polytomy_node in grouped_nodes[key_node]:
				new_tumor_map[0][polytomy_node] = {key: tumor_membership[key_node][key]}
			temp_tumor_map_list.append(copy.deepcopy(new_tumor_map))
	tumor_map_list = copy.deepcopy(temp_tumor_map_list)


edge_lists = []
edge_lists_probabilities = []
graphic_file_list = []
migration_counts = []
#i = 1
for record in tumor_map_list:
	# Generate list of tumor migration edges from tumor map
	edge_list, edge_list_probabilities, counts = generate_edge_list(initial_tree, record[0])
	edge_lists.append(edge_list)
	edge_lists_probabilities.append(edge_list_probabilities)
	migration_counts.append(counts)

# Sort list of edge lists by probability
# temp_sorted = sorted(zip(tumor_map_list, edge_lists, edge_lists_probabilities), key=lambda row: row[2], reverse=True)
# tumor_map_list, edge_lists, edge_lists_probabilities = zip(*temp_sorted)
temp_sorted = sorted(zip(tumor_map_list, edge_lists, edge_lists_probabilities, migration_counts), key=lambda row: row[2], reverse=True)
tumor_map_list, edge_lists, edge_lists_probabilities, migration_counts = zip(*temp_sorted)

topology_tree = copy.deepcopy(initial_tree)
clades = lookup_by_names(topology_tree)
for clade in clades.keys():
	clades[clade].branch_length = None
mig_file_list = []
i = 1
for (record, edge_list, edge_list_probabilities) in zip(tumor_map_list, edge_lists, edge_lists_probabilities):
	# Generate phylotree graphic for each tumor map in tumor_map_list
	#draw_tumor_map(topology_tree, record, os.path.join(args.output, os.path.splitext(os.path.basename(args.aln))[0] + "_{}_phylogeny.png".format(i)), timed_anc_seqs_tree)
	draw_tumor_map(timed_anc_seqs_tree, record, os.path.join(args.output, os.path.splitext(os.path.basename(args.aln))[0] + "_{}_phylogeny.png".format(i)), timed_anc_seqs_tree)
	graphic_file_list.append(os.path.splitext(os.path.basename(args.aln))[0] + "_{}_migration.png".format(i))
	with open(os.path.join(args.output, os.path.splitext(os.path.basename(args.aln))[0] + "_{}_Mig.txt".format(i)), 'w') as mig_file:
		mig_file_list.append(os.path.splitext(os.path.basename(args.aln))[0] + "_{}_Mig.txt".format(i))
		mig_file.write("Tumor Membership Probability: {}\n\n".format(record[1]))
		for mig_edge in zip(edge_list, edge_list_probabilities):
			#mig_file.write("{}[{}]\n".format(mig_edge[0], mig_edge[1]))
			mig_file.write("{}\n".format(mig_edge[0]))
	#seeding_graph = make_pydot_seeding_graph(edge_list, edge_list_probabilities)
	seeding_graph = make_pydot_seeding_graph(edge_list, ['-' for x in edge_list_probabilities])
	seeding_graph.set('label',"Probability: {}".format(record[1]))
	seeding_graph.set('labelloc', 't')
	try:
		seeding_graph.write_png(os.path.join(args.output, os.path.splitext(os.path.basename(args.aln))[0] + "_{}_migration.png".format(i)))
	except OSError:
		stuff = sys.exc_info()
		if stuff[1].strerror in ['"dot.exe" not found in path.', '"dot" not found in path.']:
			print("Warning: 'dot.exe' not found in path, migration paths graphic not generated, see README file...")
		else:
			raise stuff[0], stuff[1], stuff[2]
	i += 1

if read_true_paths:
	accuracy_list = []
	all_edges = set()
	for edge in true_paths:
		all_edges.add(edge)
	for edge_list in edge_lists:
		# Generate accuracy data for list of edges
		accuracy_list.append(analyze_edge_list(true_paths, edge_list))
		for edge in edge_list:
			all_edges.add(edge)
	i = 1
	for (node_membership, edge_list, edge_list_probabilities) in zip(tumor_map_list, edge_lists, edge_lists_probabilities):
		new_graph = make_proto_graph(all_edges)
		edge_probabilities = dict(zip(edge_list, edge_list_probabilities))
		for graph_edge in new_graph.get_edge_list():
			edge = graph_edge.obj_dict['attributes']['id']
			# if edge in edge_probabilities.keys():
			# 	graph_edge.obj_dict['attributes']['label'] = round(edge_probabilities[edge], 2)
			# else:
			# 	graph_edge.obj_dict['attributes']['label'] = ""
			if edge in true_paths and edge in edge_list:	# TP edge
				graph_edge.obj_dict['attributes']['style'] = ""
				if edge_probabilities[edge] > 0.9:
					graph_edge.obj_dict['attributes']['color'] = "forestgreen"
				elif edge_probabilities[edge] > 0.5:
					graph_edge.obj_dict['attributes']['color'] = "chartreuse"
				else:
					graph_edge.obj_dict['attributes']['color'] = "yellow"
			elif edge in true_paths:						# FN edge
				graph_edge.obj_dict['attributes']['style'] = ""
				graph_edge.obj_dict['attributes']['color'] = "red"
			elif edge in edge_list:							# FP edge
				graph_edge.obj_dict['attributes']['style'] = ""
				graph_edge.obj_dict['attributes']['color'] = "blue"
		try:
			new_graph.write_png(os.path.join(args.output, os.path.splitext(os.path.basename(args.aln))[0] + "_{}_migration.png".format(i)))
		except OSError:
			stuff = sys.exc_info()
			if stuff[1].strerror in ['"dot.exe" not found in path.', '"dot" not found in path.']:
				print("Warning: 'dot.exe' not found in path, migration paths graphic not generated, see README file...")
			else:
				raise stuff[0], stuff[1], stuff[2]
		i += 1


	mt_nodes_string = ';'.join(["{}_{}".format(x, ','.join(["{}({})".format(y, str(round(tumor_membership[x][y],2))) for y in tumor_membership[x].keys()])) for x in mt_nodes])

	# Generate probability/accuracy summary/graphic filename list for each record in tumor_map_list
	count_keys = ['TP(P->M)', 'FP(P->M)', 'FN(P->M)', 'TP(M->M)', 'FP(M->M)', 'FN(M->M)', 'TP(M->P)', 'FP(M->P)', 'FN(M->P)']

	header = "Tumor Map File\tProbability\t{}\n".format('\t'.join(count_keys))
	with open(os.path.join(args.output, os.path.splitext(os.path.basename(args.aln))[0] + "_mt_accuracy.txt"), 'w') as file:
		file.write(header)
		# for record in zip([x[1] for x in tumor_map_list], accuracy_list, graphic_file_list):
		for record in zip([x[1] for x in tumor_map_list], accuracy_list, graphic_file_list, migration_counts):
			# line = "{}\t{}\t{}\n".format(record[2], record[0], '\t'.join([str(record[1][x]) for x in count_keys]))
			line = "{}\t{}\t{}\t{}\n".format(record[2], record[0], '\t'.join([str(record[1][x]) for x in count_keys]), '\t'.join([str(x) for x in record[3]]))
			file.write(line)
	# Aggregate accuracy metrics (if necessary) and print to stdout
	max_prob = max([x[1] for x in tumor_map_list])
	max_prob_counts = []
	for record in zip([x[1] for x in tumor_map_list], accuracy_list):
		if record[0] == max_prob:
			max_prob_counts.append(record[1])
	max_avg_counts = {}
	for key in count_keys:
		max_avg_counts[key] = sum([x[key] for x in max_prob_counts])/float(len(max_prob_counts))
	print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(os.path.splitext(os.path.basename(args.aln))[0], max_prob, len(max_prob_counts),
										len(tumor_map_list), len(mt_nodes), '\t'.join([str(max_avg_counts[x]) for x in count_keys]),
										mt_nodes_string))

#print(timed_anc_seqs_tree)
#print(initial_tree)
nodes = lookup_by_names(timed_anc_seqs_tree)
# Draw multi tumor tree
#tumor_membership['anc_node_0'] = tumor_membership.get('anc_node_0', tumor_membership['anc_node_1'])
for node_name in nodes.keys():
	name_str = node_name
	if len(tumor_membership[node_name].keys()) == 1:
		for key in tumor_membership[node_name].keys():
			if tumor_membership[node_name][key] >= 0.995:
				name_str = key
			else:
				name_str = "{}_{}".format(key, round(tumor_membership[node_name][key],2))
			if node_name[0:8] != "anc_node":
				name_str = name_str + "({})".format(node_name)
			#name_str = "{}_{}({})".format(key, tumor_membership[node_name][key], node_name)
	else:
		for key in tumor_membership[node_name].keys():
			name_str = name_str + "\n{}: {}".format(key, round(tumor_membership[node_name][key],2))
	nodes[node_name].name = name_str


pylab.rcParams['figure.figsize'] = len(nodes)/4.0, len(nodes)/4.0
Phylo.draw(timed_anc_seqs_tree, label_func=get_label, do_show=False)
pylab.savefig(os.path.join(args.output, os.path.splitext(os.path.basename(args.aln))[0] + "_phylogeny.png"))

avg_seeding_probs = get_avg_mig_probs(mig_file_list)

temp_edge_list = list(avg_seeding_probs.keys())
temp_prob_list = [avg_seeding_probs[key] for key in temp_edge_list]

with open(os.path.join(args.output, os.path.splitext(os.path.basename(args.aln))[0] + "_Mig.txt"), 'w') as avg_mig_file:
	for record in zip(temp_edge_list, temp_prob_list):
		avg_mig_file.write("{}[{}]\n".format(record[0], round(record[1],3)))

avg_seeding_graph = make_pydot_seeding_graph(temp_edge_list, temp_prob_list)
avg_seeding_graph.set('label', "Average Migration Graph")
avg_seeding_graph.set('labelloc', 't')
try:
	avg_seeding_graph.write_png(os.path.join(args.output, os.path.splitext(os.path.basename(args.aln))[0] + "_migration.png"))
except OSError:
	exc_info = sys.exc_info()
	if exc_info[1].strerror in ['"dot.exe" not found in path.', '"dot" not found in path.']:
		print("Warning: 'dot.exe' not found in path, migration paths graphic not generated, see README file...")
	else:
		raise exc_info[0], exc_info[1], exc_info[2]


if read_true_paths:
	accuracy_list = []
	all_edges = set()
	for edge in true_paths:
		all_edges.add(edge)
	for edge in avg_seeding_probs.keys():
		all_edges.add(edge)
	new_graph = make_proto_graph(all_edges)
	edge_probabilities = dict(zip(temp_edge_list, temp_prob_list))
	with open(os.path.join(args.output, os.path.splitext(os.path.basename(args.aln))[0] + "_Mig.txt"), 'w') as avg_mig_file:
		if dump_edges: edges_file = open(args.edges_file,'a+')
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
					graph_edge.obj_dict['attributes']['color'] = "red"
					avg_mig_file.write("{}\t{}\t{}\n".format(edge, edge_probabilities[edge], "False Negative"))
			elif edge in true_paths:							# FN edge
				if dump_edges: edges_file.write("{}\t{}\t{}\t{}\n".format(dataname, edge, edge_probabilities.get(edge, 0.0), "True"))
				graph_edge.obj_dict['attributes']['style'] = ""
				graph_edge.obj_dict['attributes']['color'] = "red"
				avg_mig_file.write("{}\t{}\t{}\n".format(edge, edge_probabilities.get(edge, 0.0), "False Negative"))
			elif edge in temp_edge_list:						# FP edge
				if dump_edges: edges_file.write("{}\t{}\t{}\t{}\n".format(dataname, edge, edge_probabilities.get(edge, 0.0), "False"))
				if edge_probabilities[edge] >= avg_edge_weight_cutoff:
					graph_edge.obj_dict['attributes']['style'] = ""
					graph_edge.obj_dict['attributes']['color'] = "blue"
					avg_mig_file.write("{}\t{}\t{}\n".format(edge, edge_probabilities[edge], "False Positive"))
				else:
					graph_edge.obj_dict['attributes']['style'] = ""
					graph_edge.obj_dict['attributes']['color'] = "grey"
					avg_mig_file.write("{}\t{}\t{}\n".format(edge, edge_probabilities[edge], "True Negative"))
		if dump_edges: edges_file.close()
	new_graph.set('label', "Average Migration Graph")
	new_graph.set('labelloc', 't')
	try:
		new_graph.write_png(os.path.join(args.output, os.path.splitext(os.path.basename(args.aln))[0] + "_migration.png"))
		new_graph.write(os.path.join(args.output, os.path.splitext(os.path.basename(args.aln))[0] + "_migration.dot"))
	except OSError:
		stuff = sys.exc_info()
		if stuff[1].strerror in ['"dot.exe" not found in path.', '"dot" not found in path.']:
			print("Warning: 'dot.exe' not found in path, migration paths graphic not generated, see README file...")
		else:
			raise stuff[0], stuff[1], stuff[2]
	with open(os.path.join(args.output, os.path.splitext(os.path.basename(args.aln))[0] + "_mt_accuracy.txt"), 'a') as file:
		count_keys = ['TP(P->M)', 'FP(P->M)', 'FN(P->M)', 'TP(M->M)', 'FP(M->M)', 'FN(M->M)', 'TP(M->P)', 'FP(M->P)', 'FN(M->P)']
		avg_acc_counts = analyze_edge_list(true_paths, [x for x in avg_seeding_probs.keys() if avg_seeding_probs[x] >= avg_edge_weight_cutoff])
		avg_graphic_file = os.path.splitext(os.path.basename(args.aln))[0] + "_migration.png"
		line = "{}\t{}\t{}\n".format(avg_graphic_file, "avg", '\t'.join([str(avg_acc_counts[x]) for x in count_keys]))
		file.write(line)
	print("{}\t{}\t{}\t{}\t{}\t{}".format(os.path.splitext(os.path.basename(args.aln))[0], max_prob, len(max_prob_counts),
										len(tumor_map_list), len(mt_nodes), '\t'.join([str(avg_acc_counts[x]) for x in count_keys])))

time.sleep(0.5)
cleanup(os.path.splitext(os.path.basename(args.aln))[0])

print("runtime: {}".format(datetime.datetime.now() - start_time))


