import math
import itertools
import types
import datetime
import random
import copy


def split_list_recursive(leaves):
	n = len(leaves)
	if n == 1:
		return leaves
	splits = []
	for i in range(1, int(math.floor(float(n)/2.0))+1):
		for left_split in split_list_recursive(leaves[0:n-i]):
			for right_split in split_list_recursive(leaves[n-i:]):
				splits.append((left_split, right_split))
	return splits


def apply_labels(tree, labels):
	if isinstance(tree, str):
		return labels[tree]
	return (apply_labels(tree[0], labels), apply_labels(tree[1], labels))


def trees_are_equal(tree_l, tree_r):
	if isinstance(tree_l, str) != isinstance(tree_r, str):
		return False
	if isinstance(tree_l, str) and isinstance(tree_r, str):
		if tree_l == tree_r:
			return True
		else:
			return False
	if (trees_are_equal(tree_l[0], tree_r[1]) and trees_are_equal(tree_l[1], tree_r[0])) or (
			trees_are_equal(tree_l[1], tree_r[1]) and trees_are_equal(tree_l[0], tree_r[0])):
		return True
	return False


def permute_unique_trees(leaf_names, max_samples=None):
	# if max_samples is None:
	# 	leaf_ids = ["leaf_{}".format(i) for i in range(len(leaf_names))]
	# 	trees = {}
	# 	for topology in split_list_recursive(leaf_ids):
	# 		topo_trees = {}
	# 		for labeling in itertools.permutations(leaf_names):
	# 			labels = dict(zip(leaf_ids, labeling))
	# 			tree = apply_labels(topology, labels)
	# 			tree_is_unique = True
	# 			for key_tree in topo_trees.keys():
	# 				if trees_are_equal(key_tree, tree):
	# 					tree_is_unique = False
	# 					topo_trees[key_tree].append(tree)
	# 					break
	# 			if tree_is_unique:
	# 				topo_trees[tree] = [tree]
	# 		trees.update(topo_trees)
	# 	return trees
	leaf_ids = ["leaf_{}".format(i) for i in range(len(leaf_names))]
	trees = []
	topologies = split_list_recursive(leaf_ids)
	if math.factorial(len(leaf_names)) >= 10 * max_samples:
		labelings = []
		while len(labelings) < max_samples:
			random.shuffle(leaf_names)
			new_labeling = tuple(copy.deepcopy(leaf_names))
			if new_labeling not in labelings:
				labelings.append(new_labeling)
		labelings = [list(labeling) for labeling in labelings]
	else:
		labelings = list(itertools.permutations(leaf_names))
	# print("Topology count: {}\nLabeling count: {}".format(len(topologies), len(labelings)))
	if max_samples is None or max_samples > len(topologies) * len(labelings):
		for topology in topologies:
			for labeling in labelings:
				labels = dict(zip(leaf_ids, labeling))
				tree = apply_labels(topology, labels)
				tree_is_unique = True
				for selected_tree in trees:
					if trees_are_equal(tree, selected_tree):
						tree_is_unique = False
						break
				if tree_is_unique:
					trees.append(tree)
	else:
		i = 0
		while len(trees) < max_samples and i < max_samples * 10:
			i += 1
			labels = dict(zip(leaf_ids, random.choice(labelings)))
			tree = apply_labels(random.choice(topologies), labels)
			tree_is_unique = True
			for selected_tree in trees:
				if trees_are_equal(tree, selected_tree):
					tree_is_unique = False
					break
			if tree_is_unique:
				trees.append(tree)
	return trees


if __name__ == '__main__':
	nodes = ['anc_node_13', 'anc_node_14', 'Clo30', 'Clo32', 'anc_node_9', 'anc_node_6', 'Clo19', 'Clo17', 'Clo1', 'Clo7']
	for i in range(2, len(nodes)+1):
		start_time = datetime.datetime.now()
		# trees = permute_unique_trees(nodes[0:i])
		trees = permute_unique_trees(nodes[0:i], 1000)
		print("Polytomy with {} nodes produced {} trees and took {} to run.".format(i, len(trees), datetime.datetime.now() - start_time))
	# for key_tree in trees.keys():
	# 	print("\n\n{}".format(key_tree))
	# 	for tree in trees[key_tree]:
	# 		print("\t{}".format(tree))
