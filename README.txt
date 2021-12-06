PathFinder_v2.0
Updated July 14, 2021
======================

PathFinder is a method that infers cell migration routes accurately by using principles of Bayesian molecular phylogenetics in cancer cell genomes.


Installation
======================
PathFinder is a python script developed in a Windows 64-bit architecture, and does not need to be compiled prior to usage. To use PathFinder, the following python packages need to be installed.

	Python 3
	matplotlib
	pydot
	graphviz
	networkx
	biopython

Additionally, PathFinder uses MEGA-CC. MEGA-CC can be downloaded for Windows, Mac OS X, and Linux from (http://www.megasoftware.net).

Lastly, to generate visualizations, PathFinder relies on the graphviz dot executable, so the Graphviz software needs to be installed and the dot.exe file placed in the system PATH variable in order for PathFinder to produce visualizations.


Input files
======================

The required input files are:

	1) An alignment file in FASTA or MEGA format, specifying the sequence of each tumor clone, including one sequence labeled "Normal" specifying the normal sequence(this is the germline sequence, which represents a healthy cell without any somatic mutations).
	2) A tab delimited clone presence file, specifying the presence of each clone in each anatomical site (clones = columns, anatomical sites = rows), as in the provided example files.
	
Additionally, the primary tumor must be specified on the command line.


***IMPORTANT***

The number of tumors cannot exceed 19.



Parameters
======================
Positional arguments:
	aln			Clone sequence alignment file.
	clone_locations		Per-tumor clone presence/absence matrix file.
	primary			Name of primary tumor.

Optional arguments:
	-o, --output
		Output directory to put results in.
	-t, --true_paths
		Specifies a file listing the correct tumor migration paths causing output visualization to be accordingly color coded.
	--primary
		Name of the primary tumor.
		Default: Treat normal sequence as primary.
	--anc_tumor_threshold
		Lower limit on probability of an ancestral node belonging to a tumor, to include that possibility in the set of seeding graphs generated.
		Default: 0.15
	--mig_event_threshold
		Lower limit on the composite probability assigned to a potential migration event to consider that event to have happened.
		Default: 0.5
	--infer_primary_seq
		Tell PathFinder to infer a primary clone sequence based on Normal alleles absent from Normal and present in all other clone sequences.
	--default_normal_char
		Default character to use for Normal sequence, if Normal sequence is not specified in input alignment.
	--max_permutations
		Maximum number of tree polytomy permutations to use.
		Default: 30
	--max_graphs_per_tree
		Maximum number of node-tumor membership states to analyze per initial tree polytomy permutation.
	--use_all_weighted_outputs
		Make final weighted edge list migration graph using all outputs of all sampled permutations, instead of using hierarchical migration-count based selection.
	--use_select_weighted_outputs
		Make final weighted edge list using only outputs with minimized count-based selection, weighted by each output's composite ancestral sequences probability.
	--draw_all_outputs
		Generate migration graph and phylogeny images for all sampled permutations/configurations.
	--phy_scale_x
		Horizontal scaling factor for tree image outputs.
	--phy_scale_y
		Vertical scaling factor for tree image outputs.
	--acc_by_edge_type
		Subset accuracy counts into P->M, M->M, M->P.
	--relax_threshold
		If a node has no anatomical site with probability>threshold, fall back to selecting the most likely anatomical site(s).
	--log_ancestral_inference
		Save inputs/outputs of MEGA-CC ancestral sequence inference calculations for debugging purposes.
	--keep_ambiguous_results
		Treat nodes with fully ambiguous inferences as belonging equally to all anatomical sites, instead of discarding the ancestral inference result.


Example
======================
A set of example inputs (Example_data\input.fas, Example_data\input.txt) is provided to run PathFinder with default parameters. To run PathFinder on the example inputs, please follow commands below from the main directory.

	python pathfinder.py Example_data\input.fas Example_data\input.txt --primary P -o Example_output

To run the example with the true paths specified:

	python pathfinder.py Example_data\input.fas Example_data\input.txt --primary P -o Example_output -t Example_data\ground_truth.txt


After running PathFinder, the following output files can be found in the main directory, or in the directory specified with -o option.


Output Files
======================

All output files use the base name of the provided alignment file as their base name.

The main files found in the specified output directory are:

{basename}_Mig.png
	A text list of migration edges with their probabilities.

{basename}_migration.png
	Shows the inferred composite migration paths, displaying paths with a probability exceeding the threshold optionally provided by --mig_event_threshold (default 0.5) in black, and paths below the threshold in grey.
	If a true paths file is specified, shown paths will be color-coded according to their classification:
		True positive: Green or yellow, depending on probability.
		False positive: Red
		False Negative: Blue
		True Negative: Grey

{basename}_all_output_counts.txt
	List of count/accuracy information for each of the individual outputs generated.


Additional files can be found in the "scratch" subdirectory:

{basename}.meg:
	The alignment file used to generate the initial Newick tree.

{basename}_tumors.meg:
	The alignment file containing only tumor characters, used by MEGACC to perform the ancestral sequences inference.

{basename}_tree_{i}.nwk:
	The {i}th sampled polytomy permutation of the initial tree, used by MEGACC to perform the ancestral sequences inference.

{basename}_{i}_{j}_migration.png:
	If --draw_all_outputs is specified, shows the inferred migration path for each individual output, where i is the index of the ampled polytomy permutation of the initial tree, and j is the index of the ancestral clone tumor membership configuration within that tree.


How to run PathFinder
======================

When a dataset contains a high number of metastatic tumors or the sequence length is small, high computational time may make it infeasible to complete the run. In order to decrease run time, the –-max_permutations and/or –-max_graphs_per_tree arguments may be used to subsample the trees/graphs being analyzed, however, doing so will cause the results to be non-deterministic, so multiple runs should be used in these cases, to ensure that the result is stable.

Additionally, there may be cases where all or most result sets include a node where no anatomical site has probability greater than the threshold value, causing those result sets to be discarded and causing the run to fail. In these cases, the –-relax_threshold flag may be added to ignore the minimum threshold only for those nodes, causing the program to select the anatomical site(s) with the highest probability irrespective of the minimum threshold value.
 

How to read output files generated by PathFinder
======================

Please be mindful of the number of migration graphs generated by PathFinder. This should be considered together with the number of migration, comigration and source count and the probabilities assigned to each inferred migration graph. The probability value will also indicate the number of migration graphs generated. The migration graph with the lowest number of migration, comigration and source count and the highest probability value is the most likely migration history. For example, if we observe probability= .8, then one or a few migration graphs will be predicted. Therefore:
	1) If Probability > .5 (high), we expect the rest of the migration graphs to have lower probabilities. 
	2) If Probability < .5 (low), we expect a certain ambiguity to arise in pathways. We recommend the user to be cautious and examine any migration inferences that result into equal number of migration, comigration and source count with almost similar probabilities values. 



