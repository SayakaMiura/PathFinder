PathFinder_v1.0
Updated April 21, 2020
======================

PathFinder is a method that infers cell migration routes accurately by using principles of Bayesian molecular phylogenetics in cancer cell genomes.


Installation
======================
PathFinder is a python script developed in a Windows 64-bit architecture, and does not need to be compiled prior to usage. To use PathFinder, the following python packages need to be installed.

	Python 2
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

	1) An alignment file in FASTA format, specifying the sequence of each tumor clone, including one sequence labeled "Normal" specifying the normal sequence.
	2) A tab delimited clone presence file, specifying the presence of each clone in each tumor site (clones = columns, tumor sites = rows), as in the provided example files.


***IMPORTANT***

The number of tumors cannot exceed 19.



Parameters
======================
* -o, --output			Output Directory
Output directory to put results in.
* -t, --true_paths		True Paths File
Specifies a file listing the correct tumor migration paths causing output visualization to be accordingly color coded.
* --anc_tumor_threshold		Ancestral Tumor Threshold
Lower limit on probability of an ancestral node belonging to a tumor, to include that possibility in the set of seeding graphs generated.
Default: 0.15
* --mig_event_threshold		Migration Event Threshold
Lower limit on the composite probability assigned to a potential migration event to consider that event to have happened.
Default: 0.5
* --machina_inputs
Tell PathFinder to generate mutation tree and label files suitable for input to MACHINA.
* --infer_primary_seq
Tell PathFinder to infer a primary clone sequence based on Normal alleles absent from Normal and present in all other clone sequences.



Example
======================
A set of example inputs (Example_data\input.fas, Example_data\input.txt) is provided to run PathFinder with default parameters. To run PathFinder on the example inputs, please follow commands below from the main directory.

	python pathfinder.py Example_data\input.fas Example_data\input.txt -o Example_output


After running PathFinder, the following output files can be found in the main directory, or in the directory specified with -o option.


Output Files
======================

All output files use the base name of the provided alignment file as their base name.

The main output files are:

{basename}_phylogeny.png
	Shows the phylogency of clone sequences, annotating each node in the phylogeny with the tumor it belonged to (leaf nodes), or the tumor(s) that it was inferred to belong to, including the probability of belonging to each tumor in the case that membership is ambiguous.

{basename}_migration.png
	Shows the inferred composite migration paths, displaying paths with a probability exceeding the threshold optionally provided by --mig_event_threshold (default 0.5) in black, and paths below the threshold in grey.

{basename}_Mig.txt
	Provides the same data displayed in {basename}_migration.png in a textual format, with the first number in brackets denoting the index of the edge among edges with the same source and destination, and the second number in bracket denoting the edge's composite probability.



When tumor membership of ancestral nodes is ambiguous, each possibility is evaluated and produces its own phylogeny and migration file. For each possibility, the following outputs will be produced:

{basename}_{i}_phylogeny.png
	Shows the phylogeny of clone sequences for the {i}th most likely combination of ancestral node membership states, annotating each node in the phylogeny with the tumor it belonged to for this combination.

{basename}_{i}_migration.png
	Shows the set of migration paths inferred from the {i}th most likely combination of ancestral node membership states.

{basename}_{i}_Mig.txt
	Provides the same data displayed in {basename}_{i}_migration.png in a textual format, with the number in brackets denoting the index of the edge among edges with the same source and destination.


The output directory will also contain the following intermediate files used/produced by MEGA-CC during processing.

{basename}.meg
	The alignment of clone sequences given to MEGA-CC to generate the clone phylogeny.

{basename}_tumors_anc_seqs_in.nwk
	The clone phylogeny produced from the initial alignment of clone sequences.

{basename}_tumors.meg
	The alignment of clone sequences given to MEGA-CC, with each sequence prepended with a character representing the tumor it has membership in, to infer the tumor membership of ancestral clones.
