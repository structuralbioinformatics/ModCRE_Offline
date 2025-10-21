####                                                           ####
####      ModCRE to be deployed as a downloadable Package      ####
####                                                           ####

Installation instructions
------
1) download this repository 
2) download the pbm and pdb database folders
    a) wget http://aleph.upf.edu/modcrefiles/pbm.tgz (60 Gb)
    b) wget http://aleph.upf.edu/modcrefiles/pdb.tgz (57 Gb)
3) place these folders into the downloaded reopository
4) decompress them 
    a) tar -xvzf pbm.tgz (110 Gb)
    b) tar -xvzf pdb.tgz (88 Gb)
5) Ensure all dependencies are installed

Dependencies:

Python
----
Bio              → from Biopython (`pip install biopython`)
SBILib           → comes with package but can be installed seperately (`pip install SBILib` see https://github.com/structuralbioinformatics/SBILib for dependencies)
bs4              → BeautifulSoup (`pip install beautifulsoup4`)
bottle           → lightweight web framework (`pip install bottle`)
ihm              → Integrative Modeling library (`pip install ihm`)
matplotlib       → plotting library (`pip install matplotlib`)
numpy            → numerical computing (`pip install numpy`)
pandas           → data analysis (`pip install pandas`)
plotly           → interactive plotting (`pip install plotly`)
scipy            → scientific computing (`pip install scipy`)
seaborn          → statistical plotting (`pip install seaborn`)
sklearn          → scikit-learn (`pip install scikit-learn`)

Enviornmental
---
BLAST+ (run on 2.12.0)
CD-HIT (run on 4.8.1)
Clustal-Omega (run on 1.2.4)
ClustalW2 (run on 2.1)
EMBOSS (run on 6.6.0)
Ghostscript (run on 9.53.3)
HMMER (run on 3.3.2)
MEME (run on 5.1.1)
Modeller (run on 10.3)
Python (run on 3.8.6)
TMalign (run on Version 30012025)
dssp (run on 3.0.0)
x3dna (run on 2.5.0)


6) Change paths in ModCRElib/configure/config.ini to point to the correct location for your installations. 


----------
Running without parallel
----------


# Example 1 
Modelling a Transcription Factor (TF) from an amino acid sequence
------ 
1) We need a single or multi fasta file containing a.a. sequence of protein to be modelled (example of AHR_Example.fa is given)
2) modelling.sh contains the command to model the protein, 
    pdb="/home/pgohl/ModCRE_Package/pdb" Needs to be replaced with the location of the downloaded pdb folder on your machine
    -i uniput protein in fasta (could be multi fasta),
    -o output directory (models) 
####
We now have a folder containing the modelled Transcription factor in pdb format. 
You can view the file in chimera or online at https://www.rcsb.org/3d-view by uploading the file.
####

# Example 2
Predict TF binding specificity
------
3) bin/renumberModels.py is a script to ensure that atom and amino acid numbers of models from modelling.sh is continuous (not allways needed)
    arg 1 = folder containing the pdbs to be renumbered
    arg 2 = output folder 
    (python bin/renumberModels.py models remodels)
4) pwm.sh predicts the binding specificity
    pdb="/home/pgohl/ModCRE_Package/pdb" Needs to be replaced with the location of the downloaded pdb folder on your machine
    pbm="/home/pgohl/ModCRE_Package/pbm" Needs to be replaced with the location of the downloaded pbm folder on your machine
####
We now have a predicted binding specificity of the transcription factor in pwm and meme format in the output folder specified in the previous command. 
Remember that ModCRE is designed to use PWMs in aggregate to predict binding sites and individual PWM predictions may be more or less accurate. 
####

# Example 3
Scan dna sequence for binding sites
------
5) pwm/make_scan_ready.py is file that generates a database file of the predicted pwms from the previous step that are in the correct format for scanning.
    -line 4 references "database.txt" (the generated database file from pwm.sh) and the name of the new file to be used for scanning
6) scan_sequence.sh is the script to run a scan with.
    -i The dna sequence file in fasta format that is to be scanned
    -o the name of the output folder
    -s the species identifier (if restricting scan to TFs of a given specie)
    -ft the fimo threshold value to be used in designating hits
    --db the location of the database file to be used (the defaul option to be used is the larger pwm folder provided)
    -c Cluster complexes into connected binding sites
####
There are several ways to observe the results. The TFs binding the sequence and their binding sites can be retrieved from the file orthologs.json.
Here we can view the name of the TF experienceing the hit, the start and end index for binding along the DNA sequence and the various orthologs binding.
Alternatively, scan_sequence.sh (with parameters used above) will generate thread files of the TF and its orthologs binding the DNA sequence in the folder aux_files.
####

# Example 4
Generate a model of a TF attached to a predicted binding site along a full length of DNA
------
####
We can view the TFs binding to the full scanned DNA sequence as predicted in the previous step. For this we will need to process the files a little through.
####
7)  modelling.sh contains the option to use thread files as an input instead of the previously used multi fasta files
    -i a file containng a list of the location of the thread files to be used (Threads_list.txt is provided)
    -t indicates that threads are used instead of aa sequences
    -o the output folder location 
8) We copy the desired models (from the output of the previous step) into a folder containg all the binary interactions that we would like to use in the modelled complex
9) rename_complex_input.py is a python script that prepares scanning output file names for complex builder (eg. python rename_complex_input.py BinaryInteractions/A9YTQ3.5nj8_1A.18-29.pdb ---> BinaryInteractions/A9YTQ3.5nj8_1A.18-29:1:243_TF.pdb)
    -the name of the file must follow the following format:
        {UniprotAccession}.{PDBID}_{Chain}.{index of binding start}-{index of binding end}:{model start index}:{model end index}_{a label}.pdb
10) BuildComplex.sh is the script that will build a complex based on binary interaction files contained in a given folder
    a) /soft/system/software/x3dna/2.3/bin/fiber will generate a pdb file for a given DNA sequence
        -seq the dna sequence to be used (in the case of modelling the scanning results use the same sequence)
        -b the dna conformation to be applied followed by output location 
    b) exe/complexbuilder.py will build the complex
        -d the folder containing the binary interaction files
        -o the output folder location 
####
Now the modelled complex can be view in the output folder (Complex/fragment_1-100/dna__1-100_aa.pdb). 
####

# Example 5
Generate thread files from a modelled TF 
------
11) get_best_bindings_threads.sh produces thread files for use in modelling and retrieving scores. 
    pdb="/home/pgohl/ModCRE_Package/pdb" Needs to be replaced with the location of the downloaded pdb folder on your machine
    pbm="/home/pgohl/ModCRE_Package/pbm" Needs to be replaced with the location of the downloaded pbm folder on your machine
    -i (the pdb of the transcription factor to be used)
    -o (output directory)  
    --seq fasta sequence containing DNA sequence to be bound
    --dna nucleatide sequence of the binding site 
    --pwm meme file for the transcription factor
####
The thread file can be used to score the binding of a TF along any DNA sequence that matches the binding site length. 
previous steps need not be repeated for other DNA sequences, simply create a thread file for the relevant substituted sequence by
replaceing the DNA sequences at the bottom of the threads folder:
>dna
CAGCTGGCTGTG;0
//
>dna_fixed
CAGCTGGCTGTG;0
//

# Example 6
Generate a scoring profile of a TF-DNA interaction
------
####
12) get_best_score.sh produces a scoring profile for the transcription factor on target dna sequence.
    -i (the input is a thread file as produced by get_best_bindings_threads)
    -o (output directory)  
####
Finally we have a score profile (statistical potentials) for the TF binding along the tested DNA binding site.
####

# Example 7
Generate a scoring profile plot for a TF along a DNA sequence
------
13) get_score_profiles.sh will generate the raw scores that will be used to plot the scoring 
    -i a file containing the folders containing the models to be used
    -d a fasta format dna sequence to profile the TF binding against
    -o Output name for tables and plots (default is the name of the input FOLDER)
####
The output will be stored in each folder provided by the input file. A folder will have been generated and named 
profilerinput.txt_profiling.34_272 by default (profilerinput.txt being the input file, profiling the folder from that file).
Within can be found individual model scores (pickle files) as well as the mean tables (csv files)
####
14) bin/plotprofile.py will generate a plot from the mean score file
    arg 1 = path to the mean table to be plotted 
    arg 2 = location the plot is to be stored at
    arg 3 = column (score type) to be plotted (options will be printed out if provided isn't in file)
    (python bin/plotprofile.py profiling/profilerinput.txt_profiling.34_272/profilerinput.txt_profiling.Profiletest_1.mean.csv profiling/energy_plot.out_normal_s3dc_dd.png normal_s3dc_dd)

There is support for running jobs in parallel. In order to do this the relevant information in the config.ini Cluster field must be filled in.
After that simply run jobs with the parallel parameter
