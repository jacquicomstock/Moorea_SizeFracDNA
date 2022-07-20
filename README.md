# Moorea_SizeFracDNA
All processing &amp; analyses from raw fastq files to analyses &amp; figures are contained in the above scripts.

Processing of all raw fastq files was done on a cluster through dada2, which was run in R. After denoising with dada2, ASV taxonomy was initially assigned using SILVA v138. SAR11, SAR202, and cyanobacterial lineages were then extracted and further classified using specialized annotated databases in Phyloassigner. Taxonomic assignments from SILVA and Phyloassigner were then merged into one taxonomy file and subsequently used for all remaining analyses. 

Below is described how I installed all programs to run the above scripts. 
#### This is written step-by-step for absolute beginners.

### 0) Access your cluster
First you need to sign into your cluster. If you are from UCSB and have an account with Pod, use the following command (replace "USERNAME" with whatever you're username is) and then type in your password when prompted. If you are not using campus wifi, you will have to use the campus VPN to be able to access the cluster.

```{bash}
ssh USERNAME@pod.cnsi.ucsb.edu
```

If you are using a Mac, you can type the above command (and all following commands) into Terminal. If you are using a PC, it is easiest if you download MobaXterm (https://mobaxterm.mobatek.net/download.html) and treat that as your "terminal". From this you can have a "local" terminal window where you access your own personal computer and a terminal where you ssh into your cluster. For almost all commands you will want to be writing them in your terminal window connected to the cluster, but if you are uploading/downloading files from or onto your local computer, those commands should be run in your local terminal.

### 1) Install conda 
This install conda section is copy/pasted directly from Fabian Wittmers (https://github.com/BIOS-SCOPE/PhyloAssigner_python_UCSB)

if you have conda already installed, e.g. if you used it to install a pipline like Qiime2 through conda, you can move on to step 2).

Download the installer for linux from: https://docs.conda.io/en/latest/miniconda.html#linux-installers

How? 
One easy way of achieving this can be the tool wget. Log into your linux machine / cluster. 
Most clusters do have 'wget' pre-installed, so this line should work if you just copy-paste it: 

```{bash}
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
```

This will download a file called 'Miniconda3-py39_4.12.0-Linux-x86_64.sh' into whatever directory you are currently in. 
Next, install conda:

```{bash}
bash Miniconda3-py39_4.12.0-Linux-x86_64.sh
```

After you have finished the installation process (you can press yes whenever conda asks you something during this), you should see "(base)" at the beginning of your commandline handle. If not, log out and in again. If you still don't see the "(base)" at the beginning of the command-line, try:

```{bash}
source ~/.bashrc
conda activate base
```

Now you should have conda installed and activated on your system. Conda is called a 'package manager', it does handle all the annoying parts of installing a bioinformatics tool on your machine for you, so you don't have to worry about all of PhyloAssigners dependencies and having to make sure that they have been installed with the correct version and path, etc.

### 2) Installing R and dada2 on the cluster
An updated R version (>=4.0) is available through the conda-forge channel, so first add the channel.
```{bash}
conda config --add channels conda-forge
```

We are going to install R through conda-forge and not through the default channel, so we need to set its priority over the default channel.
NOTE: You can undo this change by setting strict priority to the default channel as described a bit further below.
```{bash}
conda config --set channel_priority strict
```

Check whether an updated R version is added in the conda search space or not.
```{bash}
conda search r-base
```

Now, it is always a good practice (recommoned here) to create a new conda environment, which will help to debug the package-specific and compatibility issues without disrupting the base environment.
You can change the name to anything you want, here I have named it after the version of R I am installing (R4.1.3)
```{bash}
conda create -n R4.1
```

Let's activate the newly create conda environment.
```{bash}
conda activate R4.1
```

And finally install the R package.
```{bash}
conda install -c conda-forge r-base
```

To reset the default conda channel, run the following
```{bash}
conda config --set channel_priority true
conda update --all
```

To confirm that the default conda channel has been restored (should say channel_priority: flexible)
```{bash}
conda config --describe channel_priority
```

And finally, to start R in the command line, simply type "R" and enter, and you should run R.

Now we need to install the necessary packages. Namely, DADA2, since that is the primary reason we are running R on a cluster. You may have to change the BiocManager version to whatever is current. This may take several minutes.

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2", version = "3.14")
```

Now you're ready to start running dada2 on the cluster!

### 3) Installing, setting up, and running Phyloassigner

This was taken directly from Fabian Wittmers, who recreated Phyloassigner in Python (https://github.com/BIOS-SCOPE/PhyloAssigner_python_UCSB)

You will need to download Phyloassigner (and accompanying databases) from Fabian Wittmer's github (https://github.com/BIOS-SCOPE/PhyloAssigner_python_UCSB). You can download everything as a single zip file by clicking the green "code" button and clicking download. 
Once you have downloaded Phyloassigner you can upload that file from your local computer to your cluster. Open up a new local terminal window (One that's NOT connected to the cluster) and type the following command:

```{bash}
scp -r FILE USERNAME@SERVERADRESSE:/DIRECTORY/
```

Here is an example of how to fill in the above script:

```{bash}
scp -r PhyloAssigner_python_UCSB-main.zip carlsonlab@pod.cnsi.ucsb.edu:/home/carlsonlab/
```
Note that this zipped file MUST be in the directory you are currently in on your computer. So if your default directory is Desktop and your zip file is in your downloads folder, this command will spit out an error message.

PhyloAssigner handles alignment, reformatting, placement and output summary all within one command and therefore relies on some dependencies. I advice to run PhyloAssigner on a linux system, idealy a cluster with a good chunk of RAM, although it can run on low RAM systems just fine (will just take longer). It installation is simple, all dependencies can be installed in 1 command by creating a new conda environment to run PhyloAssigner in.

```{bash}
conda env create -f pythonassignler_linux.yml
```

This git contains a collection of placement reference trees build by different members of the WordenLab, Luis Bolanos, and Kevin Vergin.

Currently available databases (and the paper in which they were published) comprise:

**Global 16S**

Vergin et al. 2013; "High-resolution SAR11 ecotype dynamics at the Bermuda Atlantic Time-series Study site by phylogenetic placement of pyrosequences"

**SAR11**

Bolanos et al. 2021; "Seasonality of the Microbial Community Composition in the North Atlantic"
see: https://github.com/lbolanos32/NAAMES_2020

**SAR202**

Landy et al. 2017; "SAR202 Genomes from the Dark Ocean Predict Pathways for the Oxidation of Recalcitrant Dissolved Organic Matter"
see: https://github.com/lbolanos32/NAAMES_2020

**Chrysophyceae (16S plastid)**

Worden Lab; unpublished

**Dictyochophyceae (16S plastid)**

Choi et al. 2021; "Seasonal and Geographical Transitions in Eukaryotic Phytoplankton Community Structure in the Atlantic and Pacific Oceans"

**Pelagophyceae (16S plastid)**

Choi et al. 2021; "Seasonal and Geographical Transitions in Eukaryotic Phytoplankton Community Structure in the Atlantic and Pacific Oceans"

**Stramenopiles (16S plastid)**

Choi et al. 2021; "Seasonal and Geographical Transitions in Eukaryotic Phytoplankton Community Structure in the Atlantic and Pacific Oceans"

**Cyanobacteria**

Sudek et al. 2015; "Cyanobacterial distributions along a physico-chemical gradient in the Northeastern Pacific Ocean"

**Prochlorococcus**

Worden Lab; unpublished

**Cyanobacteria + Plastid**

Sudek et al. 2015; "Cyanobacterial distributions along a physico-chemical gradient in the Northeastern Pacific Ocean"
Choi et al. 2017; "Newly discovered deep-branching marine plastid lineages are numerically rare but globally distributed"

**Viridiplantae (16S plastid)**

Worden Lab; unpublished


To activate Phyloassigner, type the following command
```{bash}
conda activate pythonassigner
```
Change directories into phyloassigner
```{bash}
cd Phyloassigner/PhyloAssigner_python_UCSB-main/
```
and print the phyloassigner help message
```{bash}
python pythonassigner_v0.9.py -h
```

Now Phyloassigner is set up and ready for use! Below is a command example written by Fabian:
#### Command Example
PhyloAssigner requires 3 reference files, that are part of each PhyloAssigner database in this git:

reference tree
reference phylogeny to place the ASVs on
reference alignment
alignment that was used to reconstruct the reference tree
reference mapping
connects each edge in the tree with its corresponding label that shall be printed in the results
Those 3 input files are provided. In addition, the user is required to provide an output directory. Optional arguments include the number of threads to use (for the alignment of query and reference sequences) and what placement algorithm to use. The default placement algorithm is PPLACER, so the user does not have to specify it.

```{bash}
python pythonassigner_v0.9.py \
    --out_dir example_output/ \
    --ref_align databases/Global_16S_refDB/ref.aln \
    --ref_tree databases/Global_16S_refDB/ref_tree.txt \
    --query_seqs your_ASVs.fasta \
    --mapping databases/Global_16S_refDB/edge.mapping \
    --threads 32 \
    --placer pplacer
```

Above is all the code that is needed to run phyloassigner. You can change the reference database names based on which you plan to use (eg global vs SAR11, etc). You should also change the name of the fasta file.
To save this chunk of code as a script, you can type the following into the command line:
```{bash}
vi runphylo.sh
```
This will create a new shell script (text file). To edit this file, press "i", paste the above code into it, then pres "Esc", followed by ":w" and Enter, then followed by ":q" and Enter. This will write the file and quite the text editor (vi).
If you type "ls" into the command line to see what files/directories are listed, you will see runphylo.sh pop up.

To run all the commands in this new shell script you just wrote, type the following:
```{bash}
bash runphylo.sh
```
