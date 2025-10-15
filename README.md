# Rmultiome
A small tool set that is being designed for 10X multiome data from brain nuclei

Rebased after learning many of the vignettes that seem to apply, don't really
for my particular data set (very high depth, brain nuclei, thus very strong 
batch effects but no need to do things like redoing peaks, etc).

#Guide to Rmultiome
##Initial Setup
1. Create a directory for the project.
2. Inside the directory create a rawdata subdirectory, and inside this directory should be all the sample directories with the cellranger output. Follows are the only files this tool will actually be using from the cellranger output:
```/projects/projectdir/sample1/atac_fragments.tsv.gz
/projects/projectdir/sample1/atac_fragments.tsv.gz.tbi
/projects/projectdir/sample1/filtered_feature_bc_matrix.h5```
3. clone the Rmultiome code repo:  
`git clone https://github.com/brianlamere/Rmultiome`
4. Rstudio won't want to go to a location outside of your home directory, because they don't like security, projects worked on by more than one person, projects that survive past one particular person, etc.  The easy way around this is to make a symbolic link in your home directory, to the project directory.  From the command line, type a command such as:
`ln -s /projects ~/projects-link`
"ln" is the link tool, the flag "-s" instructs ln to make a symbolic link.  From there it is source, then target; source is what you are linking to, target is where the link will exist.  The ~ in the example target is a shell option that simply stands for whatever your home directory is, and works on nearly all UNIX variants including OSX.
5. I won't go much into rstudio usage here, but start rstudio, select to create a new project, browse to that link you created above, then from there select the project directory you used.  In my use I have a larger "/projects" directory that has all my projects in it; you might instead have linked directly to your project.
6. Once rstudio has started, in your file browser should be the "rawdata", "references," and "Rmultiome" folders.  Open the Rmultiome/system_settings.R file, and modify the early line:  project_base_dir should be set to where your project is located (the folder that has rawdata, references, and Rmultiome).  Modify to your project path location.

##QC Phase


