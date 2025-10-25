# Guide to Using Rmultiome
## Initial Setup
1. Create a directory for the project.
2. Inside the directory create a rawdata subdirectory, and inside this directory should be all the sample directories with the cellranger output. Follows are the only files this tool will actually be using from the cellranger output:
```
/projects/projectdir/sample1/atac_fragments.tsv.gz
/projects/projectdir/sample1/atac_fragments.tsv.gz.tbi
/projects/projectdir/sample1/filtered_feature_bc_matrix.h5
```
3. clone the Rmultiome code repo:  
`git clone https://github.com/brianlamere/Rmultiome`
4. Rstudio won't want to go to a location outside of your home directory, because they don't like security, projects worked on by more than one person, projects that survive past one particular person, etc.  The easy way around this is to make a symbolic link in your home directory, to the project directory.  From the command line, type a command such as:
`ln -s /projects ~/projects-link`
"ln" is the link tool, the flag "-s" instructs ln to make a symbolic link.  From there it is source, then target; source is what you are linking to, target is where the link will exist.  The ~ in the example target is a shell option that simply stands for whatever your home directory is, and works on nearly all UNIX variants including OSX.
5. I won't go much into rstudio usage here, but start rstudio, select to create a new project, browse to that link you created above, then from there select the project directory you used.  In my use I have a larger "/projects" directory that has all my projects in it; you might instead have linked directly to your project.
6. Once rstudio has started, in your file browser should be the "rawdata" and "Rmultiome" folders.  Open the Rmultiome/system_settings.R file, and modify the early line:  project_base_dir should be set to where your project is located (the folder that has rawdata and Rmultiome).  Modify to your project path location.

## QC Phase 1: standard/1 dimensional trimming
1. Open the run_qc.R file
2. modify the very first line of run_qc.R to point to where your system_settings.R file is.  I cannot, unfortunately, make this any cleaner other than needing you to modify this line yourself
3. Run/source the lines from the top of the file, through step 1-1, which should end by listing the files in your rawdatadir directory.
4. modify the "mysample" line in step 1-2 to be the sample you're starting with, given this is a multiomic multi-sample tool.  Note that you must call the sample the exact same thing it is called in the subdirectories in your rawdatadir (which is project_data_dir/rawdata, by default).
5.  After step 1-2, run step 1-3; if this completes without error, you're set!  Note it might have a warning, especially if you had a custom reference that used underscores like mine did and didn't realize they would not be allowed (and would be renamed using hyphens).
6.  Continue on with step 1-4, which will create 8 plots for the default base object.  There will likely be warnings, since seurat is using code in their code base that is flagged as deprecated in their code base, unless you're using a newer version of Seurat than I did at time of this writing.
7.  Against my better judgement but only because the file might be confusing otherwise, I left some settings in step 1-5.  They are mostly random settings for a sample I was using when I saved this file - your tissue type might be different, your settings will almost certainly be different, and each those settings should almost without a doubt be changed, based on the plots from step 1-4.  I will have a different guide for how to use the plots to come up with values, but this guide assumes you already know how to pick QC settings based on the plots. Note I am using nCount_ATAC and nCount_RNA, but not nFeature_ATAC or nFeature_RNA.  I included a plot that should show why - the RNA count vs RNA feature plot, which should be a very clean line for which there should be an obvious and clean linear regression.  If your RNA count vs RNA feature plot does NOT show a very clear linear regression, then you may need to re-evaluate your upstream work, as it should be a clean line like this:  [[insert image with our plot]]
8.  In step 1-6, you'll be applying the trimming settings.  First you're verifying the changes to the larger trimming_settings data frame, then you're updating the larger trimming_settings data frame, then you're subsetting based on these changes.  Note this will make more sense once you realize the larger data frame has the settings for all your samples, and will be what run_pipeline1.R uses to set up all the samples.
9. In step 1-7, you'll review the plots again, now with a trimmed object.  If you need to change the trimming settings (up or down) just repeat 1-5 to 1-7, until you feel comfortable with the QC for the sample.  This guide is not a guide to what the settings should be, it is a guide to how to use this tool, but keep in mind that just because a heat map might get near or past an edge that is cut off, doesn't mean the setting is wrong; there is a max percent.MT you'll want for your sample regardless, for example, no matter what the plot looks like.
10. Once you have locked down the trimming settings you want, proceed to step 1-8, where you will save the trimming settings to a project file.

Note you may find it easier to tell what the new plots are (and thus what the new settings result in) if you clear out the plots from rstudio, which in the version being used at the time of this writing looks like a little broom on the top of the plot window frame.

## QC phase 2: KDE trimming
1. First you initialize the kde_settings, as in 2-1 in the run_qc.R
2. second, you'll change the atac_percentile and rna_percentile in step 2-2, using the plots from the last time you ran 1-7 above as guides.
3. next you'll update the main kde_settings as in 2-3
4. in step 2-4 you'll visualize what your settings would keep from what was set in 2-2.
5. in the optional 2-5, you can compare what "intersection" versus "union" does.  Intersection is everything kept by *both* atac_percentile and rna_percentile, union is everything that is in either set.
6. from here you would repeat steps 2-2 to 2-5 until you feel good with the setting
7. finally, at step 2-6 you would save the sample's kde settings to the project settings directory so it can be used in run_pipeline1.R
