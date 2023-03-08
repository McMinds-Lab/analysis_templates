Analysis templates
===
Standard lab protocols, recommendations, and a collection of scripts for basic data processing, analysis and visualization, to be modified as necessary for different projects

Organization
===
File structure (Mise en place)
___
This is something I'm actively optimizing and learning about, so my strategy might not exactly be a gold strategy, but it is extremely important to at least be thinking about!

Primary principles:

* Keep a 'raw data' directory, where files are considered completely immutable for a given project or analysis. I waver between having a system-wide raw data directory that contains multiple project sub-directories, or a raw data directory within each project directory. But the point is that there should be a place where files belong that can be considered 'original' and that have READMEs associated with them to explain where these data came from. This directory could ideally be re-generated by unzipping a version archived somewhere such as Zenodo.
    - 'intermediate data', such as processed versions of your data like filtered sequencing reads or 'cleaned' metadata files, should be generated from the 'raw' files using scripts, which are contained in a separate:
* 'Scripts' directory. This directory is ideally a project-specific Git repository that can track changes and be collaborated on through GitHub, where an archived version can be created upon manuscript submission. These scripts should ideally be cross-platform compatible, use relative filepaths, and take command line arguments, such that someone can simply download the raw data, clone this repository, specify local settings and the location of an output directory, and completely re-generate the results of your study by running your scripts as specified in a README.
    - Try to organize things around groups of 7 or less - e.g. if you have 20 individual scripts, merge some of them, tie them together with 7 wrappers, or organize them in 7 subfolders 
    - Order scripts by starting their filenames with a 2-digit number (e.g. 01_sequencing_qc.sh)
    - Have each script place all its outputs (including log files) in a folder with the same name as the script, within the overall output directory specified by the user

Resources
===
Here's a list of resources my students and I have found useful for learning the science and code:

GitHub itself:
===
https://osf.io/preprints/metaarxiv/x3p2q/

Stats
===
General principles and great examples
___
https://betanalpha.github.io/

Linear models / GLMMs
___
https://trialsjournal.biomedcentral.com/articles/10.1186/s13063-022-06097-z

Good analysis of how stats are related to causal inference
___
[Inferring Multiple Causality: The Limitations of Path Analysis](https://www.jstor.org/stable/pdf/2389934.pdf?refreqid=excelsior%3A510214cf80f9fb879aaee87301ea1e5d&ab_segments=&origin=&initiator=&acceptTC=1)

Coding
===
In R
___
https://nyu-cdsc.github.io/learningr/assets/simulation.pdf  
https://intro2r.com/  
https://www.codecademy.com/learn/learn-r  

In Stan
___
https://betanalpha.github.io/  

In Markdown
___
https://www.markdownguide.org/basic-syntax/  

