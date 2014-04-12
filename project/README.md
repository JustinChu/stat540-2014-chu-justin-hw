Stat 540 - Metagenomics Reflections
==========================
Concise summary of main role / contribution of each group member:
* Jake: Responsible for running alignments in parallel and part of the read length compensation methodology (creating simulated dataset).
* Justin: Responsible for extracting data from RapSearch2 output by writing scripts for extracting and compiling data from alignments (read counts, residues aligned, etc) into usable form and pipeline/workflow planning.
* Thuy: Responsible for obtaining differentially abundant genes by writing R scripts and part of the read length compensation methodology (R scripting and use of blocking with limma).
* Craig: Responsible for GO term analysis and performed majority of poster formatting.
Evan: Responsible for obtaining gene list from MetaPathways by writing scripts to extract genes from output and sample correlation analysis.
* Kat: Responsible for taxonomy analysis using MEGAN and much of the our poster text regarding biological interpretation.

More detailed description and reflections on your specific role / contribution:
Initially, my role consisted of helping to plan out our analysis methodology and was thus responsible for creating a pipeline that we would model our analysis on. I was later tasked with extracting information from the RAPsearch2 output. I had to write my own scripts from scratch because of the non-standard blast-like output that it produces. In addition the output was fragmented (for parallelism) and very large. I had to code to both parallelize and merge my data extraction jobs to deal the size of the data and short time frames. I also wrote code to use only our reduced gene list as our initial output was too large for R to handle. This had to be done in a very short period of time since I required all alignments for a specific library to finish before generating counts for it and the fact that downstream analysis was dependent on the output I provided.

Observations about the group dynamic:
Everyone in our group had done their sections as well as they could given the time and computational constraints that had given the size of our dataset. There were at times some miscommunications (eg. forgetting to change permissions on some hidden data files) that could have been discovered earlier if we had better communication but otherwise we worked well together. Jake was instrumental in getting the ball rolling with regards to creating methodology to evaluate read length correction. Thuy spent a lot of time implementing a lot of R code by herself. Craig spent a lot of time editing our poster.

Scientific reflections:
Processing Issues: We were unable to process all the read data due to it size. We tried to parallelize the work but it still ended up taking too much time and space. Although we eventually reduced our gene set size we should have done it earlier and with possibly fewer genes. A good way would be to not only filter by what genes had GO terms but also by the GO terms themselves (eg. remove all genes with only “putative uncharacterized protein” as their term).
Read Length issues: Though we corrected for read lengths in our linear model, and validated it, we were not able to fix it when we performed our sample correlation analysis. This meaning our correlation observed between samples in this part of our analysis may have been due to read length difference.
Differential Abundance analysis: Using methods similar to differential expression or methylation analysis, we identified many genes as being significantly different in our analysis. As there were so many, other method relating to ranking the hits may have been better, especially when analyzing the GO term analysis that we performed.




