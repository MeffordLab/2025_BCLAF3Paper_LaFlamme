# commands for running repeat analysis tools and generating waterfall plots
# load repeat-analysis-tools: https://github.com/PacificBiosciences/apps-scripts/tree/master/RepeatAnalysisTools
module load repeat-analysis-tools

# set variables
BAM="./sample.bam" # bam
REF="./hg38.analysisSet.fa" # reference fasta

# waterfall plot for BCLAF3 repeat
python extractRegion.py $BAM $REF 'chrX:19990538-19991358' | python waterfall.py -m CCG,CCA -o waterfallBCLAF3.png # +/- 385bp

# waterfall plot for FMR1 repeat
python extractRegion.py $BAM $REF 'chrX:147911666-147912495' | python waterfall.py -m CGG -o waterfallFMR1.png # +/- 385bp