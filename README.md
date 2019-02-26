SerumQC (Now deprecated and using bifrost)
=======

SerumQC is a pipeline for quality control of Whole Genome Sequencing (WGS) of bacterial samples used for surveillance. The pipeline attempts to check that samples have enough depth for future analysis and is contaminant free in a quick turn around in order to retrieve new samples if neccesary. The resulting QC report attempts to suggest if a sample is ok for further use or requires more depth or requires additional input.

Currently to run SerumQC please contact kimn@ssi.dk as an installation script has not yet been set up for handling setting up additional links and the modifications required for serum.config. SerumQC requires either the use of Slurm or Torque grid engine as well as the following list of software:
* SPAdes
* bbnorm
* kraken
* mlst
* trimmomatic
* ariba
* bwa
* elprep
* samtools
* pilon

More information will be available at a later date in regards to background and workflow.

Please contact kimn@ssi.dk if you have any questions.
