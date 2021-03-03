- Docker setup notes:
  - src/javalib: rnapeg Java code and Picard support library
  - src/perllib, util-perllib: 
    - some Perl module files from those svn packages,
    - SampleName.pm from samplename package
  - src/bin: scripts from rnapeg, plus RNApeg.sh Docker wrapper script
  - src/perllib/rnapeg
  => is there a better/cleaner way to do this?  PROBABLY!  (MNE 6/2020)

- directory listing from PC Docker setup:

/home/medmonso/rnapeg/
/home/medmonso/rnapeg/RNApeg
/home/medmonso/rnapeg/RNApeg/README.md
/home/medmonso/rnapeg/RNApeg/Dockerfile
/home/medmonso/rnapeg/RNApeg/configs
/home/medmonso/rnapeg/RNApeg/configs/genome
/home/medmonso/rnapeg/RNApeg/configs/genome/GRCh37-lite.config.txt
/home/medmonso/rnapeg/RNApeg/configs/genome/GRCh37-lite.template.txt
/home/medmonso/rnapeg/RNApeg/Dockerfile.alpine
/home/medmonso/rnapeg/RNApeg/src
/home/medmonso/rnapeg/RNApeg/src/bin
/home/medmonso/rnapeg/RNApeg/src/bin/floating_junction_fix.pl
/home/medmonso/rnapeg/RNApeg/src/bin/junction2gene_flat.pl
/home/medmonso/rnapeg/RNApeg/src/bin/junction2gene.pl
/home/medmonso/rnapeg/RNApeg/src/bin/bam_junction.pl
/home/medmonso/rnapeg/RNApeg/src/bin/bak
/home/medmonso/rnapeg/RNApeg/src/bin/bak/RNApeg.sh
/home/medmonso/rnapeg/RNApeg/src/bin/RNApeg.sh
/home/medmonso/rnapeg/RNApeg/src/bin/junction_extraction_wrapper.pl
/home/medmonso/rnapeg/RNApeg/src/bin/debug.sh
/home/medmonso/rnapeg/RNApeg/src/javalib
/home/medmonso/rnapeg/RNApeg/src/javalib/picard-2.6.0.jar
/home/medmonso/rnapeg/RNApeg/src/javalib/rnapeg-2.6.0.jar
/home/medmonso/rnapeg/RNApeg/src/perllib
/home/medmonso/rnapeg/RNApeg/src/perllib/util-perllib
/home/medmonso/rnapeg/RNApeg/src/perllib/util-perllib/ClusterLogFile.pm
/home/medmonso/rnapeg/RNApeg/src/perllib/util-perllib/Configurable.pm
/home/medmonso/rnapeg/RNApeg/src/perllib/util-perllib/WorkingFile.pm
/home/medmonso/rnapeg/RNApeg/src/perllib/util-perllib/DBTools.pm
/home/medmonso/rnapeg/RNApeg/src/perllib/util-perllib/FileUtils.pm
/home/medmonso/rnapeg/RNApeg/src/perllib/util-perllib/Cluster.pm
/home/medmonso/rnapeg/RNApeg/src/perllib/util-perllib/MethodMaker.pm
/home/medmonso/rnapeg/RNApeg/src/perllib/util-perllib/TdtConfig.pm
/home/medmonso/rnapeg/RNApeg/src/perllib/util-perllib/Reporter.pm
/home/medmonso/rnapeg/RNApeg/src/perllib/util-perllib/BucketMap.pm
/home/medmonso/rnapeg/RNApeg/src/perllib/util-perllib/DelimitedFile.pm
/home/medmonso/rnapeg/RNApeg/src/perllib/util-perllib/MiscUtils.pm
/home/medmonso/rnapeg/RNApeg/src/perllib/util-perllib/JavaRun.pm
/home/medmonso/rnapeg/RNApeg/src/perllib/util-perllib/CommandLineRebuilder.pm
/home/medmonso/rnapeg/RNApeg/src/perllib/common-perllib-genome
/home/medmonso/rnapeg/RNApeg/src/perllib/common-perllib-genome/ReferenceNameMapper.pm
/home/medmonso/rnapeg/RNApeg/src/perllib/common-perllib-genome/GenomeUtils.pm
/home/medmonso/rnapeg/RNApeg/src/perllib/common-perllib-genome/GeneAnnotationFlat.pm
/home/medmonso/rnapeg/RNApeg/src/perllib/common-perllib-genome/SampleTumorNormal.pm
/home/medmonso/rnapeg/RNApeg/src/perllib/common-perllib-genome/FAI.pm
/home/medmonso/rnapeg/RNApeg/src/perllib/common-perllib-genome/GeneAnnotation.pm
/home/medmonso/rnapeg/RNApeg/src/perllib/common-perllib-genome/SampleName.pm
/home/medmonso/rnapeg/RNApeg/src/perllib/common-perllib-genome/GeneListCollapser.pm
/home/medmonso/rnapeg/RNApeg/src/perllib/rnapeg
/home/medmonso/rnapeg/RNApeg/src/perllib/rnapeg/FloatingJunctionDetector.pm
