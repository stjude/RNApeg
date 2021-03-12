#!/usr/bin/env perl
# extract results for individual genes from rnapeg output,
# in tab-delimited and .bed format
# MNE 7/2014

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die);
# use DelimitedFile;
# use Reporter;

my %FLAGS;
my @GENES;
GetOptions(\%FLAGS,
	   "-gene=s" => \@GENES,
	   "-glob=s",
	   # glob for cross-sample-corrected junction file
	  );

die "need -gene" unless @GENES;

my @files = glob($FLAGS{glob} || die "-glob");
die "no files match" unless @files;

foreach my $jf (@files) {
  foreach my $gene (@GENES) {

    my $excerpt_file = sprintf '%s.genes.%s', basename($jf), $gene;

    unless (-s $excerpt_file) {
      my $cmd_excerpt = sprintf 'report_excerpt.pl -file %s -column genes -value %s', $jf, $gene;
      system $cmd_excerpt;
      die if $?;
    }

    my $bed_file = sprintf '%s_refgene_shade.bed', $excerpt_file;
    unless (-s $bed_file) {
      my $cmd_bed = sprintf 'junction2bed.pl -jf %s -now', $excerpt_file;
      system $cmd_bed;
      die if $?;
    }
  }
}



