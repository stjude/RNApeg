#!/bin/env perl
# create

use strict;
use warnings;

use DelimitedFile;
use Getopt::Long;
use File::Basename;

use FileUtils qw(read_simple_file);
use MiscUtils qw(dump_die);

my %FLAGS;

my @files;
GetOptions(
	   \%FLAGS,
	   "-file=s" => \@files,
	   # report file to parse
	   "-files=s",
	   # list of report files
	   "-glob=s",

	   "-min-genes=i",
	   # minimum count of annotated genes
	   # FIX ME: awkward, could be synonyms, etc. right on top of each other
	   # so not necessarily a real fusion

	   "-bidirectional=i",
	   # minimum count of reads observed on each strand

	   "-novel",
	   # novel junctions only

	   "-min-flanking=i",

	   "-min-perfect=i",
	   "-min-acceptable=i",

	   "-out=s",
	   # outfile

	  );


unless (@files) {
  if (my $l = $FLAGS{files}) {
    my $set = read_simple_file($l);
    @files = @{$set};
  } elsif (my $glob = $FLAGS{glob}) {
    @files = glob($glob);
    die "no matches" unless @files;
  } else {
    die "specify -file or -files [listfile]\n";
  }
}

my $min_genes = $FLAGS{"min-genes"};
my $bidirectional = $FLAGS{"bidirectional"};
my $novel = $FLAGS{novel};
my $min_flanking = $FLAGS{"min-flanking"};
my $min_perfect = $FLAGS{"min-perfect"};
my $min_acceptable = $FLAGS{"min-acceptable"};

die "can't use -out with multiple infiles" if $FLAGS{out} and @files > 1;

foreach my $infile (@files) {
  printf STDERR "parsing %s...\n", $infile;
  my $df = new DelimitedFile(
			     "-file" => $infile,
			     "-headers" => 1,
			    );

  my $outfile = $FLAGS{out} || sprintf '%s.excerpt.tab', basename($infile);
  my $rpt = $df->get_reporter("-file" => $outfile);

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {
    my $usable = 1;

    if ($min_genes) {
      # require minimum count of gene annotations for putative fusions
      my @g = split /,/, $row->{genes};
      $usable = 0 unless @g >= $min_genes;
    }

    if ($bidirectional) {
      # require junction to be observed on both strands with minimum read coverage
      $usable = 0 unless $row->{qc_plus} >= $bidirectional and $row->{qc_minus} >= $bidirectional;
    }

    $usable = 0 if $novel and $row->{type} ne "novel";
    $usable = 0 if $min_flanking and $row->{qc_flanking} < $min_flanking;
    $usable = 0 if $min_perfect and $row->{qc_perfect_reads} < $min_perfect;

    if ($min_acceptable) {
      my $acceptable = $row->{qc_perfect_reads} + $row->{qc_clean_reads};
      $usable = 0 if $acceptable < $min_acceptable;
    }

    $rpt->end_row($row) if $usable;
  }
  $rpt->finish();
}
