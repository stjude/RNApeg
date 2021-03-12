#!/usr/bin/env perl
# remove novel junctions present in other junction files

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die build_argv_list get_hash_option);
use DelimitedFile;
use Reporter;
use FileUtils qw(read_simple_file);
use JunctionUtils qw(parse_junction);

my %FLAGS;

my $QC_MIN_NOVEL_READS = 0;
my $QC_MIN_FLANKING = 0;
my $QC_BIDIRECTIONAL = 0;
my $QC_MIN_PERFECT = 0;
my $QC_MIN_GOOD = 0;

my @clopts = (
	      "-file=s",
	      "-files=s",
	      "-glob=s",

	      "-blacklist-files=s",
	      "-force",
	      "-verbose",

	      "-min-novel-reads=i" => \$QC_MIN_NOVEL_READS,
	      "-qc-min-flanking=i" => \$QC_MIN_FLANKING,
	      "-qc-bidirectional=i" => \$QC_BIDIRECTIONAL,
	      "-qc-min-perfect=i" => \$QC_MIN_PERFECT,
	      "-qc-min-good=i" => \$QC_MIN_GOOD
	      # minimum QC required to blacklist, i.e. only do it when
	      # given strong evidence.

	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $infiles = build_argv_list("-flags" => \%FLAGS,
			    "-single" => "file",
			    "-set" => "files",
			    "-glob" => "glob",
			   );

die "no infiles" unless @{$infiles};

my $blacklist_files = read_simple_file($FLAGS{"blacklist-files"} || die);
my %blacklist;
printf STDERR "building blacklist...\n";
my $verbose = $FLAGS{verbose};
foreach my $bf (@{$blacklist_files}) {
  printf STDERR "  %s\n", basename($bf);
  my $df = new DelimitedFile("-file" => $bf,
			     "-headers" => 1,
			     );
  while (my $row = $df->get_hash()) {
    my $blacklist = $row->{type} eq "novel" ? 1 : 0;

    if ($blacklist and $QC_BIDIRECTIONAL and
	($row->{qc_plus} == 0 or $row->{qc_minus} == 0)
       ) {
      $blacklist = 0;
    }

    if ($blacklist and $QC_MIN_NOVEL_READS and
	$row->{count} < $QC_MIN_NOVEL_READS
       ) {
      $blacklist = 0;
    }

    if ($blacklist and $QC_MIN_FLANKING and
	$row->{qc_flanking} < $QC_MIN_FLANKING
       ) {
      $blacklist = 0;
    }

    if ($blacklist and $QC_MIN_PERFECT and
	$row->{qc_perfect_reads} < $QC_MIN_PERFECT
       ) {
      $blacklist = 0;
    }

    if ($blacklist and $QC_MIN_GOOD) {
      my $count_good = $row->{qc_perfect_reads} + $row->{qc_clean_reads};
      if ($count_good < $QC_MIN_GOOD) {
	$blacklist = 0;
      }
    }

    if ($blacklist) {
      my $j = $row->{junction};
      printf STDERR "blacklist %s\n", $j if $verbose and !$blacklist{$j};
      $blacklist{$j} = 1;
    }
  }
}

printf STDERR "blacklisted junctions: %d\n", scalar keys %blacklist;

foreach my $infile (@{$infiles}) {
  my $outfile = basename($infile) . ".filtered.tab";

  if (-s $outfile and !$FLAGS{force}) {
    printf STDERR "%s exists, skipping\n", $outfile;
    next;
  }

  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
  			      "-auto_qc" => 1,
			     );
  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  my %filtered;
  while (my $row = $df->get_hash()) {
    my $j = $row->{junction} || die;

    if ($blacklist{$j}) {
      my ($start, $end) = parse_junction($j);
      $filtered{$start->{ref}}++;
    } else {
      $rpt->end_row($row);
    }
  }
  printf STDERR "dropped:\n";
  foreach my $chr (sort keys %filtered) {
    printf STDERR "  %s: %d\n", $chr, $filtered{$chr};
  }


  $rpt->finish();
}

