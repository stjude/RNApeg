#!/bin/env perl
# downgrade "known" junctions to "novel" if evidence is exclusive to
# a database we don't want to use for aberrant classification (e.g. AceView)
# Do this AFTER standard annotation/correction:
# - edge correction will be performed (important)
# - annotations will still remain if we want them
# MNE 5/2016

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use FileUtils qw(read_simple_file);
use MiscUtils qw(dump_die build_argv_list get_hash_option);
use DelimitedFile;
use Reporter;
use RefFlatFile;
use Counter;

my %FLAGS;
my @clopts = (
	      "-refflat=s",
	      "-aceview",

	      "-file=s",
	      "-files=s",
	      "-glob=s",
	     );

GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $files;
if (my $f = $FLAGS{files}) {
  $files = read_simple_file($f);
} elsif (my $one = $FLAGS{file}) {
  $files = [ $one ];
} elsif (my $p = $FLAGS{"glob"}) {
  $files = [ glob($p) ];
}
die "no files" unless $files and @{$files};

my $f_refflat = $FLAGS{refflat} || die "-refflat";
my $rf = new RefFlatFile();

my $rf_format = $FLAGS{aceview} ? "aceview" : "refgene";
$rf->parse_file(
		"-refflat" => $f_refflat,
		"-type" => $rf_format
	       );
my %blacklist;
foreach my $r (@{$rf->rows}) {
  $blacklist{$r->{name} || die} = 1;
}

my $c = new Counter($files);
foreach my $infile (@{$files}) {
  my $outfile = basename($infile) . ".downgrade.tab";

  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
  			      "-auto_qc" => 1,
			     );

  while (my $row = $df->get_hash()) {
    if ($row->{type} eq "known") {
      my %transcripts = map {$_, 1} split /,/, $row->{transcripts};
      foreach my $t (keys %transcripts) {
	delete $transcripts{$t} if $blacklist{$t};
      }
      unless (%transcripts) {
	# known call is made exclusively from blacklisted source,
	# downgrade to novel
	$row->{type} = "novel";
      }
    }
    $rpt->end_row($row);
  }
  $rpt->finish();

  $c->next($infile);
}
