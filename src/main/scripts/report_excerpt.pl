#!/usr/bin/env perl

use strict;
use warnings;

use DelimitedFile;
use Getopt::Long;
use File::Basename;

use FileUtils qw(read_simple_file);

my %FLAGS;

my $LIST_POSSIBLE = 1;
# treat entries in target column as possible lists

my @files;
GetOptions(
	   \%FLAGS,
	   "-file=s" => \@files,
	   # report file to parse
	   "-files=s",
	   # list of report files
	   "-glob=s",

	   "-column=s",
	   # column name to search
	   "-value=s",
	   # value to find
	   # FIX ME: if field is a list, etc.
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

foreach my $infile (@files) {
  printf STDERR "parsing %s...\n", $infile;
  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			    );

  my $column = $FLAGS{column};
  if ($column) {
    my $value = $FLAGS{value} || die "specify -value\n";

    my $outfile = sprintf '%s.%s.%s', basename($infile), $column, $value;

    my $rpt = $df->get_reporter("-file" => $outfile);

    # while (my $row = $df->next("-ref" => 1)) {  # headerless
    while (my $row = $df->get_hash()) {
      my $data = $row->{$column};
      $data = "" unless defined $column;
      my $usable;
      if ($LIST_POSSIBLE) {
	foreach (split /,/, $data) {
	  $usable =1 if $_ eq $value;
	}
      } else {
	$usable = 1 if $data eq $value;
	# TO DO: case sensitivity?
      }
      $rpt->end_row($row) if $usable;
    }
    $rpt->finish();
  } else {
    die sprintf "specify -column (%s)\n", join ",", @{$df->headers_raw};
  }
}
