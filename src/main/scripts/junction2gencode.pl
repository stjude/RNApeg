#!/bin/env perl
# annotate junction files re: "primary" ENSEMBL transcripts
# MNE 2/2015
#
# QUESTIONS:
# - strip version numbers from transcripts?

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die unquote build_argv_list);
use DelimitedFile;
use SJPreferredIsoform;
use TdtConfig;

my %FLAGS;

my @GTF;

GetOptions(\%FLAGS,
	   "-genome=s",

	   "-gtf=s" => \@GTF,
	   # GTF file(s) containing transcript annotations
	   # e.g.
	   # /nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/mRNA/GENCODE/gencode.v19.annotation.gtf
	   # /nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/mRNA/GENCODE/gencode.v19.long_noncoding_RNAs.gtf

	   "-sjpi=s",
	   # preferred isoform matrix
	   # e.g.
	   # /nfs_exports/apps/gnu-apps/NextGen/nextgensupport2/NCBI/RefSeq/gene_transcript_matrix.withNM.mod

	   "-file=s",
	   "-files=s",
	  );

my $files = build_argv_list("-flags" => \%FLAGS, "-single" => "file", "-set" => "files");

my $genome = $FLAGS{genome} || die "-genome";
my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";

if (@GTF) {
  print STDERR "using manually-specified GTF files\n";
} else {
  push @GTF, $config_genome->{GENCODE_FULL_GTF} || die;
  push @GTF, $config_genome->{GENCODE_NONCODING_GTF} || die;
}

my $f_sjpi = $FLAGS{sjpi} || $config_genome->{GENE_TRANSCRIPT_MATRIX} || die;
my $sjpi = new SJPreferredIsoform("-file" => $f_sjpi);
my $nm2preferred = $sjpi->nm2preferred();

my $gtf_primary;
if (0) {
  print STDERR "DEBUG: no GTF\n";
  $gtf_primary = {};
} else {
  $gtf_primary = parse_gtf();
}

foreach my $infile (@{$files}) {
  printf STDERR "processing %s...\n", $infile;
  my $outfile = basename($infile) . ".gencode.tab";
  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => [
					   qw(
					       gencode_primary
					       gencode_unknown
					       sj_primary
					       sj_unknown
					    )
					  ]
			     );

  while (my $row = $df->get_hash()) {
    my @primary;
    my @unknown;
    my @sj_primary;
    my @sj_unknown;

    if (my $ts = $row->{transcripts}) {
      foreach my $t (split /,/, $ts) {
	if ($t =~ /^ENST/) {
	  die "$t unversioned" unless $t =~ /\.\d+$/;
	  my $p = $gtf_primary->{$t};
	  if (defined $p) {
	    # data available (yea or nay)
	    push @primary, $t if $p;
#	    die "alternate $t" unless $p;
	  } else {
	    push @unknown, $t;
	  }
	} elsif ($t =~ /NM_/) {
	  my $p = $nm2preferred->{$t};
	  if (defined $p) {
	    push @sj_primary, $t if $p;
	  } else {
	    push @sj_unknown, $t;
	  }
	}
      }
    }

    $row->{gencode_primary} = join ",", @primary;
    $row->{gencode_unknown} = join ",", @unknown;
    $row->{sj_primary} = join ",", @sj_primary;
    $row->{sj_unknown} = join ",", @sj_unknown;

    $rpt->end_row($row);
  }
  $rpt->finish();
}


sub parse_gtf {
  die "specify -gtf [file] [-gtf ...]\n" unless @GTF;

  my @fields = qw(
		   seqname
		   source
		   feature
		   start
		   end
		   score
		   strand
		   frame
		   group
		);
  # again, feh: standard parser? or module

  my %primary;

  foreach my $gtf (@GTF) {
    printf STDERR "parsing %s...", $gtf;
    open(GTF, $gtf) || die;
    my $lines = 0;
    while (<GTF>) {
      next if /^#/;
      chomp;
      my @f = split /\t/, $_;
      die unless @f == @fields;
      my %r;
      @r{@fields} = @f;

      if (++$lines % 200000 == 0) {
	printf STDERR "%d%%...", tell(GTF) * 100 / -s $gtf;
      }

      if ($r{feature} eq "transcript") {
	my %attr = map {split /\s+/, $_} split /;\s*/, $r{group};
	foreach (values %attr) {
	  $_ = unquote($_);
	}
	$r{attributes} = \%attr;

	my $rname = $attr{transcript_name} || die;
	my $primary = $rname =~ /\-001/ ? 1 : 0;

	my $tid = $attr{transcript_id} || die;

	if (exists $primary{$tid}) {
	  my $old = $primary{$tid};
	  die "conflict for $tid" unless $old == $primary;
	}
	$primary{$tid} = $primary;
      }
    }
    print STDERR "done\n";
  }

  return \%primary;
}
