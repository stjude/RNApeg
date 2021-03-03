#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Basename;

use FileUtils qw(read_simple_file);
use MiscUtils qw(dump_die);
use JavaRun;
#use SJSupport;
use TdtConfig;
use List::Util qw(min max);

use Cluster;

use constant DEFAULT_GENOME => "GRCh37-lite";

my %FLAGS;
my $FORCE = 1;
my @REGIONS;
GetOptions(
	   \%FLAGS,
	   "-bam=s",
	   "-bams=s",
	   "-now",

	   "-genome=s",
	   "-type=s",

	   "-refflat-type=s",
	   "-refflat=s",
	   # override

	   "-bed",

	   "-annotate",
	   "-debug",
	   "-rgb-known=s",
	   "-rgb-novel=s",
	   "-force=i" => \$FORCE,
	   "-primary-only",
	   "-no-duplicates",

	   "-min-span=i",
	   # minimum span to report a junction

	   "-out=s",
	   # outfile for single bam only

	   "-region=s" => \@REGIONS,
	   "-2bit=s",
	   "-fasta=s",
	   "-strand=s",
	   "-ignore-incompatible",
	   "-no-config",
	  );

my $NO_CONFIG_MODE = $FLAGS{"no-config"};
# primarily used when running in Docker

#
# reference genome sequence:
#
my $fasta = $FLAGS{"fasta"};
# preferred
my $twobit = $FLAGS{"2bit"};
my $refflat_file = $FLAGS{refflat};

my ($genome, $config_genome);
if ($NO_CONFIG_MODE) {
  die "need -fasta or -twobit" unless $fasta or $twobit;
  die "-refflat" unless $refflat_file;
} else {
  $genome = $FLAGS{genome} || DEFAULT_GENOME;
  $config_genome = TdtConfig::readConfig('genome', $genome) || die "no genome config for $genome";
  $fasta = $config_genome->{FASTA} unless $fasta;
  $twobit = $config_genome->{TWOBIT} unless $twobit;
}

#
# gene models:
#
unless ($refflat_file) {
  my $refflat_type = $FLAGS{"refflat-type"} || die "-refflat-type";
  unless ($refflat_type eq "all") {
    die "-refflat-type must be \"all\", subtypes no longer supported";
    # can be added if needed, but will require config file changes
  }

  #my $sjs = new SJSupport("-genome" => $genome);
  #my $refflat_file = $sjs->get_refflat_file("-type" => $refflat_type);
  $refflat_file = $config_genome->{COMBINED_REFFLAT} || die "no COMBINED_REFFLAT in configs for genome $genome";
}

#my $t2g_file = $sjs->get_refflat_t2g_helper("-type" => $refflat_type);
my $t2g_file;
# disabled since "all" is now the only option

my $type = $FLAGS{type};
die "specify -type [all|novel|reference]" unless $type and ($type eq "all" or $type eq "novel" or $type eq "reference");

collapse_regions();
# combine overlapping query regions. This must be done or else
# duplicate junctions will be reported in the output.
# This may still happen if regions are very close (i.e. within
# the range of a mapped read).

my $files;
if (my $bam = $FLAGS{bam}) {
  $files = [ $bam ];
} else {
  $files = read_simple_file($FLAGS{bams} || die "specify -bam [file] or -bams [listfile]");
}

my %outfiles;
foreach my $bam (@{$files}) {
  die "where is $bam" unless -s $bam;
  my $outfile;
  if ($outfile = $FLAGS{out}) {
    die "-out only works with a single BAM" unless @{$files} == 1;
  } else {
    $outfile = sprintf '%s.junctions.tab', basename($bam);
    die "duplicate outfile $outfile for $bam" if $outfiles{$outfile};
    # e.g. if duplicate basename, already handled if using
    # junction_extraction_wrapper.pl
    $outfiles{$outfile} = $bam;
  }

  my $jr = new JavaRun();
  $jr->ram("5g");
  $jr->classname("org.stjude.compbio.rnapeg.SplicedReadReporter");

  my $cl = sprintf "-bam %s -of %s", $bam, $outfile;
  if ($type eq "novel" or
      $type eq "reference" or
      $FLAGS{annotate} or
      $FLAGS{bed}
     ) {
    $cl .= sprintf " -refflat %s", $refflat_file;
    $cl .= sprintf " -t2g %s", $t2g_file if $t2g_file;
  }

  if ($fasta) {
    $cl .= sprintf ' -fasta %s', $fasta;
  } elsif ($twobit) {
    $cl .= sprintf ' -2bit %s', $twobit;
  } else {
    die "ERROR: need -fasta or -2bit reference sequence";
  }
  $cl .= " -reference-only" if $type eq "reference";
  $cl .= " -annotate" if $FLAGS{annotate};
  $cl .= sprintf ' -rgb-known "%s"', $FLAGS{"rgb-known"} if $FLAGS{"rgb-known"};
  $cl .= sprintf ' -rgb-novel "%s"', $FLAGS{"rgb-novel"} if $FLAGS{"rgb-novel"};
  if ($FLAGS{bed}) {
    $cl .= " -bed";
    $cl .= " -annotate" unless $type eq "novel" or $type eq "reference";
    # if junctions file used, code assumes we want to filter
    # unless annotate mode is enabled
  }

  if ($FLAGS{"primary-only"}) {
    $cl .= " -primary-only";
  }

  $cl .= " -no-duplicates" if $FLAGS{"no-duplicates"};

  $cl .= sprintf " -min-span %d", $FLAGS{"min-span"} if $FLAGS{"min-span"};
  $cl .= sprintf " -strand %s", $FLAGS{strand} if $FLAGS{strand};

  foreach my $r (@REGIONS) {
    $cl .= sprintf ' -region %s', $r;
  }
  $cl .= " -ignore-incompatible" if $FLAGS{"ignore-incompatible"};

  my $cmd = $jr->run(
		     "-command" => $cl,
		     "-return" => 1
		    );

  if ($FLAGS{now}) {
    if (-s $outfile and !$FORCE) {
      print STDERR "skipping, $outfile\n";
    } else {
      print STDERR "$cmd\n";
      system $cmd;
    }
  } elsif ($FLAGS{debug}) {
    print STDERR "not running $cmd\n";
  } else {
    my $c = new Cluster();

    $c->node_class("");
    # disable idataplex requirement

    $c->memory_reserve_mb(8000);
    $c->memory_limit_mb(8000);
    $c->outfile($outfile);
    $c->project("PCGP");
    $c->command($cmd);
    $c->run();
  }
}


sub collapse_regions {
  if (@REGIONS) {

    while (1) {
      # this is hella inefficient, something like Algorithm::Combinatorics
      # would be better (should it be installed)
      my @remove;
      my $add;
    SCAN:
      for (my $i = 0; $i < @REGIONS; $i++) {
	for (my $j = 0; $j < @REGIONS; $j++) {
	  next if $i == $j;

	  my $r1 = $REGIONS[$i];
	  my $r2 = $REGIONS[$j];

	  my ($r1_ref, $r1_start, $r1_end) = parse_region($r1);
	  my ($r2_ref, $r2_start, $r2_end) = parse_region($r2);

	  if ($r1_ref eq $r2_ref) {
	    if ($r1_start > $r2_end or
		$r2_start > $r1_end or
		$r1_end < $r2_start or
		$r2_end < $r1_start) {
	      # no overlap
	    } else {
	      printf STDERR "overlap %s %s\n", $r1, $r2;
	      $add = sprintf '%s:%d-%d',
		$r1_ref,
		  min($r1_start, $r2_start),
		    max($r1_end, $r2_end);
	      push @remove, $r1, $r2;
	      last SCAN;
	    }
	  }
	}
      }

      if ($add) {
	# overlap found
	foreach my $r (@remove) {
	  @REGIONS = grep {$_ ne $r} @REGIONS;
	}
	push @REGIONS, $add;
	# continue
      } else {
	last;
      }
    }
  }
}

sub parse_region {
  my ($thing) = @_;
  my @f = split /:/, $thing;
  die unless @f == 2;
  my ($chr, $range) = @f;
  @f = split /\-/, $range;
  die unless @f == 2;
  my ($start, $end) = @f;
  return ($chr, $start, $end);
}
