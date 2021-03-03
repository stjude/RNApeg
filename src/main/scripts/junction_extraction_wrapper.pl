#!/usr/bin/env perl
#
# wrapper to extract RNA junctions, correct, and annotate
#
# MNE 2/2013
#

use strict;
use warnings;
use Carp qw(confess);

use 5.10.1;
# not required for this script, but required for junction2gene step
# (DBD modules?); add here so outer script halts immediately

use Getopt::Long;
use MiscUtils qw(dump_die build_argv_list);
use File::Basename;

use GeneAnnotation;

use TdtConfig;
use Cluster;
use CommandLineRebuilder;

use constant DEFAULT_GENOME => "GRCh37-lite";
use constant MAINTAINER => "Michael Edmonson <michael.edmonson\@stjude.org";

my %FLAGS;
my $FORCE = 1;
#my $CLUSTER_RAM = 8000;
#my $CLUSTER_RAM = 10000;
# 11/1/2013: Yuxin's jobs killed by cluster with apparent memory usage
# of 9 gb, huh?  Java -Xmx5g in effect.  ???
# Seems more reliable with 10g.
my $CLUSTER_RAM = 12000;
# 5/31/2016: moar

my @TAGS;
my @REGIONS;
my @GENES;

my @options = (
	       "-bam=s",
	       "-bams=s",

	       "-genome=s",

	       "-refflat=s",
	       # alternate combined refFlat file to use,
	       # e.g.
	       # /nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/mRNA/Combined/gencode/all_refFlats_gencode.txt

	       "-min-novel-reads=i",

	       "-now",
	       "-force=i" => \$FORCE,
	       "-primary-only",
	       "-no-duplicates",

	       "-perl-test",
	       "-ram=i" => \$CLUSTER_RAM,

	       "-refgene=s",
	       # refGene file for final gene annotation

	       "-min-span=i",

	       "-tag=s" => \@TAGS,

	       "-region=s" => \@REGIONS,
	       "-gene=s" => \@GENES,


	       # manually-specified reference sequence:
	       "-2bit=s",
	       "-fasta-dir=s",
	       # manually-specified .2bit file
	       "-fasta=s",

	       "-strand=s",
	       # specify + or - for sense/antisense-only BAMs

	       "-ignore-incompatible",
	       "-config-refflat",
	       "-config-combined",

	       "-gencode",

	       "-out-dir=s",
	       "-no-config",
	      );

GetOptions(\%FLAGS, @options);

#perl_sanity_check();
# causing problems w/cloud implementation of Cicero so outlived its usefulness.
# now insist users read the fine manual instead.

reference_param_sanity_check();

if ($FLAGS{"perl-test"}) {
  perl_test();
  exit(0);
}

my $NO_CONFIG_MODE = $FLAGS{"no-config"};
# primarily used when running in Docker

my ($genome, $config_genome, $ger);
unless ($NO_CONFIG_MODE) {
  $genome = $FLAGS{genome} || ($NO_CONFIG_MODE ? undef : DEFAULT_GENOME);
  $config_genome = TdtConfig::readConfig('genome', $genome) || die;
  $ger = $config_genome->{GENE_EXON_REGION_DIR};
}

my $refflat = $FLAGS{refflat};
# refGene database(s) used in correction (override)
my $refgene_ff = $FLAGS{refgene};
# refGene database used for final gene annotation (override)

if ($FLAGS{"config-refflat"}) {
  # use refFlat from config file rather than combined db
  $refflat = $config_genome->{REFSEQ_REFFLAT} || die "no config entry for REFSEQ_REFFLAT";
  $refgene_ff = $refflat;
  # also use this source for gene annotation
} elsif ($FLAGS{"config-combined"}) {
  $refflat = $refgene_ff = $config_genome->{COMBINED_REFFLAT} || die "no COMBINED_REFFLAT config";
}

if ($FLAGS{gencode}) {
  $refflat = $config_genome->{COMBINED_REFFLAT_GENCODE} || die;
}

#
#  init gene annotation source:
#
if ($FLAGS{strand} or not($ger)) {
  $refgene_ff = $config_genome->{REFSEQ_REFFLAT} unless $refgene_ff;
  # backup if (a) user didn't specify and (b) GENE_EXON_REGION not available
  if ($refgene_ff) {
    die "where is $refgene_ff" unless -s $refgene_ff;
  }
  if ($FLAGS{strand}) {
    # filtering by strand requires a gene annotations source with strand info.
    # refgene will work, GENE_EXON_REGION will not.
    printf STDERR "strand specified: taking gene annotations from %s\n", $refgene_ff;
  }

}
die "need GENE_EXON_REGION_DIR/REFSEQ_REFFLAT in genome configuration, or -refgene [refgene flatfile]" unless $ger or $refgene_ff;

init_regions() if @GENES;
# init targeted regions, if specified

my $files = build_argv_list("-flags" => \%FLAGS,
			    "-single" => "bam",
			    "-set" => "bams");

#
#  verify that outfiles won't collide:
#
my %all_out;
foreach my $bam (@{$files}) {
  my $out_base = get_outfile_base($bam);
  if (my $last = $all_out{$out_base}) {
    die sprintf "outfile collision between %s and %s: see -tag option in documentation", $bam, $last;
  }
  $all_out{$out_base} = $bam;
}

if ($FLAGS{now}) {
  #
  #  run steps as needed on this node
  #

  my $strand = $FLAGS{strand};
  if ($strand) {
    die "-strand must be + or -" unless $strand eq "+" or $strand eq "-";
  }

  foreach my $bam (@{$files}) {
    my $out_junctions = get_outfile_base($bam) . ".junctions.tab";

    my $cmd;
    #
    #  extract raw junction counts:
    #
#    $cmd = sprintf 'bam_junction.pl -type all -refflat-type all -bam %s -annotate -now -genome %s -force %d -out %s', $bam, $genome, $FORCE, $out_junctions;
#    $cmd = sprintf 'bam_junction.pl -type all -bam %s -annotate -now -genome %s -force %d -out %s', $bam, $genome, $FORCE, $out_junctions;
    $cmd = sprintf 'bam_junction.pl -type all -bam %s -annotate -now -force %d -out %s', $bam, $FORCE, $out_junctions;
    if ($genome) {
      $cmd .= sprintf " -genome %s", $genome;
    } elsif ($NO_CONFIG_MODE) {
      $cmd .= " -no-config";
    } else {
      die;
    }

    if ($refflat) {
      $cmd .= " -refflat " . $refflat;
    } else {
      $cmd .= " -refflat-type all";
    }

    $cmd .= " -primary-only" if $FLAGS{"primary-only"};
    $cmd .= " -no-duplicates" if $FLAGS{"no-duplicates"};
    $cmd .= sprintf " -min-span %d", $FLAGS{"min-span"} if $FLAGS{"min-span"};
    foreach my $region (@REGIONS) {
      $cmd .= sprintf ' -region %s', $region;
    }
    $cmd .= sprintf ' -fasta %s', $FLAGS{"fasta"} if $FLAGS{"fasta"};
    $cmd .= sprintf " -2bit %s", $FLAGS{"2bit"} if $FLAGS{"2bit"};
    $cmd .= sprintf " -strand %s", $strand if $strand;
    $cmd .= " -ignore-incompatible" if $FLAGS{"ignore-incompatible"};

    run_cmd($cmd, $out_junctions);

    #
    # floating junction correction:
    #
    my $out_shift = $out_junctions . ".shifted.tab";
#    $cmd = sprintf 'floating_junction_fix.pl -jf-shift %s -all-db-combined -now -genome %s -force %d', $out_junctions, $genome, $FORCE;
#    $cmd = sprintf 'floating_junction_fix.pl -jf-shift %s -now -genome %s -force %d', $out_junctions, $genome, $FORCE;
    $cmd = sprintf 'floating_junction_fix.pl -jf-shift %s -now -force %d', $out_junctions, $FORCE;

    if ($genome) {
      $cmd .= sprintf " -genome %s", $genome;
    } elsif ($NO_CONFIG_MODE) {
      $cmd .= " -no-config";
    } else {
      die;
    }

    if ($refflat) {
      $cmd .= " -jf-db " . $refflat;
    } else {
      $cmd .= " -all-db-combined";
    }
    $cmd .= sprintf ' -min-novel-reads %d', $FLAGS{"min-novel-reads"} if $FLAGS{"min-novel-reads"};
    $cmd .= sprintf ' -fasta-dir %s', $FLAGS{"fasta-dir"} if $FLAGS{"fasta-dir"};
    $cmd .= sprintf ' -fasta %s', $FLAGS{"fasta"} if $FLAGS{"fasta"};

    if (1) {
      # write both .tab and .bed
      $cmd .= " -write-both";
    } else {
      $cmd .= " -write -annotate";
    }
    $cmd .= sprintf " -strand %s", $strand if $strand;
    $cmd .= " -fq" if $FLAGS{"out-dir"};

    run_cmd($cmd, $out_shift);

    #
    # final gene annotation:
    #
    my $out_annotated = $out_shift . ".annotated.tab";
      #    $cmd = sprintf 'junction2gene -junctions %s -preserve -genome %s -force %d', $out_shift, $genome, $FORCE;
    $cmd = sprintf 'junction2gene.pl -junctions %s -preserve -force %d', $out_shift, $FORCE;

    if ($genome) {
      $cmd .= sprintf " -genome %s", $genome;
    } elsif ($NO_CONFIG_MODE) {
      $cmd .= " -genome GRCh37-lite";
      # genome is just used for species chromosome count,
      # config system is not involved.  Assume human in this case.
    } else {
      die;
    }

    if ($refgene_ff) {
      $cmd .= sprintf " -refgene %s", $refgene_ff if $refgene_ff;
    } else {
      $cmd .= sprintf " -gene-exon-region-dir %s", $ger if $ger;
    }
    $cmd .= sprintf " -strand %s", $strand if $strand;

    run_cmd($cmd, $out_annotated);

  }

} else {
  #
  # submit a job for each file to process each step
  #
  foreach my $bam (@{$files}) {
    my $outfile = get_outfile_base($bam) . ".junctions.tab.shifted.tab.annotated.tab";

    my $cmd;
    if (1) {
      my $clrb = new CommandLineRebuilder(
					  "-flags" => \%FLAGS,
					  "-parameters" => \@options
					 );
      $clrb->exclude_parameter("-bam");
      $clrb->exclude_parameter("-bams");
      $clrb->include_parameter("-now");
      $cmd = $clrb->get_command_line("-bam" => $bam);
    } else {
      $cmd = sprintf '%s -now -force %d -genome %s -bam %s', $0, $FORCE, $genome, $bam;
      $cmd .= sprintf ' -min-novel-reads %d', $FLAGS{"min-novel-reads"} if $FLAGS{"min-novel-reads"};
      $cmd .= " -primary-only" if $FLAGS{"primary-only"};
      $cmd .= " -no-duplicates" if $FLAGS{"no-duplicates"};
      $cmd .= sprintf " -refgene %s", $refgene_ff if $refgene_ff;
      $cmd .= sprintf ' -ram %d', $CLUSTER_RAM;
      $cmd .= sprintf ' -min-span %d', $FLAGS{"min-span"} if $FLAGS{"min-span"};
      foreach (@TAGS) {
	$cmd .= sprintf ' -tag %s', $_;
      }
      foreach my $region (@REGIONS) {
	$cmd .= sprintf ' -region %s', $region;
      }
      $cmd .= sprintf ' -out-dir %s', $FLAGS{"out-dir"} if $FLAGS{"out-dir"};
    }

    my $c = new Cluster();
    $c->force($FORCE);
    # by default, re-run even if output files already exist

#    $c->app("perl-5.10.1");
    # fail on clinical

    $c->node_class("");
    $c->memory_reserve_mb($CLUSTER_RAM);
    $c->memory_limit_mb($CLUSTER_RAM);
    $c->outfile($outfile);
    $c->project("PCGP");
    $c->command($cmd);
    $c->run();
  }
}

sub run_cmd {
  my ($cmd, $outfile) = @_;
  printf STDERR "%s: running: %s\n", scalar(localtime), $cmd;
  system $cmd;
  confess "error calling $cmd, exit $?" if $?;
  die "where is $outfile" unless -s $outfile;
}

sub perl_test {
  my $perl_binary = $^X;
  if ($perl_binary eq "/usr/bin/perl") {
    printf STDERR "ERROR: perl binary is %s! use \"module load\" instead\n", $perl_binary;
  }


  my @cmds;
  push @cmds, "/bin/env perl -v";
  push @cmds, '/bin/env perl -e \'print "perl binary=$^X\n"\'';

  foreach my $script (qw(
			  bam_junction.pl
			  floating_junction_fix.pl
			  junction2gene.pl
		       )) {
    push @cmds, sprintf '/bin/env perl -cw `which %s`', $script;
  }

  foreach my $cmd (@cmds) {
    print "$cmd\n";
    system $cmd;
  }
}

sub perl_sanity_check {
  my $perl_binary = $^X;
  if ($perl_binary eq "/usr/bin/perl") {
    die sprintf "ERROR: perl binary is \"%s\".  Use \"module load perl/5.10.1\" instead!\n", $perl_binary;
  }
}

sub get_outfile_base {
  my ($bam) = @_;
  my $outfile;
  my $out_dir = $FLAGS{"out-dir"};
  if (@TAGS) {
    die "tags not implemented w/out-dir yet" if $out_dir;
    my ($hit) = grep {$bam =~ /$_/} @TAGS;
    die "$bam doesn't match any tags" unless $hit;
    $outfile = sprintf '%s.%s', $hit, basename($bam);
  } elsif ($out_dir) {
    die "$out_dir not a directory" unless -d $out_dir;
    $outfile = sprintf '%s/%s', $out_dir, basename($bam);
  } else {
    $outfile = basename($bam);
  }
  return $outfile;
}

sub init_regions {
  my $ga;
  if ($ger) {
    $ga = new GeneAnnotation(
			     "-style" => "gene_exon_region",
			     "-gene_exon_region_dir" => $ger,
			     "-ignore_non_coding" => 0,
			     "-genome" => $genome
			    );
  } elsif ($refgene_ff) {
    die "refgene flatfile mode not yet implemented, contact " . MAINTAINER;
  }

  printf STDERR "finding gene intervals, this may take a little while...\n" if @GENES;
  foreach my $gene (@GENES) {
    my $hits = $ga->find_gene($gene) || die "no records found for $gene";
    # might be some gene symbol incompatibility
    foreach my $hit (@{$hits}) {
      my $region = sprintf '%s:%d-%d', @{$hit}{qw(chrom txStart txEnd)};
      printf STDERR "region=%s\n", $region;
      push @REGIONS, $region;
      # NOTE: it's quite possible these may overlap!
      # collapse these later to prevent duplicate output.

      # FIX ME: include any flanking sequence?
    }
  }
}

sub reference_param_sanity_check {
  my @p = qw(2bit fasta-dir);
  my $count;
  foreach my $p (@p) {
    $count++ if $FLAGS{$p};
  }

  if ($count and $count != 2) {
    die sprintf "ERROR: if any of %s are specified, both must be specified\n", join "/", @p;
  }
}

