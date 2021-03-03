#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use List::Util qw(sum);
##use lib "/hpcf/apps/pcgp/samplename/lib/perl/";
#use lib "/hpcf/apps/compbio/samplename/lib/perl/";
#use SampleName;

use Cluster;
use FileUtils qw(read_simple_file);
use Reporter;
use DelimitedFile;
#use misc;
use Reporter;
use TdtConfig;

my %FLAGS;
my @JUNCTION_FILES;

my @QC_FILES;

#my $MIN_READS_FOR_NOVEL_JUNCTION = 5;
#my $MIN_READS_FOR_NOVEL_JUNCTION = 0;
# per gang, 2/15/13 ... too noisy
my $MIN_READS_FOR_NOVEL_JUNCTION = 3;
# Gang Wu:
# Based on an evaluation on the novel junction reads in M7 projects,
# we noted that most "novel" junctions (~95%) had read count less than
# 3 across a whole cohort (14 samples), suggesting that these are just
# sporadic or random events. Therefore, we only keep "novel" junctions
# when they are supported by at least 3 supporting reads.

my $FORCE = 1;
my $MIN_CROSS_RAM = 8000;

GetOptions(
	   \%FLAGS,
	   "-jf-shift=s",
	   "-jf-shift-list=s",
	   "-glob-annotated",

	   "-jf-db=s",
	   "-all-db-combined",
	   "-all-db-individual",
	   "-refgene",

	   "-max-shift=i",
	   "-genome=s",

	   "-fasta-dir=s",
	   # legacy
	   "-fasta=s",
	   # preferred

	   "-write",
	   "-write-both",
	   "-annotate",
	   "-bed",

	   "-now",

	   "-summary",

	   "-min-novel-reads=i" => \$MIN_READS_FOR_NOVEL_JUNCTION,

	   "-qc=s" => \@QC_FILES,
	   "-qc-list=s",
	   "-force=i" => \$FORCE,

	   "-cross-sample",
	   "-ram=i",
	   "-strand=s",
	   "-fq",
	   "-no-config",
	  );

use FloatingJunctionDetector;

if ($FLAGS{summary}) {
  summary_report();
  exit(0);
} elsif (@QC_FILES or $FLAGS{"qc-list"}) {
  if (my $qcl = $FLAGS{"qc-list"}) {
    my $list = read_simple_file($qcl);
    @QC_FILES = @{$list};
  }
  my %all_novel;
  foreach (@QC_FILES) {
    qc_check($_, \%all_novel);
  }
  printf STDERR "unique novel: %d\n", scalar keys %all_novel;
  exit(0);
}

my $genome;
unless ($FLAGS{"no-config"}) {
  $genome = $FLAGS{genome} || die "-genome";
}


my $shift_files = get_shift_files();

unless ($FORCE) {
  my @filtered;
  foreach my $jf (@{$shift_files}) {
    my $shifted = get_outfile($jf);
    push @filtered, $jf unless -s $shifted;
  }
  $shift_files = \@filtered;
}

unless (@{$shift_files}) {
  printf STDERR "no files need processing\n";
  exit(0);
}

my $fjd = new FloatingJunctionDetector();

my $fasta = $FLAGS{"fasta"};
my $fasta_dir = $FLAGS{"fasta-dir"};

if ($fasta) {
  # preferred
  $fjd->fasta($fasta);
} elsif ($fasta_dir) {
  # legacy/deprecated
  printf STDERR "WARNING: -fasta-dir DEPRECATED, use -fasta\n";
  $fjd->fasta_dir($fasta_dir);
} else {
  # try config lookup
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "no genome config for $genome";
  if ($fasta = $config_genome->{FASTA}) {
    # preferred
    $fjd->fasta($fasta);
  } elsif ($fasta_dir = $config_genome->{FASTA_CHR}) {
    # legacy/deprecated
    printf STDERR "WARNING: FASTA_CHR deprecated, please migrate to FASTA\n";
    $fjd->fasta_dir($fasta_dir);
  } else {
    die "can't get reference sequence";
  }
}

$fjd->max_shift($FLAGS{"max-shift"}) if $FLAGS{"max-shift"};
$fjd->strand($FLAGS{strand}) if $FLAGS{strand};
database_setup($fjd) if $FLAGS{now};
# only load database if actually running

if ($FLAGS{"cross-sample"}) {
  cross_sample_correction();
  exit(0);
}

foreach my $jf_shift (@{$shift_files}) {

  if ($FLAGS{now}) {
    #
    #  pass 1: combine vs. reference database(s)
    #
    $fjd->find_shifts(
		      "-jf-shift" => $jf_shift,
#		      "-jf-db" => $FLAGS{"jf-db"},
		     );
    # detect shifts vs. reference database
    $fjd->combine_junctions();
    # combine junctions/counts based on results

    #
    #  pass 2: self-compare novel junctions only to see if any can be combined.
    #
    # run this step separately: some (or even many) novel junctions
    # will be merged with reference junctions and disappear.  Only try
    # to correct/combine those novel junctions which still remain
    # after this process.
    print STDERR "pass 2: detect shifts...\n";
    $fjd->find_shifts("-novel-mode" => 1);
    print STDERR "pass 2: combine...\n";
    $fjd->combine_junctions("-novel-mode" => 1);

    #
    #  pass 3: compare novel junctions vs. known exon boundaries.
    #  Useful e.g. for novel skips of known exons.
    #  Some of these will be caught in pass #2 but this requires
    #  at least one observation of the version matching the exon
    #  boundary version, and so will not work e.g. for singletons.
    #  This pass requires a double-anchored match
    #  (i.e. from the end of one known junction to the start of another)
    #
    if (1) {
      print STDERR "pass 3: detect shifts...\n";
      $fjd->find_shifts("-novel-exon-mode" => 1);
      print STDERR "pass 3: combine...\n";
      $fjd->combine_junctions("-novel-mode" => 1);
    }

    #
    #  pass 4: attempt to shift novel junctions vs. single exon starts/ends
    #
    if (0) {
      print STDERR "DEBUG, skipping pass 4\n";
    } else {
      print STDERR "pass 4: detect shifts...\n";
      $fjd->find_shifts("-novel-exon-mode" => 2);
      print STDERR "pass 4: combine...\n";
      $fjd->combine_junctions("-novel-mode" => 1);
    }

    #
    #  done, report:
    #
    $fjd->min_reads_for_novel_junction($MIN_READS_FOR_NOVEL_JUNCTION);
    if ($FLAGS{write}) {
      # write either .tab or .bed
      my $shifted = get_outfile($jf_shift);
      $fjd->write_junctions(
			    "-file" => $shifted,
			    "-annotate" => $FLAGS{annotate},
			    "-bed" => $FLAGS{bed},
			   );
    } elsif ($FLAGS{"write-both"}) {
      # write both .tab and .bed
      $FLAGS{bed} = 0; $FLAGS{annotate} = 1;
      $fjd->write_junctions(
			    "-file" => get_outfile($jf_shift),
			    "-annotate" => $FLAGS{annotate},
			    "-bed" => $FLAGS{bed},
			   );

      $FLAGS{bed} = 1; $FLAGS{annotate} = 0;
      $fjd->write_junctions(
			    "-file" => get_outfile($jf_shift),
			    "-annotate" => $FLAGS{annotate},
			    "-bed" => $FLAGS{bed},
			   );
    } else {
      printf STDERR "not writing output (specify -write or -write-both)\n";
    }
  } else {
    my @args = $0;
    push @args, "-now";
    push @args, "-jf-shift" => $jf_shift;

    if ($FLAGS{"all-db-combined"}) {
      push @args, "-all-db-combined";
    } elsif ($FLAGS{"all-db-individual"}) {
      push @args, "-all-db-individual";
    } elsif (my $db = $FLAGS{"jf-db"}) {
      push @args, "-jf-db" => $db;
    } else {
      die "need jf-db / -all-db-combined / -all-db-individual";
    }
    push @args, "-write" if $FLAGS{write};
    push @args, "-annotate" if $FLAGS{annotate};
    push @args, "-bed" if $FLAGS{bed};
    push @args, "-devel-path" if $FLAGS{"devel-path"};
    push @args, "-force $FORCE";

    my $cmd = join " ", @args;

    my $c = new Cluster();
#    $c->app("perl-5.10.1");
    $c->node_class("");
    $c->memory_reserve_mb(8000);
    $c->memory_limit_mb(8000);
    $c->outfile(get_outfile($jf_shift));
    $c->project("PCGP");
    $c->command($cmd);
    $c->run();
  }
}

sub get_outfile {
  my ($jf_shift) = @_;

  my $bn = $FLAGS{fq} ? $jf_shift : basename($jf_shift);

  if ($FLAGS{"cross-sample"}) {
    return sprintf "%s.cross_sample_corrected.tab",
      $bn;
  } else {
    return sprintf "%s.shifted.%s",
      $bn,
	$FLAGS{bed} ? "bed" : "tab";
  }
}

sub summary_report {
  my @raw = glob("*.junctions.tab");
  die unless @raw;

  my $rpt = new Reporter(
			 "-file" => "summary_tab.txt",
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   file
					   raw_junction_count
					   combined_junction_count
					   removed_count
					   decrease_percent

					   raw_novel_junctions
					   combined_novel_junctions
					   novel_decrease_percent

					   raw_novel_median_reads
					   combined_novel_median_reads

					   final_percent_novel_junctions
					   final_percent_novel_junction_reads
					)
				      ]
			);

  foreach my $fn (@raw) {
    my $raw_info = get_summary_info($fn);
#    my $shifted_info = get_summary_info($fn . ".shifted.tab");
    my $shifted_info = get_summary_info($fn . ".shifted.tab.annotated.tab.cross_sample_corrected.tab");

    my $raw_count = $raw_info->{total};
    my $shifted_count = $shifted_info->{total};

    my $raw_novel_count = $raw_info->{novel};
    my $shifted_novel_count = $shifted_info->{novel};

    my %r;
    $r{file} = $fn;
    $r{raw_junction_count} = $raw_count;
    $r{combined_junction_count} = $shifted_count;
    $r{removed_count} = $raw_count - $shifted_count;
    $r{decrease_percent} = sprintf "%.2f", 100 * (1 - ($shifted_count / $raw_count));

    $r{raw_novel_junctions} = $raw_novel_count;
    $r{combined_novel_junctions} = $shifted_novel_count;
    $r{novel_decrease_percent} = sprintf "%.2f", 100 * (1 - ($shifted_novel_count / $raw_novel_count));

    $r{raw_novel_median_reads} = $raw_info->{novel_median_reads};
    $r{combined_novel_median_reads} = $shifted_info->{novel_median_reads};

    $r{final_percent_novel_junctions} = sprintf "%.2f", $shifted_novel_count * 100 / $shifted_count;
    $r{final_percent_novel_junction_reads} = sprintf "%.2f", ($shifted_info->{reads_novel} * 100) / ($shifted_info->{reads_novel} + $shifted_info->{reads_known});

    $rpt->end_row(\%r);
  }

  $rpt->finish();
}

sub get_summary_info {
  my ($fn) = @_;
  printf STDERR "loading %s...\n", $fn;
  my %info;
  my $df = new DelimitedFile("-file" => $fn,
			     "-headers" => 1,
			     );

  my %counts;

  while (my $row = $df->get_hash()) {
    my $type = $row->{type} || die;
    $info{$type}++;
    $info{total}++;

    my $count = $row->{count} || die;
    my $key = "reads_" . $type;
    $info{$key} += $count;

    push @{$counts{$type}}, $count;
  }

  foreach my $type (qw(novel known)) {
    my $key = $type . "_median_reads";
    my $counts = $counts{$type};
    $info{$key} = median($counts);
  }

  return \%info;
}

sub count_headerless_lines {
  my ($fn) = @_;
  open(C, $fn) || die;
  my $count = 0;
  while (<C>) {
    $count++;
  }
  return $count - 1;
}

sub qc_check {
  my ($fn, $all_novel) = @_;
  open(IN, $fn) || die;
  my $total_junctions = 0;
  my $total_reads = 0;
  my %j;
  my %novel;
  while (<IN>) {
    chomp;
    next if /^junction/;
    my @f = split /\t/, $_;
    die unless @f >= 2;
    my ($j, $count, $type) = @f;

    my @j = split /,/, $j || die;
    die unless @j == 2;
    my ($from, $to) = @j;
    my ($start_chr, $start_base) = FloatingJunctionDetector::parse_junction($j[0]);
    my ($end_chr, $end_base) = FloatingJunctionDetector::parse_junction($j[1]);
    die "urgh" if $start_base == $end_base;

    $total_junctions++;
    $total_reads += $count;
    die "ambiguous junction $j" if $j{$j};
    $j{$j} = 1;
    if ($type eq "novel") {
      $novel{$j} = 1;
      $all_novel->{$j} = 1;
    }
  }
  printf "junctions:%d novel:%d reads:%d file:%s\n", 
    $total_junctions,
      scalar(keys %novel),
	$total_reads,
	  $fn;
}

sub database_setup {
  my ($fjd) = @_;
  if (my $db = $FLAGS{"jf-db"}) {
    # single manually-specified database
    $fjd->add_database("-file" => $db);
  } elsif ($FLAGS{"all-db-individual"}) {
    # load all databases, individually to add preference/priority
    my $sjs = new SJSupport();

    my $priority = 1;
    foreach my $type (qw(
			  refgene
			  ensembl
			  refflat
			  aceview
		       )) {
      my $fn = $sjs->get_refflat_file(
				       "-genome" => $genome,
				       "-type" => $type
				      ) || die;
      my $class = $priority <= 2 ? "core" : "extended";
      $fjd->add_database(
			 "-file" => $fn,
			 "-db-name" => $type,
			 "-db-priority" => $priority++,
			 "-db-class" => $class
			);
    }

    # summarize junction db:
    my $db = $fjd->junctions_db();
    my %by_db;
    foreach my $ref (sort keys %{$db}) {
      foreach my $start (sort {$a <=> $b} keys %{$db->{$ref}}) {
	foreach my $end (sort {$a <=> $b} keys %{$db->{$ref}{$start}}) {
	  my $entry = $db->{$ref}{$start}{$end};
	  $by_db{$entry->{db_name}}++;
	}
      }
    }
    print STDERR "junction summary:\n";
    foreach my $db (sort keys %by_db) {
      printf STDERR "  %s: %6d\n", $db, $by_db{$db};
    }
    printf STDERR "    total: %d\n", sum values %by_db;
  } elsif ($FLAGS{"all-db-combined"}) {
    # single file of all combined databases
    my $config_genome = TdtConfig::readConfig('genome', $genome) || die;
    my $fn = $config_genome->{COMBINED_REFFLAT} || die;
    $fjd->add_database(
		       "-file" => $fn
		      );
  } elsif ($FLAGS{"refgene"}) {
    # refGene only
    my $sjs = new SJSupport();
    my $fn = $sjs->get_refflat_file(
				     "-genome" => $genome,
				     "-type" => "refgene",
				    ) || die;
    $fjd->add_database(
		       "-file" => $fn
		      );

  } else {
    die "specify -jf-db [database]  / -all-db-individual / -all-db-combined / -refgene\n";
  }
}

sub cross_sample_correction {
  #
  #  special final pass to correct novel junctions across samples.
  #  Novel junctions which have not been anchored to a reference exon
  #  (single or double end) may be assigned to arbitrary positions
  #  in individual sample files.  These can only be detected/corrected
  #  by analyzing data from multiple samples.
  #
  my $files = get_shift_files();

  if ($FLAGS{now}) {
    #
    #  find all novel junctions in all input files:
    #
    $fjd->junctions_to_shift({});
    foreach my $jf (@{$files}) {
      my $df = new DelimitedFile(
				 "-file" => $jf,
				 "-headers" => 1,
				);
      while (my $row = $df->get_hash()) {
	my $type = $row->{type};
	if ($type eq "novel") {
	  $fjd->import_junction(
				"-row" => $row,
				"-hash" => $fjd->junctions_to_shift()
			       );

	} elsif ($type ne "known") {
	  die;
	}
      }

      my $hash = $fjd->junctions_to_shift();
      my $count = 0;
      foreach my $chr (keys %{$hash}) {
	foreach my $start (keys %{$hash->{$chr}}) {
	  foreach my $end (keys %{$hash->{$chr}{$start}}) {
	    $count++;
	  }
	}
      }
      printf STDERR "after %s: %d\n", $jf, $count;
    }

    #
    #  run -novel-mode shift, which will identify which junctions
    #  need to be moved, and where:
    #
    $fjd->find_shifts("-novel-mode" => 1);
    $fjd->combine_junctions("-novel-mode" => 1);

    #
    #  generate output, checking whether novel exons have been
    #  shifted:
    #
    my $shifted = $fjd->junctions_to_shift();
    foreach my $jf (@{$files}) {
      printf STDERR "processing $jf...\n";
      # unfortunately, the increased junction population makes it
      # likely additional collapsing possibilities will be detected.
      # So a simple old->new mapping won't always work, and
      # some sites may disappear altogether.

      my $df = new DelimitedFile(
				 "-file" => $jf,
				 "-headers" => 1,
				);
      my %novel_counts;
      #
      # collect all novel counts into tracker hash:
      #
      while (my $row = $df->get_hash()) {
	my $type = $row->{type};
	if ($type eq "novel") {
	  my $j = $row->{junction} || die;
	  die if $novel_counts{$j};
	  $novel_counts{$j} = $row;
	}
      }

      #
      #  check each entry for obsolescence.  If obsolete,
      #  flag and migrate read counts to final destination:
      #
      foreach my $j (keys %novel_counts) {
	my @j = split /,/, $j || die;
	die unless @j == 2;
	my ($from, $to) = @j;
	my ($start_chr, $start_base) = FloatingJunctionDetector::parse_junction($j[0]);
	my ($end_chr, $end_base) = FloatingJunctionDetector::parse_junction($j[1]);
	my $entry = FloatingJunctionDetector::safe_get($shifted, $start_chr, $start_base, $end_base) || die "no entry";
	if ($entry->{obsolete}) {
	  # site has moved elsewhere or been merged
	  my $target = $entry;
	  my ($new_start, $new_end);
	  while ($target->{obsolete}) {
	    # follow trail until final site found
	    ($new_start, $new_end) = @{$target->{moved_to}};
#	    print STDERR " trail: $j => $new_start $new_end\n";
	    $target = FloatingJunctionDetector::safe_get($shifted, $start_chr, $new_start, $new_end) || die "no entry2";
	  }

	  my $j_new = FloatingJunctionDetector::format_junction($start_chr, $new_start, $new_end);
	  printf STDERR "shifting %s => %s\n", $j, $j_new;
	  if (exists $novel_counts{$j_new}) {
	    # entry already exists: combine counts
	    $fjd->pool_info($novel_counts{$j}, $novel_counts{$j_new});
	  } else {
	    # create entry for new site
	    my %r = %{$novel_counts{$j}};
	    $r{junction} = $j_new;
	    # copy old data, renaming junction site
	    $novel_counts{$j_new} = \%r;
	  }
	  delete $novel_counts{$j};
	} else {
#	  printf STDERR "junction %s ok\n", $j;
	}
      }

      #
      # - parse through ordered file, reporting preserved entries only,
      #   delete entries from counts hash
      #
      my $outfile = get_outfile($jf);
      $df = new DelimitedFile(
				 "-file" => $jf,
				 "-headers" => 1,
				);

      my $rpt = $df->get_reporter(
				  "-file" => $outfile
				 );

      while (my $row = $df->get_hash()) {
	my $type = $row->{type};
	if ($type eq "novel") {
	  my $j = $row->{junction} || die;
	  if (my $cache_row = $novel_counts{$j}) {
	    # this site's position hasn't changed.
	    # the count however may have if other junctions shifted here.
	    $rpt->end_row($cache_row);
	    delete $novel_counts{$j};
	  }
	} elsif ($type eq "known") {
	  # known junction: never changes
	  $rpt->end_row($row);
	} else {
	  die;
	}
      }

      #
      # append "leftover" counts (i.e. new coordinates):
      #
      foreach my $j (sort keys %novel_counts) {
	printf STDERR "leftover: %s\n", $j;
	$rpt->end_row($novel_counts{$j});
      }
      $rpt->finish();

      system sprintf "junction2bed.pl -jf %s", $outfile;
      # generate cross-corrected version of .bed
      die "error calling junction2bed.pl" if $?;

    }
  } else {
    my @args = $0;
    push @args, "-now";
    push @args, "-cross-sample";
    push @args, "-genome" => $genome;
    if (my $listfile = $FLAGS{"jf-shift-list"}) {
      push @args, "-jf-shift-list" => $listfile;
    } elsif ($FLAGS{"glob-annotated"}) {
      push @args, "-glob-annotated";
    } else {
      die;
    }
    push @args, "-write" if $FLAGS{write};
    push @args, "-annotate" if $FLAGS{annotate};
    push @args, "-bed" if $FLAGS{bed};
    push @args, "-devel-path" if $FLAGS{"devel-path"};
    push @args, "-force $FORCE";
    if ($FLAGS{"all-db-combined"}) {
      push @args, "-all-db-combined";
    } elsif ($FLAGS{"all-db-individual"}) {
      push @args, "-all-db-individual";
    } elsif (my $db = $FLAGS{"jf-db"}) {
      push @args, "-jf-db" => $db;
    } else {
      die "need jf-db / -all-db-combined / -all-db-individual";
    }

    my $cmd = join " ", @args;

    my $ram = $FLAGS{ram};
    unless ($ram) {
      $ram = int(@{$files} * 16);
    # 12/10/2014: attempt to tailor RAM to large jobs,
    # eventual cluster kill-off attempting to process 931 files in
    # /nfs_exports/genomes/1/projects/RNASEQ/PCGP/BucketIntermediate/JunctionReadCounts/counts1
      $ram = $MIN_CROSS_RAM if $ram < $MIN_CROSS_RAM;
    }

    my $c = new Cluster();
#    $c->app("perl-5.10.1");
    $c->node_class("");
    $c->memory_reserve_mb($ram);
    $c->memory_limit_mb($ram);
    $c->outfile(get_outfile($files->[0]));
    # hack: just for one of the set!
    $c->project("PCGP");
    $c->command($cmd);
    $c->run();
  }
}

sub get_shift_files {
  my @shift_files;
  if (my $jf = $FLAGS{"jf-shift"}) {
    push @shift_files, $jf;
  } elsif (my $list = $FLAGS{"jf-shift-list"}) {
    my $files = read_simple_file($list);
    @shift_files = @{$files};
  } elsif ($FLAGS{"glob-annotated"}) {
    @shift_files = glob("*.annotated.tab");
  } else {
    die "specify -jf-shift [junction file] or -jf-shift-list [listfile]\n";
  }
  die "no files to correct!\n" unless @shift_files;
  return \@shift_files;
}


