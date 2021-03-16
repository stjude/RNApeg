#!/usr/bin/env perl
# simple reporting of strongest evidence of aberrant splicing by gene
#
# TO DO:
# - option to use count of reads *passing flanking QC* rather than the
#   raw count?
# - require minimum + and - read counts?

use strict;
use warnings;
use Carp qw(confess);

use Getopt::Long;
use File::Basename;

use DelimitedFile;
use Reporter;
use Counter;
use JunctionUtils qw(parse_junction);
use AtomicOutfile;
use UniqueRank;
use FileUtils qw(read_simple_file write_simple_file);
use List::Util qw(min);
use MiscUtils qw(dump_die average median);
use Cluster;
use CommandLineRebuilder;
use RefFlatFile;
use ChrBucketMap;

# use DelimitedFile;
# use Reporter;

##use lib "/hpcf/apps/pcgp/samplename/lib/perl/";
#use lib "/hpcf/apps/compbio/samplename/lib/perl/";
#use SampleName;
# my %info = SampleName::parse(basename($bam));
# foreach (sort keys %info) {
#   printf STDERR "%s: %s\n", $_, $info{$_};
# }
# die;

my $REFGENE_ONLY = 1;
my $MATRIX_TOP_COUNT = 50;
my $SHARE_RANK_IF_RATIO_EQUAL = 1;
my $CLUSTER_RAM = 10000;

#use constant EXCLUSIVE_ABERRANT_RATIO => 1000;
use constant EXCLUSIVE_ABERRANT_RATIO => 100;
# if a novel aberrant junction is found but the known junction
# it touches is not, special reporting.
# Not sure if these cases are very interesting (i.e. novel junction
# completely dominates/suppresses the known junction)
# or should be ignored.

my $MAX_CLOSE_ABERRANT_EDGE_DISTANCE = 100;
# delete me! replace w/absolute distance

use constant SHORT_197 => "chr21:39764366:+,chr21:39764564:+";

my @command_line_params = (
			   "-jf=s",
			   "-jf-list=s",
			   "-glob-cross",
			   "-glob-dir=s",
			   "-glob=s",

			   "-require-paired",
			   # if we find a novel junction but not its known mate,
			   # don't report it with the artificially high ratio of 1000.
			   # - pro: returns highest true comparison available
			   # - con: hides these events, which may be very interesting
			   #        and in fact the most significant event

			   "-matrix",
			   # digest results

			   "-min-novel-reads=i",
			   "-qc-min-flanking=i",
			   "-qc-min-good=i",
			   "-qc-min-perfect=i",
			   "-qc-bidirectional=i",

			   "-verbose",

			   "-short197-exclusive",
			   "-classify-short197",
			   "-generate-bam-list",

			   #
			   # subset processing parameters:
			   #
			   "-short197",
			   "-no-short197",
			   "-disease=s",
			   "-erg-del=i",

			   "-track-max-ratio=i",
			   "-min-known-reads=i",
			   "-qc-check-known",

			   "-cluster",
			   "-no-cluster",

			   "-cap-ratio=f",

			   "-use-flanking-counts",
			   # when computing the ratio of aberrant junction
			   # reads to reference reads, use the qc_flanking
			   # counts rather than the raw counts, as they
			   # are much less likely to contain marginal
			   # "fingernail" mapping reads

			   "-rank2example=s",
			   # extract summary/example evidence from
			   # a rank file
			   "-max=i",
			   # count in top X to extract

			   "-all",
			   # disable SJ cohort filtering

			   "-hack",


			   # START border check flags:
			   "-test-shared-border=i",
			   # test: how often do we see the same known edge
			   # interacted with (a) within a sample?
			   # (b) across samples?
			   "-shared-border-min-aberrant-ratio=f",
			   "-prune-low-counts",

			   "-gl-gold-cancer-ranges=s",
			   # same as used by medal_ceremony.pl
			   "-require-reviewable",
			   "-shared-min-counts=s",
			   # END border check flags


			   "-max-edge-distance=i",

			   "-gene=s",

			   "-max-one-exonic-edge",
			   "-refflat=s",
			   "-require-downstream-known",
			   "-min-aberrant-known-edge-distance=i",

			   "-digest-genes=s",
			   "-digest-shared=s",

			  );

my %FLAGS;
GetOptions(
	   \%FLAGS,
	   @command_line_params
	  );

my $clr = new CommandLineRebuilder("-parameters" => \@command_line_params,
				   "-flags" => \%FLAGS);
$clr->exclude_parameter("-glob-cross");
$clr->exclude_parameter("-glob");
$clr->exclude_parameter("-cluster");

my $MIN_NOVEL_COUNT = get_flag("min-novel-reads", 10);
my $MIN_KNOWN_READ_COUNT = get_flag("min-known-reads", 0);
# these alone aren't very useful as they can contain poor-quality alignments

my $QC_MIN_FLANKING = get_flag("qc-min-flanking", 0);
my $QC_MIN_GOOD_ALIGNMENTS = get_flag("qc-min-good", 0);
my $QC_MIN_PERFECT_ALIGNMENTS = get_flag("qc-min-perfect", 0);
my $QC_REQUIRE_BIDIRECTIONAL = get_flag("qc-bidirectional", 0);
my $TRACK_MAX_RATIO = get_flag("track-max-ratio", 1);
my $QC_CHECK_KNOWN = get_flag("qc-check-known", 0);
my $CAP_RATIO = get_flag("cap-ratio", 0);

printf STDERR "aberrant ratio: tracking maximum?: %s\n",
  $TRACK_MAX_RATIO ? "yes" : "no";

#use PrimerDB;
#my $pdb = new PrimerDB;
#my $dbi = $pdb->dbd;

#use DBTools qw(get_dbi_hg_readonly);
#my $dbi = get_dbi_hg_readonly();

if ($FLAGS{matrix}) {
  generate_matrix();
  exit(0);
} elsif ($FLAGS{"short197-exclusive"}) {
  report_short197_exclusive();
  exit(0);
} elsif ($FLAGS{"classify-short197"}) {
  classify_short197();
  exit(0);
} elsif ($FLAGS{"generate-bam-list"}) {
  generate_bam_list();
  exit(0);
} elsif ($FLAGS{"rank2example"}) {
  generate_examples_from_ranks();
  exit(0);
} elsif ($FLAGS{hack}) {
  die is_aberrant_edge_near_known("chr21:39764366:+,chr21:39764564:+", "chr21:39764366:+,chr21:39774479:+");
} elsif ($FLAGS{"test-shared-border"}) {
  test_shared_borders();
  exit(0);
} elsif ($FLAGS{"summarize-aberrant"}) {
  # create summary reports of aberrant junctions across samples in a cohort
  summarize_aberrant();
  exit(0);
} elsif ($FLAGS{"known-downgrade"}) {
  known_downgrade();
  exit(0);
} elsif ($FLAGS{"digest-genes"}) {
  # summarize junctions across samples for a set of genes
  digest_genes();
  exit(0);
} elsif ($FLAGS{"digest-shared"}) {
  digest_shared();
  exit(0);
}

my @junction_files;
if (my $jf = $FLAGS{jf}) {
  push @junction_files, $jf;
} elsif (my $l = $FLAGS{"jf-list"}) {
  my $list = read_simple_file($l);
  push @junction_files, @{$list};
} elsif ($FLAGS{"glob-cross"}) {
  my $glob = sprintf "%s/*.cross_sample_corrected.tab",
    $FLAGS{"glob-dir"} || ".";
  push @junction_files, glob($glob);
} elsif ($FLAGS{"glob"}) {
  my $glob = sprintf "%s/%s", ($FLAGS{"glob-dir"} || "."), $FLAGS{glob};
  push @junction_files, glob($glob);
  die unless @junction_files;
} else {
  die "specify -jf [junction file] | -jf-list [listfile] | -glob-cross\n";
}

if ($FLAGS{cluster} and not($FLAGS{"no-cluster"})) {
  die "no files" unless @junction_files;
  foreach my $jf (@junction_files) {
    my $outfile = $jf . ".fun.tab.aberrant.tab";
    my $cmd = $clr->get_command_line("-jf" => $jf);;

    my $c = new Cluster();
    $c->node_class("");
    $c->memory_reserve_mb($CLUSTER_RAM);
    $c->memory_limit_mb($CLUSTER_RAM);
    $c->outfile($outfile);
#    $c->queue("pcgp_heavy_io");
    $c->project("aberrant_splicing");
    $c->command($cmd);
    $c->run();

  }
  exit(0);
}

my $cbm_exons;
if ($FLAGS{"max-one-exonic-edge"}) {
  # bucket refFlats by exon for border checks
  my $f_refflat = $FLAGS{refflat} || die "-refflat";
  my $rff = new RefFlatFile();
  $rff->missing_genes_ok(1);
  # combined file contains data from multiple formats
  $rff->canonical_references_only(1);
  # ?
  printf STDERR "parsing %s...\n", $f_refflat;
  $rff->parse_file(
		   "-type" => "refflat",
		   "-refflat" => $f_refflat,
		  );
  $cbm_exons = new ChrBucketMap();
  $cbm_exons->f_chr("chrom");
  $cbm_exons->f_start("exon_start");
  $cbm_exons->f_end("exon_end");

  printf STDERR "building exon db...\n";
  foreach my $row (@{$rff->rows}) {
    foreach my $exon (@{$row->{exons}}) {
      my %r;
      $r{name} = $row->{name};
      $r{strand} = $row->{strand};
      $r{chrom} = $row->{chrom} || die;
      $r{exon_start} = $exon->{start} || die;
      $r{exon_end} = $exon->{end} || die;
      $cbm_exons->add_row(
			  "-row" => \%r,
			 );
    }
  }
}

my $c = new Counter(\@junction_files);

foreach my $jf (@junction_files) {
  my $fun = $jf . ".fun.tab";
  die "where is $fun" unless -s $fun;

  printf STDERR "processing %s...\n", $fun;

  my $df = new DelimitedFile("-file" => $fun,
			     "-headers" => 1,
			     );

  my @interesting;

  #
  #  find aberrant junctions of interest for this sample:
  #
  my %touched;

  while (my $row = $df->get_hash()) {
    my $event = $row->{event} || die;
    next unless $row->{count} >= $MIN_NOVEL_COUNT;
    my $usable;
    if ($event eq "matches_known_exon_start" or
	$event eq "matches_known_exon_end") {
      $usable = 1;
    } elsif ($event eq "novel_skip_of_known_exons") {
      # skip
    } else {
      die $event;
    }

    $usable = 0 unless junction_qc($row);

    if ($usable) {
      push @interesting, $row;
      foreach my $j (parse_junction($row->{junction})) {
	push @{$touched{$j->{ref}}{$j->{base}}}, $row;
      }
    }
  }

  #
  #  find reference junctions sharing a boundary with our target junctions:
  #
  $df = new DelimitedFile(
			  "-file" => $jf,
			  "-headers" => 1,
			 );
  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {
    next unless $row->{type} eq "known";

    if ($REFGENE_ONLY) {
      my @rg = get_refgene_transcripts($row);
      next unless @rg;
#      printf STDERR "refgenes: %s\n", join ",", @rg;
    }

    if ($QC_CHECK_KNOWN and !junction_qc($row)) {
      # require known junction evidence pass minimum QC
#      dump_die($row, "hey now, high coverage and no flanking") if $row->{count} >= 15 and not $row->{qc_flanking};
      next;
    }

#    dump_die($row);

#    foreach (sort keys %{$row}) {
#      printf "%s: %s\n", $_, $row->{$_};
#    }

    my @hits;
    foreach my $j (parse_junction($row->{junction})) {
      if ($touched{$j->{ref}}{$j->{base}}) {
	foreach my $novel_junction (@{$touched{$j->{ref}}{$j->{base}}}) {
	  push @{$novel_junction->{known_hits}}, $row;
	}
      }
    }
  }

  #
  #  each novel aberrant junction should now be associated with
  #  the known junction it touches.
  #

#  my @exclusive_aberrant;

  my $outfile_all = basename($fun) . ".aberrant_all.tab";
  my $rpt_all = new Reporter(
			 "-file" => $outfile_all,
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   gene
					   aberrant_to_known_ratio
					   aberrant_found_but_known_not
					   aberrant_junction
					   aberrant_junction_count
					   known_junction
					   known_junction_count
					   known_junction_transcripts
					   abs_aberrant_edge_distance
					)
				      ]
			 );



  my %max_aberrant;

  my $use_flanking_counts = $FLAGS{"use-flanking-counts"};

  foreach my $row (@interesting) {
    my $known_hits = $row->{known_hits};
    if ($known_hits) {
      if (0) {
	printf STDERR "%s:\n", $row->{junction};
	foreach my $k (@{$known_hits}) {
	  printf STDERR "  %s: %s\n", $k->{junction}, $k->{transcripts};
	}
      }

      my $novel_count = $row->{count};
      my %ratios;

      my @pending_all;

      my $aberrant_junction = $row->{junction};

      foreach my $known (@{$known_hits}) {
	my $known_count = $known->{count};
	next if $known_count < $MIN_KNOWN_READ_COUNT;

	my $usable = 1;

	my ($known_edge, $aberrant_edge) = find_shared_known_and_aberrant_edge($known->{junction}, $aberrant_junction);

	if ($FLAGS{"max-one-exonic-edge"}) {
	  #
	  # require the aberrant edge NOT touch an exon, i.e. possible
	  # novel exon.
	  #

	  if (my $exon_hits = $cbm_exons->find(
				       "-chr" => $aberrant_edge->{ref},
				       "-start" => $aberrant_edge->{base},
				       "-end" => $aberrant_edge->{base}
					      )) {
	    dump_die($aberrant_edge, "aberrant exonic $aberrant_junction " . $exon_hits->[0]->{name}, 1) if $FLAGS{verbose};
	    $usable = 0;
	  } else {
#	    dump_die($aberrant_edge, "aberrant not exonic $aberrant_junction", 1);
	  }
	}

	if ($FLAGS{"require-downstream-known"}) {
	  #
	  # require downstream edge to be known
	  #
	  my $known_pos = $known_edge->{base};
	  my $aberrant_pos = $aberrant_edge->{base};

	  my $exon_hits = $cbm_exons->find(
					   "-chr" => $aberrant_edge->{ref},
					   "-start" => $known_pos,
					   "-end" => $known_pos
					   ) || die "can't find known edge";
	  my %strands = map {$_->{strand}, 1} @{$exon_hits};

	  my $pass;
	  if ($strands{"+"}) {
	    if ($known_pos > $aberrant_pos) {
	      $pass = 1;
#	      die sprintf "+ pass known:%s aberrant:%s edge_known:%d edge_aberrant:%d\n", $known->{junction}, $aberrant_junction, $known_pos, $aberrant_pos;
	    } else {
	      printf STDERR "+ fail known:%s aberrant:%s edge_known:%d edge_aberrant:%d\n", $known->{junction}, $aberrant_junction, $known_pos, $aberrant_pos if $FLAGS{verbose};
	      # OK (spec)
	    }
	  }
	  if ($strands{"-"}) {
	    if ($known_pos < $aberrant_pos) {
#	      die sprintf "- pass known:%s aberrant:%s edge_known:%d edge_aberrant:%d\n", $known->{junction}, $aberrant_junction, $known_pos, $aberrant_pos;
	      $pass = 1;
	    } else {
	      printf STDERR "- fail known:%s aberrant:%s edge_known:%d edge_aberrant:%d\n", $known->{junction}, $aberrant_junction, $known_pos, $aberrant_pos if $FLAGS{verbose};
	    }
	  }

	  $usable = 0 unless $pass;
	}

	if (my $min_known_dist = $FLAGS{"min-aberrant-known-edge-distance"}) {

	  my $edge_distance = find_known_and_aberrant_edge(
							   $known->{junction},
							   $aberrant_junction,
							   "-distance" => 1
							  );
	  die if $edge_distance <= 0;
	  if ($edge_distance < $min_known_dist) {
	    printf STDERR "edge distance fail known:%s aberrant:%s\n", $known->{junction}, $aberrant_junction;
	    $usable = 0;
	  }
	}

	next unless $usable;

	my $novel_to_known_ratio;
	if ($use_flanking_counts) {
	  $novel_to_known_ratio = $row->{qc_flanking} / $known->{qc_flanking};
	} else {
	  $novel_to_known_ratio = $novel_count / $known_count;
	}

	$novel_to_known_ratio = $CAP_RATIO if $CAP_RATIO and $novel_to_known_ratio > $CAP_RATIO;

	my %saw;
	foreach my $gene (split /,/, $row->{genes}) {
	  next if $saw{$gene};
	  $saw{$gene} = 1;

	  if ($TRACK_MAX_RATIO) {
	    # track the highest ratio of aberrant gene count to known
	    # gene count.
	    # Does this really make sense though if there are multiple
	    # possible matching known junctions?
	    my $max = $max_aberrant{$gene}{ratio} || 0;
	    if ($novel_to_known_ratio > $max) {
	      $max_aberrant{$gene}{ratio} = $novel_to_known_ratio;
	      $max_aberrant{$gene}{aberrant} = $row;
	      $max_aberrant{$gene}{known} = $known;
	    }
	  } else {
	    push @{$ratios{$gene}{$aberrant_junction}}, [ $novel_to_known_ratio, $known ];
	  }

	  my %r;
	  $r{gene} = $gene;
	  $r{aberrant_junction} = $aberrant_junction;
	  $r{aberrant_junction_count} = $row->{count};
	  $r{aberrant_found_but_known_not} = 0;
	  $r{known_junction} = $known->{junction};
	  $r{known_junction_count} = $known->{count};
	  $r{known_junction_transcripts} = $known->{transcripts};

	  if ($TRACK_MAX_RATIO) {
	    $r{aberrant_to_known_ratio} = $novel_to_known_ratio;
	    # NAIVE: shouldn't be used because there may be > 1
	    # possible known junction to match to!
	    $rpt_all->end_row(\%r);
	  } else {
	    # wait until we can calculate the most-conservative ratio
	    push @pending_all, \%r;
	  }
	}
      }

      if (!$TRACK_MAX_RATIO) {
	my %max_aj;

	foreach my $gene (keys %ratios) {
	  foreach my $aj (keys %{$ratios{$gene}}) {
	    my @sorted = sort {$a->[0] <=> $b->[0]} @{$ratios{$gene}{$aj}};

	    my ($wanted, $known) = @{$sorted[0]};
	    # take the most conservative ratio of aberrant junction to
	    # reference junction (reasonable in case there is more than
	    # one possible matching known junction).

	    # GLOBAL tracking by gene:
	    my $max = $max_aberrant{$gene}{ratio} || 0;
	    if ($wanted > $max) {
	      $max_aberrant{$gene}{ratio} = $wanted;
	      $max_aberrant{$gene}{aberrant} = $row;
	      $max_aberrant{$gene}{known} = $known;
	    }

	    # LOCAL tracking by aberrant junction:
	    $max = $max_aj{$aj} || 0;
	    $max_aj{$aj} = $wanted if $wanted > $max;
	  }
	}

	foreach my $r (@pending_all) {
	  # report the most-conservative ratio vs. known junctions,
	  # i.e. against the known junction having the highest coverage.
	  $r->{aberrant_to_known_ratio} = $max_aj{$r->{aberrant_junction}};
	  $r->{abs_aberrant_edge_distance} = find_known_and_aberrant_edge(
									  $r->{known_junction},
									  $r->{aberrant_junction},
									  "-distance" => 1
									 );

	  $rpt_all->end_row($r);
	}

      }
    } elsif ($FLAGS{"require-paired"}) {
      printf STDERR "skipping observation of aberrant only junction %s for %s\n", $row->{junction}, $row->{genes};
      # can happen if:
      # - no evidence for known junction
      # - known junction doesn't meet required coverage/QC
    } else {
      # no observation of the known junction this novel junction touches
      # (or poor evidence that's been disqualified)
      my %saw;
      foreach my $gene (split /,/, $row->{genes}) {
	next if $saw{$gene};
	$saw{$gene} = 1;

	my $level = $CAP_RATIO || EXCLUSIVE_ABERRANT_RATIO;

	$max_aberrant{$gene}{ratio} = $level;
	$max_aberrant{$gene}{aberrant} = $row;
	$max_aberrant{$gene}{known} = undef;

	my %r;
	$r{gene} = $gene;
	$r{aberrant_to_known_ratio} = $level;
	$r{aberrant_junction} = $row->{junction};
	$r{aberrant_junction_count} = $row->{count};
	$r{aberrant_found_but_known_not} = 1;
	$r{known_junction} = "";
	$r{known_junction_count} = -1;
	$r{known_junction_transcripts} = "";
	$r{abs_aberrant_edge_distance} = "";
	$rpt_all->end_row(\%r);
      }
    }
  }
  $rpt_all->finish();

  #
  #  report results:
  #
  my $outfile = basename($fun) . ".aberrant.tab";
  my $rpt = new Reporter(
			 "-file" => $outfile,
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   gene
					   aberrant_to_known_ratio
					   aberrant_found_but_known_not
					   aberrant_junction
					   aberrant_junction_count
					   known_junction
					   known_junction_count
					   known_junction_transcripts
					   abs_aberrant_edge_distance
					)
				      ]
			);

  foreach my $gene (sort keys %max_aberrant) {
    my %r;
    $r{gene} = $gene;
    $r{aberrant_to_known_ratio} = $max_aberrant{$gene}{ratio};
    my $aberrant = $max_aberrant{$gene}{aberrant} || die;
    @r{qw(aberrant_junction aberrant_junction_count)} =
      @{$aberrant}{qw(junction count)};
    my $known = $max_aberrant{$gene}{known};
    if ($known) {
      $r{aberrant_found_but_known_not} = 0;
      @r{qw(known_junction known_junction_count known_junction_transcripts)} =
	@{$known}{qw(junction count transcripts)};
      if ($REFGENE_ONLY) {
	my @rg = get_refgene_transcripts($known);
	die unless @rg;
	$r{known_junction_transcripts} = join ",", @rg;
      }
#      $r{aberrant_edge_near_known} = is_aberrant_edge_near_known($r{known_junction}, $r{aberrant_junction});

      $r{abs_aberrant_edge_distance} = find_known_and_aberrant_edge(
								    $r{known_junction},
								    $r{aberrant_junction},
								    "-distance" => 1
								   );
    } else {
      $r{aberrant_found_but_known_not} = 1;
      foreach my $f (
		     qw(
			 known_junction
			 known_junction_count
			 known_junction_transcripts
			 abs_aberrant_edge_distance
		      )
		    ) {
	$r{$f} = "";
      }
    }

    $rpt->end_row(\%r);
  }

  $rpt->finish();

  $c->next($jf);
}


sub get_refgene_transcripts {
  my ($row) = @_;

  my @t = split /,/, $row->{transcripts};
  die unless @t;
  my @rg;
  foreach (@t) {
    push @rg, $_ if /^N[MR]_/;
  }
  return @rg;
}

sub generate_matrix {
  my $glob = sprintf '%s/*.aberrant.tab', $FLAGS{"glob-dir"} || ".";
  my @raw = glob($glob);
  die unless @raw;

  my $ab_files;
  if ($FLAGS{all}) {
    $ab_files = \@raw;
  } else {
    $ab_files = filter_files(\@raw);
  }

  my $ao = new AtomicOutfile("-base" => "set",
			     "-suffix" => "tab");
  $ao->add(get_report_name() || die);
  write_simple_file($ab_files, $ao->get_outfile());

  #
  #  compute a rank for each gene within each sample:
  #
  my %ratios;
  my %all_genes;
  my $c = new Counter($ab_files);
  my $max_edge_distance = $FLAGS{"max-edge-distance"};
  foreach my $af (@{$ab_files}) {
    my $top = $ratios{$af} = {};
    my $df = new DelimitedFile("-file" => $af,
			       "-headers" => 1,
			      );
    while (my $row = $df->get_hash()) {
      my $gene = $row->{gene} || die;
      my $ratio = $row->{aberrant_to_known_ratio} || 0;
      $top->{$gene}{ratio} = $ratio;
      $all_genes{$gene} = 1;
    }

    # assign a rank based on aberrant-to-known frequency:
    my $rank = 0;
    if ($SHARE_RANK_IF_RATIO_EQUAL) {
      my $last_ratio;
      foreach my $gene (sort { $top->{$b}{ratio} <=> $top->{$a}{ratio} }
			keys %{$top}) {
	# sort by ratio (descending), then gene name
	my $this_ratio = $top->{$gene}{ratio};
	if (not(defined $last_ratio) or $this_ratio != $last_ratio) {
	  $rank++;
	} else {
#	  printf STDERR "hey now: same rank $rank ratio=$this_ratio last_ratio=$last_ratio\n";
	}
	$top->{$gene}{rank} = $rank;
#	printf STDERR "gene=%s rank=%d ratio=%f\n", $gene, $rank, $this_ratio;
	$last_ratio = $this_ratio;
      }

    } else {
      foreach my $gene (sort { $top->{$b}{ratio} <=> $top->{$a}{ratio} ||
				 $a cmp $b }
			keys %{$top}) {
	# sort by ratio (descending), then gene name
	$top->{$gene}{rank} = ++$rank;
      }
    }
    $c->next($af);
  }

  #
  #  for each gene, calculate average and median rank across samples:
  #
  my (%avg_rank, %median_rank);

  my $null_rank = (scalar keys %all_genes);
  # if event not observed for a gene/sample, assign worst rank
#  my $null_rank = (scalar keys %all_genes) * 2;
  # if we assign a higher penalty for not being observed,
  # will more frequently-mutated genes rise to top?

  foreach my $gene (sort keys %all_genes) {
    # my @ranks = grep {defined $_} map {$ratios{$_}{$gene}{rank}} @ab_files;
    # FAIL: falsely high reporting if only a few examples
    my @ranks;
    foreach my $af (@{$ab_files}) {
      my $r = $ratios{$af}{$gene}{rank};
      push @ranks, defined $r ? $r : $null_rank;
    }
    my $avg_rank = average(\@ranks);
    my $median_rank = median(\@ranks);
    $avg_rank{$gene} = $avg_rank;
    $median_rank{$gene} = $median_rank;
    printf STDERR "gene=%s avg_rank=%s median_rank=%d rank_count=%d raw_ranks=%s\n", $gene, $avg_rank, $median_rank, scalar(@ranks), join ",", sort {$a <=> $b} @ranks if $FLAGS{verbose};
  }


  foreach my $is_average (0, 1) {
    my $field = $is_average ? "mean_rank" : "median_rank";
    my $hash = $is_average ? \%avg_rank : \%median_rank;

    my $ao = new AtomicOutfile("-base" => "ranks",
			       "-suffix" => "tab");
    $ao->add(get_report_name() || die);
    $ao->add($is_average ? "mean" : "median");
    my $outfile = $ao->get_outfile();

    my $rpt = new Reporter(
			   "-file" => $outfile,
			   "-delimiter" => "\t",
			   "-labels" => [
					 "meta_rank",
					 "order",
					 "gene",
					 $field,
					]
			  );

    my @genes_by_rank = sort {$hash->{$a} <=> $hash->{$b}} keys %{$hash};
    my $ur = new UniqueRank();

    foreach my $gene (@genes_by_rank) {
      my %r;
      my $value = $hash->{$gene};
      $r{gene} = $gene;
      $r{meta_rank} = $ur->get_rank($value);
      $r{order} = $ur->get_order();
      $r{$field} = $value;
      $rpt->end_row(\%r);
    }
    $rpt->finish();

    my @matrix_genes = @genes_by_rank[0 .. ($MATRIX_TOP_COUNT - 1)];

    my @labels = (
		  "gene",
		  "meta_rank",
		  "order",
		  $field
		  );
    push @labels, map {get_label($_)} @{$ab_files};

    $ao = new AtomicOutfile("-base" => "matrix",
			    "-suffix" => "tab");
    $ao->add(get_report_name() || die);
    $ao->add($is_average ? "mean" : "median");
    $rpt = new Reporter(
			"-file" => $ao->get_outfile(),
			"-delimiter" => "\t",
			"-labels" => \@labels
		       );

    $ur = new UniqueRank();
    foreach my $gene (@matrix_genes) {
      my %r;
      my $value = $hash->{$gene};
      $r{gene} = $gene;
      $r{meta_rank} = $ur->get_rank($value);
      $r{order} = $ur->get_order();
      $r{$field} = $value;
      foreach my $af (@{$ab_files}) {
	$r{get_label($af)} = $ratios{$af}{$gene}{rank} || "-1";
      }
      $rpt->end_row(\%r);
    }
    $rpt->finish();
  }

}

sub get_label {
  my ($label) = @_;
  my $new = $label;
  $new = basename($label);
  $new =~ s/\.bam.*/\.bam/;
#  printf STDERR "%s => %s\n", $label, $new;
  return $new;
}

sub report_short197_exclusive {
  my $ab_short = get_short197("-short197" => 1, "-all" => 1);
  my $ab_not = get_short197("-no-short197" => 1, "-all" => 1);

  my ($short, $df_short) = load_junctions($ab_short, "short197");
  my ($not, $df_not) = load_junctions($ab_not, "non-short197");

  my $rpt = new Reporter(
			 "-file" => "exclusive_short.tab",
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   gene
					   aberrant_junction
					   known_junction
					   known_junction_transcripts
					)
				      ]
			);

  my $short_exclusive = 0;
  my $both = 0;
  foreach my $j (keys %{$short}) {
    if ($not->{$j}) {
      # also in non-short, skip
      $both++;
    } else {
      # unique to short
      my %r = %{$short->{$j}};
      $rpt->end_row(\%r);
      $short_exclusive++;
    }
  }
  $rpt->finish();

  printf STDERR "short:%d not:%d overlap:%d short_exclusive:%d\n",
    scalar(keys %{$short}),
      scalar(keys %{$not}),
	$both,
	  $short_exclusive;

}

sub get_short197 {
  # duplicate code, big hurry, so sorry
  my (%options) = @_;

  my $short197;
  my $no_short197;
  if ($options{"-short197"}) {
    $short197 = 1;
  } elsif ($options{"-no-short197"}) {
    $no_short197 = 1;
  } else {
    die "huh? " . join ",", keys %options;
  }

  my $pattern = $options{"-all"} ? "aberrant_all.tab" : "aberrant.tab";

  my $glob = sprintf '%s/*.%s', ($FLAGS{"glob-dir"} || "."), $pattern;
  my @ab_files = glob($glob);
  die unless @ab_files;

  my %wanted;
  open(SHORT, "short197_status.txt") || die "can't open short197_status.txt: $!";
  while (<SHORT>) {
    chomp;
    s/\r$//;
    my ($sample, @stuff) = split /\s+/;
    my %status;
    foreach my $thing (@stuff) {
      if ($thing =~ /\w/) {
	#	  printf STDERR "set key '%s'\n", $thing;
	$status{$thing} = 1;
      }
    }
    #      printf STDERR "line=$_ keys=%d\n", scalar keys %status;
    delete $status{unknown};
    if (%status) {
      die "ambiguous " . join ",", keys %status if scalar keys %status > 1;
      my ($type) = keys %status;
      die unless $type eq "short197" or $type eq "noshort197";
      if (
	  ($short197 and $type eq "short197") or
	  ($no_short197 and $type eq "noshort197")
	 ) {
	$wanted{$sample} = 1;
      }
    }
  }

  #
  #  filter report list to those matching desired samples:
  #
  my @filtered;
  foreach my $f (@ab_files) {
    basename($f) =~ /^(SJ[A-Z]+\d+)/ || die;
    my $sample = $1;
    push @filtered, $f if $wanted{$sample};
  }
  @ab_files = @filtered;
  printf STDERR "filtered to %d files, %s\n", scalar @ab_files, join ",", @ab_files;

  return \@ab_files;
}

sub load_junctions {
  my ($list, $label) = @_;
  printf STDERR "loading %s novel junctions...\n", $label;
  my $c = new Counter($list);
  my %j;
  my $df;
  foreach my $af (@{$list}) {
    $df = new DelimitedFile("-file" => $af,
			       "-headers" => 1,
			      );
    while (my $row = $df->get_hash()) {
      $j{$row->{aberrant_junction} || die} = $row;
    }
    $c->next($af);
  }
  return (\%j, $df);
}

sub classify_short197 {
  my @corr = glob("*cross_sample_corrected.tab");
  die unless @corr;
  my %call;
  my $c = new Counter(\@corr);
  foreach my $corr (@corr) {
    basename($corr) =~ /^(SJ[A-Z]+\d+)/ || die "fail for $corr";
    my $subject = $1;

    my $df = new DelimitedFile("-file" => $corr,
			       "-headers" => 1,
			      );

    my $is_short197 = "noshort197";

    while (my $row = $df->get_hash()) {
      my $j = $row->{junction};
      if ($j eq SHORT_197) {
	$is_short197 = "short197";
	last;
      }
    }
    die "duplicate" if $call{$subject};
    printf STDERR "%s: %s\n", basename($corr), $is_short197;
    $call{$subject} = $is_short197;
    $c->next($corr);
  }

  my $wf = new WorkingFile("short_status_junction.txt");
  my $fh = $wf->output_filehandle();
  foreach my $s (sort keys %call) {
    printf $fh "%s\n", join "\t", $s, $call{$s};
  }
  $wf->finish();
}

sub generate_bam_list {
  my $rows = get_annotation_rows();
  my @bams;
  foreach my $row (@{$rows}) {
    my $sample = $row->{sample} || die;
    $sample =~ /^(SJ[A-Z]+)/ || die;
    my $glob = sprintf '/nfs_exports/genomes/1/projects/RNASEQ/PCGP/BucketRaw/%s/%s*bam', $1, $sample;
    my @files = glob($glob);
    die "can't find unique hit for $glob" unless @files == 1;
    push @bams, $files[0];
  }
  write_simple_file(\@bams, "bams.txt");
}

sub get_annotation_rows {
  my $df = new DelimitedFile(
			     "-file" => "SJERG_isoform_proportion_v8.txt",
			     "-headers" => 1,
			     "-skip_until" => "sample",
			     );
  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  my @rows;
  while (my $row = $df->get_hash()) {
    push @rows, $row;
  }
  return \@rows;
}

sub get_report_name {
  my @things;
  if (my $disease = $FLAGS{disease}) {
    $disease = uc($disease);
    push @things, $disease;

    my $del_mode = $FLAGS{"erg-del"};
    if (defined $del_mode) {
      die "only ERG" unless $disease eq "ERG";
      push @things, sprintf "deletion-%s", $del_mode ? "positive" : "negative";
    }

  } else {
    push @things, "cohort";
    # entire list
  }

  if ($FLAGS{short197}) {
    push @things, "short197";
  } elsif ($FLAGS{"no-short197"}) {
    push @things, "no_short197";
  }

  return join "_", @things;
}

sub old_short197 {
  my $short197 = $FLAGS{"short197"};
  my $no_short197 = $FLAGS{"no-short197"};
  my %wanted;
  my @ab_files;
  if ($short197 or $no_short197) {
    open(SHORT, "short197_status.txt") || die $!;
    while (<SHORT>) {
      chomp;
      s/\r$//;
      my ($sample, @stuff) = split /\s+/;
      my %status;
      foreach my $thing (@stuff) {
	if ($thing =~ /\w/) {
#	  printf STDERR "set key '%s'\n", $thing;
	  $status{$thing} = 1;
	}
      }
#      printf STDERR "line=$_ keys=%d\n", scalar keys %status;
      delete $status{unknown};
      if (%status) {
	die "ambiguous " . join ",", keys %status if scalar keys %status > 1;
	my ($type) = keys %status;
	die unless $type eq "short197" or $type eq "noshort197";
	if (
	    ($short197 and $type eq "short197") or
	    ($no_short197 and $type eq "noshort197")
	   ) {
	  $wanted{$sample} = 1;
	}
      }
    }

    #
    #  filter report list to those matching desired samples:
    #
    my @filtered;
    foreach my $f (@ab_files) {
      basename($f) =~ /^(SJ[A-Z]+\d+)/ || die;
      my $sample = $1;
      push @filtered, $f if $wanted{$sample};
    }
    @ab_files = @filtered;
    printf STDERR "filtered to %d files, %s\n", scalar @ab_files, join ",", @ab_files;
  }
}

sub filter_files {
  my ($raw_list) = @_;

  my $short197 = $FLAGS{"short197"};
  my $no_short197 = $FLAGS{"no-short197"};
  my $disease = $FLAGS{disease};
  $disease = uc($disease) if $disease;
  my $del_mode = $FLAGS{"erg-del"};
  if (defined $del_mode) {
    die "only ERG" unless $disease and $disease eq "ERG";
  }

  my $ar = get_annotation_rows();
  my %wanted;
  my %all;
  foreach my $row (@{$ar}) {
    if (0) {
      foreach (sort keys %{$row}) {
	printf "%s: %s\n", $_, $row->{$_};
      }
    }

    my $usable = 1;
    my $short197_ratio = $row->{e6_alt_ab};

    $usable = 0 if $short197 and !$short197_ratio;
    $usable = 0 if $no_short197 and $short197_ratio;

    if ($disease) {
      my $d = $row->{disease2} || die "no disease2";
      # duplicate columns, disease2 is 2nd instance (appended by parser)
      $usable = 0 unless $d eq $disease;

      if ($usable and defined $del_mode) {
	my $status = $row->{ERG_Del} || die;
	die "unknown status $status for " . $row->{sample} unless $status eq "Yes" or $status eq "No";
	$usable = 0 if $del_mode == 0 and $status ne "No";
	$usable = 0 if $del_mode == 1 and $status ne "Yes";
      }
    }

    $all{$row->{sample}} = 1;
    if ($usable) {
      my $sample = $row->{sample} || die;
      $wanted{$sample} = 1;
    }
  }

  my @filtered;
  my %saw;
  foreach my $f (@{$raw_list}) {
    basename($f) =~ /^(SJ[^\-]+)/ || die "no SJ for $f";
    my $sample = $1;
    die unless $all{$sample};
    # sanity check to ensure we can link our reports to .xls annotations

    die "duplicate $sample" if $saw{$sample};
    $saw{$sample} = 1;
    push @filtered, $f if $wanted{$sample};
  }
  confess "WTF: no usable rows!" unless @filtered;

  return \@filtered;
}

sub junction_qc {
  my ($row) = @_;

  my @fail;

  unless ($row->{qc_flanking} >= $QC_MIN_FLANKING) {
    push @fail, "qc_flanking";
  }

  if ($QC_REQUIRE_BIDIRECTIONAL) {
    unless ($row->{qc_plus} and $row->{qc_minus}) {
      push @fail, "qc_bidirectional";
    }
  }
  if ($QC_MIN_GOOD_ALIGNMENTS) {
    my $good = $row->{qc_perfect_reads} + $row->{qc_clean_reads};
    unless ($good >= $QC_MIN_GOOD_ALIGNMENTS) {
      push @fail, "qc_min_good";
    }
  }

  if ($QC_MIN_PERFECT_ALIGNMENTS) {
    unless ($row->{qc_perfect_reads} >= $QC_MIN_PERFECT_ALIGNMENTS) {
      push @fail, "qc_perfect";
    }
  }

  if ($FLAGS{verbose} and @fail) {
    printf STDERR "QC fail %s: %s\n", $row->{junction}, join ",", @fail;
  }

  return @fail ? 0 : 1;
}

sub get_flag {
  my ($name, $default) = @_;
  my $value;
  if (exists $FLAGS{$name}) {
    $value = $FLAGS{$name};
  } else {
    $value = $default;
  }
  return $value;
}

sub generate_examples_from_ranks {

  my @genes;
  my %gene2rank;
  # optional: won't know if gene is manually specified (future)

  #
  #  build list of genes to process:
  #
  my $single_gene = $FLAGS{gene};
  if (my $rank_file = $FLAGS{"rank2example"}) {
    my $max = $FLAGS{max} || 20;
    my $df = new DelimitedFile("-file" => $rank_file,
			       "-headers" => 1,
			      );
    my $count;
    while (my $row = $df->get_hash()) {
      my $gene = $row->{gene} || die;
      my $rank = $row->{meta_rank} || die;
      if ($single_gene) {
	next unless $gene eq $single_gene;
      } else {
	last if $rank > $max;
      }
      $gene2rank{$gene} = $rank;
      push @genes, $gene;
    }
  } else {
    die "TO DO: manual gene list, etc.";
  }

  my @aberrant_files = glob("*aberrant.tab");
  my %wanted_genes = map {$_, 1} @genes;

  my %track;
  printf STDERR "scanning aberrant files...";
  my %all_j;
  foreach my $af (@aberrant_files) {
    #
    # scan aberrant files for genes of interest
    #
    print STDERR ".";
    my $df = new DelimitedFile("-file" => $af,
			       "-headers" => 1,
			      );
    while (my $row = $df->get_hash()) {
      next if $row->{aberrant_found_but_known_not};
      # hack
      my $gene = $row->{gene} || dump_die($row, $af);
      next unless $wanted_genes{$gene};
      push @{$track{$gene}{$af}}, $row;

      foreach my $field (qw(known_junction aberrant_junction)) {
	$all_j{$row->{$field} || die} = 1;
      }
    }
  }
  print STDERR "\n";

  #
  #  cache raw junction rows for all junctions of interest.
  #  Save repeated (expensive) parsing of all files for each gene.
  #  This itself is very slow but total time is faster when generating
  #  larger rank sets.
  #
  my %junction_info;
  printf STDERR "start: %s\n", scalar localtime;
  printf STDERR "scanning junction files (this will take a while)...\n";
  my $df_j;
  my $c = new Counter(\@aberrant_files);
  foreach my $af (@aberrant_files) {
    my $infile = $af;
    $infile =~ s/\.fun.tab.*// || die;
    die "where is $infile" unless -s $infile;

    $df_j = new DelimitedFile(
			      "-file" => $infile,
			      "-headers" => 1,
			     );
    while (my $row = $df_j->get_hash()) {
      my $j = $row->{junction};
      if ($all_j{$j}) {
	die "duplicate" if $junction_info{$af}{$j};
	$junction_info{$af}{$j} = $row;
      }
    }
    $c->next($af);
  }
  printf STDERR "end: %s\n", scalar localtime;

  #
  #  summarize:
  #
  my $rpt = new Reporter(
			 "-file" => "aberrant_summary.tab",
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   gene
					   rank
					   aberrant_samples
					   unique_aberrant_junctions
					   aberrant_junctions
					   aberrant_junction_counts
					   aberrant_junction_distance
					)
				      ],
			 "-auto_qc" => 1,
			);

  my ($rpt_all_examples, $rpt_all_examples_aberrant, $rpt_all_examples_aberrant_top, $rpt_all_examples_aberrant_top_near);

  printf STDERR "start: %s\n", scalar localtime;

  my $total_samples = scalar @aberrant_files;

  foreach my $gene (@genes) {
    my $by_sample = $track{$gene};
    my $aberrant_samples = scalar keys %{$by_sample};
    my @rows;
    foreach my $sample (keys %{$by_sample}) {
      my $set = $by_sample->{$sample};
      die if @{$set} > 1;
      # unique by gene in these reports
      my $row = $set->[0];
      $row->{sample} = $sample;
      push @rows, $row;
    }

    my %aberrant;
    my %aberrant2row;
    foreach my $r (@rows) {
      $aberrant{$r->{aberrant_junction}}++;
      $aberrant2row{$r->{aberrant_junction}} = $r;
      # for singleton fields only
    }

    my @by_count = sort {$aberrant{$b} <=> $aberrant{$a}} keys %aberrant;

    my $rank = $gene2rank{$gene} || "";
    my $junction_number = 0;
    foreach my $aj (@by_count) {
#      die $aj;
      my @set = grep {$_->{aberrant_junction} eq $aj} @rows;

      my @sorted = sort {$b->{aberrant_to_known_ratio} <=> $a->{aberrant_to_known_ratio}} @set;

      my $best = $sorted[0];
      # example with strongest signal

      #
      #  create excerpt file for aberrant and reference junction:
      #
      my $ao = new AtomicOutfile("-base" => "aj",
				 "-suffix" => "tab");
      $ao->add(sprintf 'rank%03d', $rank) if $rank;
      $ao->add($gene);

      $ao->add(sprintf 'junction%d', ++$junction_number);
      $ao->add(sprintf 'count%d', scalar @set);

      my $sample = $best->{sample};
      $ao->add($1) if $sample =~ /(SJ.*bam)/;
      my $outfile = $ao->get_outfile();

      my %wanted_j = map {$_, 1} @{$best}{qw(known_junction aberrant_junction)};
      die unless scalar keys %wanted_j == 2;

      my $rpt_j = $df_j->get_reporter("-file" => $outfile);
      unless ($rpt_all_examples) {
	my @extra = "sample";
	push @extra, "rank_gene" if $rank;
	push @extra, (
		      "rank_junction",
		      # junction # within gene
		     );
	push @extra, "top_ranked_samples";
	# note this is NOT the number of samples containing the junction,
	# but rather the number of samples where this junction was the
	# most-aberrant junction
	push @extra, "aberrant_edge_near_known";

	$rpt_all_examples = $df_j->get_reporter(
						"-file" => "aberrant_all.tab",
						"-extra" => \@extra
					       );
	$rpt_all_examples_aberrant = $df_j->get_reporter(
						"-file" => "aberrant_only.tab",
						"-extra" => \@extra
					       );
	$rpt_all_examples_aberrant_top = $df_j->get_reporter(
						"-file" => "aberrant_only_top.tab",
						"-extra" => \@extra
					       );
	$rpt_all_examples_aberrant_top_near = $df_j->get_reporter(
						"-file" => "aberrant_only_top_near.tab",
						"-extra" => \@extra
					       );
      }

      foreach my $j (keys %wanted_j) {
	my $r = $junction_info{$sample}{$j} || die "no info for $sample $j";
	$rpt_j->end_row($r);

	$r->{sample} = $sample;
	$r->{rank_gene} = $rank if $rank;
	$r->{rank_junction} = $junction_number;
	$r->{top_ranked_samples} = scalar @set;

	my $is_aberrant_edge_near_known = is_aberrant_edge_near_known(keys %wanted_j);
	$r->{aberrant_edge_near_known} = $is_aberrant_edge_near_known;

	$rpt_all_examples->end_row($r);
	if ($r->{type} eq "novel") {
	  $rpt_all_examples_aberrant->end_row($r);
	  if ($junction_number == 1) {
	    $rpt_all_examples_aberrant_top->end_row($r);
	    $rpt_all_examples_aberrant_top_near->end_row($r) if $is_aberrant_edge_near_known;
	  }
	}
      }
      $rpt_j->finish();

      junction2bed($outfile);
    }

    my %r;
    $r{gene} = $gene;
    $r{rank} = $rank;
    $r{aberrant_samples} = $aberrant_samples;
    $r{unique_aberrant_junctions} = scalar keys %aberrant;
    $r{aberrant_junctions} = join ";", @by_count;
    $r{aberrant_junction_counts} = join ",", map {$aberrant{$_}} @by_count;
    $r{aberrant_junction_distance} = join ",", map {$aberrant2row{$_}{abs_aberrant_edge_distance}} @by_count;

    $rpt->end_row(\%r);
  }
  printf STDERR "end: %s\n", scalar localtime;

  foreach my $r ($rpt,
		 $rpt_all_examples,
		 $rpt_all_examples_aberrant,
		 $rpt_all_examples_aberrant_top,
		 $rpt_all_examples_aberrant_top_near
		) {
    $r->finish();
  }
}

sub is_aberrant_edge_near_known {
  my ($j1, $j2) = @_;
  # one known and one related aberrant junction (order irrelevant)

  my %edge2count;
  foreach my $j_raw ($j1, $j2) {
    foreach my $edge (parse_junction($j_raw)) {
      $edge2count{$edge->{base}}++;
    }
  }

  my @compare = grep {$edge2count{$_} == 1} keys %edge2count;
  # aberrant site comparison: one aberrant edge, one known edge
  # known edge will have count 2 (skipped)
  die unless @compare == 2;

  my $distance = abs($compare[0] - $compare[1]);

  return $distance <= $MAX_CLOSE_ABERRANT_EDGE_DISTANCE ? 1 : 0;
}

sub test_shared_borders {
  #
  #  TO DO: lower-level .bed files with read counts?
  #
  my $type = $FLAGS{"test-shared-border"};
  my @af_raw;
  printf STDERR "*** run type: ";
  if ($type == 1) {
    # top-ranked junction per gene only
    printf STDERR "top-ranked aberrant junctions only\n";
    @af_raw = glob("*aberrant.tab");
  } elsif ($type == 2) {
    printf STDERR "ALL aberrant junctions\n";
    my $pattern = $FLAGS{glob} || '*.aberrant_all.tab';
    @af_raw = glob(sprintf "%s/%s", ($FLAGS{"glob-dir"} || "."), $pattern);
  } else {
    die;
  }

  my $gl_ranges = $FLAGS{"gl-gold-cancer-ranges"} || die;
  open(RNG, $gl_ranges) || die;
  my %gl_reviewable_genes;
  while (<RNG>) {
    chomp;
    my ($gene, $range) = split /\t/, $_;
    $gl_reviewable_genes{$gene} = 1;
  }
  # FIX ME:
  # THIS CAN EASILY BREAK IF SYMBOL CHANGES BETWEEN THIS ANNOTATION
  # SOURCE AND THE ONE USED TO ANNOTATE THE JUNCTION CALLS!

  my $smc = $FLAGS{"shared-min-counts"};
  my $min_passed_junctions = 1;
  my $min_samples_to_pass = 1;

  if ($smc) {
    my @f = split /,/, $smc;
    die unless @f == 2;
    ($min_passed_junctions, $min_samples_to_pass) = @f;
  }

  my $tag;
  my $ab_files;
  if ($FLAGS{all}) {
    $tag = "all";
    $ab_files = \@af_raw;
  } else {
    $tag = get_report_name();
    $ab_files = filter_files(\@af_raw);
  }

  my $set_file = sprintf 'run_set_%s.txt', $tag;
  my $outfile = sprintf 'aberrant_shared_known_edges_%s.tab', $tag;

  printf STDERR "final input files: %d\n", scalar @{$ab_files};
  write_simple_file($ab_files, $set_file);

  my %shared;

  my $min_aberrant_ratio = $FLAGS{"shared-border-min-aberrant-ratio"} || 0;
  my $max_edge_distance = $FLAGS{"max-edge-distance"};

  foreach my $af (@{$ab_files}) {
    printf STDERR "parsing %s...\n", $af;
    my $df = new DelimitedFile(
			      "-file" => $af,
			      "-headers" => 1,
			     );
    while (my $row = $df->get_hash()) {
      my $j_known = $row->{known_junction} || die;
      my $j_aberrant = $row->{aberrant_junction} || die;

      my $ar = $row->{aberrant_to_known_ratio};
      die unless defined $ar;
      unless ($ar >= $min_aberrant_ratio) {
#	printf STDERR "toss %s, ratio=%f\n", $j_aberrant, $ar;
	next;
      }

      if ($max_edge_distance) {
	my $aberrant_distance = $row->{abs_aberrant_edge_distance} || die;
	next unless $aberrant_distance <= $max_edge_distance;
      }

      if (0) {
	print STDERR "DEBUG!!!\n";
	next unless $row->{gene} eq "KDM1A";
      }

#      my ($known_edge, $aberrant_edge) = find_known_and_aberrant_edge($j_known, $j_aberrant);
      my ($known_edge, $aberrant_edge) = find_shared_known_and_aberrant_edge($j_known, $j_aberrant);

      $row->{sample} = $af;

#      dump_die($row, sprintf("known=%d aberrant=%d", $known_edge->{base}, $aberrant_edge->{base}), 1);
#      printf STDERR "known:%s aberrant:%s known_edge:%d aberrant_edge:%d\n", $j_known, $j_aberrant, $known_edge->{base}, $aberrant_edge->{base};

      push @{$shared{$known_edge->{ref}}{$known_edge->{base}}{$aberrant_edge->{base}}}, $row;

    }
  }

  my $rpt_shared = new Reporter(
			 "-file" => $outfile,
			 "-delimiter" => "\t",
			 "-labels" => [
				       "Chr",
				       "known_edge",
				       "aberrant_edge_total",
				       "genes",
				       "gl_reviewable_gene",
				       "sample_count",
				       "aberrant_edge_counts",
				       "aberrant_edges",
				       "aberrant_samples",
				      ],
			 "-auto_qc" => 1,
			);

  my $ao = new AtomicOutfile("-base" => "sj",
			     "-suffix" => "tab",
			     "-auto_increment" => 1,
			    );

  foreach my $ref (sort keys %shared) {
    foreach my $known_edge (sort {$a <=> $b} keys %{$shared{$ref}}) {
      my @aberrant_edges = keys %{$shared{$ref}{$known_edge}};

      my %genes;
      my %edge2sample;
      my %samples;
      my %strands;

      foreach my $ae (@aberrant_edges) {
	my $rows = $shared{$ref}{$known_edge}{$ae};
	foreach my $r (@{$rows}) {
	  my $gene = $r->{gene};
	  die "comma in gene $gene" if $gene =~ /,/;
	  $genes{$r->{gene}} = 1;
	  $edge2sample{$ae}{$r->{sample}} = 1;
	  $samples{$r->{sample}} = 1;
	}
      }

      if ($FLAGS{"prune-low-counts"}) {
	# require all reported junctions to meet the minimum sample count.
	# by default all junctions are reported as long as at least one
	# of them meets the minimum.  This option declutters the output.
	my %ae = map {$_, 1} @aberrant_edges;
	foreach my $ae (keys %edge2sample) {
	  my $count = scalar keys %{$edge2sample{$ae}};
	  if ($count < $min_samples_to_pass) {
	    delete $ae{$ae};
	  }
	}

	# rebuild index info as above:
	%genes = ();
	%edge2sample = ();
	%samples = ();
	@aberrant_edges = keys %ae;
	foreach my $ae (@aberrant_edges) {
	  # rebuild index using passed edges only
	  my $rows = $shared{$ref}{$known_edge}{$ae};
	  foreach my $r (@{$rows}) {
	    my $gene = $r->{gene};
	    die "comma in gene $gene" if $gene =~ /,/;
	    $genes{$r->{gene}} = 1;
	    $edge2sample{$ae}{$r->{sample}} = 1;
	    $samples{$r->{sample}} = 1;
	  }
	}
      }

      my $gl_reviewable = 0;
      foreach (keys %genes) {
	$gl_reviewable = 1 if $gl_reviewable_genes{$_};
      }

#      if (@aberrant_edges > 1) {
      if (@aberrant_edges) {
	my @e2s_sorted = sort {scalar keys %{$edge2sample{$b}} <=> scalar keys %{$edge2sample{$a}}} keys %edge2sample;

	my $thing = join ";", map {$_ . "=" . scalar keys %{$edge2sample{$_}}} @e2s_sorted;

	my @aberrant_edge_counts = map {scalar keys %{$edge2sample{$_}}} @e2s_sorted;

	my @aberrant_samples;
	foreach my $e (@e2s_sorted) {
	  my @samples;
	  foreach my $s (sort keys %{$edge2sample{$e}}) {
	    my $sample = $s;
	    $sample = basename($sample);
	    $sample =~ s/\.bam.*//;
	    push @samples, $sample;
	  }
	  push @aberrant_samples, join ",", @samples;
	}

	my $passed = 0;
	foreach my $count (@aberrant_edge_counts) {
	  $passed++ if $count >= $min_samples_to_pass;
	}
	$passed = 0 unless $passed >= $min_passed_junctions;
	next unless $passed;

	my %r;
	$r{Chr} = $ref;
	$r{known_edge} = $known_edge;
	$r{aberrant_edge_total} = scalar @aberrant_edges;
	$r{genes} = join ",", sort keys %genes;
	$r{sample_count} = scalar keys %samples;
	$r{aberrant_edge_counts} = join ",", @aberrant_edge_counts;
	$r{aberrant_edges} = join ",", @e2s_sorted;
	$r{aberrant_samples} = join ";", @aberrant_samples;
	$r{gl_reviewable_gene} = $gl_reviewable;

	next if $FLAGS{"require-reviewable"} and !$gl_reviewable;

	$rpt_shared->end_row(\%r);

	#
	#  build example files:
	#
	$ao->reset();
	$ao->add($tag);
	$ao->add(sort keys %genes);
	my $j_outfile = $ao->get_outfile();

	my $rpt_j = new Reporter(
				 "-file" => $j_outfile,
				 "-delimiter" => "\t",
				 "-labels" => [
					       qw(
						   junction
						   count
						   type
						   genes
						   transcripts
						)
					      ]
				);

	my %known;
	my @aj;

	my %all_known_samples;

	foreach my $ae (@aberrant_edges) {
	  my $rows = $shared{$ref}{$known_edge}{$ae};
#	  my $example = $rows->[0];
#	  $known{$example->{known_junction}} = $example;
	  my $example = $rows->[0];
	  foreach my $r (@{$rows}) {
	    $known{$r->{known_junction}} = $r;
	  }

	  my @s = keys %{$edge2sample{$ae}};
	  foreach my $s (@s) {
	    $all_known_samples{$s} = 1;
	  }

	  my %r;
	  $r{junction} = $example->{aberrant_junction};
	  $r{type} = "novel";
	  $r{count} = scalar @s;
	  $r{genes} = $example->{gene};
	  $r{transcripts} = "";
	  push @aj, \%r;
	}

	if (scalar keys %known > 1) {
	  printf STDERR "NOTE: multiple known junctions in %s! %s\n", $j_outfile, join "  ", keys %known;
	  # - multiple novel junctions
	  # - each novel junction pairs with a different known junction
	  # - however both of the known junctions share a common edge,
	  #   which is why reuslts are bucketed together.
	  # - PRO: capture of shared reference site even if samples
	  #        use different isoforms
	  # - CON: different isoforms
	}

	foreach my $j (keys %known) {
	  # known junctions first
	  my %r;
	  $r{junction} = $j;
	  $r{type} = "known";
	  $r{count} = scalar keys %all_known_samples;
	  my $kj = $known{$j};
	  $r{genes} = $kj->{gene};
	  $r{transcripts} = $kj->{known_junction_transcripts};

	  $rpt_j->end_row(\%r);
	}

	# then aberrant:
	foreach my $r (sort {$b->{count} <=> $a->{count}} @aj) {
	  $rpt_j->end_row($r);
	}
	$rpt_j->finish();

	junction2bed($j_outfile);
	# convert to .bed

      }
    }
  }

  $rpt_shared->finish();

}


sub find_known_and_aberrant_edge {
  # given a known and aberrant junction, return the known and aberrant
  # edges (ignoring the known edge shared by both junctions)
  my ($j_known, $j_aberrant, %options) = @_;

  my %edges;
  foreach my $j_raw ($j_known, $j_aberrant) {
    foreach my $edge (parse_junction($j_raw)) {
      $edges{$edge->{base}}{$j_raw} = $edge;
    }
  }

#  printf STDERR "%s\n", join "\n", $j_known, $j_aberrant, keys %edges;

  my ($edge_known, $edge_aberrant);

  foreach my $base (keys %edges) {
    my @entries = keys %{$edges{$base}};
    if (@entries == 1) {
      # count of 2 = the known edge shared in both the known
      # and aberrant junction, ignore
      my $j = $entries[0];
      my $edge = $edges{$base}{$j};

      if ($j eq $j_known) {
	$edge_known = $edge;
      } elsif ($j eq $j_aberrant) {
	$edge_aberrant = $edge;
      } else {
	die;
      }
    }
  }

  if ($options{"-distance"}) {
    return(abs($edge_known->{base} - $edge_aberrant->{base}));
  } else {
    return ($edge_known, $edge_aberrant);
  }
}

sub find_shared_known_and_aberrant_edge {
  # given a known and aberrant junction, return the SHARED known edge
  # and the aberrant edge
  my ($j_known, $j_aberrant, %options) = @_;

  my %edges;
  foreach my $j_raw ($j_known, $j_aberrant) {
    foreach my $edge (parse_junction($j_raw)) {
      push @{$edges{$edge->{base}}}, $edge;
    }
  }

  my $known_pos;
  my $aberrant_pos;

  foreach my $base (keys %edges) {
    $known_pos = $base if scalar @{$edges{$base}} == 2;
  }
  die unless $known_pos;

  foreach my $edge (parse_junction($j_aberrant)) {
    my $pos = $edge->{base};
    $aberrant_pos = $pos unless $pos == $known_pos;
  }
  die unless $aberrant_pos;

  return (
	  ($edges{$known_pos}->[0] || die),
	  ($edges{$aberrant_pos}->[0] || die)
	 );
}

sub junction2bed {
  my ($jf) = @_;
  my $cmd = sprintf 'junction2bed.pl -jf %s -now', $jf;
  system $cmd;
  die $cmd if $?;
  # create .bed version
}

sub summarize_aberrant {
  die;
}

sub digest_genes {
  my %wanted = map {$_, 1} split /,/, $FLAGS{"digest-genes"};
  my $glob_pattern = $FLAGS{glob} || die "-glob";
  my $glob_dir = $FLAGS{"glob-dir"} || ".";
  my @raw = glob(sprintf '%s/%s', $glob_dir, $glob_pattern);
  die unless @raw;
  my $infiles = filter_files(\@raw);

  if (0) {
    print STDERR "DEBUG: limit input\n";
    $infiles = [ @{$infiles}[0,1] ];
  }

  my %track;
  foreach my $infile (@{$infiles}) {
    printf STDERR "%s\n", $infile;
    my $df = new DelimitedFile("-file" => $infile,
			       "-headers" => 1,
			     );
    while (my $row = $df->get_hash()) {
      my $genes = $row->{genes} || next;
      my %found;
      foreach my $g (split /,/, $genes) {
	$found{$g} = 1 if $wanted{$g};
      }
      if (%found) {
	die "multiple genes!!" if scalar keys %found > 1;
	my ($gene) = keys %found;
	my $j = $row->{junction} || die;
	$row->{file} = $infile;
	push @{$track{$gene}{$j}}, $row;
      }
    }
  }

  my @fields = (
		"count",
		"qc_flanking",
	       );

  my %mini = (
	      "count" => "reads",
	      "qc_flanking" => "flnk"
	     );


  foreach my $gene (keys %track) {

    my $df = new DelimitedFile("-file" => $infiles->[0],
			       "-headers" => 1,
			     );
    my $outfile = sprintf 'digest_%s.tab', $gene;
    my $rpt = $df->get_reporter(
				"-file" => $outfile,
				"-extra" => [ "bed_name" ],
			       );


    foreach my $j (keys %{$track{$gene}}) {
      my $rows = $track{$gene}{$j};

      my %files = map {$_->{file}, 1} @{$rows};

      my @things;
      push @things, sprintf 's=%d', scalar keys %files;

      foreach my $f (@fields) {
	my @counts = sort {$b <=> $a} map {$_->{$f}} @{$rows};
	if (@counts > 10) {
	  @counts = @counts[0 .. 8];
	  push @counts, "...";
	}
	push @things, sprintf '%s=%s', ($mini{$f} || die), join ",", @counts;
      }

      my %r = %{$rows->[0]};
      $r{count} = scalar keys %files;
      # number of samples junction observed in
      $r{"bed_name"} = join ";", @things;
      $rpt->end_row(\%r);
    }

    $rpt->finish();
    junction2bed($outfile);

  }
}

sub digest_shared {

  my $DIGEST_BUFFER = 3000;

  my $infile = $FLAGS{"digest-shared"} || die;
  my $disease = $FLAGS{disease} || die "-disease";
  my $outfile = basename($infile) . ".digest.tab";
  my $f_refflat = $FLAGS{refflat} || die "-refflat";
  my $rff = new RefFlatFile();
  $rff->missing_genes_ok(1);
  # combined file contains data from multiple formats
  $rff->canonical_references_only(1);
  # ?
  printf STDERR "parsing %s...\n", $f_refflat;
  $rff->parse_file(
		   "-type" => "refflat",
		   "-refflat" => $f_refflat,
		  );
  my $cbm = new ChrBucketMap();
  $cbm->f_chr("chrom");
  $cbm->f_start("txStart");
  $cbm->f_end("txEnd");

  printf STDERR "building transcript db...\n";
  foreach my $row (@{$rff->rows}) {
    $cbm->add_row(
		  "-row" => $row
		 );
  }


  my $rpt = new Reporter(
			 "-file" => $outfile,
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   Chr
					   genes
					   strand
					   known_edge
					   aberrant_edge
					   aberrant_edge_buffered
					   aberrant_edge_buffered_interval
					   sample_count
					   samples
					)
				      ],
			 "-auto_qc" => 1,
			);

  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );

  while (my $row = $df->get_hash()) {

    my $hits = $cbm->find(
			  "-chr" => $row->{Chr} || die,
			  "-start" => $row->{known_edge},
			  "-end" => $row->{known_edge}
			 );
    dump_die($row, "WTF: no transcript matches") unless $hits;
    my %strands;
    my $query_edge = $row->{known_edge} || die;

    foreach my $h (@{$hits}) {
      my $usable;
      foreach my $exon (@{$h->{exons}}) {
	$usable = 1 if $exon->{start} == $query_edge
	  or $exon->{end} == $query_edge;
      }
      $strands{$h->{strand}} = 1 if $usable;
    }
    dump_die($row, "WTF: no strand info") unless %strands;
    dump_die($row, "WTF, multiple strands") if keys %strands > 1;
    my ($strand) = keys %strands;

    my @ae = split /,/, $row->{aberrant_edges};
    my @aec = split /,/, $row->{aberrant_edge_counts};
    my @aes = split /;/, $row->{aberrant_samples};
    die unless @ae == @aec;
    for (my $i=0; $i < @ae; $i++) {
      my %r = %{$row};
      my $ae = $ae[$i];
      my $aes = $aes[$i];
      $r{aberrant_edge} = $ae;
      $r{sample_count} = $aec[$i];
      $r{samples} = $aes;
      $r{strand} = $strand;

      die unless split(/,/, $aes) == $aec[$i];

      my $ae_buffered;
      if ($strand eq "+") {
	$ae_buffered = $ae - $DIGEST_BUFFER;
      } elsif ($strand eq "-") {
	$ae_buffered = $ae + $DIGEST_BUFFER;
      } else {
	die;
      }
      my @pos = sort {$a <=> $b} ($ae, $ae_buffered);
      $r{aberrant_edge_buffered} = $ae_buffered;
      $r{aberrant_edge_buffered_interval} = sprintf '%s:%d-%d', $row->{Chr}, @pos;

      $rpt->end_row(\%r);
    }
  }
  $rpt->finish();

}
