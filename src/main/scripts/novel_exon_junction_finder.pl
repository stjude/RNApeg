#!/usr/bin/env perl
# in corrected junction counts, find:
# - novel skips of known exons
# - novel exons (1 end bounding known exon)
# - events spanning transcripts ("cross-pollination" slides)
# - recurrent events across samples
# MNE 2/2013

use strict;
use warnings;

use SampleName;

use Getopt::Long;
use FloatingJunctionDetector;
use DelimitedFile;
use JunctionUtils qw(parse_junction);
use Reporter;
use FileUtils qw(read_simple_file);
use Counter;
use MiscUtils qw(dump_die);

use File::Basename;
use TdtConfig;


my %FLAGS;

my $MIN_READS = 0;
my $FIELD_JUNCTION = "junction";
my $FIELD_TYPE = "type";

GetOptions(\%FLAGS,
	   "-jf=s",
	   # junction file
	   "-jf-list=s",
	   # list of junction files
	   "-glob-cross",
	   "-glob=s",
	   # glob pattern for junction files

	   "-genome=s",
	   "-min-reads=i" => \$MIN_READS,
	   "-sample=s",

	   "-basename",
	   "-force",

	   "-field-junction=s" => \$FIELD_JUNCTION,
	   "-field-type=s" => \$FIELD_TYPE,

	   "-passthrough",
	   # report all entries

	   "-refflat=s",
	   "-refgene",
	  );

my $isoform_fn = $FLAGS{refflat};

unless ($isoform_fn) {
#  my $genome = $FLAGS{genome} || "GRCh37-lite";
  my $genome = $FLAGS{genome} || die "specify -genome";
  if ($genome) {
    my $config = &TdtConfig::readConfig("genome", $genome);
    if ($FLAGS{refgene}) {
      $isoform_fn = $config->{REFSEQ_NM_REFFLAT};
    } else {
      $isoform_fn = $config->{COMBINED_REFFLAT};
    }
    die unless $isoform_fn;
  }
}

die "no isoform file" unless $isoform_fn;

my $fjd = new FloatingJunctionDetector();
$fjd->add_database(
		   "-file" => $isoform_fn
		  );
my ($exon_db, $exon_starts, $exon_ends) = $fjd->get_exon_database();
# recycle

#die join "\n", sort {$a <=> $b} keys %{$exon_starts->{chr10}};

#die $exon_starts->{chr10}{100008748};

my @junction_files;
if (my $jf = $FLAGS{jf}) {
  push @junction_files, $jf;
} elsif (my $l = $FLAGS{"jf-list"}) {
  my $list = read_simple_file($l);
  push @junction_files, @{$list};
} elsif ($FLAGS{"glob-cross"}) {
  push @junction_files, glob("*.cross_sample_corrected.tab");
} elsif (my $pattern = $FLAGS{"glob"}) {
  push @junction_files, glob($pattern);
  die "no files for $pattern" unless @junction_files;
} else {
  die "specify -jf [junction file] | -jf-list [listfile] | -glob-cross\n";
}
die "no junction files" unless @junction_files;

my %track;
# for recurrence
my $recurrent_mode = @junction_files > 1;

my $sample_counter = 0;
my $c = new Counter(\@junction_files);

my $passthrough_mode = $FLAGS{passthrough};

foreach my $jf (@junction_files) {
  $c->next($jf);

  my $sample = $FLAGS{sample};
  unless ($sample) {
    my $bn = basename($jf);
    if ($bn =~ /^SJ/) {
      my %info = SampleName::parse($bn, "WARN");
      if (%info) {
	$sample = $info{sample} || die;
      }
    }
    unless ($sample) {
      $sample = "bogus_sample_" . ++$sample_counter;
    }
  }

  my $df = new DelimitedFile(
			     "-file" => $jf,
			     "-headers" => 1,
			    );
  # while (my $row = $df->next("-ref" => 1)) {  # headerless

#  my $outfile = basename($jf) . ".fun.tab";
  my $outfile = ($FLAGS{basename} ? basename($jf) : $jf) . ".fun.tab";

  if (-s $outfile and !$FLAGS{force}) {
    printf STDERR "%s exists, skipping\n", $outfile;
    next;
  }

  my $rpt = $df->get_reporter("-file" => $outfile,
			      "-extra" => [
					   "event",
					   "transcripts_start",
					   "transcripts_end",
					   "intra_transcript",
					  ]
			     );

  # my $rpt = new Reporter(
  # 			 "-file" => $outfile,
  # 			 "-delimiter" => "\t",
  # 			 "-labels" => [
  # 				       qw(
  # 					   event
  # 					   junction
  # 					   count
  # 					   genes
  # 					   transcripts_start
  # 					   transcripts_end
  # 					   intra_transcript
  # 					)
  # 				      ]
  # 			);


  while (my $row = $df->get_hash()) {
    my $type = $row->{$FIELD_TYPE} || die "no type field $FIELD_TYPE";
#    next unless $type eq "novel";

    if ($MIN_READS) {
      die "no count field" unless exists $row->{count};
      next unless $row->{count} >= $MIN_READS;
    }

    # only reporting on novel events
    #    my ($js, $je) = parse_junction($row->{junction} || die, "-ucsc" => 1);

    my $junction = $row->{$FIELD_JUNCTION} || die "no junction field $FIELD_JUNCTION";

    my $event;

    if ($type eq "novel") {
      my ($js, $je) = parse_junction($junction, "-ucsc" => 1);
#      dump_die($js, $junction);

      my $found_start = $exon_starts->{$js->{ref}}{$js->{base}};
      my $found_end = $exon_ends->{$je->{ref}}{$je->{base}};
      if ($found_start and $found_end) {
	$event = "novel_skip_of_known_exons";
      } elsif ($found_start) {
	$event = "matches_known_exon_end";
	# novel exon, starts at known exon's end
	# (i.e. start of junction skip)
      } elsif ($found_end) {
	$event = "matches_known_exon_start";
      }
      $row->{transcripts_start} = $found_start ? join(",", sort keys %{$found_start->{transcript2gene}}) : "";
      $row->{transcripts_end} = $found_end ? join(",", sort keys %{$found_end->{transcript2gene}}) : "";

      my @intra;
      # cases where the start and end appear within the same transcript
      if ($found_start and $found_end) {
	foreach my $t (keys %{$found_start->{transcript2gene}}) {
	  push @intra, $t if $found_end->{transcript2gene}{$t};
	}
      }
      $row->{intra_transcript} = join ",", sort @intra;
    } else {
      foreach my $f (qw(
			 transcripts_start
			 transcripts_end
			 intra_transcript
		      )) {
	$row->{$f} = "";
      }

 #     my ($js, $je) = parse_junction($junction, "-strip" => 1);
      # my ($js, $je) = parse_junction($junction, "-ucsc" => 1);
      # my $found_start = $exon_starts->{$js->{ref}}{$js->{base}};
      # my $found_end = $exon_ends->{$je->{ref}}{$je->{base}};

      # printf STDERR "DEBUG: j=%s ref=%s start=%d end=%s fs=%s fe=%s\n",
      # 	$junction,
      # 	  $js->{ref},
      # 	    $js->{base},
      # 	      $je->{base},
      # 		($found_start || "n"),
      # 		  ($found_end || "n");
      # debug/sanity check

    }

    if ($passthrough_mode ? 1 : $event) {
      #
      # per-sample report:
      #
      my %r = %{$row};
      $r{event} = $event || "";
      $rpt->end_row(\%r);

      if ($recurrent_mode) {
	# files for multiple samples being processed:
	# track events to find recurrency
	die "duplicate sample hit" if $track{$junction}{$sample};
	$row->{event} = $event;
	$track{$junction}{$sample} = $row;
      }
    }
  }
  $rpt->finish();
}

if ($recurrent_mode) {
  #
  #  look for recurrent events across samples
  #
  my $rpt = new Reporter(
			 "-file" => "recurrency.tab",
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   event
					   sample_count
					   samples
					   junction
					   counts
					   genes
					)
				      ]
			);


  foreach my $j (keys %track) {
    my @samples = sort keys %{$track{$j}};
    if (@samples > 1) {
      my $one_row = $track{$j}{$samples[0]};
      my %r = %{$one_row};
      # global
      $r{sample_count} = scalar @samples;
      $r{samples} = join ",", @samples;
      $r{counts} = join ",", map {$track{$j}{$_}{count}} @samples;
      $rpt->end_row(\%r);
    }
  }
  $rpt->finish();
}
