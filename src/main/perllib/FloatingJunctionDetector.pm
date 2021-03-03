package FloatingJunctionDetector;

use strict;
use warnings;

use File::Basename;
use Carp qw(confess cluck);
use MiscUtils qw(dump_die);

use Configurable;
use WorkingFile;
use FAI;

@FloatingJunctionDetector::ISA = qw(Configurable Exporter);
@FloatingJunctionDetector::EXPORT_OK = qw(
parse_junction
safe_get
);

use constant MATCHES_EXON_BOUNDARIES => "matches_exon_boundaries";
use constant MATCHES_SINGLE_EXON_BOUNDARY => "matches_single_exon_boundary";

use constant QC_FIELDS_NUMBER => qw(
qc_perfect_reads
qc_clean_reads
qc_flanking
qc_plus
qc_minus
    );

use constant QC_FIELDS_FLAG => qw(
    );

use MethodMaker qw(
                   fasta
                   fasta_dir

	           max_shift
	           bed_type

junctions_to_shift
junctions_db
database_files

jf_shift

min_reads_for_novel_junction
has_flanking_qc
strand
		  );

my $RGB_KNOWN = "0,128,0";
my $RGB_NOVEL = "192,0,0";

my $VERBOSE_EXONDB = 0;
my $VERBOSE_FIND = 0;
my $VERBOSE_COMBINE = 0;
my $VERBOSE_LOCK = 0;
my $VERBOSE_GENERATE = 0;
my $VERBOSE_TOSS_LOW_COVERAGE = 0;

use constant QC_FIELDS => qw(
			   qc_flanking
			   qc_plus
			   qc_minus
			   qc_perfect_reads
			   qc_clean_reads
			);

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->database_files([]);
  $self->junctions_db({});
#  $self->max_shift(15);
  $self->max_shift(20);
#  $self->max_shift(30);

#  $self->flank_check(20);
#  $self->flank_check(40);
#  $self->flank_check(50);

  $self->configure(%options);
  return $self;
}

sub load_junctions {
  my ($self, %options) = @_;
  my $fn = $options{"-file"} || die "-file";
  my $db_name = $options{"-db-name"};
#  confess keys %options;
  my $db_priority = $options{"-db-priority"};
  my $db_class = $options{"-db-class"};
  print STDERR "parse $fn\n";
  open(J, $fn) || die "can't open $fn: $!";
  my $j = $options{"-hash"} || {};
  my $first = 1;
  my $is_tcga = 0;
  my @headers;
  my $has_flanking_qc;
  my $strand_filter = $self->strand();

  while (<J>) {
    chomp;
    my @f = split /\t/, $_;

    if ($first) {
      # first line: disambiguate .bed from 5-field delimited
#      cluck "first line: $_";
      $first = 0;

      if ($f[0] eq "junction" and $f[1] eq "count") {
	$is_tcga = 1;
	@headers = @f;
	$has_flanking_qc = grep {$_ eq "qc_flanking"} @headers;
	$self->has_flanking_qc(1) if $has_flanking_qc;
	next;
      }
    }

    my ($start_chr, $end_chr, $start_base, $end_base, $count);
    my $is_refgene;
    if (@f == 2 or $is_tcga) {
      # TCGA junction file
      $count = $f[1];
      next unless $count;
      # count of supporting reads: often 0 in cufflinks output (?)
      my @j = split /,/, $f[0];
      die unless @j == 2;
      my ($from, $to) = @j;
      ($start_chr, $start_base) = parse_junction($j[0]);
      ($end_chr, $end_base) = parse_junction($j[1]);
      die if $start_base > $end_base;
      # ever happens?
      die if $start_chr ne $end_chr;
      #    die "duplicate $start_chr $start_base $end_base" if $j{$start_chr}{$start_base}{$end_base};
      # TopHat seems to have a separate entry for +/+ and -/- (w/same counts??)

      if ($has_flanking_qc) {
	my %r;
	@r{@headers} = @f;

	foreach my $f (QC_FIELDS) {
	  $j->{$start_chr}{$start_base}{$end_base}{$f} = $r{$f};
	}
      }

    } elsif (@f == 3) {
      # combined junctions file
      # e.g. /nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/mRNA/Combined/all_junct2gene.txt
      next if $f[0] eq "Junction";
      $is_refgene = 1;
      # sorta
      my @j = split /,/, $f[0];
      die unless @j == 2;
      my ($from, $to) = @j;
      ($start_chr, $start_base) = parse_junction($j[0]);
      ($end_chr, $end_base) = parse_junction($j[1]);
      die if $start_base > $end_base;
      # ever happens?
      die if $start_chr ne $end_chr;

      my $gene = $f[1];
      my @accs = split /,/, $f[2];

      foreach my $acc (@accs) {
	$j->{$start_chr}{$start_base}{$end_base}{transcript2gene}{$acc}{$gene} = 1;
      }
    } elsif (@f == 5) {
      # BED
      my $bed_type = $self->{bed_type} || die "for BED input, need -bed_type [star-1-based-intron]";
      if ($bed_type eq "star-1-based-intron") {
	# STAR .bed format:
	# 1. reports 1-based base numbers of intron/skip bounds
	# 2. use of 1-based start position appears to violate BED spec
	($start_chr, $start_base, $end_base, $count) = @f;
	$start_base--;
	# convert 1st base of intron to last base of previous exon
	$end_base++;
	# convert last base of intron to 1st base of next exon
	$end_chr = undef;
#	$start_chr =~ s/^chr//;
      } else {
	die "unknown bed type $bed_type";
      }
    } elsif (@f == 16 or @f == 13 or @f == 11) {
      # parse exon annotations from
      #   16: refGene.txt
      #   13: knownGene_refFlat.txt
      #   11: acembly.txt
      $is_refgene = 1;
      my ($nm, $start_chr, $strand, $starts, $ends, $gene) = @f[1, 2, 3, 9, 10, 12];
      if (@f == 11) {
	# acembly: parse gene sym from transcript ID
	$gene = (split(/\./, $nm))[0];
      }
#      $start_chr =~ s/^chr//;

      my $usable = 1;
      if ($strand_filter) {
	$usable = 0 unless $strand_filter eq $strand;
      }

      if ($usable) {
	my @starts = split /,/, $starts;
	my @ends = split /,/, $ends;
	die unless @starts == @ends;
	for (my $i = 0; $i < @starts - 1; $i++) {
	  # convert exon locations to junctions
	  my $this_end = $ends[$i];
	  # ends are 1-based
	  my $next_start = $starts[$i + 1] + 1;
	  # starts are 0-based

	  if ($db_name) {
	    # databases are being loaded separately with names
	    # and priorities.
	    if (my $ref = $j->{$start_chr}{$this_end}{$next_start}) {
	      # entry already exists
	      my $p = $ref->{db_priority} || die "no priority for existing entry";
	      die unless $db_priority;
#	    printf STDERR "priority: existing:%d new:%d\n", $p, $db_priority;
	      next unless $db_priority < $p;
	      # respect existing entry unless new entry is higher priority
	      # (lower values are higher priority)
	    }

#	  print STDERR "DEBUG $start_chr $this_end $next_start\n";
	    $j->{$start_chr}{$this_end}{$next_start}{db_name} = $db_name;
	    $j->{$start_chr}{$this_end}{$next_start}{db_priority} = $db_priority;
	    $j->{$start_chr}{$this_end}{$next_start}{db_class} = $db_class;
	  }

	  $j->{$start_chr}{$this_end}{$next_start}{transcript2gene}{$nm}{$gene} = 1;
	  # same junction can appear in multiple transcripts/genes

#	printf STDERR "generating %s: %d-%d\n", $gene, $this_end, $next_start;
	}
      }
    } else {
      die sprintf "unknown annotation file format, %d columns: ", scalar(@f), $_;
    }

    unless ($is_refgene) {
      $j->{$start_chr}{$start_base}{$end_base}{count} = $count;
    }
  }
  close J;
  return $j;
}

sub parse_junction {
  # STATIC
  my ($j) = @_;
  my @f = split /:/, $j;
  die unless @f == 3;
  my ($ref, $base) = @f;
#  $ref =~ s/^chr//;
  # standardize
  return ($ref, $base);
}

sub add_database {
  # load database of known junctions to compare against
  # (may be called multiple times)
  my ($self, %options) = @_;
  my $jf_db = $options{"-file"} || die "-file";
  push @{$self->database_files}, $jf_db;
  $self->load_junctions(
			%options,
			"-hash" => $self->junctions_db || die,
		       );
}

sub find_shifts {
  #
  #  identify shiftable junctions in specified set
  #
  my ($self, %options) = @_;
  my $novel_mode = $options{"-novel-mode"};
  my $novel_exon_mode = $options{"-novel-exon-mode"};

  my ($SELF_COMPARE_MODE, $j1, $j2);
  my ($exon_db, $exon_starts, $exon_ends);

  if ($novel_mode) {
    # compare novel-only junctions vs. themselves (2nd pass)
    $SELF_COMPARE_MODE = 1;
    $j1 = $self->junctions_to_shift() || die "shift junctions must be loaded already";
    $j2 = $self->generate_novel_junction_db();
    $self->flag_junctions_matching_exon_boundaries();
  } elsif ($novel_exon_mode) {
    $SELF_COMPARE_MODE = 0;
    $j1 = $self->junctions_to_shift() || die "shift junctions must be loaded already";
    # should be already flagged too (novel mode run already)
    $j2 = $self->junctions_db() || die "no junction db";
    # to skip junctions matching db already
    ($exon_db, $exon_starts, $exon_ends) = $self->get_exon_database();
    # for actual checking
  } else {
    my $jf_shift = $options{"-jf-shift"} || die "-jf-shift";
    # junctions to shift/correct
    $self->jf_shift($jf_shift);
    my $db_files = $self->database_files();
    $SELF_COMPARE_MODE = (@{$db_files} == 1 and $jf_shift eq $db_files->[0]);
    $j1 = $self->load_junctions("-file" => $jf_shift);
    $self->junctions_to_shift($j1);
    $j2 = $self->junctions_db() || die "no junction database";
  }
  print STDERR "self-compare mode\n" if $SELF_COMPARE_MODE;

  my $max_shift = $self->max_shift || die "max_shift";

#  my $FLANK_CHECK = $self->flank_check || die "flank_check";
  my $FLANK_CHECK = $max_shift + 12;
  $FLANK_CHECK = 20 if $FLANK_CHECK < 20;
  # adapt rather than hardcode:
  # MUST be larger than max_shift to avoid substr-related crashes below.
  # extent to which it's larger is the # of bases beyond the junction
  # for $max_shift checks

  my @deltas = grep {$_ != 0} (- $max_shift .. $max_shift);

  my $total_junctions = 0;
  my $shiftable_junctions = 0;

#  my %wanted_ref = map {$_, 1} (1 .. 22, qw(X Y M MT));

  my $fasta = $self->fasta();

  # legacy:
  my $fa_dir;
  unless ($fasta) {
    $fa_dir = $self->fasta_dir() || die "fasta_dir";
    printf STDERR "FASTA dir: %s\n", $fa_dir;
    die "no genome FASTA" unless $fa_dir and -d $fa_dir;
  }

  foreach my $ref (sort keys %{$j1}) {
#    unless ($wanted_ref{$ref}) {
#      printf STDERR "ignoring reference %s\n", $ref;
#      next;
#    } else {
    print STDERR "process reference $ref...\n";
#    }

    my $reference_seq = "";

    if ($fasta) {
      # preferred
      my $fai = new FAI("-fasta" => $fasta);
      # hack: new instance with every query to avoid sequence caching
      if ($fai->find_name($ref)) {
	my $seq_ref = $fai->get_sequence("-id" => $ref);
	$reference_seq = $$seq_ref;
	# bleah
      } else {
	printf STDERR "can't find FASTA entry for %s, skipping\n", $ref;
	next;
      }
    } else {
      # legacy/deprecated
      my @fa_try;

      my $thing = $ref;
      $thing =~ s/^chr//;
      push @fa_try, sprintf '%s/chr%s.fa', $fa_dir, $thing;
      push @fa_try, sprintf '%s/%s.fa', $fa_dir, $thing;
      # accept either plain or UCSC filenames
      if ($thing eq "M") {
	# standardization from SplicedReadReporter
	push @fa_try, sprintf '%s/chrMT.fa', $fa_dir;
	push @fa_try, sprintf '%s/MT.fa', $fa_dir;
      }

      my $fa;
      foreach (@fa_try) {
	if (-s $_) {
	  $fa = $_;
	  last;
	}
      }

      unless ($fa and -s $fa) {
	printf STDERR "can't find FASTA file for %s, tried %s: skipping\n", $ref, join ",", @fa_try;
	next;
      }

      print STDERR "loading $fa...";
      open(FA, $fa) || die;
      my $first = <FA>;
      while (<FA>) {
	chomp;
	$reference_seq .= uc($_);
      }
      print STDERR "\n";
    }

    my $reference_len = length($reference_seq);
    my $max_check_len = $reference_len - $FLANK_CHECK;

    my %saw;

    foreach my $start (sort {$a <=> $b} keys %{$j1->{$ref}}) {
      foreach my $end (sort {$a <=> $b} keys %{$j1->{$ref}{$start}}) {
	$total_junctions++;
	my $total_matches = 0;
	printf STDERR "junction:%s:%d-%d\n", $ref, $start, $end if $VERBOSE_FIND;

	bounds_check($end, $reference_len);

	unless ($SELF_COMPARE_MODE) {
#	  next if $j2->{$ref}{$start}{$end};
	  next if safe_get($j2, $ref, $start, $end);
	  # If a junction in the shift set is already in the reference
	  # set, respect this and never attempt to shift it.
	  # (unless of course we're looking for this sort of thing,
	  # [all_junct2gene.txt, I'm looking in your direction])
	}

	if ($j1->{$ref}{$start}{$end}{obsolete}) {
	  if ($novel_mode or $novel_exon_mode) {
	    next;
	  } else {
	    die "say what? obsolete j1" ;
	  }
	}

	if (($novel_mode or $novel_exon_mode) and 
	    $j1->{$ref}{$start}{$end}{MATCHES_EXON_BOUNDARIES()}) {
	  # never attempt to shift a novel junction that already matches
	  # known exon boundaries.  Other junctions ambiguous with 
	  # this site should be shifted here instead (should happen
	  # in later iteration)
	  printf STDERR "considering %s:%d-%d locked (double), skipping\n", $ref, $start, $end if $VERBOSE_COMBINE;
	  next;
	}

	if (($novel_mode or $novel_exon_mode) and 
	    $j1->{$ref}{$start}{$end}{MATCHES_SINGLE_EXON_BOUNDARY()}) {
	  printf STDERR "considering %s:%d-%d locked (single), skipping\n", $ref, $start, $end if $VERBOSE_COMBINE;
	  next;
	}
	    

	foreach my $delta (@deltas) {
	  my $s2 = $start + $delta;
	  my $e2 = $end + $delta;

	  my $possible;

	  if ($novel_exon_mode) {
	    if ($novel_exon_mode == 2) {
	      # a start or end matches (i.e. novel exon)
	      $possible = ($exon_starts->{$ref}{$s2} or
			   $exon_ends->{$ref}{$e2}) ? 1 : 0;
	    } else {
	      # both start and end match (i.e. novel skip of known exon)
	      $possible = ($exon_starts->{$ref}{$s2} and
			   $exon_ends->{$ref}{$e2}) ? 1 : 0;
	    }
	  } else {
#	    $possible = $j2->{$ref}{$s2}{$e2};
	    $possible = ($j2->{$ref}{$s2} and $j2->{$ref}{$s2}{$e2});
	    # don't create blank $s2 entries
	  }

	  if ($possible) {
	    #
	    #  possible shift:
	    #
	    printf STDERR "possible:%d-%d vs %d-%d\n", $start, $end, $s2, $e2 if $VERBOSE_FIND;

	    next unless bounds_check($e2, $reference_len,
				     ($novel_exon_mode and $novel_exon_mode == 2));

	    my $compare_key = join "_", sort
		sprintf("%s>%s", $start, $end),
		sprintf("%s>%s", $s2, $e2);
	    if ($SELF_COMPARE_MODE) {
	      next if $saw{$compare_key};
	      # already processed
	      $saw{$compare_key} = 1;
	    }

	    print STDERR "check $delta $ref:$start-$end $ref:$s2-$e2\n" if $VERBOSE_FIND;

	    if ($start < $FLANK_CHECK or $s2 < $FLANK_CHECK) {
	      printf STDERR "can't compare, too close to ref start\n" if $VERBOSE_FIND;
	      next;
	    }

	    if ($end > $max_check_len or $e2 > $max_check_len) {
	      printf STDERR "can't compare, too close to ref end\n" if $VERBOSE_FIND;
	      next;
	    }

	    my $before = substr($reference_seq, $start - $FLANK_CHECK, $FLANK_CHECK);
	    my $after = substr($reference_seq, $end - 1, $FLANK_CHECK);

	    my $b2 = substr($reference_seq, $s2 - $FLANK_CHECK, $FLANK_CHECK);
	    my $a2 = substr($reference_seq, $e2 - 1, $FLANK_CHECK);
	    printf STDERR "debug, %d\n%s-%s\n%s-%s\n", $delta, $before, $after, $b2, $a2 if $VERBOSE_FIND;

	    my $second_v2;
	    if ($delta > 0) {
	      my $b3 = substr($reference_seq, $s2 - $FLANK_CHECK - $delta,
			      $FLANK_CHECK + $delta);
	      my $a3 = substr($reference_seq, $e2 - 1, $FLANK_CHECK - $delta);

	      $second_v2 = $b3 . $a3;
	    } else {
	      my $b3 = substr($reference_seq, ($s2 - $FLANK_CHECK) + abs($delta),
			      $FLANK_CHECK - abs($delta));
	      my $a3 = substr($reference_seq, $e2 - 1, $FLANK_CHECK + abs($delta));

	      $second_v2 = $b3 . $a3;
	    }


	    my $first = $before . $after;
	    # first junction sequence (site and flanking)

	    my $second = $b2 . $a2;
	    # second junction sequence (site and flanking)

	    my $match;

	    if (1) {
	      $match = $first eq $second_v2;
	      my $insane = $FLANK_CHECK * 3;
	      die "WTF1" if length($first) > $insane;
	      die sprintf "WTF2 %s", join ",", length($second_v2), $delta, $s2, $e2, length($reference_seq) if length($second_v2) > $insane;
	      die "epic fail" unless length($first) eq length($second_v2);
	    } else {
	      my ($search_seq, $lookup);
	      if ($delta < 0) {
		# shift 2nd sequence right to look up in 1st sequence
		$lookup = substr($second, abs($delta));
		$search_seq = $first;
	      } else {
		# shift 1st sequence left to look up in 2nd
		$lookup = substr($first, $delta);
		$search_seq = $second;
	      }
	      $match = index($search_seq, $lookup) == 0;
	    }

	    # alternate/simpler approach that specifically looks for
	    # an ambiguous sequence on either side of the junction?
	    # -> look for cases where methods disagree??
	    # - case 1: maps to 1st, skips second
	    # - case 2: skips 1st, maps to second
	    # => that's essentially what above code does,
	    #    it just checks for more

	    # additional checks?:
	    # check aligned read after change, compare mismatch counts?

	    if ($match) {
	      printf STDERR "HIT %s:%d-%d vs %d-%d delta=%d\n",
		$ref,
		  $start, $end,
		    $s2, $e2,
		      $delta if $VERBOSE_FIND;
	      $total_matches++;

	      my $proceed = 1;
	      if ($novel_mode) {
		#
		#  assign the shift to the junction with the higher
		#  read count.  This may mean deferring the shift
		#  until the inverse comparison is made, to avoid
		#  gunking up the code below (particularly "delta")
                #
		my $count_here = $j1->{$ref}{$start}{$end}{count} || die;
		my $count_there = $j1->{$ref}{$s2}{$e2}{count} || die;

		printf STDERR "check shift of %s:%d-%d (%d) to %d-%d (%d)\n",
		$ref, $start, $end, $count_here, $s2, $e2, $count_there if $VERBOSE_FIND;

		if ($count_here > $count_there) {
		  # more reads at other site: potentially defer
		  if ($j1->{$ref}{$s2}{$e2}{MATCHES_EXON_BOUNDARIES()}) {
		    # always merge to target if target matches exon bounds
		    printf STDERR "prefer reads at target, but target matches double exon boundaries, so merging anyway\n" if $VERBOSE_FIND;
		  } elsif ($j1->{$ref}{$s2}{$e2}{MATCHES_SINGLE_EXON_BOUNDARY()}) {
		    printf STDERR "prefer reads at target, but target matches single exon boundary, so merging anyway\n" if $VERBOSE_FIND;

		  } else {
		    # defer
		    printf STDERR "deferring shift of %s:%d-%d (%d) to %d-%d (%d)\n",
		    $ref, $start, $end, $count_here, $s2, $e2, $count_there if $VERBOSE_FIND;

#		    die "gurgle1" if $j1->{$ref}{$start}{$end}{obsolete};
#		    die "gurgle2" if $j1->{$ref}{$s2}{$e2}{obsolete};
		    
		    delete $saw{$compare_key};
		    # remove lock for this comparison
		    $proceed = 0;
		    die "say what now?" unless $s2 > $start;
		  }
		}
	      }

	      if ($proceed) {
		push @{$j1->{$ref}{$start}{$end}{matches}}, {
		  start => $s2,
		  end => $e2,
		  delta => $delta
		};
	      }

	    } else {
#	      printf STDERR "MISS $delta: $search_seq $lookup\n";
	      printf STDERR "MISS $delta: $first $second_v2\n" if $VERBOSE_FIND;
	    }
	  }
	}
	$shiftable_junctions++ if $total_matches;
	printf STDERR "WARNING: %d plausible shifts for %s:%d-%d!\n", $total_matches, $ref, $start, $end if $total_matches > 1;
      }
    }
  }
  printf STDERR "total junctions:%d  shiftable:%d\n", $total_junctions, $shiftable_junctions;
  die "QC error: no junctions, really??" unless $total_junctions;

}

sub get_exon_database {
  my ($self) = @_;
  my %eb;
  my %starts;
  my %ends;
  my $db = $self->junctions_db() || die "no junction database";
  confess "empty junction database" unless %{$db};
  dead_key_check($db);
  foreach my $ref (keys %{$db}) {
    foreach my $start (sort {$a <=> $b} keys %{$db->{$ref}}) {
#      printf STDERR "exondb %s site:%d\n", $ref, $start;
      my @ends = keys %{$db->{$ref}{$start}};
      die "WTF, no ends for $ref $start" unless @ends;
      # make sure bogus keys haven't been created by hash queries
#      $eb{$ref}{$start} = 1;
#      $starts{$ref}{$start} = 1;

      foreach my $end (sort {$a <=> $b} @ends) {
#	$eb{$ref}{$end} = 1;
#	$ends{$ref}{$end} = 1;
	printf STDERR "exondb %s %d-%d site:%d site:%d %s\n", $ref, $start, $end, $start, $end, join ",", sort keys %{$db->{$ref}{$start}{$end}{transcript2gene}} if $VERBOSE_EXONDB;
	my $ts = $db->{$ref}{$start}{$end}{transcript2gene};
	foreach my $t (keys %{$ts}) {
	  foreach my $g (keys %{$ts->{$t}}) {

	    $eb{$ref}{$start}{transcript2gene}{$t}{$g} = 1;
	    $starts{$ref}{$start}{transcript2gene}{$t}{$g} = 1;

	    $eb{$ref}{$end}{transcript2gene}{$t}{$g} = 1;
	    $ends{$ref}{$end}{transcript2gene}{$t}{$g} = 1;

#	    printf STDERR "t=%s g=%s\n", $t, $g;
	  }
	}
      }
    }
  }
  die unless wantarray();
  return (\%eb, \%starts, \%ends);
}

sub flag_junctions_matching_exon_boundaries {
  my ($self) = @_;
  #
  #  map of exon boundaries in junction database:
  #
  my ($eb, $exon_starts, $exon_ends) = $self->get_exon_database();

  #
  #  flag junctions matching these boundaries:
  #
  my $js = $self->junctions_to_shift() || die "shift junctions must be loaded already";
  
  if (0) {
    foreach my $base (qw(2055352 2152762)) {
      printf STDERR "WTF %s: %s\n", $base, $eb->{11}{$base};
    }
    die "FIX ME";
  }

  foreach my $ref (keys %{$js}) {
    foreach my $start (keys %{$js->{$ref}}) {
      foreach my $end (keys %{$js->{$ref}{$start}}) {
	my $flag = ($exon_starts->{$ref}{$start} and
		    $exon_ends->{$ref}{$end}) ? 1 : 0;
	if ($flag) {
	  printf STDERR "locking %s:%d-%d\n", $ref, $start, $end if $VERBOSE_LOCK;
	  die "gah" unless (
	    ($exon_starts->{$ref}{$start} or $exon_ends->{$ref}{$start}) and
	    ($exon_starts->{$ref}{$end} or $exon_ends->{$ref}{$end})
	      );
	  my $entry = $js->{$ref}{$start}{$end};
	  $entry->{MATCHES_EXON_BOUNDARIES()} = 1;
	} else {
	  my $single = ($exon_starts->{$ref}{$start} or $exon_ends->{$ref}{$end}) ? 1 : 0;
	  if ($single) {
	    my $entry = $js->{$ref}{$start}{$end} || die;
	    $entry->{MATCHES_SINGLE_EXON_BOUNDARY()} = 1;
	  }
	}

      }
    }
  }

}

sub combine_junctions {
  #
  # combine junctions based on analysis
  #
  my ($self, %options) = @_;

  my $js = $self->junctions_to_shift || die;
  my $db = $self->junctions_db || die;
  my $novel_mode = $options{"-novel-mode"};

#  dead_key_check($js); dead_key_check($db);

  foreach my $ref (sort keys %{$js}) {
#    dead_key_check($js); dead_key_check($db);
    foreach my $start (sort {$a <=> $b} keys %{$js->{$ref}}) {
      foreach my $end (sort {$a <=> $b} keys %{$js->{$ref}{$start}}) {
	my $entry = $js->{$ref}{$start}{$end};
	if ($entry->{matches}) {
	  # junction can be shifted
	  if (@{$entry->{matches}} > 1) {
	    printf STDERR "WARNING: %d plausible matches for %s\n", scalar @{$entry->{matches}}, format_junction($ref, $start, $end);
	    foreach my $m (@{$entry->{matches}}) {
	      my ($s2, $e2) = @{$m}{qw(start end)};
	      printf STDERR "  matches %d-%d", $s2, $e2;
	      my $db_entry = safe_get($db, $ref, $s2, $e2);
	      printf STDERR " %s %s", $db_entry->{db_name}, $db_entry->{db_priority} if $db_entry and $db_entry->{db_name};
	      # db entry may not be available in novel junction mode
	      print STDERR "\n";
	    }
	  }

	  my $match = $entry->{matches}->[0];

	  my @match_ref_exons;
	  foreach my $e (@{$entry->{matches}}) {
	    printf STDERR "try %s\n", join " ", $ref, $e->{start}, $e->{end} if $VERBOSE_COMBINE;

#	    if ($js->{$ref}{$e->{start}}{$e->{end}}{MATCHES_EXON_BOUNDARIES()}) {
	    my $m = safe_get($js, $ref, $e->{start}, $e->{end});
	    if ($m and $m->{MATCHES_EXON_BOUNDARIES()}) {
	      push @match_ref_exons, $e;
	    }
	  }

	  if (@match_ref_exons) {
	    die "ambig exon match possibilities" if @match_ref_exons > 1;
	    $match = $match_ref_exons[0];
	    printf STDERR "found match using reference exon boundaries: %d-%d\n", $match->{start}, $match->{end};
	  } else {
	    # choose the target with the highest junction count
	    # in the source list
	    my $best_count = 0;
	    foreach my $e (@{$entry->{matches}}) {
#	      my $t = $js->{$ref}{$e->{start}}{$e->{end}};
	      my $t = safe_get($js, $ref, $e->{start}, $e->{end});
	      if ($t) {
		# target exists in source list
		my $count = $t->{count};
		if ($count and $count > $best_count) {
		  printf STDERR "upgrade %d -> %d\n", $best_count, $count if $best_count > 0;
		  $match = $e;
		  $best_count = $count;
		}
	      }
	    }
	  }

	  my ($s2, $e2) = @{$match}{qw(start end)};
	  printf STDERR "shift %s:%d-%d by %d => %d-%d, count=%d\n", $ref, $start, $end, $match->{delta}, $s2, $e2, $entry->{count} if $VERBOSE_COMBINE;

	  my $entry_target = $js->{$ref}{$s2}{$e2};
	  if ($entry_target) {
	    printf STDERR "existing entry for target %s:%d-%d, count=%d, adding %d\n", $ref, $s2, $e2, $entry_target->{count}, $entry->{count} if $VERBOSE_COMBINE;
	    $entry_target->{combine_count}++;
	  } else {
	    # no junction entry for this site, create
	    $entry_target = $js->{$ref}{$s2}{$e2} = {};
	    $entry_target->{created} = 1;
	  }

	  if ($entry_target->{obsolete}) {
	    if ($novel_mode) {
	      # in novel junction mode, shifts may step on each other,
	      # despite our attempts to avoid this by preferring to
	      # merge ambiguous junctions with the highest-read
	      # alternative site available.
	      my $levels = 0;
	      printf STDERR "WARNING: target for $ref $start-$end ($s2-$e2) is obsolete; following breadcrumb trail...\n";
	      while ($entry_target->{obsolete}) {
		printf STDERR "CHECK ME: following move from %d-%d to %d-%d, depth=%d\n", $s2, $e2, @{$entry_target->{moved_to}}, ++$levels;
		($s2, $e2) = @{$entry_target->{moved_to}};
		$entry_target = $js->{$ref}{$s2}{$e2} || die;
	      }
	    } else {
	      # shouldn't happen in traditional mode vs. database
	      die "WTF, target for $ref $start-$end ($s2-$e2) is obsolete!";
	    }
	  }

	  $self->pool_info($entry, $entry_target);

	  push @{$entry_target->{combine_info}}, sprintf "%d-%d:%d", $start, $end, $entry->{count};

	  die "WTF2, target for $ref $start-$end ($s2-$e2) is obsolete AGAIN!" if $entry_target->{obsolete};
	  die sprintf "WTF, %d-%d is obsolete already", $start, $end if $entry->{obsolete};
	  $entry->{obsolete} = 1;
	  $entry->{moved_to} = [$s2, $e2];
	  die "multi sets of matches" if $entry->{matches_old};
	  $entry->{matches_old} = $entry->{matches};
	  delete $entry->{matches};

	}
      }
    }
  }
}

sub write_junctions {
  my ($self, %options) = @_;
  my $outfile = $options{"-file"} || die "-file";
  my $annotate = $options{"-annotate"};
  my $bed = $options{"-bed"};
  my $js = $self->junctions_to_shift() || die;
  my $js_db = $self->junctions_db() || die;
  my $wf = new WorkingFile($outfile);
  my $fh = $wf->output_filehandle();

  my $min_novel_reads = $self->min_reads_for_novel_junction() || 0;

  if ($annotate) {
    my @labels = qw(
                   junction
	           count
		   type
		   genes
		   transcripts
		   );
    if ($self->has_flanking_qc()) {
      push @labels, QC_FIELDS;
    }
    printf $fh "%s\n", join "\t", @labels;
  } elsif ($bed) {
    my $track_name = "junctions";
    my $track_desc = "BAM junctions";

    if (my $fn = $self->jf_shift) {
      $fn =~ s/\.bam.*//;
      $track_name = $track_desc = $fn;
    }

#    print $fh "track name=junctions description=\"BAM junctions\" visibility=3 itemRgb=\"On\"\n";
    printf $fh "track name=%s description=\"%s\" visibility=3 itemRgb=\"On\"\n", $track_name, $track_desc;
  }

  foreach my $ref (sort keys %{$js}) {
    foreach my $start (sort {$a <=> $b} keys %{$js->{$ref}}) {
      foreach my $end (sort {$a <=> $b} keys %{$js->{$ref}{$start}}) {
	my $entry = $js->{$ref}{$start}{$end};
	next if $entry->{obsolete};
	# junction moved elsewhere
#	my $j = sprintf "chr%s:%d:+,chr%s:%d,+", $ref, $start, $ref, $end;
	# FAIL

#	my $j = sprintf "chr%s:%d:+,chr%s:%d:+", $ref, $start, $ref, $end;
	my $j = format_junction($ref, $start, $end);

	unless ($js_db->{$ref}{$start}{$end}) {
	  #
	  #  filters for NOVEL junctions only:
	  #
	  if ($entry->{count} < $min_novel_reads) {
	    # doesn't meet minimum read support
	    printf STDERR "tossing novel %s: only %d reads\n", $j, $entry->{count} if $VERBOSE_TOSS_LOW_COVERAGE;
	    next;
	  }

	  foreach my $f (QC_FIELDS) {
	    die "missing QC field $f" unless exists $entry->{$f};
	  }

	  unless ($entry->{qc_flanking}) {
	    # doesn't meet minimum flanking sequence
	    printf STDERR "tossing novel %s: fails flanking QC, raw count: %d\n", $j, $entry->{count};
	    # separate field for count so sortable / QC
	    next;
	  }

	  my $total_decent = $entry->{qc_perfect_reads} +
	      $entry->{qc_clean_reads};
	  # make sure we pool these for some checks, since 
	  # qc_clean_reads does NOT include perfect reads

	  if (($entry->{qc_plus} and $entry->{qc_minus}) or
	      $entry->{qc_perfect_reads} >= 1 or
	      $total_decent >= 2) {
	    # acceptable
	  } else {
	    printf STDERR "tossing novel %s: not bi-directional, fails read counts\n", $j;
	    next;
	  }
	}

	if ($annotate) {
	  my @f;
	  push @f, $j;
	  push @f, $entry->{count};
	  my (@genes, @transcripts);
	  my $type;
	  if (my $known = $js_db->{$ref}{$start}{$end}) {
	    $type = "known";
	    die "no transcript/gene info" unless $known->{transcript2gene};
	    foreach my $nm (sort keys %{$known->{transcript2gene}}) {
	      foreach my $gene (sort keys %{$known->{transcript2gene}{$nm}}) {
		push @transcripts, $nm;
		push @genes, $gene;
	      }
	    }
	  } else {
	    $type = "novel";
	  }

	  push @f, $type;
	  push @f, join ",", @genes;
	  push @f, join ",", @transcripts;
	  if ($self->has_flanking_qc()) {
	    push @f, map {$entry->{$_}} QC_FIELDS;
	  }
	  printf $fh "%s\n", join "\t", @f;
	} elsif ($bed) {
	  # UCSC .bed format
	  my @f;

	  my $ref_clean = $ref;
#	  $ref_clean = "chr" . $ref unless $ref =~ /^chr/;
	  # .bed requires chrX format
	  # update: since we're no longer munging formatting
	  # (e.g. for zebrafish) leave this alone

	  push @f, $ref_clean;
	  # 1. reference/chrom name
	  push @f, $start - 1;
	  # 2. start base number (0-based)
	  push @f, $end;
	  # 3. end base number: base is not counted, so don't need to adjust
	  # :/
	  push @f, $entry->{count};
	  # 4. count
	  push @f, $entry->{count};
	  # 5. score (count, again)
	  push @f, "+";
	  # 6. strand
	  push @f, $f[1];
	  push @f, $f[2];
	  # 7. thickStart
	  # 8. thickEnd

	  my $color;
	  if ($js_db->{$ref}{$start}{$end}) {
	    $color = $RGB_KNOWN;
	  } else {
	    $color = $RGB_NOVEL;
	  }

	  push @f, $color;
	  # 9. itemRGB
	  printf $fh "%s\n", join "\t", @f;
	} else {
	  printf $fh "%s\n", join "\t", $j, $entry->{count};
	}
      }
    }
  }
  $wf->finish();
}

sub format_junction {
  # STATIC
  my ($ref, $start, $end) = @_;
#  return sprintf "chr%s:%d:+,chr%s:%d:+", $ref, $start, $ref, $end;
  return sprintf "%s:%d:+,%s:%d:+", $ref, $start, $ref, $end;
}

sub generate_novel_junction_db {
  my ($self) = @_;
  my $js = $self->junctions_to_shift() || die "shift junctions must be loaded already";
  my $db = $self->junctions_db() || die "need junction db";
  my %j;

  dead_key_check($js);
  dead_key_check($db);

  foreach my $ref (keys %{$js}) {
    foreach my $start (sort {$a <=> $b} keys %{$js->{$ref}}) {
      foreach my $end (sort {$a <=> $b} keys %{$js->{$ref}{$start}}) {
	my $entry = $js->{$ref}{$start}{$end};
	next if $entry->{obsolete};
	# junction has already been moved elsewhere
	my $known = safe_get($db, $ref, $start, $end);
	if ($known) {
	  # known junction: don't use
	  die "where is transcript2gene for known" unless $known->{transcript2gene};
	} else {
	  # novel junction: usable
	  $j{$ref}{$start}{$end} = $entry;
	  printf STDERR "generated $ref:$start-$end\n" if $VERBOSE_GENERATE;
	}
      }
    }
  }

  dead_key_check($js);
  dead_key_check($db);

  return \%j;
}

sub safe_get {
  # STATIC
  # query hash without creating keys and filling it up with junk
  my ($hash, $ref, $start, $end) = @_;
  my $result;
  if (exists $hash->{$ref} and
      exists $hash->{$ref}{$start} and
      exists $hash->{$ref}{$start}{$end}) {
    $result = $hash->{$ref}{$start}{$end};
  }
  return $result;
}

sub dead_key_check {
  # STATIC
#  cluck("dead key check...");
  my ($js) = @_;
  my $error;
  foreach my $ref (sort keys %{$js}) {
    foreach my $start (keys %{$js->{$ref}}) {
      my @ends = keys %{$js->{$ref}{$start}};
      unless (@ends) {
	cluck "FIX ME: hash corrupted, empty key for $ref $start\n" unless @ends;
	$error = 1;
      }
    }
  }
  die "gurgle" if $error;
}



sub get_transcripts_and_genes {
# 	my $t2g = $db->{$ref}{$start}{$end}{transcript2gene} || die "no transcript2gene";
# 	my %genes;
# 	my %transcripts;
# 	foreach my $nm (keys %{$t2g}) {
# 	  $transcripts{$nm} = 1;
# 	  foreach my $gene (keys %{$t2g->{$nm}}) {
# 	    $genes{$gene} = 1;
# 	  }
# 	}
}

sub import_junction {
  my ($self, %options) = @_;
  my $row = $options{"-row"} || die "-row";
  my $hash = $options{"-hash"} || die "-hash";
  my $count = $row->{count} || die;
  my @j = split /,/, $row->{junction} || die;
  die unless @j == 2;
  my ($from, $to) = @j;
  my ($start_chr, $start_base) = parse_junction($j[0]);
  my ($end_chr, $end_base) = parse_junction($j[1]);
  $hash->{$start_chr}{$start_base}{$end_base}{count} += $count;
}

sub bounds_check {
  my ($end, $reference_len, $warn) = @_;
  my $result = 1;
  if ($end > $reference_len) {
    if ($warn) {
      # if attempting to anchor vs. a single start or end only,
      # we may inadvertently attempt to shift out of bounds.
      # No harm done in this case (vs. the much more serious
      # issue of an invalid junction in the BAM data)
      $result = 0;
    } else {
      die "FATAL: junction end $end beyond ref len of $reference_len";
    }
  }
  return $result;
}

sub pool_flag {
  # STATIC
  my ($entry, $entry_target, $field) = @_;
  my $combined = ($entry->{$field} or $entry_target->{$field}) ? 1 : 0;
#  printf STDERR "%s from:%d target:%d final:%d\n", $field, $entry->{$field}, $entry_target->{$field}, $combined;
  $entry_target->{$field} = $combined;
}

sub pool_number {
  # STATIC
  my ($entry, $entry_target, $field) = @_;
  my $combined = ($entry->{$field} || 0) + ($entry_target->{$field} || 0);
  # entry won't exist if newly-created
#  printf STDERR "%s from:%d target:%d final:%d\n", $field, $entry->{$field}, $entry_target->{$field}, $combined;
  $entry_target->{$field} = $combined;
}

sub pool_info {
  my ($self, $entry, $entry_target) = @_;

  $entry_target->{count} += $entry->{count};
#	  if (exists $entry_target->{qc_flanking}) {
  # FAILS for newly-created
  if ($self->has_flanking_qc() or exists $entry->{qc_flanking}) {
    # pool_flag($entry, $entry_target, "qc_flanking");
    # pool_flag($entry, $entry_target, "qc_plus");
    # pool_flag($entry, $entry_target, "qc_minus");
    # pool_number($entry, $entry_target, "qc_perfect_reads");
    # pool_number($entry, $entry_target, "qc_clean_reads");
    foreach my $f (QC_FIELDS_NUMBER) {
      pool_number($entry, $entry_target, $f);
    }
    foreach my $f (QC_FIELDS_FLAG) {
      pool_flag($entry, $entry_target, $f);
    }
  }
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/               
