#!/usr/bin/env perl
# add annotations to a junction file from a GTF file (e.g. GenCode)
#
# MNE 2/2015

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use File::Basename;
use DelimitedFile;
use FloatingJunctionDetector qw(parse_junction safe_get);
use RefFlatFile;

use MiscUtils qw(dump_die build_argv_list);
# use DelimitedFile;
# use Reporter;

my $F_TRANSCRIPT = "annot_transcript";
my $F_TRANSCRIPT_PRIMARY = "annot_transcript_primary";
my $PROGRESS_TRANSCRIPTS = 5000;

printf STDERR "*** WARNING: use only as a last resort, it may be much better to use the desired isoforms during junction correction instead!\n";
# important safety tip

my %FLAGS;
GetOptions(\%FLAGS,
	   "-file=s",
	   "-files=s",
	   # junction file(s)

	   "-gtf=s",
	   # GTF to annotate against, e.g.
	   # /nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/mRNA/GENCODE/gencode_junctions.txt

	   "-refgene=s",
	   # refgene to annotate against, e.g.
	   # /nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/mRNA/GENCODE/wgEncodeGencodeCompV19.txt

	   "-verbose",
	   "-progress",

    );

my $VERBOSE = $FLAGS{verbose};
my $files = build_argv_list("-flags" => \%FLAGS, "-single" => "file", "-set" => "files");

my $DB;
if ($FLAGS{gtf}) {
  $DB = parse_gtf();
} elsif ($FLAGS{refgene}) {
  $DB = parse_refgene();
} else {
  die "specify -refgene [file] or -gtf [file]\n";
}


foreach my $infile (@{$files}) {
  my $outfile = basename($infile) . ".annotated.tab";

  my $df = new DelimitedFile(
			     "-file" => $infile,
			     "-headers" => 1,
			     );


  my @extra = ( $F_TRANSCRIPT );
  push @extra, $F_TRANSCRIPT_PRIMARY if $FLAGS{gtf};
  # only available in GTF data

  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => \@extra
			     );

  while (my $row = $df->get_hash()) {
    my $j = $row->{junction} || die;
    my @j = split /,/, $j || die;
    die unless @j == 2;
    my ($from, $to) = @j;
    my ($start_chr, $start_base) = parse_junction($j[0]);
    my ($end_chr, $end_base) = parse_junction($j[1]);
    my $c = cook_chrom_name($start_chr);

    my $hits = safe_get($DB, $c, $start_base, $end_base);
    my $transcript = "";
    my $transcript_primary = "";
    if ($hits) {
      my @t;
      my @p;
#      printf STDERR "hits: %d\n", scalar @{$hits};
      foreach my $h (@{$hits}) {
	push @t, $h->{transcript};
	push @p, $h->{transcript} if $h->{primary};
      }
      $transcript = join ",", @t;
      $transcript_primary = join ",", @p;
    }
    $row->{$F_TRANSCRIPT} = $transcript;
    $row->{$F_TRANSCRIPT_PRIMARY} = $transcript_primary;

    $rpt->end_row($row);
  }

  $rpt->finish();

}


sub parse_gtf {
  #
  # would try Bio::Tools::GFF, but in a hurry, & simple enough format
  #
  # no doubt a standard parser, but in a hurry
  #
  # http://genome.ucsc.edu/FAQ/FAQformat.html#format4
  my $gtf = $FLAGS{gtf} || die "-gtf";
  open(IN, $gtf) || die;

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

  my @f;

  my %f_ignore = map {$_ => 1} qw(
				   gene
				   CDS
				   start_codon
				   stop_codon
				   UTR
				   Selenocysteine
				);

  my @exons;
  my %junctions;
  my $transcripts = 0;
  my $progress = $FLAGS{progress};

  while (<IN>) {
    chomp;
    @f = split /\t/, $_;
    die unless @f == 9;
    my %row;
    @row{@fields} = @f;

    my $feature = $row{feature};
    next if $f_ignore{$feature};

    # will break if attribute names repeat
#    my @pairs = split /;\s*/, $row{group};
#    die join "\n", @pairs;

#    printf "%s\n", $row{group};
    my %attr = map {split /\s+/, $_} split /;\s*/, $row{group};
    $row{attributes} = \%attr;

    if ($feature eq "transcript") {
      flush_transcript("-exons" => \@exons, "-junctions" => \%junctions);
      @exons = ();
      # reset
      if ($progress and ++$transcripts % $PROGRESS_TRANSCRIPTS == 0) {
	my $bytes = tell(IN);
	printf STDERR "  progress: %d transcripts, %.2f%%...\n",
	  $transcripts,
	    ($bytes * 100 /-s $gtf);
      }
    } elsif ($feature eq "exon") {
      push @exons, \%row;
    } else {
      die "unhandled feature $feature";
    }
  }
  flush_transcript("-exons" => \@exons, "-junctions" => \%junctions);
  close IN;

  return \%junctions;
}

sub flush_transcript {
  my %options = @_;
  my $exons = $options{"-exons"} || die "-exons";
  my $junctions = $options{"-junctions"} || die "-junctions";

  if (@{$exons}) {
    # something to process

    my @sorted = sort {$a->{start} <=> $b->{start}} @{$exons};

#    foreach my $e (@sorted) {
#      dump_die($e, "debug" ,1);
#    }

    my $last = @sorted - 2;

    my $attr = $sorted[0]->{attributes};
    my $transcript = $attr->{transcript_id} || die;
    my $ref_name = cook_chrom_name($sorted[0]->{seqname} || die);

    my %transcript_info;
    $transcript_info{transcript} = $transcript;
    my $tname = $attr->{transcript_name} || die;

    my $primary = $tname =~ /\-001$/ ? 1 : 0;
    $transcript_info{primary} = $primary;

    for (my $i = 0; $i <= $last; $i++) {
      my $js = $sorted[$i]->{end} || die;
      my $je = $sorted[$i + 1] ->{start} || die;

      push @{$junctions->{$ref_name}{$js}{$je}}, \%transcript_info;

      printf STDERR "%s %s:%d:+,%s:%d:+\n",
	$transcript,
	  $ref_name, $js,
	    $ref_name, $je if $VERBOSE;
    }
  }
}

sub cook_chrom_name {
  my ($raw) = @_;
  $raw =~ s/^chr//i;
  return $raw;
}

sub parse_refgene {
  my $fn = 
  my $rf = new RefFlatFile();
  $rf->parse_file(
		  "-refflat" => ($FLAGS{refgene} || die),
		  "-type" => "refgene",
		 );

  my %junctions;

  foreach my $r (@{$rf->rows()}) {
    my $ref_name = cook_chrom_name($r->{chrom} || die);
    my $transcript = $r->{name} || die;

    my %transcript_info;
    $transcript_info{transcript} = $transcript;

    foreach my $j (@{$r->{junctions}}) {
      my $js = $j->{start} || die;
      my $je = $j->{end} || die;

      push @{$junctions{$ref_name}{$js}{$je}}, \%transcript_info;

      printf STDERR "%s %s:%d:+,%s:%d:+\n",
	$transcript,
	  $ref_name, $js,
	    $ref_name, $je if $VERBOSE;
    }
  }

  return \%junctions;
}
