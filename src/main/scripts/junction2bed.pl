#!/usr/bin/env perl
# convert tab-delimited annotated junction file to .bed format
# (w/various optional filters)
#
# MNE 4/2013

use strict;
use warnings;

use Getopt::Long;

use DelimitedFile;
use JunctionUtils qw(parse_junction);
use AtomicOutfile;
use Counter;
use CommandLineRebuilder;
use Cluster;
use MiscUtils qw(dump_die);

my $RGB_KNOWN = "0,128,0";
#my $RGB_KNOWN_LQ = "0,80,0";
#my $RGB_KNOWN_LQ = "45,90,45";
my $RGB_KNOWN_LQ = "128,255,128";
my $RGB_NOVEL = "192,0,0";

my %FLAGS;

my $REFGENE_SHADE = 1;
my @params = (
	      "-jf=s",
	      # junction file
	      "-glob-cross",

	      "-min-novel-reads=i",
	      # minimum required reads to include a novel junction

	      "-refgene-only",
	      # only include known junctions for refGene isoforms

	      "-gene=s",
	      # restrict to a gene

	      "-refgene-shade=i" => \$REFGENE_SHADE,
	      # for known junctions, use a brighter color for refGene
	      # and a duller color for others

	      "-devel-path",

	      "-now",
	    );

GetOptions(\%FLAGS,@params);

my $clr = new CommandLineRebuilder("-parameters" => \@params, "-flags" => \%FLAGS);
$clr->exclude_parameter("-glob-cross");
$clr->include_parameter("-now");

my @jf;
if ($FLAGS{jf}) {
  @jf = $FLAGS{jf};
} elsif ($FLAGS{"glob-cross"}) {
  @jf = glob("*.cross_sample_corrected.tab");
} else {
  die "specify -jf | -glob-cross\n";
}
die "no junction files" unless @jf;

my $min_novel_reads = $FLAGS{"min-novel-reads"} || 0;
my $refgene_only = $FLAGS{"refgene-only"};

my $GENE_RESTRICT;
my %gene_restrict;

if (my $gene = $FLAGS{gene}) {
  $gene_restrict{$gene} = 1;
  $GENE_RESTRICT = 1;
}

my $c = new Counter(\@jf);
foreach my $jf (@jf) {
  $c->next($jf);

  my $df = new DelimitedFile("-file" => $jf,
			     "-headers" => 1,
			     );

  my $ao = new AtomicOutfile("-base" => $jf, "-suffix" => "bed");
  $ao->add("refgene_only") if $refgene_only;
  $ao->add("refgene_shade") if $REFGENE_SHADE;
  $ao->add("min_novel", $min_novel_reads) if $min_novel_reads;
  $ao->add("gene", sort keys %gene_restrict) if $GENE_RESTRICT;
#  my $outfile = $jf . ".bed";
  my $outfile = $ao->get_outfile();

  unless ($FLAGS{now}) {
    #
    #  submit job to cluster
    #
    my $cmd = $clr->get_command_line("-jf" => $jf);
    my $c = new Cluster();
#    $c->app("perl-5.10.1");
    $c->node_class("");
    $c->memory_reserve_mb(1000);
    $c->memory_limit_mb(1000);
    $c->outfile($outfile);
    $c->project("PCGP");
    $c->command($cmd);
    $c->run();

    next;
  }


  my $wf = new WorkingFile($outfile);
  my $fh = $wf->output_filehandle();

  my $track_name = $jf;
  $track_name =~ s/\.bam.*//;

  printf $fh "track name=%s description=\"%s\" visibility=3 itemRgb=\"On\"\n", $track_name, $track_name;

  while (my $row = $df->get_hash()) {
    if (0) {
      foreach (sort keys %{$row}) {
	printf "%s: %s\n", $_, $row->{$_};
      }
    }

    my $j = parse_junction($row->{junction} || die);
    die unless $j->[0]->{ref} eq $j->[1]->{ref};

    my $ref = $j->[0]->{ref};
    my $start = $j->[0]->{base};
    my $end = $j->[1]->{base};

    my $count = $row->{count} || dump_die($row, "no count");
    my $name = $row->{bed_name} || $count;

    my $type = $row->{type} || die;
    die unless $type eq "known" or $type eq "novel";

    next if $min_novel_reads and
      $type eq "novel" and
	$count < $min_novel_reads;

    my $is_refgene;
    if ($type eq "known") {
      my @t = split /,/, $row->{transcripts} || die "no transcripts";
      $is_refgene = (grep {/N[MR]_/} @t) ? 1 : 0;
      next if $refgene_only and not($is_refgene);
    }

    if ($GENE_RESTRICT) {
      if ($row->{genes}) {
	my @g = split /,/, $row->{genes} || die "no genes";
	next unless grep {$gene_restrict{$_}} @g;
      } else {
	next;
      }
    }

    my @f;

    my $ref_clean = $ref;
    $ref_clean = "chr" . $ref unless $ref =~ /^chr/;
    # .bed requires chrX format

    push @f, $ref_clean;
    # 1. reference/chrom name
    push @f, $start - 1;
    # 2. start base number (0-based)
    push @f, $end;
    # 3. end base number: base is not counted, so don't need to adjust
    # :/
    push @f, $name;
    # 4. count (or custom)
    push @f, $count;
    # 5. score (count)
    push @f, "+";
    # 6. strand
    push @f, $f[1];
    push @f, $f[2];
    # 7. thickStart
    # 8. thickEnd

    my $color;
    if ($type eq "known") {
      if ($REFGENE_SHADE) {
	$color = $is_refgene ? $RGB_KNOWN : $RGB_KNOWN_LQ;
      } else {
	$color = $RGB_KNOWN;
      }
    } else {
      $color = $RGB_NOVEL;
    }

    push @f, $color;
    # 9. itemRGB
    printf $fh "%s\n", join "\t", @f;

  }
  $wf->finish();


}
