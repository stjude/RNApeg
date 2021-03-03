package ChrBucketMap;
# bucket data by chromosome and (rough) position

use strict;
use Exporter;

use Configurable;
use BucketMap;
use GenomeUtils qw(cook_chromosome_name);
use MiscUtils qw(dump_die);

@ChrBucketMap::ISA = qw(Configurable Exporter);
@ChrBucketMap::EXPORT_OK = qw();

use constant DEFAULT_CHUNK_SIZE => 100000;

use MethodMaker qw(
		    chunk_size
		    f_chr
		    f_start
		    f_end
		    maps
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->maps({});
  $self->chunk_size(DEFAULT_CHUNK_SIZE);
  $self->configure(%options);
  return $self;
}

sub add_row {
  my ($self, %options) = @_;
  my $row = $options{"-row"} || die;
  my $maps = $self->maps || die;

  my $f_chr = $self->f_chr || die;
  my $f_start = $self->f_start;
  my $f_end = $self->f_end;

  my $chr = cook_chromosome_name(
				 ($row->{$f_chr} || dump_die($row, "no $f_chr")),
#				 "-genome" => $self->genome,
				 "-return-unknown" => 1) || die;

  my $map = $maps->{$chr};
  unless ($map) {
    $map = $maps->{$chr} = new BucketMap("-chunk" => $self->chunk_size);
  }

  my $start = $options{"-start"} || $row->{$f_start};
  my $end = $options{"-end"} || $row->{$f_end};

  die unless defined($start) and defined($end);


  $map->add_range(
		  "-start" => $start,
		  "-end" => $end,
		  "-value" => $row
		 );
}

sub find {
  my ($self, %options) = @_;
  my $chr_raw = $options{"-chr"} || die;
  my $start = $options{"-start"} || die;
  my $end =  $options{"-end"} || $start;

  my $chr = cook_chromosome_name(
				 $chr_raw,
#				 "-genome" => $self->genome,
				 "-return-unknown" => 1) || die;

  my $map = $self->maps->{$chr};
  my $results;
  if ($map) {
    my @hits;
    if ($start == $end) {
      push @hits, @{$map->fuzzy_find("-site" => $start)};
    } else {
      die unless $end > $start;
      push @hits, @{$map->fuzzy_find("-start" => $start, "-end" => $end)};
    }

    my $f_start = $self->f_start;
    my $f_end = $self->f_end;
    foreach my $h (@hits) {
      next if $h->{$f_end} < $start;
      next if $h->{$f_start} > $end;

      $results = [] unless $results;
      push @{$results}, $h;
    }
#    printf STDERR "raw hits:%d final:%d\n", scalar(@hits), $results ? scalar(@{$results}) : 0;

  }
  return $results;
}


1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
