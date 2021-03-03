package UniqueRank;
# assign a simple rank based on sorted input and unique values

use strict;
use Configurable;

@UniqueRank::ISA = qw(Configurable Exporter);
@UniqueRank::EXPORT_OK = qw();

use MethodMaker qw(
	last_value
counter_rank
counter_order

		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->counter_rank(0);
  $self->counter_order(0);
  $self->configure(%options);
  return $self;
}

sub get_rank {
  my ($self, $value) = @_;
  my $last_value = $self->last_value();
  my $counter = $self->counter_rank();
  $self->counter_order($self->counter_order + 1);
  if (!defined($last_value) or $value ne $last_value) {
    $self->counter_rank(++$counter);
  }
  $self->last_value($value);
  return $counter;
}

sub get_order {
  # must be called AFTER get_rank()!
  my ($self) = @_;
  return $self->counter_order;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/               
