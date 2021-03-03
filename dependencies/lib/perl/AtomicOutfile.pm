package AtomicOutfile;

use strict;
use Configurable;

@AtomicOutfile::ISA = qw(Configurable Exporter);
@AtomicOutfile::EXPORT_OK = qw();

use MethodMaker qw(
	base
        suffix
elements
auto_increment
counter
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->counter({});
  $self->reset();
  $self->configure(%options);
  return $self;
}

sub add {
  my ($self, @stuff) = @_;
  push @{$self->elements}, @stuff;
}

sub get_outfile {
  my ($self, %options) = @_;

  my $elements_raw = $self->elements();
  my $counter = $self->counter;
  
  my $base = join("_", $self->base(), @{$self->elements});

  if ($self->auto_increment()) {
    my $saw_count = $counter->{$base} || 0;
    $counter->{$base} = ++$saw_count;
    $base .= "_" . $saw_count if $saw_count > 1;
  }

  my $file = sprintf '%s.%s', $base, $self->suffix();

  return $file;
}

sub reset {
  my ($self) = @_;
  $self->elements([]);
}


1;
