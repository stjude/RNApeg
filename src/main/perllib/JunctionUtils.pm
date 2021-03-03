package JunctionUtils;

use strict;
use Configurable;

@JunctionUtils::ISA = qw(Configurable Exporter);
@JunctionUtils::EXPORT_OK = qw(
parse_junction
parse_junction_entry
);

use MethodMaker qw(
	
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  return $self;
}

sub parse_junction {
  my ($j, %options) = @_;
  my @j = split /,/, $j;
  die "format error for $j" unless @j == 2;
  my @results = (
    parse_junction_entry($j[0], %options),
    parse_junction_entry($j[1], %options)
      );
  return wantarray ? (@results) : \@results;
}

sub parse_junction_entry {
  my ($j, %options) = @_;
  my @f = split /:/, $j;
  die unless @f == 3;
  my ($ref, $base, $strand) = @f;
  if ($options{"-ucsc"}) {
    $ref =~ s/^chr//i;
    # standardize to UCSC format
    $ref = "chr" . $ref;
  } elsif ($options{"-strip"}) {
    # strip leading chr if present
    $ref =~ s/^chr//i;
  }
  my %j;
  $j{ref} = $ref;
  $j{base} = $base;
  $j{strand} = $strand;
  return \%j;
}


1;
