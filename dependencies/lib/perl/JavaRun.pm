package JavaRun;

use strict;
use Configurable;

use Carp qw(confess);

@JavaRun::ISA = qw(Configurable Exporter);
@JavaRun::EXPORT_OK = qw();

use MethodMaker qw(

		   jvm
		   ram
		   classname

		   classpath
                   jar
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
#  $self->jvm("/h1/edmonsom/local/Linux_x86_64/jre1.6.0_15/bin/java");
#  $self->jvm("/usr/bin/env java");
  $self->ram(512);
  $self->jvm("/usr/bin/env java");
#  $self->classpath("/home/medmonso/lib/bambino.jar");
  $self->configure(%options);
  return $self;
}

sub av_development {
  confess "av_development: obsolete method!";
}

sub av_manual {
  confess "av_manual: obsolete method!";
}

sub blob_jar_mode {
  confess "blob_jar_mode: obsolete method!";
}

sub run {
  my ($self, %options) = @_;
  my $cmdline = $options{"-command"} || die "-command";

  my $cmd = $self->jvm();

  if (my $ram = $self->ram()) {
    if ($ram =~ /^\d+$/) {
      # no units specified
      $cmd .= sprintf ' -Xmx%dm', $ram;
    } else {
      $cmd .= sprintf ' -Xmx%s', $ram;
    }
  }

  if ($self->jar()) {
    # .jar mode: no classname needed (contains default)
    $cmd .= sprintf " -jar %s", $self->jar;
  } else {
    $cmd .= sprintf " -cp %s", $self->classpath if $self->classpath;
    $cmd .= sprintf " %s", $self->classname || die "-classname";
  }
  $cmd .= " $cmdline";

  if (my $cmdfile = $options{"-save-command-file"}) {
    # save command line to text file
    open(CMDSAVE, ">$cmdfile") || die "can't write to $cmdfile";
    printf CMDSAVE "%s\n", $cmd;
    close CMDSAVE;
  }

  if (my $outfile = $options{"-save-out"}) {
    $cmd .= " > " . $outfile;
  }

  if (my $errfile = $options{"-save-err"}) {
    $cmd .= " > " . $errfile;
  }

  my $result;
  if ($options{"-capture"}) {
    print STDERR "$cmd\n";
    $result = `$cmd`;
    chomp $result;
  } elsif ($options{"-return"}) {
    # just return command line
    $result = $cmd;
  } elsif ($options{"-debug"}) {
    printf STDERR "NOT running: %s\n", $cmd;
    $result = -1;
  } elsif ($options{"-pipe"}) {
    $result = new FileHandle();
    $result->open($cmd . "|");
  } else {
    print STDERR "running: $cmd\n";
    system $cmd;
    $result = $?;
  }
  return $result;
  
}

1;
