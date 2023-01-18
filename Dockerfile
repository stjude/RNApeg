FROM maven:3.6-openjdk-11 as builder

COPY pom.xml /tmp/build/pom.xml
COPY src /tmp/build/src
COPY dependencies /tmp/build/dependencies

RUN cd /tmp/build && mvn install

FROM alpine:latest
# build from Alpine for much smaller final image

# Set working dir as /RNApeg
WORKDIR /RNApeg

# Tell the OS that this is a non-interactive frontend only for build
ARG DEBIAN_FRONTEND=noninteractive

#
#  Alpine-based install:
#
RUN apk add --no-cache openjdk8-jre
# don't need full JDK, just runtime
RUN apk add --no-cache curl
RUN apk add --no-cache make
RUN apk add --no-cache gcc
RUN apk add --no-cache perl
RUN apk add --no-cache perl-utils
# for "cpan" command-line utility
RUN apk add --no-cache perl-dev
RUN apk add --no-cache musl-dev
# these two for headers required to build perl DBI module
RUN apk add --no-cache perl-doc
# required for perl "use diagnostics", used by SampleName.pm
RUN apk add --no-cache db-dev
RUN apk add --no-cache expat-dev
# required for dependencies of BioPerl
RUN apk add --no-cache openssl
RUN apk add --no-cache openssl-dev
RUN apk add --no-cache zlib
RUN apk add --no-cache zlib-dev
# required for WDL workflows to use RNApeg
RUN apk add --no-cache bash

# Install perl modules
RUN cpan App:cpanminus
RUN cpanm --no-wget Data::Compare && chown -R root:root /root/.cpanm
RUN cpanm --no-wget DB_File XML::Parser::PerlSAX XML::Twig XML::DOM
RUN cpanm --no-wget Bio::Tools::CodonTable
RUN cpanm --no-wget DBI
RUN cpanm --no-wget Set::IntSpan

# Put code into place
COPY src/main/scripts /RNApeg/src/bin
COPY src/main/java /RNApeg/src/javalib
COPY src/main/perllib /RNApeg/src/perllib
COPY dependencies/bin /RNApeg/src/bin
COPY dependencies/lib/java /RNApeg/src/javalib
COPY dependencies/lib/perl /RNApeg/src/perllib
COPY src/main/docker/RNApeg.sh /RNApeg/src/bin

COPY --from=builder /tmp/build/target /RNApeg/src/javalib

# Change environment variables
ENV PATH="/RNApeg/src/bin:${PATH}"
ENV PERL5LIB="/RNApeg/src/perllib:${PERL5LIB}"
ENV CLASSPATH=/RNApeg/src/javalib/*

WORKDIR /results
ENTRYPOINT ["/RNApeg/src/bin/RNApeg.sh"]
CMD ["-h"]
