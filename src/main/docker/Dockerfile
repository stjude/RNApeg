# Build off of ubuntu:16.04
#FROM ubuntu:16.04

FROM alpine:latest
# build from Alpine for much smaller final image

# Set working dir as /RNApeg
WORKDIR /RNApeg

# Tell the OS that this is a non-interactive frontend only for build
ARG DEBIAN_FRONTEND=noninteractive

#
#  Alpine-based install:
#
#RUN apk add --no-cache parallel
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

#
# Ubuntu-based install:
#
# Install some basic things
#RUN apt-get update && apt-get install --assume-yes apt-utils
#RUN apt-get -y upgrade
#RUN apt-get -y install parallel
# needed??
#RUN apt-get -y install default-jdk
#RUN apt-get -y install default-jre
# don't need full JDK, just runtime
#RUN apt-get -y install curl
#RUN apt-get -y install make
#RUN apt-get -y install gcc
# gcc needed for cpanm

# Install perl modules
RUN cpan App:cpanminus
# turn off temporarily for dnanexus debugging
#RUN cpanm --no-wget Data::Compare
RUN cpanm --no-wget Data::Compare && chown -R root:root /root/.cpanm
# - is this used by any of the main scripts, or just optional/supplemental?
# - this layer doesn't seem to work in DNAnexus for some reason

#RUN cpanm --no-wget DBI
# turn this off temporarily for dnanexus debugging

# Put code into place
COPY src/bin /RNApeg/src/bin
COPY src/javalib /RNApeg/src/javalib
COPY src/perllib /RNApeg/src/perllib
#COPY configs /RNApeg/configs

# hack for locally-built Perl modules:
#COPY src/cpanm.local /cpanm.local
#ENV PERL5LIB="/cpanm.local/lib/perl5:${PERL5LIB}"
#ENV PERL5LIB="/cpanm.local/lib/perl5/x86_64-linux-gnu-thread-multi/:${PERL5LIB}"
#
# FAIL -- this doesn't work anyway for whatever reason:
#
#  DB<1> use DBI
#DBI.c: loadable library and perl binaries are mismatched (got handshake key 0xdb00080, needed 0xde00080)
# 

# Change environment variables
ENV PATH="/RNApeg/src/bin:${PATH}"
ENV PERL5LIB="/RNApeg/src/perllib/common-perllib-genome:${PERL5LIB}"
ENV PERL5LIB="/RNApeg/src/perllib/rnapeg:${PERL5LIB}"
ENV PERL5LIB="/RNApeg/src/perllib/util-perllib:${PERL5LIB}"
ENV CLASSPATH=/RNApeg/src/javalib/*
