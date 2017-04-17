#!/usr/bin/perl

# parse the parameter from the perl standin
use strict;
use FileHandle;
use File::Temp;

my $spieces        = $ARGV[0]
my $bamFilePath    = $ARGV[1]
my $tempFileName   = get_temp_filename()

sub get_temp_filename {
    my $fh = File::Temp->new(
             TEMPLATE => 'tempXXXXX',
             DIR      => 'mydir',
             SUFFIX   => '.dat',
    );
    return $fh->filename;
}

my $filename = get_temp_filename();



