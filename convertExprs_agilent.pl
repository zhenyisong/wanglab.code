#!/usr/bin/perl

use strict;
use FileHandle;

my $file_dir           = '/home/zhenyisong/data/wanglilab/GSE76720';
chdir($file_dir);
my $annotion_file_name = 'GPL13497-9755.txt';
my $matrix_file_name   = 'GSE76720_series_matrix.txt';
my $input_fh           = FileHandle->new($matrix_file_name);
my $annotation         = readAnnotationTable($annotion_file_name);
my $result_dir         = '/home/zhenyisong/data/wanglilab/non_standard_results';
my $output_file        = $matrix_file_name;
chdir($result_dir);
my $output_fh          = FileHandle->new(">$output_file");

while (my $line = $input_fh->getline) {
    chomp $line;
    if ($line =~ /^!/) {
        next;
    }
    $line =~ s/\"//g;
    if ($line =~ /^ID/) {
        $output_fh->print($line);
        $output_fh->print("\n");
    }
    my @buffer = split(/\t/,$line);
    if (exists $annotation->{$buffer[0]}) {
        $buffer[0] = $annotation->{$buffer[0]};
    }
    else {
        next;
    }
    my $last_element = pop(@buffer);
    foreach my $item (@buffer) {
       $output_fh->print($item);
       $output_fh->print("\t"); 
    }
    $output_fh->print($last_element);
    $output_fh->print("\n");
}

$output_fh->close;
$input_fh->close;



sub readAnnotationTable {
    my $file   = shift;
    my $fh     = FileHandle->new($file);
    my $result = undef;
    while(my $line = $fh->getline) {
        if ($line =~ /^A_/) {
            my @buffer = split(/\t/,$line);
            if ($buffer[6] =~ /\w+/) {
                $result->{$buffer[0]} = $buffer[6];
            }
        }
    }
    $fh->close;
    return $result;
}
