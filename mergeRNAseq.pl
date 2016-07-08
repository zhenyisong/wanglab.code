#!/usr/bin/perl

use strict;
use FileHandle;

# move all RNA-seq parsed results to
# the current directory
# cp -r SRP*/*.counts.txt RNA_seq_results/

my $file_dir          = '/home/zhenyisong/data/wanglilab/RNA_seq_results/';
#my $file_dir          = '/home/zhenyisong/data/wanglilab/temp/';
chdir($file_dir);
my $output_filename   = 'final_rna_seq.cos';
my $output_fh         = FileHandle->new(">$output_filename");
my @rna_results       = `ls *.counts.txt`;
my $exprs_results     = undef;
my $initial_file      = pop @rna_results;
my $initial_result    = readRNAseqOutput($initial_file);
$exprs_results->{$initial_file} = $initial_result;
my @intersection                = keys %{$initial_result->{'genelist'}};
foreach my $file (@rna_results) {
    my $result               = readRNAseqOutput($file);
    $exprs_results->{$file}  = $result;
    my @genenames            = keys %{$result->{'genelist'}};
    my %intersection         = map{$_ =>1} @intersection;
    @intersection            = grep( $intersection{$_}, @genenames);
}


# print the file header

my @file_keys  = keys %{$exprs_results};
my $last_file  = $file_keys[-1];

foreach my $file (@file_keys[0..$#file_keys - 1]) {
    my $header = $exprs_results->{$file}->{'header_@@'};
    foreach my $srr_id (@{$header}) {
        $output_fh->print($srr_id);
        $output_fh->print("\t");
    }
}
my $header  = $exprs_results->{$last_file}->{'header_@@'};
my @buffer  = @{$header};
my $last_id = pop @buffer;
foreach my $srr_id (@buffer) {
    $output_fh->print($srr_id);
    $output_fh->print("\t");
}
$output_fh->print($last_id);
$output_fh->print("\n");


# print out the gene expression rows

foreach my $gene (@intersection) {

    $output_fh->print($gene);
    $output_fh->print("\t");

    my $last_file = $file_keys[-1];
    foreach my $file (@file_keys[0..$#file_keys - 1]) {
        my $exprs_record = $exprs_results->{$file}->{'genelist'}->{$gene};
        foreach my $value (@{$exprs_record}) {
           $output_fh->print($value); 
           $output_fh->print("\t");
        }
    }

    my $exprs_record = $exprs_results->{$last_file}->{'genelist'}->{$gene};
    my @buffer       = @{$exprs_record};
    my $last_value   = pop @buffer;
    foreach my $value (@buffer) {
       $output_fh->print($value); 
       $output_fh->print("\t"); 
    }
    $output_fh->print($last_value); 
    $output_fh->print("\n");
    
}

$output_fh->close;



# the END
sub readRNAseqOutput {
    my $file        = shift;
    my $fh          = FileHandle->new($file);
    my $result_hash = undef;
    my $count       = 0;
    while(my $line = $fh->getline) {
        $count++;
        chomp $line;
        if($count == 1) {
            my @headers = ($line =~ /(SRR\d+)/g);
            $result_hash->{'header_@@'} = [@headers];
        }
        else {
           my @buffer   = split(/\t/,$line);
           my $genename = uc $buffer[0];
           #print @buffer;
           $result_hash->{'genelist'}->{$genename} = [@buffer[1..$#buffer]];
            
        }
        
    }
    $fh->close;
    return $result_hash;
}