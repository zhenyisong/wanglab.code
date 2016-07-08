#!/usr/bin/perl

use strict;
use FileHandle;

#my $file_dir          = '/home/zhenyisong/data/wanglilab/non_standard_results/';
my $file_dir          = '/home/zhenyisong/data/wanglilab/temp/';
chdir($file_dir);

my $output_filename   = 'final_nons.cos';
my $output_fh         = FileHandle->new(">$output_filename");

my @nons_files        = `ls *discard`;
my $nons_results      = undef;

# to cacluate the median of the exprs values

foreach my $file (@nons_files) {
    $nons_results->{$file} = readNonSfile($file); 
    foreach my $gene (keys %{$nons_results->{$file}->{'genelist'}}) {
        my $value_exprs = $nons_results->{$file}->{'genelist'}->{$gene};
        my $result      = undef;
        my $first       = $value_exprs->[0];
        my $array_num   = scalar @{$value_exprs};
        my $exprs_num   = scalar @{$first};

        foreach my $i (0..$exprs_num - 1) {
            my @buffer = ();
            foreach my $j (0..$array_num - 1) {
                push @buffer, $value_exprs->[$j]->[$i];
            }
            my $temp = median(@buffer);
            push @{$result},$temp;
        }

        $nons_results->{$file}->{'genelist'}->{$gene} = $result;   
    } 
}

#print @nons_files;
#print @{$nons_results->{"GSE76720_series_matrix.txt\n"}->{'genelist'}->{'ZCCHC11'}},"\n";

#exit(0);

# merge the name of nons_results data
my $initial_file      = $nons_files[-1];
my $initial_result    = $nons_results->{$initial_file};
my @intersection      = keys %{$initial_result->{'genelist'}};
foreach my $file (@nons_files[0..$#nons_files - 1]) {
    my $result               = $nons_results->{$file};;
    my @genenames            = keys %{$result->{'genelist'}};
    my %intersection         = map{$_ =>1} @intersection;
    @intersection            = grep( $intersection{$_}, @genenames);
    print $file,"\t",scalar @genenames,"\t",scalar @intersection,"\n";
}

# output the file header
my $last_file  = $nons_files[-1];
foreach my $file (@nons_files[0..$#nons_files - 1]) {
    my $header     = $nons_results->{$file}->{'header_@@'};
    my @buffer     = @{$header};
    foreach my $item (@buffer) {
        $output_fh->print($item);
        $output_fh->print("\t");     
    }  
}
my $header     = $nons_results->{$last_file}->{'header_@@'};
my @buffer     = @{$header};
my $last_id    = pop @buffer;
foreach my $srr_id (@buffer) {
    $output_fh->print($srr_id);
    $output_fh->print("\t");
}
$output_fh->print($last_id);
$output_fh->print("\n");


# output the gene rows

#print @nons_files;

foreach my $gene (@intersection) {

    $output_fh->print($gene);
    $output_fh->print("\t");
    my $count = 0;

    my $last_file = $nons_files[-1];
    foreach my $file (@nons_files[0..$#nons_files - 1]) {
        my $exprs_record = $nons_results->{$file}->{'genelist'}->{$gene};
        foreach my $value (@{$exprs_record}) {

           $output_fh->print($value); 
           $output_fh->print("\t");
        }
    }

    my $exprs_record = $nons_results->{$last_file}->{'genelist'}->{$gene};
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


sub readNonSfile {
   my $file   = shift;
   my $fh     = FileHandle->new($file);
   my $result = undef;
   my $count  = 0;
   my $header_len = 0;
   while(my $line = $fh->getline) {
       $count++;
       chomp $line;
       
       if($count == 1) {
           my @headers = split(/\t/,$line);
           $result->{'header_@@'} = [@headers[1..$#headers]];
           $header_len = scalar @headers - 1;
       }
       else {
           my @buffer   = split(/\t/,$line);
           my $genename = uc $buffer[0];
           $genename =~ s/^\.(\W+)//g;
           $genename =~ s/\W+//g;
           #print $genename,"__________________\n";
           my @array = ();
           foreach my $i (1..$#buffer) {
               $buffer[$i] =~ s/^\.(\W+)//g;
           }
           if($header_len != $#buffer) {
               my $mean = 0;
               foreach my $i (@buffer[1..$#buffer]) {
                   $mean += $i;
               }
               $mean = $mean/$#buffer;
               foreach my $i (1..($header_len - $#buffer)) {
                   push @array, $mean;
               }
               push @{$result->{'genelist'}->{$genename}}, [@buffer[1..$#buffer],@array];
           } 
           else {
               push @{$result->{'genelist'}->{$genename}}, [@buffer[1..$#buffer]]; 
           }
       }
   }
   $fh->close;
   return $result;
}

#
# this code is adapted from 
# http://www.perlmonks.org/index.pl?node_id=474564
#
sub median {
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}