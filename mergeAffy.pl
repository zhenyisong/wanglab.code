#!/usr/bin/perl

use strict;
use FileHandle;


my $file_dir          = '/home/zhenyisong/data/wanglilab/affy_results/';
#my $file_dir          = '/home/zhenyisong/data/wanglilab/temp/';
#my $file_dir          = 'D:\\wangli_data\\read_results\\';
chdir($file_dir);

my $output_filename   = 'final_affy.cos';
my $output_fh         = FileHandle->new(">$output_filename");

my @affy_files     = `ls *.txt`;
my $affy_results   = undef;

# read all input of the affy data

foreach my $file (@affy_files) {
    $affy_results->{$file} = readAffyfile($file);    
}

# median gene epxrs value;
foreach my $file (@affy_files) { 
    foreach my $gene (keys %{$affy_results->{$file}->{'genelist'}}) {
        my $value_exprs = $affy_results->{$file}->{'genelist'}->{$gene};
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
        $affy_results->{$file}->{'genelist'}->{$gene} = $result;   
    } 
}


# merge the name of affy_results data
my $initial_file      = $affy_files[-1];
my $initial_result    = $affy_results->{$initial_file};
my @intersection      = keys %{$initial_result->{'genelist'}};

foreach my $file (@affy_files[0..$#affy_files - 1]) {
    my $result               = $affy_results->{$file};;
    my @genenames            = keys %{$result->{'genelist'}};
    my %intersection         = map{$_ =>1} @intersection;
    @intersection            = grep( $intersection{$_}, @genenames);
}

# output the file header
my $last_file  = $affy_files[-1];
foreach my $file (@affy_files[0..$#affy_files - 1]) {
    my $header     = $affy_results->{$file}->{'header_@@'};
    my @buffer     = @{$header};
    foreach my $item (@buffer) {
        $output_fh->print($item);
        $output_fh->print("\t");     
    }  
}
my $header     = $affy_results->{$last_file}->{'header_@@'};
my @buffer     = @{$header};
my $last_id    = pop @buffer;
foreach my $srr_id (@buffer) {
    $output_fh->print($srr_id);
    $output_fh->print("\t");
}
$output_fh->print($last_id);
$output_fh->print("\n");

#print @intersection;
# output the gene rows

foreach my $gene (@intersection) {
    #print $gene,"\n";
    $output_fh->print($gene);
    $output_fh->print("\t");

    my $last_file = $affy_files[-1];
    foreach my $file (@affy_files[0..$#affy_files - 1]) {
        my $exprs_record = $affy_results->{$file}->{'genelist'}->{$gene};
        foreach my $value (@{$exprs_record}) {
           $output_fh->print($value); 
           $output_fh->print("\t");
        }
    }

    my $exprs_record = $affy_results->{$last_file}->{'genelist'}->{$gene};
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


sub readAffyfile {

   my $file       = shift;
   my $fh         = FileHandle->new($file);
   my $result     = undef;
   my $count      = 0;
   my $file_type  = 0;
   my $column_num = 0;
   if($file =~ /^pd/) {
       $file_type = 1;
   }
   while(my $line = $fh->getline) {
       $count++;
       chomp $line;
       if($count == 1) {
           my @headers = split(/\t/,$line);
           foreach my $i (1..$#headers) {
               $headers[$i] =~ s/^\.(\W+)//g;
           }
           $result->{'header_@@'} = [@headers[1..$#headers]];
           $column_num = scalar @headers;
           $column_num = $column_num - 1;
       }
       else {
           my @buffer   = split(/\t/,$line);
           my $genename = uc $buffer[0];
           if($file_type == 1) {
               my @string = split/\/\//,$genename;
               if(scalar @buffer > 1) {
                   $genename = $string[1];
                   $genename =~ s/^\.(\W+)//g;
                   $genename =~ s/\W+//g;
                   #print $genename,"____________","\n";
                   #print $genename,"\n";
                   foreach my $i (-$column_num..-1) {
                       $buffer[$i] =~ s/^\.(\W+)//g;
                   }
                   push @{$result->{'genelist'}->{$genename}} , [@buffer[-$column_num..-1]];
               }
           }
           else {
               $genename =~ s/^\.(\W+)//g;
               

               foreach my $i (1..$#buffer) {
                   $buffer[$i] =~ s/^\.(\W+)//g;
               }
               push @{$result->{'genelist'}->{$genename}} , [@buffer[1..$#buffer]];
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