#!/usr/bin/perl

# test
# how-to-use
# perl RNASeq_QC.pl hg19 temp.bam.txt 1 <option: result.txt>
# parse the parameters from the perl standin

# @author  Yisong Zhen
# @since   2017-02-22
# @update  2017-04-05
# @parameter
#     @spieces:
#            mm10, hg38, hg19
#     @bam_file_name
#     @strand_specificity
#     @outputfileName
#          use default, random generated unique file name instead.

# this script needs these packages as the surppotives
# @picard version: 2.9.0
#     http://broadinstitute.github.io/picard/
# @RSeQC version : RSeQC v2.6.4
#    http://rseqc.sourceforge.net/
# @perl v5.16.3
#     and dependant module
# @java java -version > 1.8.x
#
# picard.script mm10 *.bam none false zhen34.temp

# @how to import the script generated data into R
# read.csv(file= outputFileName_with_path, header = FALSE, sep="\t");
#
#     
# see more help on picard
# http://broadinstitute.github.io/picard/command-line-overview.html
#
# how to generate the refflat file
# for mm10
# please download the corresponding file from here
# http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/
# http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/refFlat.txt.gz
# http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
# then I renamed to refFlat_mm10.txt
# mv refFlat.txt refFlat_mm10.txt
# mv refFlat.txt refFlat_hg38.txt

# how to generate mm10_ribosome_interval_list.txt"
# see: https://www.biostars.org/p/120145/
# see:http://seqanswers.com/forums/showthread.php?p=136425
# You can find the intervals using the UCSC Table browser. 
# For this, you go to 
# http://genome.ucsc.edu/cgi-bin/hgTables
# and then set group:all tables, table:rmsk, 
# and filter to "repClass (does match) rRNA" 
# then output it as a GTF file.
# 

# CollectRnaSeqMetrics: Options choice
# Please see the CollectRnaSeqMetrics definitions for a complete description of the metrics produced by this tool.
# Option	Description
# @REF_FLAT (File)	Gene annotations in refFlat form. 
# Format described here: http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat Required.
# @RIBOSOMAL_INTERVALS (File)	
# Location of rRNA sequences in genome, in interval_list format. 
# If not specified no bases will be identified as being ribosomal. 
# Format described here: Default value: null.
# @STRAND_SPECIFICITY (StrandSpecificity)	
# For strand-specific library prep. 
# For unpaired reads, use FIRST_READ_TRANSCRIPTION_STRAND 
# if the reads are expected to be on the transcription strand. Required. 
# Possible values: {NONE, FIRST_READ_TRANSCRIPTION_STRAND, SECOND_READ_TRANSCRIPTION_STRAND}
# @MINIMUM_LENGTH (Integer)	When calculating coverage based values (e.g. CV of coverage) 
# only use transcripts of this length or greater. Default value: 500. 
# This option can be set to 'null' to clear the default value.
# @CHART_OUTPUT (File)	The PDF file to write out a plot of normalized position vs. coverage. Default value: null.
# @IGNORE_SEQUENCE (String)	If a read maps to a sequence specified with this option, all the 
# bases in the read are counted as ignored bases. These reads are not counted as Default value: null. 
# This option may be specified 0 or more times.
# RRNA_FRAGMENT_PERCENTAGE (Double)	
# This percentage of the length of a fragment must overlap one of the ribosomal 
# intervals for a read or read pair to be considered rRNA. 
# Default value: 0.8. This option can be set to 'null' to clear the default value.
# METRIC_ACCUMULATION_LEVEL (MetricAccumulationLevel)	
# The level(s) at which to accumulate metrics. Default value: [ALL_READS]. 
# This option can be set to 'null' to clear the default value. 
# Possible values: {ALL_READS, SAMPLE, LIBRARY, READ_GROUP} 
# This option may be specified 0 or more times. 
# This option can be set to 'null' to clear the default list.
# @INPUT (File)	Input SAM or BAM file. Required.
# @OUTPUT (File)	File to write the output to. Required.
# @ASSUME_SORTED (Boolean)	If true (default), then the sort order in the header file will be ignored. 
# Default value: true. This option can be set to 'null' to clear the default value. Possible values: {true, false}
# @STOP_AFTER (Long)	Stop after processing N reads, mainly for debugging. 
# Default value: 0. This option can be set to 'null' to clear the default value.
#---end
# I mkdir the picardlib in the /bioware directory


# RSeQC parameters
# I downloaded all the required data from Liguo's website
# these data is saved in the wanglab computer
# /bioware/RSeQClib
#---
use Modern::Perl;
use autodie;
use FileHandle;
use File::Temp qw(tempfile);
use Switch;
use File::Basename;
use Cwd;


# http://stackoverflow.com/questions/84932/how-do-i-get-the-full-path-to-a-perl-script-that-is-executing
# this method is the better way.
#---
# my $WORKING_DIR = dirname(__FILE__);
# the above discarded! 
# 
my $WORKING_DIR = getcwd();
chdir $WORKING_DIR;

=head
test linux shell script to down-sample the original bam file

count=0
bamFile='S256_08B_CHG010228-HX1-1_L005_R1.fastq.bam'
FLOAT=0.01
while [[ $count -le 3 ]]; do
   echo $count
   count=$((count+1))
   
   seed=$(echo "$count + $FLOAT" | bc -l )
   samtools view -s $seed -b $bamFile > bamfile_${count}_downSample.bam
done
ls *downSample.bam > temp.bam.txt
echo 'finished, please start debuging next step'

the above script generate the 10 random bam files for the debuging purpose;
now you can try 
perl RNASeq_QC.pl hg19 temp.bam.txt 1
=cut

#
# define the required parameters for PICARD
#
my $JAR           = '/bioware/picardlib/picard.jar';
my $REFFLAT_mm10  = '/bioware/picardlib/refFlat_mm10.txt';
my $REFFLAT_hg19  = '/bioware/picardlib/refFlat_hg19.txt';
my $REFFLAT_hg38  = '/bioware/picardlib/refFlat_hg38.txt';
my $RIBO_INTERVAL_LIST_mm10 = '/bioware/picardlib/mm10_ribosome_interval_list.txt';
my $RIBO_INTERVAL_LIST_hg38 = '/bioware/picardlib/hg38_ribosome_interval_list.txt';
my $RIBO_INTERVAL_LIST_hg19 = '/bioware/picardlib/hg19_ribosome_interval_list.txt';

# define the required parameters for RSeQC
# you may change to use CONSTANT t odefine these outside parameters
#---
my $mm10_houseKeeping    = '/bioware/RSeQClib/mm10.HouseKeepingGene.bed';
my $mm9_houseKeeping     = '/bioware/RSeQClib/mm9.HouseKeepingGene.bed';
my $mm10_RefSeq          = '/bioware/RSeQClib/mm10_RefSeq.bed';
my $mm9_RefSeq           = '/bioware/RSeQClib/mm9_RefSeq.bed';
my $hg19_houseKeeping    = '/bioware/RSeQClib/hg19.HouseKeepingGenes.bed';
my $hg38_houseKeeping    = '/bioware/RSeQClib/hg38.HouseKeepingGenes.bed';
my $hg19_RefSeq          = '/bioware/RSeQClib/hg19_RefSeq.bed';
my $hg38_RefSeq          = '/bioware/RSeQClib/hg38_RefSeq.bed';


my $REFFLAT            = undef;
my $RIBO_INTERVAL_LIST = undef;
my $HOUSE_KEEPING      = undef;
my $REF_SEQ            = undef;


my $spieces        = $ARGV[0];
my $bamFilePath    = $ARGV[1];
my $stranded       = $ARGV[2];
my $outputFileName = $ARGV[3];
$outputFileName    ||= get_RNAQC_filename();

# this is the unique-name check, if exist, then throw-out error!

if( -e $outputFileName) {
    print "unlukily, this file already exited. Please choose another file name to save your result\n";
    exit 6;
}

# this step is essential as the outputFile use ">>' as the directed output pipeline
# without this step will result in mixed result;
unlink "$outputFileName" if -e $outputFileName;

# According to the above table, exit codes 1 - 2, 126 - 165, and 255 [1] 
# have special meanings, and should therefore be avoided for user-specified 
# exit parameters. Ending a script with exit 127 would certainly cause confusion 
# when troubleshooting (is the error code a "command not found" or a user-defined one?). 
# However, many scripts use an exit 1 as a general bailout-upon-error. 
# Since exit code 1 signifies so many possible errors, it is not particularly useful in debugging.

switch ($spieces) {
		case 'mm10'    { $REFFLAT = $REFFLAT_mm10; 
                         $RIBO_INTERVAL_LIST = $RIBO_INTERVAL_LIST_mm10;
                         $HOUSE_KEEPING = $mm10_houseKeeping;
                         $REF_SEQ = $mm10_RefSeq;}
		case 'hg19'    { $REFFLAT = $REFFLAT_hg19; 
                         $RIBO_INTERVAL_LIST = $RIBO_INTERVAL_LIST_hg19;
                         $HOUSE_KEEPING = $hg19_houseKeeping;
                         $REF_SEQ = $hg19_RefSeq; }
		case 'hg38'    { $REFFLAT = $REFFLAT_hg38; 
                         $RIBO_INTERVAL_LIST = $RIBO_INTERVAL_LIST_hg38;
                         $HOUSE_KEEPING = $hg38_houseKeeping;
                         $REF_SEQ = $hg38_RefSeq;}
        case '?'       { helpMeNow(); exit 0;}
        case '-help'   { helpMeNow(); exit 0;}
        case '--help'  { helpMeNow(); exit 0;}
		else		   { print "input parameter error!\n\n"; helpMeNow(); exit 3;}
}                      

switch ($stranded) {
		case 1  { $stranded = 'NONE';}
		case 2  { $stranded = 'FIRST_READ_TRANSCRIPTION_STRAND'; }
		case 3  { $stranded = 'SECOND_READ_TRANSCRIPTION_STRAND';}
		else    { print "input parameter error!\n"; exit 3;}
}


my @allBamFiles = getAllBamFiles($bamFilePath);

while(my $singleBam = shift @allBamFiles) {
    my ($bamID)        = ($singleBam =~ /([0-9A-Za-z\-\_]+)\..*?$/);
    my $tempFileName   = get_temp_filename();
    my $tempFileName1  = $tempFileName .'.bam';
    my $tempFileName2  = $tempFileName .'.txt';
    my $tempFileName3  = $tempFileName .'.strand';
    my $tempFileName4  = $tempFileName .'.read';
    my $sortBamInput   = `samtools sort -m 1000000000 $singleBam $tempFileName`;
    my $indexBamInput  = `samtools index $tempFileName1`;
    my $samHeader      = `samtools view -H $tempFileName1 > $tempFileName2`;
    my $RIBO_LIST      = `cut -s -f 1,4,5,7,9  $RIBO_INTERVAL_LIST >> $tempFileName2`;
    my $PICARD         = `java -Xmx2g -jar $JAR CollectRnaSeqMetrics REF_FLAT=$REFFLAT RIBOSOMAL_INTERVALS=$tempFileName2 STRAND_SPECIFICITY=$stranded CHART_OUTPUT=null METRIC_ACCUMULATION_LEVEL=ALL_READS INPUT=$tempFileName1 OUTPUT=$tempFileName2  ASSUME_SORTED=true 2> /dev/null`;
    my $picardResult   = getPicardMetrics($tempFileName2);
    my $outputString   = "$bamID\tPCT_ROBOSOME\t$picardResult";
    my $writeResult    = `echo "$outputString" >> $outputFileName`;
    my $geneBody       = `geneBody_coverage.py -i $tempFileName1 -r $REF_SEQ -o $tempFileName 2> /dev/null`;
    my $tempFileName5  = $tempFileName.'.geneBodyCoverage.txt';
    $outputString      = readGeneBody($tempFileName5);
    $outputString      = "$bamID\tgeneBodyCoverage\t$outputString";
    $writeResult       = `echo "$outputString"  >>  $outputFileName`;
    my $exonNumber     = `read_distribution.py -i $tempFileName1 -r $REF_SEQ > $tempFileName4 2> /dev/null`;
    $outputString      = readsIntronExonRatios($tempFileName4);
    $outputString      = "$bamID\tintronExonRatio\t$outputString";
    $writeResult       = `echo "$outputString"  >>  $outputFileName`;
    my $strandSpec     = `infer_experiment.py -r $REF_SEQ -i $tempFileName1 > $tempFileName3 2> /dev/null`;
    my @stranded       = readStrandedInfo($tempFileName3);
    $outputString      = "$bamID\tStrandedInfo\t$stranded[0]\t$stranded[1]\t$stranded[2]";
    $writeResult       = `echo "$outputString" >>  $outputFileName`;
    my $GCbiase        = `read_GC.py -i $tempFileName1 -o $tempFileName 2> /dev/null`;
    $tempFileName5     = $tempFileName.'.GC.xls';
    $outputString      = readGCbiase($tempFileName5);
    $outputString      = "$bamID\tGCcontent\t$outputString";
    $writeResult       = `echo "$outputString" >>  $outputFileName`;
    # remove the rubbish
    unlink <$tempFileName*>;
    
}


# @para
#    filename: the input is the file name which saves the 
#              bam file pathes
# @return
#   the array reference which contains the bam file path
#---

sub getAllBamFiles {
    my $bamFilePath  = shift;
    my $fh           = FileHandle->new("$bamFilePath");
    my @fileArray    = ();
    while(my $line = $fh->getline()) {
        chomp $line;
        if($line =~ /\.bam$|\.BAM$/) {
            push @fileArray, $line;
        }
    }
    $fh->close;
    return (@fileArray);
}

# @parameters
#    no
# @return
#    the random file name with a unique name
#---

sub get_temp_filename {
    my $filename;
    (undef, $filename) = tempfile('tempXXXXXX', SUFFIX   => '.txt', OPEN => 0);
    return $filename;
}

# @parameters
#    no
# @return
#   the output file (.txt) name for the QC results
#----

sub get_RNAQC_filename {

    my $filename;
    (undef, $filename) = tempfile('RNAQCXXXXXX', SUFFIX   => '.txt', OPEN => 0);
    return $filename;
}

# @parameters
#    the file name from the picard analysis program
# @return
# the specified column, in this case, is the PCT_RibosomeRNA
#----

sub getPicardMetrics {
    my $fileName = shift;
    my $fh       = FileHandle->new($fileName);
    my @metricsHeader = ();
    my @metricsNumber = ();
    while(my $line = $fh->getline()) {
        chomp $line;
        if( $line =~ /^PF_BASES/) {
            
            @metricsHeader = split $line,"\t";
            $line = $fh->getline();
            chomp $line;
            @metricsNumber = ($line =~ /\S+/g);
            last;
        }
    }
    #print @metricsNumber,"\n";
    my $result = $metricsNumber[15];
    $fh->close;
    return $result;
    
}


# deprecated
# not implemented
#---
sub readGCbiase {
    my $fileName  = shift;
    my $fh        = FileHandle->new("$fileName");
    my $result    = 0;
    my @array     = ();
    while( my $line = $fh->getline()) {
        chomp $line;
        next if $line =~ /GC%/;
        my ($percent, $count) = ($line =~ /\S+/g);
        #print $percent,"\t",$count,"\n";
        push(@array, [$percent, $count]);        
    }
    $fh->close;
    my $totalReads = 0;
    foreach my $count (@array) {
        $totalReads += $count->[1];
    }
    foreach my $gc (@array) {
        #print $gc->[0],"\t",$gc->[1],"\t$totalReads\n";
        $result += ($gc->[0] * $gc->[1])/$totalReads;
    }
    return $result;
}

sub readStrandedInfo {
    my $fileName = shift;
    my $fh        = FileHandle->new("$fileName");
    my $FR        = 0.000001;
    my $RF        = 0.000001;
    my $result    = 0;
    while(my $line = $fh->getline()) {
        chomp $line;
        if($line =~ /1\+\+,1\-\-,2\+\-,2\-\+\":\s+(\S+)$/ || $line =~ /\+\+,\-\-\":\s+(S+)$/) {
             $RF = $1;
        }
        elsif($line =~ /1\+\-,1\-\+,2\+\+,2\-\-\":\s+(\S+)$/ || $line =~ /\+\-,\-\+\":\s+(S+)$/) {
             $FR = $1;
        }
        else {
            next;
        }
    }
    $fh->close;
    return ($RF/$FR, $RF, $FR);
}

# @param
#    the output file name from 

sub readsIntronExonRatios {
    my $fileName = shift;
    my $fh       = FileHandle->new($fileName);
    my $exons    = 0;
    my $introns  = 0;
    while(my $line = $fh->getline()) {
        if($line =~ /5\'UTR_Exons\W+(\d+)\W+(\d+)/) {
            $exons += $2;
        } elsif($line =~ /CDS_Exons\W+(\d+)\W+(\d+)/) {
            $exons += $2;
        } elsif($line =~ /3\'UTR_Exons\W+(\d+)\W+(\d+)/) {
            $exons += $2;
        } elsif($line =~ /Introns\W+(\d+)\W+(\d+)/) {
            $introns = $2;
        } else {
            next;
        }
        
    }
    $fh->close;
    return $introns/$exons; 
}

sub readGeneBody {
    my $fileName   = shift;
    my $fh         = FileHandle->new($fileName);
    my @percentile = ();
    my @reads      = ();
    while(my $line = $fh->getline()) {
        chomp $line;
        if($line =~ /Percentile/) {
            @percentile = ($line =~ /\d+/g);
            next;
        }
        @reads = ($line =~ /\S+/g);        
    }
    $fh->close;
    my $sumPrefix = 0;
    map { $sumPrefix += $_ } @reads[11..30];
    my $sumMiddle = 0;
    map { $sumMiddle += $_ } @reads[31..70];
    my $sumSuffix = 0;
    map { $sumSuffix += $_ } @reads[71..90];
    my $rounded1 = sprintf "%.2f", $sumMiddle/$sumPrefix;
    my $rounded2 = sprintf "%.2f", $sumSuffix/$sumPrefix;
    my $result  = '1' .':' . $rounded1 .':' . $rounded2;
    return $result;
}

sub helpMeNow {
   print "usage: RNASeq_QC.pl <spiecesName> <bamFileSavedInTxt> <Stranded> <outputFileName>\n";
   print "\n\n";
   print "predefined parameter for\n\n";
   print "<spiecesName>: hg19, hg38, mm9, mm10 (case-sensitive)";
   print "<bamFileSavedInTxt: the raw txt file (name) contains the paths to bamFiles\n";
   print "<Stranded>: 1 , 2, 3 which stands for (No stranded), FIRST_READ_TRANSCRIPTION_STRAND, SECOND_READ_TRANSCRIPTION_STRAND\n";
   print "<outputFileName>: unique file name or leave this blank, the script will generate a unique name with prefix RNAQCXXXX\n";
   print " \nThe script will have exit status: 0, , 3, 4, 6. 0 stands for OK!\n";
   print " \nComments please send to zhenyisong\@gmail.com\n";
}

if($? == 0) {
    print "we complete processing the data\n!";
    exit $outputFileName;
} else {
    print "we lost battle in the final step, please re-run the script!\n";
    exit 4;
}



