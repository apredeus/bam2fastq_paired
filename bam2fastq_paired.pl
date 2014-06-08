#!/usr/bin/perl 

use strict;
use warnings;
use IO::Handle;

sub which_read($); 
sub dec2bin($); 
sub read_line;  

## this is a special script dedicated to splitting BAM file 
## while the problem is embarrasingly easy for single-end alignment, it's quite interesting for paired-end. 
## problem of course arises when dealing with singleton alignments (one is aligned and one is not).
## what we want to achieve is the following:
## 1. for the properly paired reads, we want to correctly find which read is R1 and which is R2 using the flags. 
## 2. for the discordant reads, we want to find if this read is R1 or R2. 
## 3. for case 2, we want to make up a fake complementary read that would not be mapped under any sircumstances. 
## This way we should obtain two FASTQ files that would pass any formal paired-end QC (i.e. each line has matching names @name123456/1 and @name123456/2,etc) 
## we assume to start with a 1st COLUMN SORTED SAM FILE WITH NO HEADER. 

my $sam = shift @ARGV; 
my $tag = shift @ARGV; 


## the whole program is based on comparing the current read to the previous one. 
## this can become clunky, but saves a lot of time processing things - you only 
## have few reads in the memory and read the huge SAM file line by line. 

my $prev_name = "";
my $curr_name = ""; 
my $last_read = ""; 

my $R1_prev_name = ""; 
my $R2_prev_name = ""; 
my $R1_prev_seq  = ""; 
my $R2_prev_seq  = ""; 
my $R1_prev_qual = ""; 
my $R2_prev_qual = ""; 

my $R1_curr_name = ""; 
my $R2_curr_name = ""; 
my $R1_curr_seq  = ""; 
my $R2_curr_seq  = ""; 
my $R1_curr_qual = ""; 
my $R2_curr_qual = ""; 

die "Please provide a 1st column-sorted SAM file with no header and a common tag for output FASTQ names!\n" if ($sam eq "" || $tag eq "" ); 

open (SAM,"<",$sam) or die "$!"; 
open (R1,">","$tag.R1.fastq") or die "$!"; 
open (R2,">","$tag.R2.fastq") or die "$!"; 

my $SAM_io = IO::Handle->new_from_fd(fileno(SAM),"r");

## read the first line of SAM file, generate lines "subseq" and "subqual"
## subseq is the line of NNNNNNNNNNNNNN you put for a non-existing mate 
## subqual is the same for qualities (lowest possible quality is #, so it's ########)

read_line; 
my $subseq; 
my $subqual; 
my $count = 1; 

if ($last_read eq "R1") { 
  $subseq   = $R1_curr_seq; 
  $subqual  = $R1_curr_qual;
} elsif ($last_read eq "R2") {
  $subseq   = $R2_curr_seq; 
  $subqual  = $R2_curr_qual; 
}

$subseq  =~ tr/[A-Za-z]/N/; 
$subqual =~ s/./#/g;

$prev_name = $curr_name; 
if ($last_read eq "R1") { 
  $R1_prev_name  = $R1_curr_name; 
  $R1_prev_seq   = $R1_curr_seq; 
  $R1_prev_qual  = $R1_curr_qual; 
  $R2_prev_name  = $R2_curr_name; 
  $R2_prev_seq   = $subseq; 
  $R2_prev_qual  = $subqual; 
} elsif ($last_read eq "R2") {
  $R2_prev_name  = $R2_curr_name; 
  $R2_prev_seq   = $R2_curr_seq; 
  $R2_prev_qual  = $R2_curr_qual; 
  $R1_prev_name  = $R1_curr_name; 
  $R1_prev_seq   = $subseq; 
  $R1_prev_qual  = $subqual; 
}
    
## notice that we basically start from line 2. 
B

while ( read_line ) { 
  $count++; 
  if ($curr_name ne $prev_name) {
    printf R1 "%s\n%s\n+\n%s\n",$R1_prev_name,$R1_prev_seq,$R1_prev_qual;
    printf R2 "%s\n%s\n+\n%s\n",$R2_prev_name,$R2_prev_seq,$R2_prev_qual;

    ## now clear all writeable fields in "prev":
    $prev_name = $curr_name; 
    if ($last_read eq "R1") { 
      $R1_prev_name  = $R1_curr_name; 
      $R1_prev_seq   = $R1_curr_seq; 
      $R1_prev_qual  = $R1_curr_qual; 
      $R2_prev_name  = $R2_curr_name; 
      $R2_prev_seq   = $subseq; 
      $R2_prev_qual  = $subqual; 
    } elsif ($last_read eq "R2") {
      $R2_prev_name  = $R2_curr_name; 
      $R2_prev_seq   = $R2_curr_seq; 
      $R2_prev_qual  = $R2_curr_qual; 
      $R1_prev_name  = $R1_curr_name; 
      $R1_prev_seq   = $subseq; 
      $R1_prev_qual  = $subqual; 
    }
  } else {
    if ($last_read eq "R1") { 
      $R1_prev_seq   = $R1_curr_seq; 
      $R1_prev_qual  = $R1_curr_qual; 
    } elsif ($last_read eq "R2") {
      $R2_prev_seq   = $R2_curr_seq; 
      $R2_prev_qual  = $R2_curr_qual;
    }  
  } 
} 

## finally, we need to print that one read that is left. 

printf R1 "%s\n%s\n+\n%s\n",$R1_prev_name,$R1_prev_seq,$R1_prev_qual;
printf R2 "%s\n%s\n+\n%s\n",$R2_prev_name,$R2_prev_seq,$R2_prev_qual;

close SAM; 
close R1; 
close R2; 


## function that gives "R1" or "R2" identifier to a read depending on the flag (second column in SAM file)
sub which_read ($) { 
  my $flag = shift; 
  my $read = ""; 
  if (substr(dec2bin($flag),-7,1)) { 
    $read = "R1"; 
  } elsif (substr(dec2bin($flag),-8,1)) {
    $read = "R2"; 
  }
  return $read; 
} 

## function that decodes decimal flag into binary representation, but returns actual STRING. 
sub dec2bin ($) {
  my $str = unpack("B32", pack("N", shift));
  $str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
  return $str;
}

## funciton that reads SAM file and populates "curr" variables
sub read_line {
  my $sam_line = $SAM_io->getline();
  return undef unless defined $sam_line; # End of file?
  my @tt = split (/\t+/,$sam_line); 
  $curr_name = $tt[0]; 
  $R1_curr_name = $tt[0]."/1";  
  $R2_curr_name = $tt[0]."/2";  
  if (which_read($tt[1]) eq "R1") {
    $R1_curr_seq   = $tt[9];  
    $R1_curr_qual  = $tt[10];  
    $last_read = "R1"; 
  } elsif (which_read($tt[1]) eq "R2") {
    $R2_curr_seq   = $tt[9];  
    $R2_curr_qual  = $tt[10];  
    $last_read = "R2"; 
  }  
}
  
