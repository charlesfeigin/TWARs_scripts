#!/usr/bin/perl

use strict;
use warnings;
use Bio::AlignIO;
use Bio::SimpleAlign;

my $inFileName = $ARGV[0];

my $input = Bio::AlignIO->new(-file => "$inFileName", -format => "fasta");

my $aln = $input->next_aln();

my $len = $aln->length;

my $pwd = $aln->percentage_identity;

my $weighted = $pwd / $len;


print "$pwd\t$len\n";
