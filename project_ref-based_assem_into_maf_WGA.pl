#!/usr/bin/perl
use strict;
use warnings;
use Bio::Perl;

#----------------------------prep----------------------#

my ($in_name, $out_name, $refspec, $refassem, $suffix, $min, $prefix) = @ARGV;


=pod

This script takes a whole-genome maf alignment, a reference-based genome assembly of a species
not present in the maf alignment but referenced against a species which is present, the id of
the reference species. This script: 
1) Reads through the alignment
2) Finds subalignments where the reference species is present
3) Extracts metadata from the reference (i.e. chr/scaff number, start, end, strand etc)
4) Uses that information to extract the equivalent region from the reference based assembly

Several assumptions are made:
1) The coordinate systems of the reference genome and the reference-based assembly are identical
2) The scaffold IDs are the same

This script contains some hard-coding (esp. regarding expectations of what input genome's scaffold
ID's will look like) and is only suitable in its current form for my own data sets. Sorry! Modules
for parsing maf alignments are pretty garbage and I'm not a developer. For an upcoming project I 
will need to do something similar and will write a more generally applicable version in python 
using the same logic but without any hard-coding.

# Collect arguments
# in_name - input maf alignment
# out_name - output maf alignment
# refspec - ID of the species that acted as the reference genome for the input
#            reference-based assembly
# refassem - filename of the reference based assembly
# suffix - suffix for labelling
# minimum number of species
# prefix for labelling

=cut


# Make names for temp files
my $temp1 = "$prefix" . "_temp1.fasta";
my $temp2 = "$prefix" . "_temp2.fasta";
my $temp3 = "$prefix" . "_temp3.fasta";

# Open input and output files
open (IN, "$in_name");
open (OUT, ">$out_name");

# Get the first line (maf header)
my $first_line = <IN>;

# Print it to the output file
print OUT $first_line;

# Make a flag to determine whether the script is currently in a proper maf entry (subalignment)
# and an array to hold the current subalignment, and a counter for the number of species
my $entry_flag = 0;
my @subaln;


#---------------------------parse-----------------------#

# While parsing the lines of the input maf file...
while (my $line = <IN>) {

	# If the line is a score, indicating the first line of a subalignment...
        if ($line =~ /^a score/) {

		# Flip the flag to say we're currently making an entry
		$entry_flag = 1;

		# Add the score as the first line of the subalignment
                push (@subaln, $line);

	# If the line begins with an "s" or "i" or "e"
        } elsif (($line =~ /^s /)) { #or ($line =~ /^i/)) { #or ($line =~ /^e/)) {

		# Add it to the current subalignment
		push (@subaln, $line);

	# If the line is empty AND we are in the process of making an entry
	# then process and close the entry
        } elsif (($line =~ /^\s*$/) and ($entry_flag == 1)) {

		# First, check the subalignment for the reference species
		# If present, collect its metadata and make a samtools faidx 
		# command to extract the equivalent region from the refassem
		my @metadata = GetRef($refspec, @subaln);

		# If the reference species was not found, then check if the subaln
		# has the minimum number of species to be retained. If it does,
		# print it to the file 
		if ($metadata[0] eq 'false') {
		
			# Check species subroutine
			my $check = CheckSpeciesNumber($min, @subaln);

			# If it passes the check
			if ($check eq 'true') {
		
				# Write it to the output file
				print OUT @subaln;

				print OUT "\n\n";

			}

		# If the reference species was present
		} else { 

			# Reformat the subaln as fasta and write to temp1.fasta
			MafToFasta($temp1, @subaln);

			# Shift the faidx command off of metadata onto its own variable
			my $faidx_cmd = shift (@metadata);

			# Run the samtools faidx command and write the resulting fasta
			# to temp2.fasta
			GetSeqsFaidx($temp2, $suffix, $faidx_cmd, @metadata);

			# Run Mafft command to align the sequence extracted from the refassem
			# to the rest of the subalignment
			my $mafft_err = "$prefix" . "_mafft_err.txt";
			my $mafft_command = "mafft-linsi --add $temp2 --keeplength $temp1 > $temp3 2>$mafft_err";
			system($mafft_command);

			# Now open the new mafft alignment, extract the added, aligned sequence
			# and reformat it as MAF
			my $reformatted = FastaSeqToMaf($temp3, $suffix, @metadata);

			# Add the reformatted line to the MAF subalignment
			push (@subaln, $reformatted);

			# Check that the total number of species meets the min requirement
			my $check = CheckSpeciesNumber($min, @subaln);

			# If it passes the check
			if ($check eq 'true') {
		
				# Write it to the output file
				print OUT @subaln;

				print OUT "\n\n";

			}

		}
			

		#unlink glob "temp*.fasta";
		undef @subaln;
		$entry_flag = 0;


	}

}

unlink glob "temp*.fasta";


################################################################
sub GetRef {

	# Get name of ref species
	my $refspec = shift @_;

	# Get subalignment
	my @subaln = @_;

	my @metadata;

	# Initialize a variable to hold the samtools faidx command
	my $faidx_cmd = 'false';

	# Parse subalignment
	while (my $line = shift @subaln) {

		# If a given sequence (s ) line matches the ref species...
		if ($line =~ /s $refspec/g) {

			# Split the line and collect its relevant components for later use
			my @line = split(" ", $line);
			my $src = $line[1];	# Name of the source sequence
			my $start = $line[2];	# Start position in the source sequence (0-based)
			my $size = $line[3];	# Size of the aligned sequence
			my $strand = $line[4];	# Strand of the aligned sequence (+/-)
			my $srcSize = $line[5];	# Size of the source sequence
			my $seq = $line[6];	# The aligned sequence itself including gaps
			# Collect elements 1 through 5 as metadata
			@metadata = @line[1..5];

			# If the sequence was on the plus strand...
			if ($strand eq '+') {

				# Calculate the end position like so
				my $end = $start + ($size); #- 1);

				# Shift the start position into 1-based coord system by adding 1
				$start = $start + 1;

				$src = ParseIDs($src);
				#print "$src\n";

				# Make the samtools faidx command
				# Note: Don't output fasta in the command, need to reverse complement 
				# it later if on the minus strand
				$faidx_cmd = "samtools-1.3.1 faidx $refassem " . "$src" . ":" . "$start" . "-" . "$end";

			# If the sequence was on the minus strand...
			} elsif ($strand eq '-') {

				# Calculate the end position like so
				my $end = $srcSize - $start;

				# Calculate the start position like so
				$start = $end - ($size - 1);

				$src = ParseIDs($src);

				# Make samtools faidx command (see note above)
				$faidx_cmd = "samtools-1.3.1 faidx $refassem " . "$src" . ":" . "$start" . "-" . "$end";
				
			}

		}

	}


	# Add the samtools faidx command to the beginning of the 
	# metadata array and return it
	unshift @metadata, $faidx_cmd;	
	#print "$faidx_cmd\n";
	return (@metadata);

}

sub CheckSpeciesNumber {

	# Collect subaln and min
	my $min = shift @_;
	my @subaln = @_;

	# Initialize check value to false
	my $check = 'false';

	# Determine if the number of species meets the minimum requirement
	# Initialize a counter at zero, add 1 for each sequence line
	my $count = 0;

	# If line begins with "s" add one to the counter
	foreach (@subaln) {

		if ($_ =~ /^s/) {
	
			$count++;

		}

	}	

	# If the count is greater than or equal to min species, print it to the output
	if ($count >= $min) {

		$check = 'true';

		return ($check);

	}		

}

sub MafToFasta {

	# Collect args
	my $temp1 = shift @_;
	my @subaln = @_;

	# Open temp file
	open (TEMP1, ">$temp1");

	# Parse sequence lines from maf subaln
	foreach (@subaln) {
	
		if ($_ =~ /^s /) {

			my @split = split(" ", $_);

			my $id = $split[1];

			my $seq = $split[6];

			print TEMP1 ">$id\n$seq\n";

		}

	}	

	close TEMP1;

}

sub GetSeqsFaidx {

	my $temp2 = shift @_;

	my $suffix = shift @_;

	my $faidx_cmd = shift @_;

	my @metadata = @_;

	my @fasta = `$faidx_cmd`;

	open (TEMP2, ">$temp2");
	
	my $id = shift @fasta;
	my $seq = join('', @fasta);

	$id =~ s/\s+//g;
	$seq =~ s/\s+//g;

	my $header = ">$suffix" . "_";

	$id =~ s/\>/$header/; 

	my @split = split(':', $id);
	$id = $split[0];

	if ($metadata[3] =~ /\+/) {


		print TEMP2 "$id\n$seq";
		close TEMP2;
		return

	} elsif ($metadata[3] =~ /\-/) {

		my $seq_obj = Bio::Seq->new(-id => "$id", 
					    -seq => "$seq",
					    -alphabet => 'dna');

		#my $revcom = revcom($fasta[1]);
		my $rev_obj = $seq_obj->revcom();
		
		$seq = $rev_obj->seq();

		print TEMP2 "$id\n$seq";
		close TEMP2;		
		return

	}

	print "Something is unstranded. Adjust accordingly\n@metadata";
	exit;

}

sub ParseIDs {

	my $long = $_[0];
	my $short;

	# Devil ids will be of the form: sarHar1.chr1_GL834894_random
	# Need to reformat to GL834894.1
	if ($long =~ /sarHar/g) {

		# First split at full stops leaving : sarHar1 and chr1_GL834894_random
		my @split = split(/\./, $long);

		# Get the second element and split at underscores leaving chr1 GL834894 random
		@split = split('_', $split[1]);

		# Get the second element again (GL834894) and add .1
		$short = "$split[1]" . ".1";

	# Dog ids will be of the form: canFam3.chr18
	# Need to reformat to chr18
	} elsif ($long =~ /canFam/g) {

		if ($long =~ /chr[X0-9]+/g) {

			my @split = split(/\./, $long);

			$short = $split[1];

			#($short) = $split[1] =~ m/chr[A-Za-z0-9]+/g; 

		} elsif ($long =~ /chrUn/g) {

			my @split = split(/_/, $long);

			$short = "$split[1]" . ".1";

		} else {

			print "Don't understand this header: $long\n";

		}
				 



	}

	return ($short);
	
}

sub FastaSeqToMaf {

	# Collect args
	my $temp3 = shift @_;
	my $suffix = shift @_;
	my @metadata = @_;
	my @maf;
	my $maf;
	my $id;
	my $seq;

	# Open mafft output file
	my $seqIO_obj = Bio::SeqIO->new(-file => "$temp3", -format => 'fasta', -alphabet => 'dna');

	# Find the line matching suffix (refassem sequence)
	while (my $seq_obj = $seqIO_obj->next_seq()) {

		$id = $seq_obj->id();

		if ($id =~ /$suffix/g) {

			$seq = $seq_obj->seq();

			last;	
	
		}

	}

	# Reformat the data into mafft format using the metadata from the ref species
	push (@maf, 's');
	push (@maf, @metadata);
	push (@maf, "$seq");

	$maf = "$maf[0] $id\t\t       $maf[2]  $maf[3] $maf[4]     $maf[5] $maf[6]";

	return ($maf);

}
















