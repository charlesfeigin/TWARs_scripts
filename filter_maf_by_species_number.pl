#!/usr/bin/perl
use strict;
use warnings;
use Bio::Perl;

#----------------------------prep----------------------#

my $pause;

my ($in_name, $out_name, $min) = @ARGV;

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

		# Check the number of species in the subalignment
		my $check = CheckSpeciesNumber($min, @subaln);

		# If it passes the check
		if ($check eq 'true') {
		
			# Write it to the output file
			print OUT @subaln;

			print OUT "\n\n";

		}

		undef @subaln;
		$entry_flag = 0;

	}

}





###########################################################################
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

























