#!/usr/bin/perl
use strict;
use warnings;
use Bio::Perl;

#----------------------------prep----------------------#

my ($in_name, $out_name, $refspec, $min) = @ARGV;

# Open input and output files
open (IN, "$in_name");
open (OUT, ">$out_name");

# Get the first line (maf header)
my $first_line = <IN>;

# Print it to the output file
print OUT $first_line;

#---------------------------parse-----------------------#

# While parsing the lines of the input maf file...
while (my $line = <IN>) {

	# If the line is anything other than a sequence line
	if ($line !~ /^s /) {

		# Write it to the output file
		print OUT "$line";

	# If the line matches the reference species
	} elsif ($line =~ /$refspec/g) {

		# Write it to the output file
		print OUT "$line";

	# Otherwise, if its a sequence line that isn't the reference
	} else {

		# Check the percent of sequence composed of nucleotides
		my $check = CheckNuclComp($min, $line);

		# If it passes the check
		if ($check eq 'true') {

			# Write it to the output file
			print OUT "$line";

		}		
			
	}

}


###########################################################################
sub CheckNuclComp {

	# Collect subaln and min
	my $min = shift @_;
	my $line = shift @_;

	# Initialize check value to false
	my $check = 'false';

	# Split the line and get the sequence
	my @split = split(" ", $line);

	# Collect the sequence
	my $seq = $split[6];

	# Strip any whitespace characters
	$seq =~ s/\s+//g;
	chomp ($seq);

	# Split the sequence into characters
	@split = split('', $seq);

	# Get the total length
	my $total_len = scalar(@split);

	# Count for number of nucleotides
	my $num_nucs = 0;

	# Count nucleotides
	foreach (@split) {

		if (($_ =~ /a/i) 
                or ($_ =~ /t/i) 
		or ($_ =~ /g/i) 
		or ($_ =~ /c/i)) {

			$num_nucs++;

		}

	}

	# Get percent of length that's nucleotides
	my $percent_nucs = ($num_nucs/$total_len) * 100;

	# If its greater than or equal to the minimum
	if ($percent_nucs >= $min) {

		# It passes the check
		$check = 'true';

	}

	# Return the result of the check
	return ($check)

}
