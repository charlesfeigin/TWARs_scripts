#!/usr/bin/perl

=pod

This script will accept a bed-formatted annotation file (e.g. phastCons conserved elements intersected with ChIP elements),
determine their distances from each other going from left to right, and merge elements which
are less than a pre-specified number of bases apart. For example, if annotated p300/H3K27a
peaks have poor overall coverage by conserved elements, but are tiled by many small, relatively 
closely spaced conserved elements, merging adjacent elements which are avg bases apart
could improve the percent coverage of annotations by conserved elements. This may make sense in 
cases where the average distance is quite high but median is quite small (i.e. many small
elements are closely clustered together, but clusters themselves are much farther apart). 

This doesn't need to be done across the whole genome, only the intersection of ChIP and phastCons
data. Thus, merging can be done between conserved regions within the bounds of a larger
ChIP peak.

=cut


use strict;
use warnings;

my $pause;

my ($bed, $merged, $merge_dist) = @ARGV;

open(BED, $bed);
open(MERGED, ">$merged");

# Handle first element to initialize everything
my $first_line = <BED>;

my @first = split(' ', $first_line);

my %seen;


my $previous_chr = $first[0];
my $previous_start = $first[1];
my $previous_end = $first[2];	
my $previous_element = $first[3];
$seen{$previous_element} = 1;

while (my $line = <BED>) {

	my @split = split(' ', $line);
 
	# If we're still working on the same chromosome as last iteration...
	if ($split[0] eq $previous_chr) {
	
		# Calculate the distance from the last element
		my $dist = $split[1] - (1 + $previous_end);

		# If the distance is less than or equal to the maximum...
		if ($dist <= $merge_dist) {

			# Then replace the previous element's end point 
			# with that of the current element
			$previous_end = $split[2];

		# If the distance is too large...
		} else {

			# Complete the previous entry and write it to the output file
			my @entry = ($previous_chr, $previous_start, $previous_end, $previous_element);			
			my $entry = join("\t", @entry);
			chomp($entry);
			print MERGED "$entry\n";

			# Update "previous" variables to reflect the current segment of the 
			# current element (the first two should be redundant)
			$previous_element = $split[3];
			$previous_chr = $split[0];
			$previous_start = $split[1];
			$previous_end = $split[2];

		}

	# If the current chromosome doesn't match the previous...
	} else {

		# Complete the previous entry as before
		my @entry = ($previous_chr, $previous_start, $previous_end, $previous_element);			
		my $entry = join("\t", @entry);
		chomp($entry);
		print MERGED "$entry\n";

		# Update ALL element varaibles
		$previous_element = $split[3];
		$previous_chr = $split[0];
		$previous_start = $split[1];
		$previous_end = $split[2];

	}
}


# Complete the last entry
my @entry = ($previous_chr, $previous_start, $previous_end, $previous_element);			
my $entry = join("\t", @entry);
chomp($entry);
print MERGED "$entry\n";


exit;
