#!/usr/bin/perl

use strict;
use warnings;
use List::Util 'sum';
use Data::Dumper;

=pod

Try to figure out if there is a smarter way to do my weighted averages. 
Maybe instead of giving greater weight to elements of greater length 
(as a proportion of the sum of peak lengths).

Notes before starting:

1) Both the cds and ChIP intersection (isect) files must be sorted in a specific way.
open in Microsoft Excel/LibreOffice Calc and sort the data in ascending order by the 
following columns.

cds gff:
	sort key 1: column 9 (gff attributes) -> unique id
	sort key 2: column 6 (gff score) -> Use script to set score to exon number before intersecting
	sort key 3: column 4 (gff start) -> coordinate for start of entry

ChIP bed
	sort key 1: column 4 (bed ID)
	sort key 2: column 2 (bed start)


=cut


# Pause for debugging
my $pause;
# Flags for gff/bed specific index numbering
my $gff_flag = 'gff';
my $bed_flag = 'bed';

# Collect arguments
my ($cds_isect, $chip_isect, $target_seg_num, $max_filled_perc) = @ARGV;

# Open isect gff files
open (CDS, $cds_isect);
open (CHIP, $chip_isect);

# Collect gff file contents. These will be used to collect initial data
# then subsequently edited with each merge cycle.
my @cds_gff = <CDS>;
my @chip_bed = <CHIP>;


#------------------------------Collect initial stats----------------------------#

# Get an initial (weighted) average number of segments per exon in the cds gff file
# This data will be updated progressively with each merge cycle until the target
# average is met, or the max percent of filled gaps (see below) is exceeded 
my $avg_seg_num_exons = GetAverageSegNum(@cds_gff);

# Set the initial average gap filled percentage in ChIP isect elements to zero
# since nothing's been filled in yet
my $avg_filled_perc_chip = 0;

=pod
Make a hash that holds the lengths of each isect element. This hash will be 
initialize as the sum of the lengths of all isect elements (isect elements being
the intersection of ChIP and conservation data) for each exon (as numbered in annotations).
Each time gaps are filled in the merging algorithm, the number of gap bases filled
will be added to the corresponding element in this hash. Thus, the percent gap filled for 
each exon will be calculated as the number of filled bases in that exon divided by the total
length of the exon (original isect length plus length of filled bases). These percents will
then be used to calculated a weighted average for the percent of filled bases over the total
isect dataset and compared to the maximum allowed.
=cut
my %total_chip_lengths = GetElementLengths(@chip_bed);

# Make another hash that will hold the lengths of only the filled bases in ChIP isect elements.
# These will be the values divided by the total lengths in the total_lengths hash above
my %filled_chip_lengths;

print "Initial average segment number: $avg_seg_num_exons\n";
print "Initial average percent filled gaps: 0\n\n";


#-----------------------------Begin Merging Algorithm--------------------------#

=pod
This algorithm will iteratively merge adjacent segments in the cds and gff isect datasets
increasing the allowed distance between merging elements by 1 each cycle until one of
two sets of conditions are met.

The first (ideal) condition is that the (weighted) average number of segments per cds exon
isect is less than or equal to the average target number, and the (weighted) average percent 
length  of ChIP isects composed of filled gaps is as close tot he max percentage as possible 
while still being less than or equal to it.

The second (not as good) condition is that the average percent filled gaps in ChIP elements
reaches the maximum before the target average segment number per exon isect has been met. 
This is an allowable outcome, but it may mean the user needs to adjust their expectations 
for either the target average number of segments per exon isect or the percent of ChIP elements
that may be composed of filled gaps

After each cycle that neither of these conditionls are met, the merge distance will be increased
by 1 and the results will be re-checked
=cut

# Set the initial merge distance to 1 base (nucleotide)
my $merge_distance = 0;

# Merging algorithm loop
#while (($avg_seg_num_exons> $target_seg_num) and ($avg_filled_perc_chip < $max_filled_perc)) {
while ($avg_filled_perc_chip < $max_filled_perc) {

	# 1) Merge segments at the current distance in the exon (cds) isect dataset
	@cds_gff = MergeSegmentsGFF($merge_distance, @cds_gff);

	# Calculate a new weighted average number of segments per exon
	$avg_seg_num_exons = GetAverageSegNum(@cds_gff);
	
	# 2) Merge segments at the current distance in the ChIP isect dataset
#	@chip_bed = MergeSegmentsBed($merge_distance, @chip_bed);
	&MergeSegmentsBed($merge_distance, @chip_bed);

	# Calculate a new weighted average percent of filled gaps in ChIP peaks
	my $avg_filled_perc_chip = &GetAverageFillPerc; 

	# Perform checks on updated stats

	# If that average now meets the target...
	if (($avg_filled_perc_chip > $max_filled_perc) and ($avg_seg_num_exons <= $target_seg_num)) {

		# Take merge distance 1 step back
		$merge_distance = $merge_distance - 1;

		print "Met target average segment number per exon: $avg_seg_num_exons\n";
		print "Average percent filled gaps in ChIP peaks exceeded the maximum this iteration ($avg_filled_perc_chip)\n";
		print "Use previous merge distance: $merge_distance\n\n";
		exit;	

	} elsif ($avg_filled_perc_chip > $max_filled_perc) {

		# Take merge distance 1 step back
		$merge_distance = $merge_distance - 1;

		print "Failed to meet target average segment number per exon. Average segments per exon: $avg_seg_num_exons\n";
		print "Average percent filled gaps in ChIP peaks exceeded the maximum the this iteration ($avg_filled_perc_chip)\n"; 
		print "Use previous merge distance: $merge_distance\n\n";
		exit;

	} elsif ($avg_seg_num_exons <= $target_seg_num) {

		print "Met target average segment number per exon: $avg_seg_num_exons\n";
		print "Current average percent filled gaps in ChIP peaks: $avg_filled_perc_chip\n";
		print "Current merge distance: $merge_distance\n";		
		print "Continuing until average percent filled gaps reaches maximum\n\n";
		$merge_distance += 1;

	# If the average still doesn't meet the target...
	} else {

		print "Updated average segment number per exon: $avg_seg_num_exons\nUpdated average percent filled gaps in ChIP peaks: $avg_filled_perc_chip\nMerge distance: $merge_distance\n\n";

		# Increment the merge distance, and the loop will continue
#		print ("merge dist: $merge_distance\tavg seg num: $avg_seg_num_exons\n");
		$merge_distance += 1;	

	}

}



exit;


# Subroutine definitions
#################################################################################
sub GetAverageSegNum {

	# Get gff file from main body
	my @gff = @_;

	# Initialize variable for the total length (sum of sizes) to zero
	# Will be used to calculate weights for each exon
	my $total_length = 0;

	# Initialize variable for the weighted average
	my $weighted_avg = 0;

	# First pass over the gff file to get total length of all 
	# cds intersects (segments) in the gff file
	foreach my $segment (@gff) {

		chomp ($segment);

		# Split the line up into its columns
		my @split = split(' ', $segment);

		# Unless the line's only whitespace or begins with "#"...
		unless (($segment =~ /^\s*$/) or ($segment =~ /^#/)) {

			# Get the size of the segment
			my $size = ($split[4] - $split[3]) + 1;

			# Add it to the total length
			$total_length += $size;

		}

	}

	# In the second pass over the gff file, for each exon (delimited by both the ID in 
	# column 8 and number in column 5), get the number of segments covering the exon,
	# calculate its weight (by dividing its length by the total length), multiply the
	# number of segments by the weight, and add it to the weighted average variable above

	# Initialize variables to keep track of the previous exon, its size and number of segs
	my $previous_exon_id = 'false';
	my $previous_exon_size = 0;
	my $previous_exon_seg_num = 0;

	# Begin second pass of the gff isects
	foreach my $segment (@gff) {

		# Remove newline
		chomp ($segment);

		# Unless the line's only whitespace or begins with "#"...
		unless (($segment =~ /^\s*$/) or ($segment =~ /^#/)) {

			# Split the line up into its columns
			my @split = split(' ', $segment);

			# Get the size of the segment
			my $size = ($split[4] - $split[3]) + 1;

			# Get the segments exon-specific id
			my $id = "$split[8]" . "_$split[5]";

			# If the current id matches the previous id
			if ($id eq $previous_exon_id) {

				# Add its size to the length of the exon
				$previous_exon_size += $size;

				# Increment the number of segments for this exon
				$previous_exon_seg_num++;

			# Else, the id doesn't match, so we've finished the previous entry.
			# Calculate and add it to the weighted average, and start the next entry
			} else { 

				# First, calculate exon's weight
				my $weight = $previous_exon_size / $total_length;

				# Now multiply the number of segments by the weight
				my $weighted_seg_number = $previous_exon_seg_num * $weight;

				# Add the weighted seg number to the weighted average
				$weighted_avg += $weighted_seg_number;

				# Finally, close the previous entry by updating loop variables
		
				# Update the 'prevous' exon id to the current exon id
				$previous_exon_id = $id;

				# Update the 'prevous' exon size to the size of the current exon seg
				$previous_exon_size = $size;

				# Update the previous exon seg num to 1
				$previous_exon_seg_num = 1;

			}

		}

	}

	# Return the weighted average to the main body
	return ($weighted_avg);

}


sub GetElementLengths {

	# Get gff from main body
	my @bed = @_;

	# Initialize a hash to hold lengths of exon isect segments
	my %hash;

	# Loop through gff
	foreach my $segment (@bed) {

		# Remove newline
		chomp ($segment);

		# Unless the line's only whitespace or begins with "#"...
		unless (($segment =~ /^\s*$/) or ($segment =~ /^#/)) {

			# Split the line up into its columns
			my @split = split(' ', $segment);

			# Get the size of the segment
			my $size = ($split[2] - $split[1]) + 1;

			# Get the segments exon-specific id
			my $id = $split[3];

			# If the segment's id already exists in the hash
			if (exists $hash{$id}) {
			
				# Add its length to that of the exon
				$hash{$id} += $size;

			# If the appropriate element doesn't exist (only on first iteration)
			} else {

				# Create an entry for the current exon with the segment's length
				$hash{$id} = $size;
	
			}

		}

	}

	return (%hash);

}	

sub MergeSegmentsGFF {

	# Collect arguments
	my $merge_distance = shift @_;
	my @gff = @_;
	my @merged;

	# Handle first line of the spreadsheet to initialize element (exon/chip isect) variables
	# First collect and split the first line
	my $first_line = shift @gff;
	my @first = split(' ', $first_line);

	# Initialize variables. "Prevous" is used to indicate the working cds/exon because in each
	# Loop iteration, the current segment will be compared to the one before it (previous)
	# Thus "prevous" values are held from the last iteration for comparison to the 'current'
	my $previous_cds = $first[8]; # ID of previous cds
	my $previous_exon = $first[5]; # Number of previous exon
	my $previous_start = $first[3];
	my $previous_end = $first[4];

	# Store the first line as the previous merged segment. This is to retain all
	# gff fields for output. Before completing an entry, replace indicies 3 and 4 with
	# $prevous_start and $previous_end respectively. Update it for each new cds or exon.
	# Quick note, since we're not actually using these merged entries, their
	# strandedness doesn't matter. Coordinates are still always from left to right
	# in the intersection gff file.
	my @merged_segment = @first;

	# Iterate over segments of gff exons
	foreach my $segment (@gff) {

		# Remove newlines
		chomp($segment);

		# Unless the line's only whitespace or begins with "#"...
		unless (($segment =~ /^\s*$/) or ($segment =~ /^#/)) {

			# Split the current segment
			my @split = split(' ', $segment);

			# Collect its cds and exon IDs and see if its the same as the previous. 
			# If they both are...
			if (($split[8] eq $previous_cds) and ($split[5] == $previous_exon)) {

				# Get the segment's distance from the last segment
				my $distance = $split[3] - (1 + $previous_end);

				# If the distance is within the merge distance...
				if ($distance <= $merge_distance) {

					# Replace the prevous end point with that of the current segment
					$previous_end = $split[4];

					
				# If the distance is too large then the 'previous' merged segment is done
				# and its time to start a new one with the current segment
				} else {

					# Complete the previous merged segment's entry by replacing its 
					# start and end with the "previous" entries
					$merged_segment[3] = $previous_start;
					$merged_segment[4] = $previous_end;
				
					# Join them with tabs	
					my $merged_segment = join("\t", @merged_segment);
				
					# Add it to the current merged gff array
					push(@merged, $merged_segment);

					# Update the start and end position of the "prevous" exon
					# to the current. And the @merged_segment array. The exon itself 
					# hasn't changed yet, only the segment
					$previous_start = $split[3];
					$previous_end = $split[4];

				}


			# If the cds ID is the same but the exon is different, then we've hit a new exon.
			# the last segment of the previous exon can be completed the same way as before
			} elsif (($split[8] eq $previous_cds) and ($split[5] !~ $previous_exon)) {

				# Complete the previous merged segment's entry by replacing its 
				# start and end with the "previous" entries
				$merged_segment[3] = $previous_start;
				$merged_segment[4] = $previous_end;

				# Join with tabs
				my $merged_segment = join("\t", @merged_segment);

				# Add it to the current merged gff array
				push(@merged, $merged_segment);

				# Update the start and end position of the "prevous" cds
				# to the current. And the @merged_segment array. The cds itself 
				# hasn't changed yet, only the exon (and in turn, segment)
				# This time update the merged_segment array to reflect the new exon
				# As well as the exon number.
				@merged_segment = @split; 
				$previous_exon = $split[5];
				$previous_start = $split[3];
				$previous_end = $split[4];	


			# If the exon ID is different, then we've hit a new cds
			} else {

				# Complete the entry same way, but also update the ID of the cds
				# and the new exon number.
				$merged_segment[3] = $previous_start;
				$merged_segment[4] = $previous_end;
				
				# Join them with tabs	
				my $merged_segment = join("\t", @merged_segment);
				
				# Add it to the current merged gff array
				push(@merged, $merged_segment);

				# Update ONLY the start and end positions of the "prevous" exon
				# to the current. The cds itself hasn't changed yet, only the 
				# segment
				@merged_segment = @split;
				$previous_cds = $split[8];
				$previous_exon = $split[5];
				$previous_start = $split[3];
				$previous_end = $split[4];

			}

		}

	}	

	# Finish last entry, dont' need to update other variables
	$merged_segment[3] = $previous_start;
	$merged_segment[4] = $previous_end;
	my $merged_segment = join("\t", @merged_segment);
	push(@merged, $merged_segment);

	return(@merged);

}

sub MergeSegmentsBed {

	undef %total_chip_lengths;
	undef %filled_chip_lengths;

	# Collect arguments
	my $merge_distance = shift @_;
	my @bed = @_;
	my @merged;

	# Handle first line of the spreadsheet to initialize element (cds/exon isect) variables
	# First collect and split the first line
	my $first_line = shift @bed;
	my @first = split(' ', $first_line);

	# Initialize variables. "Prevous" is used to indicate the working peak because in each
	# Loop iteration, the current segment will be compared to the one before it (previous)
	# Thus "prevous" values are held from the last iteration for comparison to the 'current'
	my $previous_peak = $first[3]; # ID of previous peak
	my $previous_chr = $first[0];
	my $previous_start = $first[1];
	my $previous_end = $first[2];
	my $first_length = ($previous_end - $previous_start) + 1;
	$total_chip_lengths{$previous_peak} = $first_length;
=pod
# [1 2 3] 4 5 [6 7 8 9 10]
print "$first_length\n";
print Dumper(\%total_chip_lengths);
$pause = <STDIN>;
=cut
	# Store the first line as the previous merged segment. This is to retain all
	# bed fields for output. Before completing an entry, replace indicies 1 and 2 with
	# $prevous_start and $previous_end respectively. Update it for each new peak or exon.
	my @merged_segment = @first;

	# Iterate over segments of bed exons
	foreach my $segment (@bed) {

		# Remove newlines
		chomp($segment);

		# Unless the line's only whitespace or begins with "#"...
		unless (($segment =~ /^\s*$/) or ($segment =~ /^#/)) {

			# Split the current segment
			my @split = split(' ', $segment);

			# Collect its chromsome and see if its the same as the previous. 
			# If it is...
			if ($split[0] eq $previous_chr) {

				# Get the segment's distance from the last segment
				my $distance = $split[1] - (1 + $previous_end);

				# If the distance is within the merge distance...
				if ($distance <= $merge_distance) {

					# Add the distance to its element in the filled_chip_lengths hash
					# If there is already an element present (will be after 1st iteration)
					### Replace this with a subroutine later
					if (exists $filled_chip_lengths{$previous_peak}) {
			
						# Add distance to current length
						$filled_chip_lengths{$previous_peak} += $distance;


					# If the appropriate element doesn't exist (only on first iteration)
					} else {

						# Create one, set it to current distance
						$filled_chip_lengths{$previous_peak} = $distance;
	
					}

					# Add the distance AND length to its element in the total_chip_lengths hash
					# If there is already an element present (will be after 1st iteration)
					### Replace this with a subroutine later
					if (exists $total_chip_lengths{$previous_peak}) {

						# calculate its length
						my $length = ($split[2] - $split[1]) + 1;

						# Add its length and distance to total length
						$total_chip_lengths{$previous_peak} += ($length + $distance);
=pod
print "$length\n$distance\n";
print Dumper(\%total_chip_lengths);
$pause = <STDIN>;
=cut

					# If the appropriate element doesn't exist (only on first iteration)
					} else {

						# Create one, set it to the length plus distance
						my $length = ($split[2] - $split[1]) + 1;
						$total_chip_lengths{$previous_peak} = ($length + $distance);
	
					}

					# Replace the prevous end point with that of the current segment
					$previous_end = $split[2];

					
				# If the distance is too large then the 'previous' merged segment is done
				# and its time to start a new one with the current segment
				} else {

					# Complete the previous merged segment's entry by replacing its 
					# start and end with the "previous" entries
					$merged_segment[1] = $previous_start;
					$merged_segment[2] = $previous_end;
				
					# Join them with tabs	
					my $merged_segment = join("\t", @merged_segment);

					# Update the id, start and end positions of the "prevous" peak
					# to the current. 
					@merged_segment = @split;
					$previous_peak = $split[3];
					$previous_chr = $split[0];
					$previous_start = $split[1];
					$previous_end = $split[2];
					my $length = ($previous_end - $previous_start) + 1;
					$total_chip_lengths{$previous_peak} = $length;
					
				}

			# If the chromosome is different
			} else {

				# Complete the entry same way, but also update the ID of the peak
				$merged_segment[1] = $previous_start;
				$merged_segment[2] = $previous_end;
				
				# Join them with tabs	
				my $merged_segment = join("\t", @merged_segment);
				
				# Add it to the current merged bed array
				push(@merged, $merged_segment);

				# Update the id, start and end positions of the "prevous" peak
				# to the current. 
				@merged_segment = @split;
				$previous_peak = $split[3];
				$previous_chr = $split[0];
				$previous_start = $split[1];
				$previous_end = $split[2];
				my $length = ($previous_end - $previous_start) + 1;
				$total_chip_lengths{$previous_peak} = $length;

			}

		}

	}	

	# Finish last entry, dont' need to update other variables
	$merged_segment[1] = $previous_start;
	$merged_segment[2] = $previous_end;
	my $merged_segment = join("\t", @merged_segment);
	push(@merged, $merged_segment);


#print Dumper(\%filled_chip_lengths);
#$pause = <STDIN>;

	#return(@merged);

}

sub GetAverageFillPerc {

	my $weighted_average = 0;

	my @array = sort keys %total_chip_lengths;

	my $total_chip_length = sum values %total_chip_lengths;

	foreach my $key (@array) {

		if (exists $filled_chip_lengths{$key}) {

			my $filled_peak_length = $filled_chip_lengths{$key};

			my $total_peak_length = $total_chip_lengths{$key};

			my $percent_filled = 100 * ($filled_peak_length / $total_peak_length);

#			print "$percent_filled\n";
#			print "$key\n";
#			print "$filled_peak_length\n";
#			print "$total_peak_length\n\n";
#			$pause = <STDIN>;

			my $weight = $total_peak_length / $total_chip_length;

			my $weighted_length = $percent_filled * $weight;

			$weighted_average += $weighted_length;

		} else {

			my $filled_peak_length = 0;

			my $total_peak_length = $total_chip_lengths{$key};

			

		}

	}

	return($weighted_average);

}








