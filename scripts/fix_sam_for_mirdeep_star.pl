#!/usr/bin/perl

# Mirdeep* removes bowtie1 mappings with mismatches based on the NM tag. However, one mismatch should be allowed to account for SNPs. 
# This script changes a NM:i:1 tag to NM:i:0 given that the mapping a) covers the entire read and b) does not contain indels.

my $line = "";
my @fields = ();
my $read_len = "";

while($line = <STDIN>){
	chomp($line);
	# if read and 1-mismatch mapping
	if($line!~/^\@/ and $line=~/\bNM:i:1\b/){
		@fields = split(/\t/,$line);
		$read_len = length($fields[9]);
		
		# minimum read length of 20 and mapping covers entire read, then change NM:i:1 to NM:i:0
		if($read_len>=20 && $fields[5]=~/\b$read_len\M\b/){
			$line=~s/NM:i:1/NM:i:0/;
		}
	}
	print $line."\n";
}
