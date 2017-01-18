#!/usr/bin/perl

use File::Basename;

my $script = basename($0);
my $fastq = $ARGV[0] || die "$script trimmed_reads.fastq[.gz] gecko_lib.csv > trimmed_reads_rescued.fastq";
my $sgRNAs = $ARGV[1] || die "$script trimmed_reads.fastq[.gz] gecko_lib.csv > trimmed_reads_rescued.fastq";

# read sgRNAs and build hash
my %sgRNAs = ();
my $id = "";
my $seq = "";
my $subseq = "";
my @vals = ();
open(SGRNAS,"<",$sgRNAs);
while(<SGRNAS>){
	chomp;
	next if(m/^id/);
	@vals = split(/,/,$_);
	if(@vals && scalar(@vals)!=3){
		die "$script: file $sgRNAs does not have three columns!\n";
	}
	
	$seq = $vals[1];
	$id = $vals[0];

	$sgRNAs{$seq} = [$id,$seq];
	
	# rescue search: remove first and/or last bp and add the respective sequences to the hash
	$subseq = substr($seq,1,length($seq)-1);
	if($sgRNAs{$subseq} && $sgRNAs{$subseq}->[1] ne $seq){
		warn "$script: After truncation of the first bp $seq points to another different sequence - skip this!\n";	
	}else{
		$sgRNAs{$subseq} = $sgRNAs{$seq};
	}

	$subseq = substr($seq,0,length($seq)-1);
	if($sgRNAs{$subseq} && $sgRNAs{$subseq}->[1] ne $seq){
		warn "$script: After truncation of the last bp $seq points to another different sequence - skip this!\n";	
	}else{
		$sgRNAs{$subseq} = $sgRNAs{$seq};
	}
}
close(SGRNAS);

# now go through fastq, search and rescue if needed
my $line = "";
my $hit = "";
my $total = 0;
my $valid = 0;
my $original = 0;
my $rescued = 0;
my $unmapped = 0;
open(FASTQ,($fastq=~/\.gz$/ ? "zcat" : "cat")." $fastq |");
while(!eof(FASTQ)){
	chomp;
	
	#### id line ###########
	$line = <FASTQ>;
	print $line;
	$total++;
	
	#### sequence line ###########
	$line = <FASTQ>;
	chomp($line);

	$seq = $line;
	# look up sequence
	$hit = $sgRNAs{$seq} || "";

	if($hit){$original++;}
	# not hit? Remove first bp and look for hits	
	else{
		$subseq = substr($seq,1,length($seq)-1);
		$hit = $sgRNAs{$subseq} || "";
		
		# still no hit? Remove last bp and look for hits
		if(!$hit){
			$subseq = substr($seq,0,length($seq)-1);
			$hit = $sgRNAs{$subseq} || "";
		}
		if($hit){$rescued++;}
	}

	# do we have hit? Use sequence of sgRNA otherwise original sequence
	if($hit){
		$seq = $hit->[1];
		$valid++;
	}else{
		$unmapped++;
	}
	print $seq."\n";
	
	#### quality id line ###########
	$line = <FASTQ>;
	print $line;
	
	#### quality line ###########
	$line = <FASTQ>;
	chomp($line);
	
	# quality line longer than sequence line? remove last quality values
	if(length($line)>length($seq)){
		print substr($line,0,length($seq))."\n";
	}
	# quality line shorter than sequence line? fill up with last quality values
	elsif(length($line)<length($seq)){
		print $line.(substr($line,length($line)-1,1) x (length($seq) - length($line)))."\n";
	}
	# same length? do nothing
	else{
		print $line."\n";	
	}
}
close(FASTQ);

print STDERR "Total: $total, thereof valid: $valid (original: $original, rescued: $rescued) and unmapped: $unmapped\n";
