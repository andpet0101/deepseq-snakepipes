#!/usr/bin/perl

my @fields = ();

my @read_locations = ();
my $rl = "";
my %read_ids = ();
my $ri = "";
my @ri = ();
my %umis = ();
my @umis = ();
my $sum = 0;

my $header = <STDIN>;
chomp($header);
@fields = split(/\t/,$header);
$fields[11] = "reads with UMI";
$fields[12] = "number of UMI";
$fields[13] = "UMI distribution";
print join("\t",@fields)."\n";

while(<STDIN>){
	chomp;
	@fields = split(/\t/,$_);

	# mirRNA score chr strand location_precursor expression location_miRNA precursor precursor_sec_struct mature_seq_precursor_seq reads 
	
	# parse reads and count different UMIs
	@read_locations = split(/;/,$fields[11]);
	%read_ids = ();
	foreach $rl (@read_locations){
		$rl =~ s/^\[//;
		$rl =~ s/\].+//;
	
		foreach $ri (split(/,/,$rl)){
			$read_ids{$ri} = 1;
		}	
	}
	%umis = ();
	foreach $ri (keys(%read_ids)){
		@ri = split(/:/,$ri);		
		if($ri[-1]=~/[ACGT]/){
			$umis{$ri[-1]} ||= 0;
			$umis{$ri[-1]}++;
		}
	}
	
	@umis = sort(keys(%umis));
	$sum = 0;
	map{$sum+=$umis{$_}}@umis;
	
	if(@umis){
		$fields[11] = $sum;
		$fields[12] = scalar(@umis);
		$fields[13] = join(",",map{$_.":".$umis{$_}}@umis);	
	}
	else{
		$fields[11] = 0;
		$fields[12] = 0;
		$fields[13] = "";
	}
	
	print join("\t",@fields)."\n";
}
