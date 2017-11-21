#!/usr/bin/perl

use Data::Dumper 'Dumper';
use List::Util qw(min max);

# known miRNA contains all mature miRNA and must be a file like this:
#
#chr1    .       miRNA   1312502 1312521 .       -       .       ID=MIMAT0027356;Alias=MIMAT0027356;Name=hsa-miR-6727-3p;Derives_from=MI0022572
#chr1    .       miRNA   1339682 1339703 .       -       .       ID=MIMAT0027516;Alias=MIMAT0027516;Name=hsa-miR-6808-5p;Derives_from=MI0022653
#chr1    .       miRNA   1339650 1339670 .       -       .       ID=MIMAT0027517;Alias=MIMAT0027517;Name=hsa-miR-6808-3p;Derives_from=MI0022653


my $known_miRNA_file = $ARGV[0] || die "cat mirdeep_star.result | parse_mirdeep_star_results.pl known_miRNA.gff";
my @fields = ();

my @known_miRNAs = ();
my ($mirna_id,$mirna_name,$mirna_derives_id) = ("","","");
open(MIRNAS,"<",$known_miRNA_file);
while(<MIRNAS>){
	chomp;
	@fields = split(/\t/,$_);
	next if(scalar(@fields)<9);	
	
	# mature miRNA ID is in attribute Alias
	if($fields[8]=~/Alias=([^;]+)/){$mirna_id = $1;}
	else{die "Cannot find Alias attribute in gff file";}
	
	if($fields[8]=~/Name=([^;]+)/){$mirna_name = $1;}
	else{die "Cannot find ID attribute in gff file";}

	if($fields[8]=~/Derives_from=([^;]+)/){$mirna_derives_id = $1;}
	else{die "Cannot find Derives_from attribute in gff file";}
	
	push(@known_miRNAs,[$fields[0],$fields[3],$fields[4],$fields[6],$mirna_id,$mirna_name,$mirna_derives_id]);
}
close(MIRNAS);

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

my $novel_ct = 1;
my $miRNA = "";
my $i = 0;
my ($mat_chr,$mat_start,$mat_end) = ("","","");
my ($known_chr,$known_start,$known_end) = ("","","");
my $overlap = "";

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

	# finally check if miRNA is known or novel
	$mat_chr = $fields[2];
	$mat_strand = $fields[3];
	($mat_start,$mat_end) = split(/-/,$fields[6]);
	for($i=0;$i<scalar(@known_miRNAs);$i++){
		($known_chr,$known_start,$known_end,$known_strand,$mirna_id,$mirna_name,$mirna_derives_id) = @{$known_miRNAs[$i]};
		$overlap = 0;
		
		if($mat_chr eq $known_chr && $mat_strand eq $known_strand){
			$overlap = max(min($known_end,$mat_end)-max($known_start,$mat_start)+1,0);
			if($overlap>=19){
				$fields[0] = join("/",$mirna_id,$mirna_derives_id,$mirna_name);
				last;		
			}
		}
	}
	if($i==scalar(@known_miRNAs)){
		$fields[0] = "novelMiR_".$novel_ct;
		$novel_ct++;
	}
	
	print join("\t",@fields)."\n";
}


