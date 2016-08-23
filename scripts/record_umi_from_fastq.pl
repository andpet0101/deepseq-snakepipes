#!/usr/bin/perl

use Getopt::Long;
use File::Basename;
use Data::Dumper;

my $script = basename($0);
my ($fqin1,$fqin2,$fqout1,$fqout2) = ("","","","");
my ($fhin1,$fhin2,$fhout1,$fhout2) = ("","","","");
my $cutadapt_info = "";
my %adapter_info = ();
my $umi_lib_pattern = "";
my $umi_trim_spec = "";
my @vals = ();
my @umi_trim_spec = ();
my @adapter_required = ();
my @umi_len = ();
my @umi_parts = ();
my $trimmed_seq = "";
my $trimmed_qual = "";
my $umi_regexp_str = "";
my @id = ();
my $min_len = 0;
my $umi_distibution = "";
my %umi_distribution = ();
my %umi_stats = ();

GetOptions("fqin1=s" => \$fqin1,
	"fqin2:s" => \$fqin2,
	"fqout1=s" => \$fqout1,
	"fqout2:s" => \$fqout2,
	"cutadapt:s" => \$cutadapt_info,
	"umi_libs:s" => \$umi_lib_pattern,
	"umi_spec=s" => \$umi_trim_spec,
	"min_len:i" => \$min_len,
	"umi_distribution:s" => \$umi_distribution);

die "$script: Please specify an input and an output fastq file via --fqin1 and --fqout1!" if(!$fqin1 || !$fqout1);
die "$script: Could not find the first input file!" if(!-e $fqin1);
die "$script: Please specify a second output fastq file via --fqout2!" if($fqin2 && !$fqout2);
die "$script: Could not find the second input file!" if($fqin2 && !-e $fqin2);
die "$script: Please specify the position of the UMI: 'first_nt_read1,last_nt_read1[;first_nt_read2,last_nt_read2]' ; positions with '!' require to have an adapter trimmed in advance." if(!$umi_trim_spec);
die "$script: If you use '!' for the UMI, then you have to provide a cutadapt info file via --cutadapt!" if($umi_trim_spec=~/!/ && !$cutadapt_info);
die "$script: Could not find the cutadapt info file!" if($cutadapt_info && !-e $cutadapt_info);

# parse the UMI position specification
@umi_trim_spec = split(/[;]/,$umi_trim_spec);
if($umi_trim_spec[0]=~/^(\d+)(!?),(\d+)(!?)$/){
	push(@umi_len,$1,$3);
	push(@adapter_required,$2,$4);
}

if(scalar(@umi_trim_spec)==2){
	if($umi_trim_spec[1]=~/^(\d+)(!?),(\d+)(!?)$/){
		push(@umi_len,$1,$3);
		push(@adapter_required,$2,$4);
	}
}

# read cutadapt info if neccessary
if($cutadapt_info){
	open(CUT,"<","$cutadapt_info");
	while(<CUT>){
		chomp;
		@vals = split(/\t/,$_);
		
		$adapter_info{$vals[0]} = [@vals];
	}
	close(CUT);
}

# read in files
open(IN1,(($fqin1=~/\.gz$/ ? "gunzip -cd" : "cat")." $fqin1 |"));
$fhin1 = *IN1;
if($fqin2){
	open(IN2,(($fqin2=~/\.gz$/ ? "gunzip -cd" : "cat")." $fqin2 |")) if($fqin2);
	$fhin2 = *IN2;
}

open(OUT1,("| ".($fqout1=~/\.gz$/ ? "gzip -c" : "cat")." > $fqout1"));
$fhout1 = *OUT1;
if($fqout2){
	open(OUT2,("| ".($fqout2=~/\.gz$/ ? "gzip -c" : "cat")." > $fqout2")) if($fqout2);
	$fhout2 = *OUT2;
}

my %umi_stats = ("Total fragments"=>0,"Fragments with UMI"=>0,"Fragment length with UMI"=>0,"Fragments that were too short"=>0);
my $u = 0;

while(@fqs = read_fastq_records($fhin1,$fhin2)){
	# record all identified UMIs
	@umi_parts = ("","");
	
	$umi_stats{"Total fragments"}++;

	# if the option umi_libs is used: 1 - always search for UMI, 0 - search only in files whose name correspond to the pattern in umi_libs (regexp)
	if($umi_lib_pattern==1 || $fqin1=~/^$umi_lib_pattern/){
		
		# read 1 - sequence and quality
		$trimmed_seq = $fqs[0]->[1];
		$trimmed_qual = $fqs[0]->[3];
		
		if($umi_len[0]+$umi_len[1]>0 && $umi_len[0]+$umi_len[1]<length($trimmed_seq)){
			# last nts of read 1
			if(!$adapter_required[1] || 1){
				$trimmed_seq =~ s/^(\S+)(\S{$umi_len[1],$umi_len[1]})$/$1/;
				$umi_parts[1] = $2;
				$trimmed_qual =~ s/^(\S+)(\S{$umi_len[1],$umi_len[1]})$/$1/;
			}

			# first nts of read 1
			if(!$adapter_required[0] || 1){
				$trimmed_seq =~ s/^(\S{$umi_len[0],$umi_len[0]})(\S+)$/$2/;
				$umi_parts[0] = $1;
				$trimmed_qual =~ s/^(\S{$umi_len[0],$umi_len[0]})(\S+)$/$2/;
			}
		}else{
			$umi_stats{"Fragments that were too short"}++;
			next;
		}
		
		$fqs[0]->[1] = $trimmed_seq;
		$fqs[0]->[3] = $trimmed_qual;
		
		if(scalar(@fqs)==2){
			push(@umi_parts,"","");
			
			# read 2 - sequence and quality
			$trimmed_seq = $fqs[1]->[1];
			$trimmed_seq = $fqs[1]->[3];

			if($umi_len[2]+$umi_len[3]>0 && $umi_len[2]+$umi_len[3]<length($trimmed_seq)){
				# last nts of read 2
				if(!$adapter_required[3] || 1){
					$trimmed_seq =~ s/^(\S+)(\S{$umi_len[3],$umi_len[3]})$/$1/;
					$umi_parts[3] = $2;
					$trimmed_qual =~ s/^(\S+)(\S{$umi_len[3],$umi_len[3]})$/$1/;
				}

				# first nts of read 2
				if(!$adapter_required[2] || 1){
					$trimmed_seq =~ s/^(\S{$umi_len[2],$umi_len[2]})(\S+)$/$2/;
					$umi_parts[2] = $1;
					$trimmed_qual =~ s/^(\S{$umi_len[2],$umi_len[2]})(\S+)$/$2/;
				}
			}else{
				$umi_stats{"Fragments that were too short"}++;
				next;
			}
			
			$fqs[1]->[1] = $trimmed_seq;
			$fqs[1]->[3] = $trimmed_qual;
		}
		
		# filter for min len
		
		if(length($fqs[0]->[1])<$min_len || (scalar(@fqs)==2 && length($fqs[1]->[1])<$min_len)){
			$umi_stats{"Fragments that were too short"}++;
			next;
		}
		
		# add UMI separated by '+' to the sequence name
		# http://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
		
		@id = split(/\s/,$fqs[0]->[0]);
		$id[0] .= ":".join("+",@umi_parts);
		$fqs[0]->[0] = join(" ",@id); 

		if(scalar(@fqs)==2){
			@id = split(/\s/,$fqs[1]->[0]);
			$id[0] .= ":".join("+",@umi_parts);
			$fqs[1]->[0] = join(" ",@id);
		}
		
		$umi_stats{"Fragments with UMI"}++;
		$umi_stats{"Fragment length with UMI"} += length($fqs[0]->[1]);

		if(scalar(@fqs)==2){
			$umi_stats{"Fragment length with UMI"} += length($fqs[1]->[1]);
		}

		if($umi_distribution){
			$u = join("+",@umi_parts);
			$umi_distribution{$u} ||= 0;
			$umi_distribution{$u}++;
		}
	}else{
		print STDERR "$script: Not trimming UMIs in $fqin1".($fqin2 ? " and $fqin2 " : "")." because basename does not match pattern specified by 'umi_libs'!\n";
	}

	write_fastq_records(\@fqs,$fhout1,$fhout2);
}

close($fhout1);
close($fhout2) if($fqout2);

close($fhin1);
close($fhin2) if($fqin2);

if($umi_distribution){
	open(UMI,("| ".($umi_distribution=~/\.gz$/ ? "gzip -c" : "cat")." > $umi_distribution"));
	print UMI "UMI\tcount\n";
	foreach $u (sort(keys(%umi_distribution))){
		print UMI "$u\t".$umi_distribution{$u}."\n";
	}
	close(UMI) if($umi_distribution);
}

foreach my $s ("Total fragments","Fragments with UMI","Fragment length with UMI","Fragments that were too short"){
	print "$s: ".$umi_stats{$s}."\n";
}


sub read_fastq_records{
	my $fh1 = shift;
	my $fh2 = shift || "";
	
	# normal - all fastq records read
	return () if( ($fh1 && !$fh2 && eof($fh1)) || ($fh1 && $fh2 && eof($fh1) && eof($fh2) ));
	# not normal - differing number of fastq records between fh1 and fh2
	die "$script: Paired fastq files do not have the same number of lines/entries" if($fh1 && $fh2 && (eof($fh1) ^ eof($fh2)));
	
	my @fq1 = ();
	my @fq2 = ();
	my $line = "";

	foreach my $i (1..4){
		die "$script: Last record incomplete in fastq file 1" if(eof($fh1));
		$line = <$fh1>;
		chomp($line);
		push(@fq1,$line);

		if($fh2){
			die "$script: Last record incomplete in fastq file 2" if(eof($fh2));
			$line = <$fh2>;
			chomp($line);
			push(@fq2,$line);
		}
	}

	if(@fq1 && @fq2){return ([@fq1],[@fq2]);}
	elsif(@fq1){return ([@fq1]);}
	else{return ();}
}

sub write_fastq_records{
	my $fq_ref = shift;
	my $fh1 = shift;
	my $fh2 = shift || "";

	return "" if(!@{$fq_ref} || !@{$fq_ref->[0]});

	print $fh1 join("\n",@{$fq_ref->[0]})."\n";
	
	if($fh2){
		print $fh2 join("\n",@{$fq_ref->[1]})."\n";
	}

	return "";
}

