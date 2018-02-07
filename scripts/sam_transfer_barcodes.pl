#!/usr/bin/perl

while(<STDIN>){
	chomp;
	unless(m/^@/){
		if(s/:BARCODE_([ACGTN-]+)//){$_ .= "\tBC:Z:$1";}
		if(s/:CELL_([ACGTN-]+)//){$_ .= "\tXC:Z:$1";}
		if(s/:UMI_([ACGTN-]+)//){$_ .= "\tXM:Z:$1";}
	}
	print $_."\n";
}
