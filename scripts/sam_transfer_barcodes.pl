#!/usr/bin/perl

while(<STDIN>){
	chomp;
	unless(m/^@/){
		if(m/:CELL_([ACGTN-]+)/){$_ .= "\tXC:Z:$1";}
		if(m/:UMI_([ACGTN-]+)/){$_ .= "\tXM:Z:$1";}
	}
	print $_."\n";
}
