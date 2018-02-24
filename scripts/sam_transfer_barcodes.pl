#!/usr/bin/perl

my @f = ();

while(<STDIN>){
	chomp;
	unless(m/^@/){
		# UMI is at the end of the read name: leave it there and just add a MI tag for molecular identifier
		@f = split(/\t/,$_);
		if($f[0]=~/:([A-Z]+$)/){
			$_ .= "\tMI:Z:$1";
		}
	}
	print $_."\n";
}
