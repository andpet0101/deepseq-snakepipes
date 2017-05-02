#!/usr/bin/perl

my @f = ();

while(<STDIN>){
	if(m/^\@/){
		print;
		next;
	}else{
		@f = split(/\t/,$_);
		if(scalar(@f)<11){
			print;
			next;
		}

		$f[5]=~s/(\d+)N(\d+)I/$1D$2I/g;
		$f[5]=~s/(\d+)I(\d+)N/$1I$2D/g;
		
		print join("\t",@f);
		next;
	}
}
