#!/usr/bin/perl

my @fields = ();
my @extra_fields = ();
my %extra_field_keys = ();
my $is_extra_field_key = 0;
my $e = "";
my @l = ();

while(<STDIN>){
	chomp;
	@fields = split(/\t/,$_);
	if(m/^##/){
		if(m/^## Extra column keys:/){
			$is_extra_field_key = 1;
		}
		elsif($is_extra_field_key && m/^## (\S+)/){
			push(@extra_field_keys,$1);
		}
		print "$_\n";
		next;
	}
	elsif(m/^#Uploaded_variation/){
		# split location column into chr and position column
		splice(@fields,1,1,("Chr","Position"));

		# remove 'Extra' column header and insert a new header for each extra column
		pop(@fields);
		push(@fields,@extra_field_keys);
	}
	else{
		# split location column into chr and position column
		@l = split(/:/,$fields[1]);
		splice(@fields,1,1,@l);
		
		# parse and add extra columns
		#
		# convert string 'field1=value1;field2=value2' into ('field1','value1','field2','value2') and interpret as hash
		%extra_fields = split(/[;=]/,$fields[-1]);
		pop(@fields);
		foreach $e (@extra_field_keys){
			if(defined($extra_fields{$e}) && length($extra_fields{$e})>0){
				push(@fields,$extra_fields{$e});
			}
			else{
				push(@fields,"-");			
			}
		}
		
		#also fix location if neccessary
		
	}

	print join("\t",@fields)."\n";
}
