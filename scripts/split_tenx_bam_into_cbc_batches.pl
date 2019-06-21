#!/usr/bin/perl

use File::Basename;

# bam file must be sorted by CB:Z (samtools sort)
my $script = basename($0);
my $bam_file = $ARGV[0] || die "$script cbc_sorted_bam_file directory_for_output";
my $output_directory = $ARGV[1] || die "$script cbc_sorted_bam_file directory_for_output";

$output_directory =~ s/\/$//;

my $prefix = basename($bam_file);$prefix =~ s/\.bam//;

# split into 2mio records files; but keep reads with the same CB:Z in the same file
my $records_per_bucket = 2000000;

# also write names of the 2mio record files to separate files (10 record files each) 
my $buckets_per_file = 10;

# save SAM header 
my $sam_header = "";
my $line = "";
open(SAM_IN,"samtools view -H $bam_file |");
while($line = <SAM_IN>){
    $sam_header .= $line;
}
close(SAM);

# now parse 2mio reads and save in new bam file; but keep reads with the same CB:Z in the same file
my $i = 1;
my $f = 1;
my $num_records_parsed = 0;
my $last_cbc = "-";
my $cbc = "";
my $files_per_fofn = 1;

open(SAM_IN,"samtools view -@ 8 $bam_file |");
open(BAM_OUT,"| samtools sort -l 1 -@ 8 -m 2G -T /tmp/$prefix.$i -O BAM -o $output_directory/$prefix.$i.bam -");
open(FOFN,">","$output_directory/$prefix.$f.fofn");
print FOFN "$output_directory/$prefix.$i.bam\n";
print "$output_directory/$prefix.$f.fofn\n";

print BAM_OUT $sam_header;

while($line = <SAM_IN>){    
    # parse cell barcode
    if($line=~/CB:Z:(\S+)/){
        $cbc = $1;
    }else{
        $cbc = "unassigned";
        $line =~ s/\n/\tCB:Z:unassigned\n/;
    }
    
    # check if 10mio records have been reached; but keep reads with the same CB:Z in the same file
    if($num_records_parsed>=$records_per_bucket && $cbc ne $last_cbc){
        # close current file
        close(BAM_OUT);
        
        # close fofn if there are more than 10 file names in it
        if($files_per_fofn>$buckets_per_file){
            close(FOFN);
            $f++;
            open(FOFN,">","$output_directory/$prefix.$f.fofn");
            print "$output_directory/$prefix.$f.fofn\n";
            $files_per_fofn = 0;
        }
        
        # open new file and start with the header
        $i++;
        open(BAM_OUT,"| samtools sort -O BAM -o $output_directory/$prefix.$i.bam -");
        print BAM_OUT $sam_header;
        $num_records_parsed = 0;
        $files_per_fofn++;
        
        print FOFN "$output_directory/$prefix.$i.bam\n";
    }
    
    # write record and remember last cbc
    print BAM_OUT $line;
    $num_records_parsed++;
    $last_cbc = $cbc;
}

close(FOFN);
close(BAM_OUT);
close(SAM_IN);
