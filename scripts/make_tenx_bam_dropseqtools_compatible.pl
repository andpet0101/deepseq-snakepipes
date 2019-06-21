#!/usr/bin/perl

use File::Basename;

# bam file must be sorted by CB:Z (samtools sort)
my $script = basename($0);
my $tenx_bam = $ARGV[0] || die "$script tenx_bam tenx_bam_for_dropseqtools [5p|3p]";
my $tenx_bam_for_dropseqtools = $ARGV[1] || die "$script tenx_bam tenx_bam_for_dropseqtools [5p|3p]";
my $assay = $ARGV[2] || die "$script tenx_bam tenx_bam_for_dropseqtools [5p|3p]";

# assay: 5p - 5'prime, 3p - 3'prime

# save SAM header 
my $sam_header = "";
my $line = "";
open(SAM_IN,"samtools view -H $tenx_bam |");
while($line = <SAM_IN>){
    $sam_header .= $line;
}
close(SAM);

open(SAM_IN,"samtools view -@ 8 $tenx_bam |");
open(BAM_OUT,"| samtools view -b -o $tenx_bam_for_dropseqtools -");

print BAM_OUT $sam_header;

my @fields = ();
my $val = "";

while($line = <SAM_IN>){    
    chomp($line);
    
    
    # gene region/function: RE tag indicates the region type of this alignment => convert to gf tag
    if($line =~ m/RE:A:(\S+)/){
        $val = $1;
        if($val eq "E"){
            $line.="\tgf:Z:UTR";
        } elsif($val eq "N"){
            $line.="\tgf:Z:INTRONIC";
        } elsif($val eq "I"){
            $line.="\tgf:Z:INTERGENIC";
        }
    }
    
    # gene strand: TX tag should be present, then strand of read indicates gene strand
    if($line =~ m/TX:Z:/){
        @fields = split(/\t/,$line);
        if($fields[1] & 16){
            if($assay eq "5p"){
                # 5p: read is antisense to gene
                $line.="\tgs:Z:+";
            }else{
                # 3p: read is sense to gene
                $line.="\tgs:Z:-";
            }
        }else{
            if($assay eq "5p"){
                # 5p: read is antisense to gene
                $line.="\tgs:Z:-";
            }else{
                # 3p: read is sense to gene
                $line.="\tgs:Z:+";
            }
        }
    }
    
    # Add gene name tag (GX)
    if($line =~ /GX:Z:(\S+)/){
        $line.="\tgn:Z:$1";
    }
    
    print BAM_OUT $line."\n";
}

close(BAM_OUT);
close(SAM_IN);
