#!/usr/bin/perl

use FindBin;
use File::Basename;
use Getopt::Long qw(:config pass_through);
use Cwd 'abs_path';

############
# Settings #
############

# script, snakemake bin, snakemake workflow directory and configuration directory
my $script = basename($0);
my $script_dir = $FindBin::RealBin;

our $snakemake_bin = "snakemake";

our $snakemake_workflow_dir = ".";
if($ENV{'SNAKEMAKE_WORKFLOW_DIR'}){
	$snakemake_workflow_dir = $ENV{'SNAKEMAKE_WORKFLOW_DIR'};
}elsif(-e "$script_dir/../snakemake"){
	$snakemake_workflow_dir = "$script_dir/../snakemake"
}

our $snakemake_conf_dir = ".";
if($ENV{'SNAKEMAKE_CONF_DIR'}){
	$snakemake_conf_dir = $ENV{'SNAKEMAKE_CONF_DIR'};
}elsif(-e "$script_dir/../config"){
	$snakemake_conf_dir = "$script_dir/../config"
}


# time and file permissions
my $start_time = time();
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
$mon++;
$year += 1900;

umask(002);

# options and help
if(scalar(@ARGV)==0){
	help();	
	# exits automatically
}

my ($workflow,$snakefile,$configfile,$cluster_configfile,$description,$jobscript) = ("","","","","","");
my $logdir = "snakemake_logs";
my $jobs = 1;
my $dryrun = 0;
my $use_drmaa = 1;
my $use_qsub = 0;
my $show_snakemake_options = 0;
my $send_mail = "";

GetOptions('help|h' => \&help,
	'snakemake_help' => \&snakemake_help,
	'workflow=s' => \$workflow,
	'snakefile|s=s' => \$snakefile,
	'configfile=s' => \$configfile,
	'cluster-config|u=s' => \$cluster_configfile,
	'log-dir' => \$logdir,
	'use_drmaa' => \$use_drmaa,
	'use_qsub' => \$use_qsub,
	'jobs|j=i' => \$jobs,
	'dryrun' => \$dryrun,
	'description' => \$description,
	'mail=s' => \$send_mail);

my $snakemake_params = join(" ",@ARGV);

# workflows and snakemake file checks
if($workflow && !$snakefile){
	$snakefile = "$snakemake_workflow_dir/$workflow.snakemake";
}

if($snakefile){
	if(!-e $snakefile){
		print STDERR "$script: The snakemake file for the analysis ($snakefile) does not exist.\n";
		exit(1);	
	}
}else{		
	print STDERR "$script: Please select a workflow from the help or provide a snakemake file.\n";
	exit(1);
}

# config file checks
if($configfile && !-e $configfile){
	print STDERR "$script: The configuration file for the analysis ($configfile) does not exist.\n";
	exit(1);	
}

# cluster config file checks
if($cluster_configfile){
	if(!-e $cluster_configfile){
		print STDERR "$script: The configuration file for the cluster ($cluster_configfile) does not exist.\n";
		exit(1);	
	}
	if(basename($cluster_configfile) !~ /\.json$/){
		print STDERR "$script: The configuration file for the cluster must be in JSON format and the name must end on '.json'.\n";
		exit(1);
	}
}elsif(!$cluster_configfile && -e "$snakemake_conf_dir/cluster.json"){
	$cluster_configfile = "$snakemake_conf_dir/cluster.json";
}else{
	print STDERR "$script: No cluster configuration file provided or found - all jobs will use standard resources.\n";
}

# jobscript checks
$jobscript = "$snakemake_conf_dir/jobscript.sh";
if(!-e $jobscript){
	print STDERR "$script: Could not find the jobscript template script ($jobscript) - continue without!\n";
	$jobscript = "";
}

# check and create log dir
if(!-e dirname($logdir)){
	print STDERR "$script: The directory for the logs ($logdir) was not found.\n";
	exit(1);
}

$logdir = abs_path($logdir);
$logdir .= sprintf("/snakemake_%02d%02d%4d_%02d%02d",$mday,$mon,$year,$hour,$min);
if(-e $logdir && !$dryrun){system("rm -r $logdir");}	
if(!$dryrun){system("mkdir -p $logdir");}


# two methods for cluster submission: DRMAA (default) and qsub
$use_drmaa = 0 if($use_qsub);


#################################################
# Print snakemake file description upon request #
#################################################

if($snakefile && $description){
	my $line_nr = 0;
	my $in_header = 0;
	open(SNAKE,"<","$snakefile");
	while(<SNAKE>){
		chomp;
		$line_nr++;
		# if first line starts with '#', we are in the header comment		
		if($line_nr==1 && m/^#/){
			$in_header = 1;		
		}
		# one of the following lines starts with anything else then '#', we are not in the header comment anymore
		if($line_nr>1 && !m/^#/){
			$in_header = 0;
		}
		if($in_header){
			s/^#//;
			print;		
		}
	}
	exit(0);
}

###########
# Command #
###########
my $cmd = "";

$cmd .= $snakemake_bin;

# for proper colour-coding
unless(`which unbuffer`=~/^which: no unbuffer/){
	$cmd = "unbuffer $cmd";
}

$cmd = "module load apps/snakemake;".$cmd;

$cmd .= " --snakefile $snakefile";
$cmd .= " --configfile $configfile" if($configfile);
$cmd .= " --printshellcmds --dryrun" if($dryrun);
$cmd .= " --timestamp ";
$cmd .= " --latency-wait 600 ";

# if run on the submit node of biocluster3, run as cluster job 
my $host = `hostname`;
if($host=~/^login-0-0/){
	if(!$ENV{'SGE_ROOT'}){
		print STDERR "$script: Started snakemake on submit node $host but SGE_ROOT is not set.\n";
		exit(1);	
	}	
	
	# use DRMAA for job control
	if($use_drmaa){
		$cmd = "export DRMAA_LIBRARY_PATH=".$ENV{'SGE_ROOT'}."/lib/linux-x64/libdrmaa.so;".$cmd;
		$cmd .= " --drmaa \" -b n -o $logdir -e $logdir ";
	}

	# use qsub for job control, add a 'sleep 10' to give the queue a chance to update after each job
	else{
		$cmd .= " --cluster \"sleep 10;qsub -b n -o $logdir -e $logdir ";
	}
	$cmd .= " -q ngs.q ";

	# pass remaining cluster configuration
	if($cluster_configfile){
		$cmd .= " -S {cluster.shell} -q ngs.q -l h_rt={cluster.runtime} -l mem_free={cluster.memory} -l {cluster.other_resources} -pe smp {cluster.cpu} \" --cluster-config $cluster_configfile";
	}else{
		$cmd .= "\"";
	}
	
	# pass custom job script
	if($jobscript){
		$cmd .= " --jobscript $jobscript";
	}	
}

# number of jobs
$cmd .= " --jobs $jobs";

# add additional settings for snakemake
$cmd .= " $snakemake_params";

# redirect log
if($snakefile && !$dryrun){
	$cmd .= " 2>&1 | tee $logdir/snakemake.log";
}

# run snakemake command
print "Snakemake call: $cmd\n";
print "==============================================================================================================\n\n";
my $err = system("$cmd");


# post-processing
if($snakefile && !$dryrun && -e "$logdir/snakemake.log"){

	# remove color codes from snakemake log
	system('sed -i -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g" '."$logdir/snakemake.log");
	
	# send mail if requested
	if($send_mail){
		my $subject = "snakemake ".basename($snakefile)." finished";


		# snakemake can fail in advance due to syntax errors of the snakemake file
		# I do not need an email for that therefore a minimum of 60 seconds is required for sending an email
		my $end_time = time();
		if($end_time - $start_time > 60){
			system("mail -s \"$subject\" $send_mail < $logdir/snakemake.log");
		}
	}
}


sub help{
	print "==============================================================================================================\n";
	print "These are the settings for the wrapper script.\n";
	print "==============================================================================================================\n\n";
	print "--help|-h				Show available arguments for the wrapper script\n";
	print "--full_help				Show all available arguments (wrapper script and Snakemake)\n";
	print "--workflow workflow			Snakemake workflow to run - select from the list below\n";
	print "--snakefile|-s file.snakemake		Snakemake file to run\n";
	print "--description 				Print a short description about what the Snakemake workflow/file does\n";
	print "--dryrun 				Do not run the Snakemake file, just print the commands\n";
	print "--configfile file			Configuration file needed to run the snakemake file\n";
	print "--cluster-config|-u file.json 		Configuration file to specify required cluster resources\n";
	print "--cluster-log-dir			Directory where to save cluster stdout and stderr messages\n";
	print "--use_drmaa 				Use DRMAA for submitting/controlling cluster jobs\n";
	print "--use_qsub 				Use standard qsub for submitting/controlling cluster jobs (default)\n";
	print "--jobs|-j number 			Number of jobs to run in parallel (1)\n";
	print "--mail test\@test.de			Send a mail when snakemake finishes\n";
	print "\n";
	print "The following Snakemake workflows are available:\n";
	
	my @files = glob("$snakemake_workflow_dir/*snakemake");
	
	if(length(@files)>0){
		foreach my $f (@files){
			$f = basename($f);
			$f =~ s/\.snakemake$//;
			print "     $f\n";
		}
	}else{
		print "     none\n";
	}
	print "\n";

	exit(0);
}


sub snakemake_help{
	print "==============================================================================================================\n";
	print "These are the settings available for snakemake. All specified settings are passed directly to the program call\n";
	print "and will overwrite any settings made by the wrapper script.\n";
	print "==============================================================================================================\n\n";
	system("module load apps/snakemake;$snakemake_bin --help");
	exit(0);
}

__END__

