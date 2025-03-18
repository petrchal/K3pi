#!/usr/bin/perl
use strict;
use warnings;

# Check for job prefix argument
my $job_prefix = shift @ARGV or die "Usage: $0 <job_name_prefix>\n";

# Define file/directory names
my $session_file   = "$job_prefix.session.xml";
my $csh_dir        = "./csh";
my $production_dir = "./production";

# Check that the session file exists
unless (-e $session_file) {
    die "Session file '$session_file' not found!\n";
}

# Find all job scripts in ./csh matching <job_prefix>_*.csh
print "Searching for jobs $csh_dir/${job_prefix}_*.csh \n";
my @job_scripts = glob("$csh_dir/${job_prefix}_*.csh");
unless (@job_scripts) {
    die "No job scripts found in $csh_dir for prefix '$job_prefix'.\n";
}

my $num_found=@job_scripts;
print "Found $num_found scheduled jobs.\n";


# Array to store indices of failed jobs
print "Searching for results $production_dir/kaon_${job_prefix}_*.root \n";
my @failed_indices = ();

foreach my $script (@job_scripts) {
    if ($script =~ /${job_prefix}_(\d+)\.csh$/) {
        my $index = $1;
        #my $result_file = "$production_dir/kaon_${job_prefix}_${index}.root";
        #my $result_file = "$production_dir/kaon_*_${index}.root";
        my @result_files = glob("$production_dir/kaon_*_${index}.root");
        # If the expected result file does not exist, mark this job as failed
        #unless (-e $result_file) {
        unless(@result_files){
            push @failed_indices, $index;
        }
    }
}

my $num_failed = @failed_indices;
if (@failed_indices) {
    print "Detected $num_failed failed jobs \n";# for indices: \n ". join(', ', @failed_indices) . "\n";
    #print "Resubmitting failed jobs using: star-submit -r $session_file\n";
    #print "star-submit -r @failed_indices  $session_file ";
    # Execute the star-submit -r command
    my $exit_status = system("star-submit -r " . join(',', @failed_indices) . " $session_file");
    if ($exit_status != 0) {
        die "Error: star-submit -r command failed with exit code $exit_status\n";
    }
} else {
    print "All jobs appear to have finished successfully. No failed jobs found.\n";
}

