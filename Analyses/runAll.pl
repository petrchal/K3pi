#!/usr/bin/perl
use strict;
use warnings;

# Input arrays
#my @DCA = (0.5, 1., 1.5, 2,2.5);  # Example DCA values, modify as needed
#my @nhits = (10, 15, 20, 25);    # Example nhits values, modify as needed

my @DCA = (0.5, 1., 2,);  # Example DCA values, modify as needed
my @nhits = (10,  20);    # Example nhits values, modify as needed


# Header file name
my $header_file = 'inputs.h';

# Total number of combinations for progress calculation
my $total_combinations = scalar(@DCA) * scalar(@nhits) * 2; # Multiply by 2 for left and right versions
my $current_progress = 0;

# Loop over all combinations of DCA and nhits
foreach my $dca (@DCA) {
    foreach my $hit (@nhits) {
        foreach my $side ("left", "right") {
            # Update progress
            $current_progress++;
            my $progress_percentage = sprintf("%.2f", ($current_progress / $total_combinations) * 100);
            print "Processing combination $current_progress/$total_combinations ($progress_percentage%)...\n";

            # Determine the macro definition
            my $macro = $side eq "left" ? "#define _leftSIDE" : "#define _rightSIDE";

            # Create a new header file with the current DCA, nhits, and macro values
            open(my $header_fh, '>', $header_file) or die "Cannot create $header_file: $!";
            print $header_fh "$macro\n";
            print $header_fh "const float c_DCA=$dca;\n";
            print $header_fh "const int c_nhits=$hit;\n";
            close($header_fh);

            # Build the command to execute ROOT
            my $log_file = "DCA${dca}_nhits${hit}_${side}Vz.log";
            my $command = "root -b -q plotEfficiencies.C++ > $log_file 2>&1";

            # Execute the command
            print "Executing: $command\n";
            system($command) == 0 or warn "Command failed: $!";

            # Clean up the header file (optional, depending on your needs)
            unlink $header_file if -e $header_file;
        }
    }
}

print "All combinations processed successfully.\n";
