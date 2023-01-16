#!/usr/bin/perl -w

use strict;
use File::Basename;
use File::Copy;
use Getopt::Std; 
use Switch;
use threads;
use threads::shared;
use Sys::Info;
use Sys::Info::Constants qw( :device_cpu );
use Sys::MemInfo qw(totalmem freemem totalswap);
 

my ($rs,@files)=@ARGV;

# remember $# is last index ..size-1
unless ($#files>=0) {
  print "NO  files - quitting\n";
    exit;
    }

my $maxFileSize=200000000; #.2GB
my $info = Sys::Info->new;
my $cpu  = $info->device( CPU => my %options );
my $ncpu= $cpu->count;
print " found ${ncpu} cpus\n";

my $nchilds=$ncpu-1; #keep one for a system :-)
if ($nchilds<2) {$nchilds=2;}

#mem size
my @words = split ' ',  `free -b |grep Mem`;
my $totalMem=$words[1];
my $memavail=$words[3];
my $memrat=100*$memavail/$totalMem;
print "Available memory ".(${memavail}/1000000000)."GB out of ".($totalMem/1000000000)."GB ... ${memrat}% \n";
#claim at least 50% of the memory
if ($memrat<50) { $memavail=0.5*$totalMem;}

#filesizes (output size)
my $fsum = 0.1;
my @flist=@files;
my $fnum=$#flist+1;
 while (@flist) {
       my $f=pop @flist;
       chomp $f;
      $fsum += -s $f;
   }
 print "Sum file sizes=".($fsum/1000000000)."GB \n";
  

my $rescount=0; #final results

   
RunDivision($rs,0,@files);


sub RunDivision{
    my($res,$level,@file)=@_;
    chomp $res;
    my $nn=@file;
    print " entering runDivision: ${res} at level ${level} with ${nn} files\n";# \n @files \n";   

    #optimal file division per child
    my $optimumNfiles=$memavail/(2.0*$fsum/${nn}); #2...include the output file in memory sum
    print "Optimum # of files= ${optimumNfiles} for the  whole memory \n";

    my $nparse= int($optimumNfiles/$nchilds);
   if ($nparse<2) { # at end go for two files only
      $nparse=2;
      $nchilds=int($memavail/(2*$fsum/${nn}));
   }

  #too many files?
  #  if($nchilds * $nparse > ${nn}) {$nparse=int($fnum/$nchilds)+1};
  #$nparse=2; override
 print "Runnig with ${nchilds} threads of $nparse files each. \n";

  #last final file..just copy
    if ($#file==0){
        my $fl=$file[0];
        chomp $fl;
	my $cmd="cp $fl ${res}";
	print "$cmd \n";
	`$cmd`;
	return 0;
    }

    my $outName="${res}_l${level}t";
    my $count=0;
    my @childs;
    while (@file){
	my $n=0;
	my $list=" ";
	while( @file && ($n<$nparse)){
	    my $f=pop @file;
	    chomp $f;
	    my $size=-s $f;
	    #print " $f size is $size \n";
	    if($size>$maxFileSize){
	      print "$f is above max file size ... freezing \n";
	      my $cmd="cp $f _${rs}$rescount.root";
	      $rescount++;
	      runCmd($cmd);
	      next; 
	    }
	    $list="$list $f";
	    $n++;
	} #while
	my $cmd="";
	if ($n==0) {next;}
	if ($n==1){
	    $cmd="cp $list ${outName}$count.root";
	}
	if ($n>1){ 
            $cmd="hadd -k -f ${outName}$count.root $list";
       }
	my $pid = fork();
        if ($pid) {
	    # parent
	    #print "pid is $pid, parent $$\n";
	    push(@childs, $pid);
        } elsif ($pid == 0) {
	    # child
	    runCmd($cmd);
	    exit 0;
        } else {
	    die "couldnt fork: $!\n";
        }

	if (@childs == $nchilds){
	    foreach (@childs) {
		my $tmp = waitpid($_, 0);
		print "done with pid $tmp\n";
	    }
	    @childs=();
	}

	$count++;
    } #while
    
    foreach (@childs) {
	my $tmp = waitpid($_, 0);
	print "done with pid $tmp\n";
    }
    
    my @newlist=`ls ${outName}*.root`;
    if (@newlist>0) {
      RunDivision("${outName}.root",$level+1,@newlist);
      my $cmd="mv ${outName}.root ${res}";
      print "$cmd \n";
      `${cmd}`;
      `rm ${outName}*.root`;
     }
}

sub runCmd{
    my $cmd=shift;
    print "$cmd \n";
    my $tmp=`$cmd`;
    print "$tmp \n";
}
