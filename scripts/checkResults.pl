#!/usr/bin/perl -w

use strict;
use File::Basename;
use File::Copy;
use Getopt::Std;
use Switch;
#locking code

my $resultLog="checkSanity.log";

#use Cwd qw(cwd getcwd);
#my $workDir = cwd();


#my %opt;
#getopts('drtmsilcn',\%opt);
#print %opt;print "\n";
#my $usage = "usage: -drtmsilcn  dataType[Y12_UU_All, Y12_UU_Cent, Y12_UU_MinB] [outDir] [scratchDir] \n";
# read parameters

#if ($opt{l}){
#   if ($#ARGV != 2) {print "Three input parameters are needed when using -l option ($#ARGV): \n";die $usage;}
##}
#else{
#  if ($#ARGV < 1) {print "Two input parameters are needed ($#ARGV): \n";die $usage;}
#}
my $script="./check.C";
if ((-e $script)) {`rm $script`;}
open(SC,">$script") or die "$script: $!\n";
  print SC <<EOF;
#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TObjString.h"
   
 TObjArray *makeList(const char *file_list){ 

   TObjArray* FileList = new TObjArray();
  //  FileList->Add( TFile::Open("hsimple1.root") );
  //  FileList->Add( TFile::Open("hsimple2.root") );
  
    ifstream fList(file_list);
    if (!fList)
     {
       cout << "!!! Can't open file " << myList << endl;
       return NULL;
     }
  
    char lineFromFile[255];
    while(fList.getline(lineFromFile, 250)){
      cout << "Opening: " << lineFromFile << endl;
      //  
      auto str=new TObjString(lineFromFile);
      FileList->Add(str);
    }
    return FileList;
  }

void check(const char *list){
       auto fileList=makeList(list);
       assert(fileList);
       cout<<" entries="<<fileList->GetEntriesFast()<<endl;

      TIter next(fileList);
      TObjString* strO;
      int i=0;
      while ((strO = (TObjString*)next())){
        i++;
        TString str=strO->String();
        TFile* f = new TFile(str,"read");
        if (f->IsZombie()){
          cout<<str<<" ZOMBIE"<<endl;
          delete f; continue;
        }
       auto tree = (TTree*)f->Get("kaons");
       cout<<str<<" ENTRIES="<< tree->GetEntriesFast()<<endl; 
       f->Close(); delete f;
      }
     }
EOF

 
my %results;
my @fileList=@ARGV;
foreach (@fileList){
    my $file=$_;
    if ((-e $resultLog)) {`rm $resultLog`;}
    my $cmd="nohup root.exe -b -q '${script}(\"${file}\")' &> $resultLog &";
    print  "$cmd \n";
    #my $res=
    `$cmd`;
    #print "$res \n";

  }
 
 #`rm $script`;
