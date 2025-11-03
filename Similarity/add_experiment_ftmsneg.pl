#!/usr/local/bin/perl -w

use strict;

#Input files needed to fill the database
my $expfile=$ARGV[0];
my $nodefile=$ARGV[1];
my $ms2scanfile=$ARGV[2];
my $ms3scanfile=$ARGV[3];
my $ms4scanfile=$ARGV[4];
my $ms5scanfile=$ARGV[5];

#Convert files to unix format, just in case
`dos2unix $nodefile`;
if(defined($ms2scanfile)){`dos2unix $ms2scanfile`;}
if(defined($ms3scanfile)){`dos2unix $ms3scanfile`;}
if(defined($ms4scanfile)){`dos2unix $ms4scanfile`;}
if(defined($ms5scanfile)){`dos2unix $ms5scanfile`;}

#Some general parameters
my $m=0.6;
my $n=3;
my $timewindow_in_sec=20;
my $mw=0.005;
my $tw=$timewindow_in_sec/60.0;

#Read in the database and find the highest ids for expid, compid, specid, ms2id, ms3id, ms4id and ms5id
my $db_expid=0;
my $db_compid=0;
my $db_specid=0;
my $db_ms2id=0;
my $db_ms3id=0;
my $db_ms4id=0;
my $db_ms5id=0;

my @dbfiles=("CSV/experiment.csv","CSV/compound.csv","CSV/spectra.csv","CSV/MS2spectra.csv","CSV/MS3spectra.csv","CSV/MS4spectra.csv","CSV/MS5spectra.csv");
my @ids=(\$db_expid,\$db_compid,\$db_specid,\$db_ms2id,\$db_ms3id,\$db_ms4id,\$db_ms5id);
my @positions=(1,1,1,2,2,2,2);

for(my $i=0;$i<scalar(@ids);$i++){
  update($dbfiles[$i],$ids[$i],$positions[$i]);
}

#Read all retention times from nodefile (aligned spectra)
my %cmp_rtime;
my %cmp_mass;
my %cmp_ms2idlist;
my %cmp_ms3idlist;
my %cmp_ms4idlist;
my %cmp_ms5idlist;

open(IN,"<$nodefile");
#Skip header
$_=<IN>;
while(<IN>){
  my $line=$_;
  chomp($line);
  my @vals=split/\t/,$line;
  $cmp_rtime{$vals[2]}=$vals[0];
  $cmp_mass{$vals[2]}=$vals[1];
}
close IN;

#Arrays of hashes for all the MSX scans
my @msx_rtime;
my @msx_mpeak_intensity;
my @msx_mpeak_mass;
my @msx_parent_ion_mass;
my @msx_parent_spec_id;
my @msx_cmp_mass;
my @msx_dpeak_mass;
my @msx_dpeak_intensity;
my @msx_nextmsx_idlist;

#Parse all MSX scans
if(defined($ms2scanfile)){parse_MSX_scanfile($ms2scanfile,2);}
if(defined($ms3scanfile)){parse_MSX_scanfile($ms3scanfile,3);}
if(defined($ms4scanfile)){parse_MSX_scanfile($ms4scanfile,4);}
if(defined($ms5scanfile)){parse_MSX_scanfile($ms5scanfile,5);}

#Open all datafiles for appending:
open(EXP,">>".$dbfiles[0]);
open(COMP,">>".$dbfiles[1]);
open(SPEC,">>".$dbfiles[2]);
open(MS2,">>".$dbfiles[3]);
open(MS3,">>".$dbfiles[4]);
open(MS4,">>".$dbfiles[5]);
open(MS5,">>".$dbfiles[6]);

open(MATCH,">matches.txt");

#Store experiment info in the database
open(IN,"<$expfile");
$_=<IN>;
my @expvals=split/\t/,;
print EXP "$db_expid\t".join("\t",@expvals);
close(IN);

#Check until where to go when compering current compounds to add to existin ones in the database
my $prevmax_compid=$db_compid-1;

#Now recursively process all compounds that have an MS2 spectrum associated
foreach my $cmp(keys %cmp_mass){
    
  #Make the mappings from the internal ids to the database ids
  my %prog2dbidms2;
  my %prog2dbidms3;
  my %prog2dbidms4;  
  my %prog2dbidms5;  
    
  #Check number of MS2 spectra
  if(defined($cmp_ms2idlist{$cmp})){
    my $num_ms2=scalar(@{$cmp_ms2idlist{$cmp}});  
    #When multiple MS2 spectra are found, take the one closest in retention time
    my $closesttime=10000000;
    my $closestid=-1;
    if($num_ms2 >1){
      my $crtime=$cmp_rtime{$cmp};  
      for(my $i=0;$i<$num_ms2;$i++){
        my $curms2rtime=$msx_rtime[2]{${$cmp_ms2idlist{$cmp}}[$i]};  
        if(abs($crtime-$curms2rtime)<$closesttime){
          $closesttime=abs($crtime-$curms2rtime);
          $closestid=${$cmp_ms2idlist{$cmp}}[$i];    
        }
      }
      #print "Closest id = $closestid\n";
      #Keep only the closest MS2 spectrum and remove the other one
      @{$cmp_ms2idlist{$cmp}}=($closestid);
    }

    #Store this compound in the compound table
    store_compound($cmp);
    print MATCH "COMPID $db_compid [$cmp]";
    compare_compound_to_database($cmp);

    #Process the MS2 spectrum
    foreach my $ms2id(@{$cmp_ms2idlist{$cmp}}){
      $prog2dbidms2{$ms2id}=$db_ms2id;        
      print "$ms2id:\n";  
      #Check which MS3 spectra linked to this cmp can be linked to this MS2 spectrum
      my @validatedms3=();
      foreach my $ms3id(@{$cmp_ms3idlist{$cmp}}){
        print "$ms3id\t";
        my $cur_p_ion=round($msx_parent_ion_mass[3]{$ms3id});
        my $found=-1;
        for(my $i=0;$i<scalar(@{$msx_dpeak_mass[2]{$ms2id}});$i++){
          my $curdp_mass=${$msx_dpeak_mass[2]{$ms2id}}[$i];
          if(abs($curdp_mass-$cur_p_ion) < $mw){
            #Match found
            $found=$i;            
          }
        }
        if($found > -1){
          print "validated [$cur_p_ion -> ".(${$msx_dpeak_mass[2]{$ms2id}}[$found])."]\n";  
          push(@validatedms3,$ms3id);
          $prog2dbidms3{$ms3id}=$db_ms3id;
          $db_ms3id++;
          push(@{$msx_nextmsx_idlist[2]{$ms2id}},$ms3id);
          $msx_parent_spec_id[3]{$ms3id}=$ms2id;
        }
        else{
          print "$ms3id not validated\n";    
        }
      }
      #Update cmp list
      @{$cmp_ms3idlist{$cmp}}=@validatedms3;   
      
      #Add to the MS2SPECTRA table
      print MS2 "$db_compid\t$db_ms2id\t";
      #print parent ion
      print MS2 "".($msx_parent_ion_mass[2]{$ms2id})."\t";
      #print MS2 peak list
      print MS2 "".(join(",",@{$msx_dpeak_mass[2]{$ms2id}}))."\t";
      print MS2 "".(join(",",@{$msx_dpeak_intensity[2]{$ms2id}}))."\t";   
      if(scalar(@validatedms3) > 0){
        print MS2 "".(join(",",map($prog2dbidms3{$_},@validatedms3)))."\n";
      }
      else{
        print MS2 "NULL\n";
      }  
      $db_ms2id++;    
    }
    
    #Process the MS3 spectra
    foreach my $ms3id(@{$cmp_ms3idlist{$cmp}}){
      print "$ms3id:\n";  
      #Check which MS4 spectra linked to this cmp can be linked to this MS3 spectrum
      my @validatedms4=();
      foreach my $ms4id(@{$cmp_ms4idlist{$cmp}}){
        print "$ms4id\t";
        my $cur_p_ion=round($msx_parent_ion_mass[4]{$ms4id});
        my $found=-1;
        for(my $i=0;$i<scalar(@{$msx_dpeak_mass[3]{$ms3id}});$i++){
          my $curdp_mass=${$msx_dpeak_mass[3]{$ms3id}}[$i];
          if(abs($curdp_mass-$cur_p_ion) < $mw){
            #Match found
            $found=$i;            
          }
        }
        if($found > -1){
          print "validated [$cur_p_ion -> ".(${$msx_dpeak_mass[3]{$ms3id}}[$found])."]\n";  
          push(@validatedms4,$ms4id);
          $prog2dbidms4{$ms4id}=$db_ms4id;
          print "prog2dbidms4[$ms4id]=$db_ms4id\n";          
          $db_ms4id++;
          push(@{$msx_nextmsx_idlist[3]{$ms3id}},$ms4id);
          $msx_parent_spec_id[4]{$ms4id}=$ms3id;
        }
        else{
          print "$ms4id not validated\n";    
        }
      }
      #Update cmp list
      @{$cmp_ms4idlist{$cmp}}=@validatedms4;    

      #Add to the MS3SPECTRA table
      my $curid=$prog2dbidms3{$ms3id};
      print MS3 "$db_compid\t$curid\t";
      #print parent ion
      print MS3 "".($msx_parent_ion_mass[3]{$ms3id})."\t";
      #print MS3 peak list
      print MS3 "".(join(",",@{$msx_dpeak_mass[3]{$ms3id}}))."\t";
      print MS3 "".(join(",",@{$msx_dpeak_intensity[3]{$ms3id}}))."\t";   
      if(scalar(@validatedms4) > 0){
        print MS3 "".(join(",",map($prog2dbidms4{$_},@validatedms4)))."\n";
      }
      else{
        print MS3 "NULL\n";
      }
    }

    #Process the MS4 spectra
    foreach my $ms4id(@{$cmp_ms4idlist{$cmp}}){
      print "$ms4id:\n";  
      #Check which MS5 spectra linked to this cmp can be linked to this MS4 spectrum
      my @validatedms5=();
      foreach my $ms5id(@{$cmp_ms5idlist{$cmp}}){
        print "$ms5id\t";
        my $cur_p_ion=round($msx_parent_ion_mass[5]{$ms5id});
        my $found=-1;
        for(my $i=0;$i<scalar(@{$msx_dpeak_mass[4]{$ms4id}});$i++){
          my $curdp_mass=${$msx_dpeak_mass[4]{$ms4id}}[$i];
          if(abs($curdp_mass-$cur_p_ion) < $mw){
            #Match found
            $found=$i;            
          }
        }
        if($found > -1){
          print "validated [$cur_p_ion -> ".(${$msx_dpeak_mass[4]{$ms4id}}[$found])."]\n";  
          push(@validatedms5,$ms5id);
          $prog2dbidms5{$ms5id}=$db_ms5id;          
          $db_ms5id++;
          push(@{$msx_nextmsx_idlist[4]{$ms4id}},$ms5id);
          $msx_parent_spec_id[5]{$ms5id}=$ms4id;
        }
        else{
          print "$ms5id not validated\n";    
        }
      }
      #Update cmp list
      @{$cmp_ms5idlist{$cmp}}=@validatedms5;      
      
      #Check the code in comments ! Still some strange bug with $curid
      
      #Add to the MS4SPECTRA table
      #print "Here [$ms4id]\n";
      #my $curid=$prog2dbidms4{$ms4id};
      #print "".($prog2dbidms4{$ms4id})."\n";
      #print MS4 "$db_compid\t$curid\t";
      #print parent ion
      #print MS4 "".($msx_parent_ion_mass[4]{$ms4id})."\t";
      #print MS4 peak list
      #print MS4 "".(join(",",@{$msx_dpeak_mass[4]{$ms4id}}))."\t";
      #print MS4 "".(join(",",@{$msx_dpeak_intensity[4]{$ms4id}}))."\t";   
      #if(scalar(@validatedms5) > 0){
      #  print MS4 "".(join(",",map($prog2dbidms5{$_},@validatedms5)))."\n";
      #}
      #else{
      #  print MS4 "NULL\n";    
      #}
      
      #Add to the MS5SPECTRA table
      for(my $ms5num=0;$ms5num<scalar(@validatedms5);$ms5num++){
        my $curms5id=$validatedms5[$ms5num];  
        my $curid=$prog2dbidms5{$curms5id};
        print MS5 "$db_compid\t$curid\t";
        #print parent ion
        print MS5 "".($msx_parent_ion_mass[5]{$curms5id})."\t";
        #print MS5 peak list
        print MS5 "".(join(",",@{$msx_dpeak_mass[5]{$curms5id}}))."\t";
        print MS5 "".(join(",",@{$msx_dpeak_intensity[5]{$curms5id}}))."\n";   
      }
      
    }
    #Store spectra of this compound in the SPECTRA table
    my $sumarr=0;
    if(defined($cmp_ms2idlist{$cmp})){$sumarr+=scalar(@{$cmp_ms2idlist{$cmp}});}
    if(defined($cmp_ms3idlist{$cmp})){$sumarr+=scalar(@{$cmp_ms3idlist{$cmp}});}
    if(defined($cmp_ms4idlist{$cmp})){$sumarr+=scalar(@{$cmp_ms4idlist{$cmp}});}
    if(defined($cmp_ms5idlist{$cmp})){$sumarr+=scalar(@{$cmp_ms5idlist{$cmp}});}
    
    if($sumarr > 0){
        print SPEC "$db_compid\t";
        if(defined($cmp_ms2idlist{$cmp})){
        if(scalar(@{$cmp_ms2idlist{$cmp}}) > 0){
            print SPEC "".(join(",",map($prog2dbidms2{$_},@{$cmp_ms2idlist{$cmp}})))."\t";
        }
        else{
            print SPEC "NULL\t";    
        }
        }
        else{print SPEC "NULL\t";}    
    
        if(defined($cmp_ms3idlist{$cmp})){
            if(scalar(@{$cmp_ms3idlist{$cmp}}) > 0){
            print SPEC "".(join(",",map($prog2dbidms3{$_},@{$cmp_ms3idlist{$cmp}})))."\t";
        }
        else{
            print SPEC "NULL\t";    
        }
        }
        else{print SPEC "NULL\t";}    
    
        #if(defined($cmp_ms4idlist{$cmp})){
        #if(scalar(@{$cmp_ms4idlist{$cmp}}) > 0){
        #    print SPEC "".(join(",",map($prog2dbidms4{$_},@{$cmp_ms4idlist{$cmp}})))."\t";
        #}
        #else{
        #    print SPEC "NULL\t";    
        #}
        #}
        #else{print SPEC "NULL\t";}    
    
        if(defined($cmp_ms5idlist{$cmp})){
        if(scalar(@{$cmp_ms5idlist{$cmp}}) > 0){
            print SPEC "".(join(",",map($prog2dbidms5{$_},@{$cmp_ms5idlist{$cmp}})))."\n";
        }
        else{
            print SPEC "NULL\n";    
        }
        }
        else{print SPEC "NULL\n";}    
    }
      
    $db_compid++;
  }
}

#Close database file
close(EXP);
close(COMP);
close(SPEC);
close(MS2);
close(MS3);
close(MS4);
close(MS5);

close(MATCH);

sub parse_MSX_scanfile {
    (my $scanfile,my $msx)=@_;
  
    my $numscans;
    my $firstscan;
    my $lastscan;

    open(IN,"<$scanfile");
    while(<IN>){
      my $line=$_;
      if($line=~/^RunHeaderInfo.*/){
        #Look for scan numbers 
        while($line!~/^first\_scan.*/){
          $_=<IN>;   
          $line=$_;
        }
        if($line=~/^first_scan.*/){
          $line=~/^first_scan = (\d+), last_scan = (\d+),.*/;
          $firstscan=$1;
          $lastscan=$2;
          #print "$firstscan -> $lastscan\n";      
        }
        else{
          die "Could not find scan info\n";   
        }
      }
      if($line=~/ScanHeader # (\d+).*/){
        my $id=$msx."_".$1;
        my $numreads=0;
        #print "$id\n";

        #First look for mother peak retention time
        while($line!~/^start_time.*/){
          $_=<IN>;   
          $line=$_;
        }
        if($line=~/^start_time.*/){
          $line=~/^start_time = (\d+\.*\d*),.*/;
          my $time=$1;
          $msx_rtime[$msx]{$id}=$time;
        }
        else{
          die "Could not find retention time info\n";   
        }

        while($line!~/^num_readings.*/){
          $_=<IN>;   
          $line=$_;
        }
        if($line=~/^num_readings.*/){
          $line=~/^num_readings = (\d+),.*/;
          $numreads=$1;
        }
        else{
          die "Could not find scan info\n";   
        }
                
        #Next look for the mother peak mass and intensity
        while($line!~/^uScanCount.*/){
          $_=<IN>;   
          $line=$_;
        }
        if($line=~/^uScanCount.*/){
          $line=~/^uScanCount = (\d+), PeakIntensity = (\d+\.*\d*), PeakMass = (\d+\.*\d*).*/;
          my $intensity=$2;
          $msx_mpeak_intensity[$msx]{$id}=$intensity;
          my $mass=$3;
          $msx_mpeak_mass[$msx]{$id}=$mass;
          #print "".($mpeak_time{$id})."\t$mass -> $intensity\n";      
        }
        else{
          die "Could not find scan info\n";   
        }
        
        while($line!~/^Precursor.*/){
          $_=<IN>;   
          $line=$_;
        }
        if($line=~/^Precursor.*/){
          my $mass;
          my $cmass;          
          if($msx==2){  
            $line=~/^Precursor Mass  (\d+\.*\d*).*/;
            $mass=$1;
            $cmass=$1;
            $msx_parent_ion_mass[$msx]{$id}=$cmass;            
          }
          elsif($msx==3){
            $line=~/^Precursor Mass  (\d+\.*\d*)  (\d+\.*\d*).*/;
            $mass=$1;
            $cmass=$2;
            $msx_parent_ion_mass[$msx]{$id}=$cmass;
          }
          elsif($msx==4){
            $line=~/^Precursor Mass  (\d+\.*\d*)  (\d+\.*\d*)  (\d+\.*\d*).*/;
            $mass=$1;
            $cmass=$3;
            $msx_parent_ion_mass[$msx]{$id}=$cmass;
          }
          elsif($msx==5){
            $line=~/^Precursor Mass  (\d+\.*\d*)  (\d+\.*\d*)  (\d+\.*\d*)  (\d+\.*\d*).*/;
            $mass=$1;
            $cmass=$4;
            $msx_parent_ion_mass[$msx]{$id}=$cmass;
          }
          $msx_cmp_mass[$msx]{$id}=$mass;
        }
        else{
          die "Could not find scan info\n";   
        }

        #Start reading all the daughter peaks from the ion trap  
        while($line!~/^DataPeaks.*/){
          $_=<IN>;   
          $line=$_;
        }    
        
        my $numtrueread=0;
        
        #Now read $numreads readings
        for(my $i=0;$i<$numreads;$i++){
          #First an empty line
          $_=<IN>;
          #Next the actual packet info
          $_=<IN>;
          my $line=$_;
          $line=~/.*intensity = (\d+\.*\d*), mass\/position = (\d+\.*\d*).*/;
          #my $mass=sprintf("%.2f",$2);
          my $mass=round($2);
          my $intensity=$1;
          if(($intensity >= 100) and ($mass < $msx_cmp_mass[$msx]{$id})){
            push(@{$msx_dpeak_mass[$msx]{$id}},$mass);
            push(@{$msx_dpeak_intensity[$msx]{$id}},$intensity);
            #print "\t$mass -> $intensity\n";
            $numtrueread++;
          }
          #Next line is not used for the moment
          $_=<IN>;
       }
       
       my $foundcompound=find_compound($msx,$id);
       
       if(($numtrueread > 0) and ($foundcompound ne "NULL")){
          #Add this MS2ID to the cmp MS2IDLIST
          if($msx==2){
            push(@{$cmp_ms2idlist{$foundcompound}},$id);
          }
          if($msx==3){
            push(@{$cmp_ms3idlist{$foundcompound}},$id);
          }
          if($msx==4){
            push(@{$cmp_ms4idlist{$foundcompound}},$id);
          }
          if($msx==5){
            push(@{$cmp_ms5idlist{$foundcompound}},$id);
          }
          #print "Found $msx\n";
       }
       else{
          delete $msx_rtime[$msx]{$id};
          delete $msx_mpeak_intensity[$msx]{$id};
          delete $msx_mpeak_mass[$msx]{$id};
          delete $msx_parent_ion_mass[$msx]{$id};
          delete $msx_cmp_mass[$msx]{$id};
          delete $msx_dpeak_mass[$msx]{$id};
          delete $msx_dpeak_intensity[$msx]{$id};
       }
      }
    }
    close IN;

    #Eliminate duplicates that arose through rounding, keep the highest intensity peak
    foreach my $key(keys %{$msx_dpeak_mass[$msx]}){
       my %tmphash;
       for(my $i=0;$i<scalar(@{$msx_dpeak_mass[$msx]{$key}});$i++){
         my $mass=${$msx_dpeak_mass[$msx]{$key}}[$i];
         my $intensity=${$msx_dpeak_intensity[$msx]{$key}}[$i];
         if(defined($tmphash{$mass})){
            if($intensity > $tmphash{$mass}){
               $tmphash{$mass}=$intensity;
            }
         }
         else{
           $tmphash{$mass}=$intensity; 
         }
       }
       @{$msx_dpeak_mass[$msx]{$key}}=();
       @{$msx_dpeak_intensity[$msx]{$key}}=();
       foreach my $key2(keys %tmphash){
         push(@{$msx_dpeak_mass[$msx]{$key}},$key2);
         push(@{$msx_dpeak_intensity[$msx]{$key}},$tmphash{$key2});
       }
    }
}

sub round {
    my($number) = shift;
    return int($number + .5);
}

sub update {
   (my $file, my $idx, my $pos) = @_;
   open(IN,"<$file");
   #Skip header
   $_=<IN>;
   #Read rest of file
   my $max=0;
   while(<IN>){
     my $line=$_;
     chomp($line);
     my @vals=split/\t/,$line;
     if($vals[$pos-1] > $max){
       $max=$vals[$pos-1];   
     }
   }
   close(IN);
   $$idx=$max+1;
}

sub find_compound {
  (my $msx,my $id)=@_;

  my $ptime=$msx_rtime[$msx]{$id};
  my $pmass=$msx_cmp_mass[$msx]{$id};
  my $bestname="NULL";
  my $bestdiff=100000;      
      
  foreach my $cmp(keys %cmp_rtime){
    my $ctime=$cmp_rtime{$cmp};
    my $cmass=$cmp_mass{$cmp};
 
    my $diff=abs($ptime-$ctime);
    if(($diff < $tw) and ($diff < $bestdiff) and (abs($pmass-$cmass) < $mw) ){
      $bestname=$cmp;
      $bestdiff=$diff;
    }
  }      
  return($bestname);
}

sub store_compound {
  (my $cmp)=@_;
  my @rowvals=($db_compid,$cmp,"NULL","NULL",$cmp_mass{$cmp},"NULL","NULL",$cmp_rtime{$cmp},$db_expid,"NULL","NULL","NULL","NULL","NULL","NULL","NULL");  
  my $str=join("\t",@rowvals);
  print COMP "$str\n";
}

sub compare_compound_to_database{
  (my $cmp)=@_;
  
  my $bestmatch="NULL";
  my $bestcommon=0;
  my $bestgdp=0;
  
  #Get compound mpeak, masslist and peaklist
  my $mpeak1=$cmp_mass{$cmp};
  my $ms2id=${$cmp_ms2idlist{$cmp}}[0];
  my @masslist1=@{$msx_dpeak_mass[2]{$ms2id}};
  my @peaklist1=@{$msx_dpeak_intensity[2]{$ms2id}};
  
  #print "[$mpeak1]:";
  #foreach my $key (sort(@masslist1)){
  #  print "\t$key";
  #}
  #print "\n";
  
  open(IN,"<".($dbfiles[3]));
  #Skip header
  $_=<IN>;
  while(<IN>){
    my $line=$_;
    chomp($line);
    my @vals=split/\t/,$line;
    if($vals[0]<=$prevmax_compid){
      #Calculate global dot product and global common peaks
      my $mpeak2=$vals[2];
      my @masslist2=split/,/,$vals[3];
      my @peaklist2=split/,/,$vals[4];
      my @result=symmetric_dotproduct_combined($mpeak1, $mpeak2, \@masslist1, \@masslist2, \@peaklist1, \@peaklist2);
      my $curncommon=$result[0];
      my $curgdp=$result[1];
      my $minpeaks=$result[2];
      #print "Common=$curncommon out of $minpeaks, Global dot product=$curgdp\n";
      if($curgdp > 0.8){
        if($curncommon > $bestcommon){
          $bestcommon=$curncommon;
          $bestgdp=$curgdp;
          $bestmatch=$vals[0];  #Compid        
        }
        elsif($curncommon == $bestcommon){
          if($curgdp > $bestgdp){
            $bestgdp=$curgdp;
            $bestmatch=$vals[0];
          }
        }
      }  
    }
    else{
      last;   
    }
  }
  close IN;
  if($bestmatch eq "NULL"){
    print MATCH "\tNo match\n";
  }
  else{
    print MATCH "\tBest match [$bestmatch]: Globalcommon=$bestcommon Global similarity=$bestgdp\n";
  }
  #die;
}


sub symmetric_dotproduct_combined{
  (my $mpeak1,my $mpeak2, my $mslist1, my $mslist2, my $pklist1, my $pklist2)=@_;   
  
  my @masslist1=@{$mslist1};
  my @masslist2=@{$mslist2};
  my @peaklist1=@{$pklist1};
  my @peaklist2=@{$pklist2};
  
  my $minpeaks=scalar(@masslist1);
  if(scalar(@masslist2)<$minpeaks){
    $minpeaks=scalar(@masslist2);   
  }

  #@masslist1=(200,400,370);
  #@masslist2=(200,400,370);
  #@peaklist1=(100,500,1000);
  #@peaklist2=(1000,5000,10000);

  #my @peaks1=sort {$a<=>$b} @masslist1;
  #my @peaks2=sort {$a<=>$b} @masslist2;
  #print "Orig 1:\t";
  #for(my $i=0;$i<scalar(@peaks1);$i++){
  #   print " ".($peaks1[$i]);
  #}
  #print "\n";
  #print "Orig 2:\t";
  #for(my $i=0;$i<scalar(@peaks2);$i++){
  #   print " ".($peaks2[$i]);
  #}
  #print "\n";

  #Next print the same, but for neutral losses
  #print "Neutral loss 1:\t";
  #for(my $i=0;$i<scalar(@peaks1);$i++){
  #   #print " ".(round($ldpeak_mass{$matches{$p1}})-$peaks1[$i]);
  #   print " ".(round($mpeak1)-$peaks1[$i]);
  #}
  #print "\n";
  #print "Neutral loss 2:\t";
  #for(my $i=0;$i<scalar(@peaks2);$i++){
  #   #print " ".(round($ldpeak_mass{$matches{$p2}})-$peaks2[$i]);
  #   print " ".(round($mpeak2)-$peaks2[$i]);
  #}
  #print "\n";

    
  #Rescale peaks
  my $Fd_denom1=0;
  my $Fd_denom2=0;  
  for(my $i=0;$i<scalar(@masslist1);$i++){
     $peaklist1[$i]=($peaklist1[$i]**$n)*($masslist1[$i]**$m);
     $Fd_denom1+=($peaklist1[$i])*($peaklist1[$i]);
  }
  for(my $i=0;$i<scalar(@masslist2);$i++){
     $peaklist2[$i]=($peaklist2[$i]**$n)*($masslist2[$i]**$m);
     $Fd_denom2+=($peaklist2[$i])*($peaklist2[$i]);     
  }

#Calculate global common peaks, once using normal peaks, once using neutral losses

 #Calculate traditional part

  #Do the same for the neutral losses
  my @masslist1nl;
  my @masslist2nl;
  foreach my $key(@masslist1){
    push(@masslist1nl,round($mpeak1)-$key);   
  }
  foreach my $key(@masslist2){
    push(@masslist2nl,round($mpeak2)-$key);   
  }
  
  my @intersection=();
  my %isect=();
  foreach my $key(@masslist1nl){
    if(!defined($isect{$key})){
      $isect{$key}=1
    };    
  }
  foreach my $key(@masslist2nl){
    if(defined($isect{$key})){
      $isect{$key}=2;
    }
  }
  foreach my $key(keys %isect){
    if($isect{$key}>1){
      push(@intersection,$key);
      #print "NL $key\n";
    }
  }
  my $globalcommon=scalar(@intersection);

  #Get intersection
  @intersection=();
  %isect=();
  foreach my $key(@masslist1){
    if(!defined($isect{$key})){
      $isect{$key}=1
    };    
  }
  foreach my $key(@masslist2){
    if(defined($isect{$key})){
      $isect{$key}=2;
    }
  }
  foreach my $key(keys %isect){
    if($isect{$key}>1){
      push(@intersection,$key);
      #print "$key\n";
    }
  }
  my $tmpcommon=scalar(@intersection);
  $globalcommon+=$tmpcommon;

  my %inter;
  foreach my $key(@intersection){
    $inter{$key}=1;   
  }
   
  my %mp1=();
  for(my $i=0;$i<scalar(@masslist1);$i++){
    if(defined($inter{$masslist1[$i]})){
       $mp1{$masslist1[$i]}=$peaklist1[$i];
    }
  }
  my %mp2=();
  for(my $i=0;$i<scalar(@masslist2);$i++){
    if(defined($inter{$masslist2[$i]})){
       $mp2{$masslist2[$i]}=$peaklist2[$i];
    }
  }
  if(scalar(keys %mp1) != scalar(keys %mp2)){
     print "MP1:";
     foreach my $key(keys %mp1){
       print " $key";
     }       
     print "\nMP2:";
     foreach my $key(keys %mp2){
       print " $key";
     }            
     print "\n";
     die "Scalars incompatible\n";   
  }

  
  #Now calculate the similarity
  my $FD=0;
  my $Fd_nom=0;
  foreach my $key(keys %mp1){
    $Fd_nom+=$mp1{$key}*$mp2{$key};   
  }
  
  #Remove the common mass-peak pairs
  foreach my $key(keys %mp1){
    for(my $i=0;$i<scalar(@masslist1);$i++){
      if($key==$masslist1[$i]){
        splice(@masslist1,$i,1);
        splice(@peaklist1,$i,1);
        last;
      }
    }
  }
  foreach my $key(keys %mp2){
    for(my $i=0;$i<scalar(@masslist2);$i++){
      if($key==$masslist2[$i]){
        splice(@masslist2,$i,1);
        splice(@peaklist2,$i,1);
        last;
      }
    }
  }

#Calculate neutral losses: only for noncommon peaks of the first part
  my @masslist1b=();
  my @masslist2b=();
  for(my $i=0;$i<scalar(@masslist1);$i++){
    $masslist1b[$i]=round($mpeak1)-$masslist1[$i];
  }
  for(my $i=0;$i<scalar(@masslist2);$i++){
    $masslist2b[$i]=round($mpeak2)-$masslist2[$i];      
  }
  
  #Get intersection
  @intersection=();
  %isect=();

  foreach my $key(@masslist1b){
    if(!defined($isect{$key})){
      $isect{$key}=1
    };    
  }
  foreach my $key(@masslist2b){
    if(defined($isect{$key})){
      $isect{$key}=2;
    }
  }
  foreach my $key(keys %isect){
    if($isect{$key}>1){
      push(@intersection,$key);
      #print "NL2 $key\n";
    }
  }
  $tmpcommon+=scalar(@intersection);
  #print "$tmpcommon peaks in common\n";
  #$cptotal=$tmpcommon;
  %inter=();
  foreach my $key(@intersection){
    $inter{$key}=1;   
  }
  
  %mp1=();
  for(my $i=0;$i<scalar(@masslist1);$i++){
    if(defined($inter{$masslist1b[$i]})){
       $mp1{$masslist1b[$i]}=$peaklist1[$i];
    }
  }
  %mp2=();
  for(my $i=0;$i<scalar(@masslist2);$i++){
    if(defined($inter{$masslist2b[$i]})){
       $mp2{$masslist2b[$i]}=$peaklist2[$i];
    }
  }
  if(scalar(keys %mp1) != scalar(keys %mp2)){
     print "MP1:";
     foreach my $key(keys %mp1){
       print " $key";
     }       
     print "\nMP2:";
     foreach my $key(keys %mp2){
       print " $key";
     }            
     print "\n";
     die "Scalars incompatible\n";   
  }

  
  #Now calculate the similarity
  foreach my $key(keys %mp1){
    $Fd_nom+=$mp1{$key}*$mp2{$key};   
  }
  #print "".($Fd_nom*$Fd_nom)." ".($Fd_denom1*$Fd_denom2)."\n";
  $FD=($Fd_nom*$Fd_nom)/($Fd_denom1*$Fd_denom2);
  
  my @result=($globalcommon,$FD,$minpeaks);
  return(@result);
}

