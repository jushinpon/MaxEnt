@files = <00MC*.data>;

#foreach my $i (@files){
#unlink "$i";
#print "Kill $i file\n";
#}

unlink "./maxent.exe";
$temp = system("gfortran -fopenmp -O3 -o maxent.exe 00MAXENT_main.f90");
#$temp = system("gfortran -O3 -o maxent.exe 00MAXENT_main.f90");
#$temp = system("gfortran -g -fcheck=all -Wall 00MAXENT_main.f90");
 die "Compiling failed" if($temp);
#sleep(1);
$temp = './maxent.exe';#.' > 00printout.txt';
print $temp;
system($temp);
