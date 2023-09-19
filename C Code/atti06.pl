#!/opt/perl/bin/perl -w

$sys = 0 ;

$DATA_DIR   = 'D:/StarTracker-20230706T212348Z-001/StarTracker/01_compile_ver3/01_compile_ver3' ; 
$SOURCE_DIR = $DATA_DIR . '/att' ;
$IDL_DIR = $SOURCE_DIR ;

$IN_DIR = $DATA_DIR . '/data' ;
$OU_DIR = $DATA_DIR . '/output' ;
unless (-e $OU_DIR || -d $OU_DIR || -w $OU_DIR) {
     system("mkdir $OU_DIR") ; }

$BIN_DIR = $OU_DIR ;

$ts = 0.0 ;
$tf = 5670.0 ;

print $DATA_DIR, "\n", "INPUT: ", $IN_DIR, "\n", "OUPUT: ", $OU_DIR, "\n" ;
system("gcc -g ${SOURCE_DIR}/atti06.c -lm -o${BIN_DIR}/atti06.exe") ;
$sys = system("${BIN_DIR}/atti06.exe $DATA_DIR $IN_DIR $OU_DIR > ${OU_DIR}/atti.rul") ;
# exit(100) ;

print $DATA_DIR, 'quat_true_six.c', "\n" ;
system("gcc -g ${SOURCE_DIR}/quat_true_six.c -lm -o${BIN_DIR}/quat_true_six.exe") ;
$sys = system("${BIN_DIR}/quat_true_six.exe $IN_DIR $OU_DIR > ${OU_DIR}/quat_true_six.rul") ;
# exit(100) ;

print $DATA_DIR, 'quat_true_quest.c', "\n" ;
system("gcc -g ${SOURCE_DIR}/quat_true_quest.c -lm -o${BIN_DIR}/quat_true_quest.exe") ;
$sys = system("${BIN_DIR}/quat_true_quest.exe $IN_DIR $OU_DIR > ${OU_DIR}/quat_true_quest.rul") ;
# exit(100) ;

# Plot results using IDL
$pid = open(WRITEME, '>', 'C:/Users/timot/Desktop/EKF_c/test.png') or die "Couldn't fork: $!\n";
print WRITEME "!path = '$IDL_DIR:' + !path \n" ;
print WRITEME "quat_true_quest, $ts, $tf, '$OU_DIR'\n" ;  
print WRITEME "quat_true_six, $ts, $tf, '$OU_DIR'\n" ;  
close WRITEME ;
print "PERL> pid = $pid\n" ;

