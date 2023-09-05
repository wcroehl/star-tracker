#!/opt/perl/bin/perl -w

$sys = 0 ;

$DATA_DIR   = 'C:/01_compile_ver3';
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
$sys = system("${BIN_DIR}/atti06.exe $DATA_DIR $IN_DIR $OU_DIR > ${OU_DIR}/atti.rul") ;
# exit(100) ;

print $DATA_DIR, 'quat_true_six.c', "\n" ;
$sys = system("${BIN_DIR}/quat_true_six.exe $IN_DIR $OU_DIR > ${OU_DIR}/quat_true_six.rul") ;
# exit(100) ;

print $DATA_DIR, 'quat_true_quest.c', "\n" ;
$sys = system("${BIN_DIR}/quat_true_quest.exe $IN_DIR $OU_DIR > ${OU_DIR}/quat_true_quest.rul") ;
# exit(100) ;

# Plot results using IDL
$pid = open(WRITEME, '>', 'C:/Users/timot/Desktop/EKF_c/test.png') or die "Couldn't fork: $!\n";
print WRITEME "!path = '$IDL_DIR:' + !path \n" ;
print WRITEME "quat_true_quest, $ts, $tf, '$OU_DIR'\n" ;  
print WRITEME "quat_true_six, $ts, $tf, '$OU_DIR'\n" ;  
close WRITEME ;
print "PERL> pid = $pid\n" ;
