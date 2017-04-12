#!/usr/bin/perl -w

use strict;

use Cwd qw(getcwd);
use File::Basename;
use File::Copy qw(copy);
use List::Util qw[min max];

#change directory to script location
chdir(dirname(__FILE__));
use lib dirname(__FILE__);

use gmxFile::Path qw(remove_tree);
use List::Util qw(sum);

#disable quotes as they could screw up pattern matching
$ENV{GMX_NO_QUOTES}='NO';
#disable backups because otherwise tests fail after being run 99x
if(!defined $ENV{GMX_MAXBACKUP}) {
    $ENV{GMX_MAXBACKUP}=-1
}

my $tmpi_ranks = 0;
my $omp_threads = 0;
my $npme_ranks = -1;
my $mpi_ranks = 0;
my $double   = 0;
my $crosscompiling = 0;
my $bluegene = 0;
my $verbose  = 5;
my $xml      = 0;
# energy file comparision tolerance (potentials, not virials or pressure)
my $etol_rel = 0.001;
my $etol_abs = 0.05;
# topology comparision tolerance
my $ttol_rel = 0.0001;
my $ttol_abs = 0.001;
my $suffix   = '';
my $autosuffix   = '1';
my $prefix   = '';
# virial - this tests the shifted force contribution.
# However, it is a sum of very many large terms, so it is
# often numerically imprecise.
my $virtol_rel   = 0.01;
my $virtol_abs   = 0.1;

# force tolerance is measured by calculating scalar products.
# tolerance 0.001 means the scalar product between the two
# compared forces should be at least 0.999.
my $ftol_rel     = 0.001;
my $ftol_abs     = 0.05;
my $ftol_sprod   = 0.001;

# global variables to flag whether to explain some situations to the user
my $addversionnote = 0;
my $only_subdir = qr/.*/;
my $tolerance_factor = 1;

# trickery for program and reference file names, command line options, etc.
my $mdprefix  = sub {''};
my $mdparams  = '';
my $gpu_id    = '';
my $ref       = '';
my $mpirun    = 'mpirun';
my $parse_cmd = '';
my $gmx_cmd   = "gmx";
my $mdrun_cmd = "";
my %progs = ( 'grompp'   => '',
              'mdrun'    => '',
              'pdb2gmx'  => '',
              'check'    => '',
              'editconf' => '',
              'make_edi' => '',
              'nmeig'    => '' );

# Names for filenames used by individual test cases to limit the
# amount of parallelism that mdrun can attempt
my $max_openmp_threads_filename = 'max-openmp-threads';
my $max_mpi_ranks_filename = "max-mpi-ranks";

# List of all the generic subdirectories of tests; pdb2gmx and essentialdynamics are treated
# separately.
my @all_dirs = ('simple', 'complex', 'kernel', 'freeenergy', 'rotation', 'extra');

sub add_gpu_id
{
    my ($gpuid_string, $pp_ranks) = @_;

    # we may or may not be using GPUS and/or separate PME ranks
    # depending on the gmxtest.pl command line and the test involved,
    # so the number of PP ranks has to be used to choose a sensible
    # subset of the gpuid string.
    return ($gpuid_string && $pp_ranks) ? "-gpu_id " . substr($gpuid_string, 0, $pp_ranks) : "";
}

sub specify_number_of_pme_ranks {
    my ($num_ranks, $npme_ranks, $grompp_mdp, $pp_ranks_ref) = @_;
    # Only try -npme if using some kind of PME, with enough ranks and
    # the user asked for it. Note that mdrun -npme > 0 is supported
    # only with 3+ total ranks.
    if ((find_in_file("(coulomb[-_]?type|vdw[-_]?type)\\s*=\\s*(pme|PME)", $grompp_mdp) > 0) &&
        (($npme_ranks == 0) ||
         (($npme_ranks > 0) && ($num_ranks >= 3))))
    {
        # Specify number of PME-only ranks
        $$pp_ranks_ref = $num_ranks - $npme_ranks;
        return "-npme $npme_ranks";
    }
    else {
        # Use default behaviour for MPMD with PME
        $$pp_ranks_ref = $num_ranks;
        return "";
    }
}

sub setup_vars()
{
    # We assume that the name of executables match the pattern 
    # ${prefix}mdrun[_mpi][_d]${suffix} where ${prefix} and ${suffix} are 
    # as defined above (or over-ridden on the command line), "_d" indicates 
    # a double-precision version, and (only in the case of mdrun) "_mpi" 
    # indicates a parallel version compiled with MPI.
    if ( $mpi_ranks > 0 ) {
	if ($autosuffix) {
            $gmx_cmd .= "_mpi";
	}
	if ( $bluegene > 0 )
	{
	    # edit the next line if you need to customize the call to mpirun
	    $mdprefix = sub { "$mpirun -n $_[0] --cwd " . getcwd() . " :" };
	} elsif ( $mpirun =~ /(ap|s)run/ ) {
	    $mdprefix = sub { "$mpirun -n $_[0]" };
	} else {
	    # edit the next line if you need to customize the call to mpirun
	    $mdprefix = sub { "$mpirun -np $_[0] -wdir " . getcwd() };
	}
    }
    if ($autosuffix && ( $double > 0)) {
        $gmx_cmd .= "_d";
    }
    foreach my $prog ( keys %progs ) {
        $progs{$prog} = "$gmx_cmd$suffix $prog";
    }
    if ( "$mdrun_cmd" ne "" ) {
        $progs{"mdrun"} = $mdrun_cmd;
    }
    $ref = 'reference_' . ($double > 0 ? 'd' : 's');
    
    # now do -tight or -relaxed stuff
    foreach my $var ( $etol_rel, $etol_abs, $ttol_rel, $ttol_abs, $virtol_rel, $virtol_abs, $ftol_rel, $ftol_abs, $ftol_sprod ) {
	$var *= $tolerance_factor;
    }
}

# Wrapper function to call system(), and then perhaps a callback based on the
# value of the return code. When no callback exists, a generic error is
# displayed
sub do_system
{
    my $command = shift;
    my $normalreturn = shift;
    my $callback = shift;
    $normalreturn = 0 unless(defined $normalreturn);

    if ( $verbose > 2 ) {
	print "$command\n";
    }
    my $returnvalue = system($command) >> 8;
    if ($normalreturn != $returnvalue)
    {
	if (defined $callback)
	{
	    $returnvalue = &$callback($returnvalue);
	}
        if($normalreturn != $returnvalue)
	{
	    print "\nAbnormal return value for '$command' was $returnvalue\n";
	}
    }
    return $returnvalue;
}

# Built-in replacement for some of grep
# Returns number of matches for pattern (1st arg) in file (2nd arg)
# If the file is grompp.mdp the search is case-insensitive.
# Optionally writes matching lines to file named in 3rd arg, which
#  simulates "grep pattern file > redirectfile"

sub find_in_file {
    my $pattern = shift;
    my $filename_to_search = shift;
    my $filename_for_redirect = shift;
    my $return=0;
    my $case_sensitive=1;

    defined($filename_to_search) || die "find_in_file: Missing argument\n";
    my $do_redirect = defined $filename_for_redirect;
    if ($do_redirect) {
        open(REDIRECT, ">$filename_for_redirect") || die "Could not open redirect file '$filename_for_redirect'\n";
    }

    # If the file doesn't exist, then the pattern doesn't
    # match. Fortunately, the clients of find_in_file that care about
    # the number of matches found use the name of a file that always
    # exists.
    open(FILE,$filename_to_search) || return 0;
    if ($filename_to_search =~ "grompp.mdp\$") {
        $case_sensitive=0;
    }

    while(<FILE>) {
        if ($case_sensitive ? /$pattern/ : /$pattern/i) {
            $return++;
            if ($do_redirect) {
                print REDIRECT $_;
            }
        }
    }
    close(FILE) || die "Could not close file '$filename_to_search'\n";
    if ($do_redirect) {
        close(REDIRECT) || die "Could not close file '$filename_for_redirect'\n";
    }
    return $return;
}

sub check_force($$)
{
    my $traj = shift;
    my $reftrr = shift;
    my $cfor = "checkforce.out";
    my $cfor2 = "checkforce.err";
    my $nerr_force = 0;
    do_system("$progs{'check'} -f $reftrr -f2 $traj -tol $ftol_rel -abstol $ftol_abs >$cfor 2>$cfor2", 0,
	      sub { print "\ngmx check failed on the .edr file while checking the forces\n"; $nerr_force = 1; });
    
    open(FIN,"$cfor");
    while(my $line=<FIN>)
    {
	if ($line =~ /^b([^ ]*) .*TRUE/) {
	    print("Existence of $1 doesn't match!\n");
	    $nerr_force++;
	    next;
	}elsif ($line =~ /^End of file on/) {
	    print("Different number of frames!\n");
	    $nerr_force++;
	    next;
	}elsif ($line =~ /^\r?$/ || $line =~ /^[xv]\[/ || $line =~ /Both files read correctly/ ) {
	    next;
	}elsif (!($line =~ /^f\[/)) {
	    print("Unknown Error: $line!\n");
	    $nerr_force++;
	    next;
	}
	my @ll=split("[()]",$line);
	my @f1=split(" ",$ll[1]);
	my @f2=split(" ",$ll[3]);
	
	my $l1 = sqrt($f1[0]*$f1[0]+$f1[1]*$f1[1]+$f1[2]*$f1[2]);
	my $l2 = sqrt($f2[0]*$f2[0]+$f2[1]*$f2[1]+$f2[2]*$f2[2]);
        my $denominator = $l1 * $l2;
        if (0.0 == $denominator &&
            $l1 != $l2) {
            # Hopefully this means there was an error somewhere,
            # rather than an atomic force that was (correctly) binary
            # identical to zero...
            $nerr_force++;
        } else {
            my $sprod = ($f1[0]*$f2[0]+$f1[1]*$f2[1]+$f1[2]*$f2[2]) / $denominator;

            if( $sprod < (1.0-$ftol_sprod))
            {
                $nerr_force = $nerr_force + 1;
            }
        }
    }     
    close(FIN);
    if ($nerr_force == 0) {
      unlink($cfor,$cfor2);
    }
    return $nerr_force;
}

sub check_virial($)
{
    my $refedr = shift;
    my $cvir = "checkvir.out";
    my $cvir2 = "checkvir.err";
    my $nerr_vir = 0;

    do_system("$progs{'check'} -e $refedr -e2 ener -tol $virtol_rel -abstol $virtol_abs -lastener Vir-ZZ >$cvir 2>$cvir2", 0,
	sub { print "\ngmx check failed on the .edr file while checking the virial\n"; $nerr_vir = 1; });
    
    open(VIN,"$cvir");
    while(my $line=<VIN>)
    {
	next unless $line =~ /Vir-/;
	my @v1=split(" ",substr($line,26,14)); #TODO replace substr with split to make more reliable
	my @v2=split(" ",substr($line,52,13)); #if again reactiving check_virial
	
	my $diff = abs($v1[0]-$v2[0]);
	
	my $norm = abs($v1[0])+abs($v2[0]);
	
	if((2*$diff > $virtol_rel *$norm) && ($diff>$virtol_abs))
	{
	    $nerr_vir = $nerr_vir + 1;
	}
    }     
    close(VIN);
    if ($nerr_vir == 0) {
      unlink($cvir,$cvir2);
    }
    return $nerr_vir;
}

sub check_xvg {
    my $refx = shift;
    my $kkk  = shift;
    my $ndx = shift;
    my $pdb2gmx_test_names = shift;

    my $nerr = 0;
    if ((-f $refx) && (-f $kkk)) {
	open(REF,"$refx") || die "Could not open file '$refx'\n";
	open(KKK,"$kkk") || die "Could not open file '$kkk'\n";
	my $n = 0;
	my $header = 0;
	while (my $line = <REF>) {
	    my $line2=<KKK>;
            if (not defined($line2)){#REF has more lines
                $nerr++;
                next;
            }
	    next if $line =~ /^[@#]/;
	    chomp($line);
	    chomp($line2);
	    my @tmp = split(' ',$line);
	    my @tmp2 = split(' ',$line2);
	    my $x1 = $tmp[$ndx];
	    my $x2 = $tmp2[$ndx];
            my $error;
            my $hasPdb2gmx_test_name = defined $pdb2gmx_test_names && defined $$pdb2gmx_test_names[$n];
            print XML "<testcase name=\"$$pdb2gmx_test_names[$n]\">\n" if ($xml && $hasPdb2gmx_test_name);
            my $tol;
            if ($x1+$x2==0) { 
              $error = abs($x1-$x2);
              $tol = $etol_abs;
            } else {
              $error = abs(($x1-$x2)/($x1+$x2));
              $tol = $etol_rel;
            }
            if ($error > $tol) {
                $nerr++;
                if (!$header) {
            	$header = 1;
            	print("Here follows a list of the lines in $refx and $kkk which did not\npass the comparison test within tolerance $etol_rel\nIndex  Reference   This test       Error  Description\n");
                }
                printf("%4d  %10g  %10g  %10g  %s\n",$n+1,$x1,$x2, $error, 
            	   $hasPdb2gmx_test_name ? $$pdb2gmx_test_names[$n] : 'unknown');
                printf(XML "<error message=\"Reference: %g Result: %g Error: %g\"/>\n",$x1,$x2, $error) 
            	if ($xml && $hasPdb2gmx_test_name);
                
            }
            print XML "</testcase>\n" if ($xml && $hasPdb2gmx_test_name);
            $n++;
	}
        while (my $line2=<KKK>) {#KKK has more lines
          $nerr++
        }
	close REF;
	close KKK;
    }
    return $nerr;
}

# Callback used only when mdrun has returned an error, so we can work
# out how to try to call it so it works
sub how_should_we_rerun_mdrun {
    my ($tmpi_ranks_ref, $omp_threads_ref, $mpi_ranks_ref, $gpu_id_ref) = @_;

    # Is it because we are using too many cores, or trying to use -nt
    # with a reference build, or running a test that does not run in
    # parallel?
    open(MDOUT,"mdrun.out");
    my(@lines) = <MDOUT>;
    close(MDOUT);
    my $rerun = -1;
    foreach my $line (@lines) {
        if ($line =~ /There is no domain decomposition for/) {
            my $new_mpi_ranks = 8;
            if ($$mpi_ranks_ref > 0) {
                $$mpi_ranks_ref = $new_mpi_ranks;
            } else {
                $$tmpi_ranks_ref = ${new_mpi_ranks};
            }
            print ("Mdrun cannot use the requested (or automatic) number of ranks, retrying with $new_mpi_ranks.\n");
            $rerun = 1;
            last;
        }
        elsif ($line =~ /Setting the number of thread-MPI .* is only supported/) {
            printf ("Mdrun was compiled without thread-MPI or MPI support. Retrying with only 1 rank.\n");
            $$tmpi_ranks_ref = 0;
            $rerun = 1;
            last;
        }
        elsif ($line =~ /OpenMP threads were requested/) {
            $$omp_threads_ref = 8;
            print ("Mdrun cannot use the requested (or automatic) number of OpenMP threads, retrying with $$omp_threads_ref.\n");
            $rerun = 1;
            last;
        }
        elsif ($line =~ /Domain decomposition does not support simple neighbor searching/) {
            my $new_mpi_ranks = 1;
            print ("Mdrun cannot use the requested (or automatic) number of MPI ranks, retrying with ${new_mpi_ranks}.\n");
            if ($$mpi_ranks_ref > 0) {
                $$mpi_ranks_ref = $new_mpi_ranks;
            } else {
                $$tmpi_ranks_ref = ${new_mpi_ranks};
            }
            $$gpu_id_ref = substr($$gpu_id_ref, 0, $new_mpi_ranks);
            $rerun = 1;
            last;
        }
        elsif ($line =~ " was started with .* PP .*, but .* (.*) GPU") {
            # The above will match the strings generated by mdrun for
            # the default behaviour, and when -gpu_id was specified.
            my $new_mpi_ranks = $1;
            if ($$mpi_ranks_ref > 0) {
                $$mpi_ranks_ref = $new_mpi_ranks;
            } else {
                $$tmpi_ranks_ref = ${new_mpi_ranks};
            }
            $$gpu_id_ref = substr($$gpu_id_ref, 0, $new_mpi_ranks);
            $rerun = 1;
            last;
        }
        elsif ($line =~ "Your choice of .* MPI rank.* and the use of .* total threads .* leads to the use of .* OpenMP threads, whereas we expect the optimum to be with more MPI ranks" ||
               $line =~ "Your choice of number of MPI ranks and amount of resources results in using" ) {
            # On large nodes this error needs to be handled.  It can
            # be converted to a warning with setting the environment
            # variable GMX_BYPASS_EFFICIENCY_CHECK, but we don't want
            # anybody testing GROMACS to have to do that, whether with
            # "make check" or stand-alone.
            my $new_omp_threads = 6; # matches nthreads_omp_mpi_target_max in the code
            if ($$omp_threads_ref > 6)
            {
                $$omp_threads_ref = $new_omp_threads;
            }
            else
            {
                # Do something that will work!
                $$omp_threads_ref = 1;
            }
            $rerun = 1;
            last;
        }
    }
    return $rerun;
}

sub run_mdrun {
    # Copy all parameters by value, which is useful so we can modify
    # them if we need to, and have the changes local to this test
    my ($tmpi_ranks, $omp_threads, $mpi_ranks, $npme_ranks, $gpu_id, $mdprefix, $mdparams, $grompp_mdp) = @_;
    # Only one of tmpi_ranks or mpi_ranks may be greater than zero, but
    # this is checked for sanity after parsing user input.

    # Set up and enforce the maximum number of OpenMP threads to
    # try for this test case
    if ( -f $max_openmp_threads_filename )
    {
        open my $fh, '<', $max_openmp_threads_filename or die "error opening $max_openmp_threads_filename: $!";
        $omp_threads = do { local $/; <$fh> };
        chomp $omp_threads;
    }

    # Set up and enforce the maximum number of MPI ranks to try
    # for this test case
    if ( -f $max_mpi_ranks_filename ) {
        open my $fh, '<', $max_mpi_ranks_filename or die "error opening $max_mpi_ranks_filename: $!";
        my $max_ranks = do { local $/; <$fh> };
        chomp $max_ranks;
        if ($mpi_ranks > 0) {
            $mpi_ranks = ($mpi_ranks < $max_ranks) ? $mpi_ranks : $max_ranks;
        } elsif ($tmpi_ranks > 0) {
            $tmpi_ranks = ($tmpi_ranks < $max_ranks) ? $tmpi_ranks : $max_ranks;
        } else {
            # The user specified nothing, but the default must
            # still honour this maximum!
            #
            # If mdrun is serial, this code will trigger a second run
            # where -ntmpi will not be set. This is not ideal, but
            # only a few test cases specify the maximum number of
            # ranks, and pretty much only Jenkins should compile
            # the serial version of mdrun. In the absence of an
            # explicit serial mode for this script, it should lead
            # to a net improvement in Jenkins throughput.
            $tmpi_ranks = $max_ranks;
        }
    }

    # There is no general way to make all the tests pass by default on
    # all possible hardware, because there's not yet a good way for
    # the harness to find out about the hardware and mdrun
    # configuration in time to moderate the user (or default) choice
    # for various parallelism settings. Some tests can't run with more
    # than a few domains, yet machines might have many cores and/or
    # inconvenient ratios of cores to GPUs, or mdrun might be
    # deliberately compiled without OpenMP, etc. So we're forced to
    # inspect the output of mdrun when it fails, and try to react
    # in an appropriate way.
    #
    # First, we try running with whatever number of whatever the
    # user/default wants, then adapt to the error messages, for
    # at most three attempts.
    for my $attempt (0 .. 2)
    {
        my $ntmpi_opt = '';
        my $ntomp_opt = '';
        if (0 < $tmpi_ranks) {
            $ntmpi_opt = "-ntmpi $tmpi_ranks";
        }
        if (!find_in_file("cutoff[-_]scheme.*=.*group", $grompp_mdp) > 0 && $omp_threads > 0) {
            $ntomp_opt = "-ntomp $omp_threads";
        }

        my $pp_ranks = undef;
        my $num_ranks = $tmpi_ranks < $mpi_ranks ? $mpi_ranks : $tmpi_ranks;
        my $npme_opt = specify_number_of_pme_ranks($num_ranks, $npme_ranks, $grompp_mdp, \$pp_ranks);
        my $gpuid_opt = add_gpu_id($gpu_id, $pp_ranks);
        my $command = $mdprefix->($mpi_ranks)
            . " $progs{'mdrun'} $ntmpi_opt $npme_opt $gpuid_opt $ntomp_opt $mdparams >mdrun.out 2>&1";

        #do_system using the special callback will return:
        #1 for known error: rerun
        #0 exit value was 0 (and thus this callback wasn't called)
        #-1 unknown error
        my $nerror = do_system($command, 0,
                               sub { how_should_we_rerun_mdrun(\$tmpi_ranks, \$omp_threads, \$mpi_ranks, \$gpu_id) });
        if ($nerror > 0) {
            print "Retrying mdrun with better settings...\n";
        } else {
            return $nerror;
        }
    }
    # mdrun always failed, so pass the error upstream
    return 1;
}

# Parameters:
#  dir        : The directory in which the test will run, and output files will be created/found
#  input_dir  : The relative path to the directory in which the input and reference files can be found (normally ".")
#  test_name  : The name of the test case to report to the user
#
# Returns : 1 for success, 0 for failure
sub test_case {
    my ($dir, $input_dir, $test_name) = @_;
    my $success = 0;

    my $cwd = getcwd();
    chdir($dir);
    if ($verbose > 1) {
        print "Testing $test_name . . . ";
    }

    my $nerror = 0;
    my $ndx = "";
    if ( -f "$input_dir/index.ndx" ) {
        $ndx = "-n $input_dir/index";
    }
    my $grompp_mdp = "$input_dir/grompp.mdp";
    my $grompp_out = "grompp.out";
    my $grompp_err = "grompp.err";
    $nerror = do_system("$progs{'grompp'} -f $grompp_mdp -c $input_dir/conf -p $input_dir/topol -maxwarn 10 $ndx >$grompp_out 2>$grompp_err");

    my @error_detail;
    if (! -f "topol.tpr") {
        print ("No topol.tpr file in $dir. grompp failed\n");
        $nerror = 1;
    }
    if ($nerror == 0) {
        my $reftpr = "$input_dir/${ref}.tpr";
        if (! -f $reftpr) {
            print ("No $reftpr file for $test_name\n");
            print ("This means you are not really testing $test_name\n");
            copy('topol.tpr', $reftpr);
        } else {
            my $tprout="checktpr.out";
            my $tprerr="checktpr.err";
            do_system("$progs{'check'} -s1 $reftpr -s2 topol.tpr -tol $ttol_rel -abstol $ttol_abs >$tprout 2>$tprerr", 0,
                      sub { print "Comparison of input .tpr files failed!\n"; $nerror = 1; });
            $nerror |= find_in_file("^(?!comparing|WARNING)","$tprout");
            if ($nerror > 0) {
                push(@error_detail, ("checktpr.out", "checktpr.err"));
                print "topol.tpr file different from $reftpr. Check files in $dir for $test_name\n";
            }
            if (find_in_file ('reading tpx file (reference_[sd].tpr) version .* with version .* program',"$tprout") > 0) {
                print "\nThe GROMACS version being tested may be older than the reference version.\nPlease see the note at end of this output.\n";
                $addversionnote = 1;
            }
            if ($nerror == 0) {
                unlink($tprout,$tprerr);
            }
        }
    } else {
        push(@error_detail, ($grompp_err, $grompp_out));
    }

    if ($nerror == 0) {
        my $grompp_warn = "grompp.warn";
        open(GROMPP,$grompp_err) || die "Could not open file '$grompp_err'\n";
        open(WARN,"> $grompp_warn") || die "Could not open file '$grompp_warn'\n";
        my $p=0;
        while(<GROMPP>) {
            $p=1 if /^WARNING/;
            print WARN if ($p);
            $p=0 if /^\r?$/; # match a blank line, even with Windows line ending
        }
        close(GROMPP) || die "Could not close file '$grompp_err'\n";
        close(WARN) || die "Could not close file '$grompp_warn'\n";
        my $ref_warn = "$input_dir/${ref}.warn";
        if (! -f $ref_warn) {
            print("No $ref_warn file for $test_name\n");
            print ("This means you are not really testing $test_name\n");
            copy($grompp_warn, $ref_warn);
        } else {
            open(WARN1,$grompp_warn) || die "Could not open file '$grompp_warn'\n";
            open(WARN2,"$ref_warn") || die "Could not open file '$grompp_warn'\n";
            while (my $line1=<WARN1>) {
                my $line2=<WARN2>;
                if (not defined($line2)){#FILE1 has more lines
                    print("$grompp_warn has more lines than $ref_warn\n");
                    $nerror++;
                    next;
                }
                $line1 =~ s/(e[-+])0([0-9][0-9])/$1$2/g; #hack on windows X.Xe-00X -> X.Xe-0X (posix)

                # Hack to avoid issues based only on a difference in
                # relative path to grompp.mdp, which might happen in a
                # non-problematic way when the same test case is run
                # in a different output directory
                $line1 =~ s/file .*grompp.mdp/file grompp.mdp/;
                $line2 =~ s/file .*grompp.mdp/file grompp.mdp/;

                if ("$line2" ne "$line1") {
                    $nerror++;
                    # Note that the next line uses the fact that
                    # $line1 and $line2 still have their end-of-line
                    # termination characters
                    print("$grompp_warn had line\n$line1 which did not match line\n$line2 from $ref_warn\n");
                }
            }
            while (my $line2=<WARN2>) {#FILE2 has more lines
                print("$grompp_warn has fewer lines than $ref_warn\n");
                $nerror++
            }
            if ($nerror>0) {
                print("Different warnings in $ref_warn and $grompp_warn\n");
                push(@error_detail, ($grompp_err, $grompp_out));
            } else {
                unlink($grompp_warn);
            }
        }
    }
    if ($nerror == 0) {
        # Do the mdrun at last!

        # With tunepme Coul-Sr/Recip isn't reproducible
        my $local_mdparams = $mdparams . " -notunepme";
        if ($dir =~ /nb_kernel.*CSTab/) {
            # All group-kernel test cases with cubic splines need
            # tables for their interactions.
            $local_mdparams .= " -table $input_dir/../table -tablep $input_dir/../tablep";
        }
        my $part = "";
        if ( -f "continue.cpt" ) {
            $local_mdparams .= " -cpi continue -noappend";
            $part = ".part0002";
        }
        if ($test_name =~ /-cpu-only/) {
            $local_mdparams .= " -nb cpu";
        }
        $nerror = run_mdrun($tmpi_ranks, $omp_threads, $mpi_ranks, $npme_ranks, $gpu_id, $mdprefix, $local_mdparams, $grompp_mdp);
        if ($nerror != 0) {
            if ($parse_cmd eq '') {
                push(@error_detail, ("mdrun.out", "md.log"));
            } else {
                do_system("$parse_cmd <mdrun.out >mdrun_parsed.out");
                push(@error_detail, ("mdrun_parsed.out", "md.log"));
            }
        }

        my $ener = "ener${part}.edr";
        my $traj = "traj${part}.trr";
        my $log = "md${part}.log";
        if ($nerror!=0) {
            $nerror=1;
        } elsif ((-f "$ener" ) && (-f "$traj")) {  # Check whether we have any output
            # Now check whether we have any reference files
            my $refedr = "$input_dir/${ref}.edr";
            if (! -f  $refedr) {
                print ("No $refedr file for $test_name.\n");
                print ("This means you are not really testing $test_name\n");
                copy("$ener", $refedr);
            } else {
                my $potout="checkpot.out";
                my $poterr="checkpot.err";
                # Now do the real tests
                do_system("$progs{'check'} -e $refedr -e2 $ener -tol $etol_rel -abstol $etol_abs -lastener Potential >$potout 2>$poterr", 0,
                          sub {
                              if($nerror != 0) {
                                  print "\ngmx check failed on the .edr file, probably because mdrun also failed";
                              } else {
                                  print "\ngmx check FAILED on the .edr file";
                                  $nerror = 1;
                              }
                          });
                my $nerr_pot = find_in_file("step","$potout");
                push(@error_detail, "$potout ($nerr_pot errors)") if ($nerr_pot > 0);

                my $nerr_vir   = 0; #TODO: check_virial();
                push(@error_detail, "checkvir.out ($nerr_vir errors)") if ($nerr_vir > 0);

                $nerror |= $nerr_pot | $nerr_vir;
                unlink($potout,$poterr) unless $nerr_pot;
            }
            my $reftrr = "$input_dir/${ref}.trr";
            if (! -f $reftrr ) {
                print ("No $reftrr file for $test_name.\n");
                print ("This means you are not really testing $test_name\n");
                copy("$traj", $reftrr);
            } else {
                # Now do the real tests
                my $nerr_force = check_force($traj, $reftrr);
                push(@error_detail, "checkforce.out ($nerr_force errors)") if ($nerr_force > 0);
                $nerror |= $nerr_force;
            }
            my $reflog = "${ref}.log";
            if (! -f $reflog ) {
                copy($log, $reflog);
            }

            # This bit below is only relevant for free energy tests
            my $refxvg = "$input_dir/${ref}.xvg";
            my $nerr_xvg = check_xvg($refxvg,'dgdl.xvg',1);
            push(@error_detail, "$refxvg ($nerr_xvg errors)") if ($nerr_xvg > 0);
            $nerror |= $nerr_xvg;
        }
        else {
            print "mdrun output files ${ener} and/or ${traj} were not found.\n";
            $nerror = 1;
        }
    }
    my $error_detail = join(', ', @error_detail) . ' ';
    print XML "<testcase name=\"$test_name\">\n" if ($xml);
    if ($nerror > 0) {
        my $error_detail = join(', ', @error_detail) . ' ';
        print "FAILED. Check ${error_detail}file(s) in $dir for $test_name\n";
        if ($xml) {
            print XML "<error message=\"Errors in ${error_detail}\">\n";
            print XML "<![CDATA[\n";
            foreach my $err (@error_detail) {
                my @err = split(/ /, $err);
                my $errfn = $err[0];
                print XML "$errfn:\n";
                if (!open FH, $errfn) {
                    print XML "failed to open $errfn";
                } else {
                    while(my $line=<FH>) {
                        $line=~s/\x00//g; #remove invalid XML characters
                        print XML $line;
                    }
                }
                print XML "\n--------------------------------\n";
                close FH;
            }
            print XML "]]>\n";
            print XML "</error>";
        }
    }
    else {
        if ($verbose > 0) {
            if (find_in_file(".", $grompp_mdp) < 50) {
                # if the input .mdp file is trivially short, then
                # the diff test below will always fail, however this
                # is normal and expected for the usefully-short
                # kernel test .mdp files, so we don't compare the
                # .mdp files in this case
                print "PASSED\n";
            }
            else {
                my $mdp_result = 0;
                foreach my $reference_mdp ( $grompp_mdp ) {
                    if (-f $reference_mdp) {
                        open(FILE1,$reference_mdp) || die "Could not open file '$reference_mdp'\n";
                        open(FILE2,"mdout.mdp") || die "Could not open file 'mdout.mdp'\n";
                        my $diff=0;
                        while (my $line1=<FILE1>) {
                            my $line2=<FILE2>;
                            next if $line1 =~ /(data|host|user|generated)/;
                            next if $line2 =~ /(data|host|user|generated)/;
                            if (not defined($line2)){#FILE1 has more lines
                                $diff++;
                                next;
                            }
                            $diff++ unless ("$line2" eq "$line1");
                        }
                        while (my $line2=<FILE2>) {#FILE2 has more lines
                            $diff++
                        }
                        $mdp_result++ if $diff > 2;
                        close(FILE1) || die "Could not close file '$reference_mdp'\n";
                        close(FILE2) || die "Could not close file 'mdout.mdp'\n";
                    }
                }
                if ($mdp_result > 0) {
                    print("PASSED but check mdp file differences\n");
                }
                else {
                    print "PASSED\n";
                    unlink("mdout.mdp");
                }
            }
        }
        $success = 1;
    }
    print XML "</testcase>\n" if ($xml);
    chdir($cwd);
    return $success;
}

sub test_systems {
    my ($npassed, $nn, @subdirs) = @_;
    $$npassed = 0;
    $$nn = 0;
    foreach my $dir ( @subdirs ) {
        my $test_name = $dir;
        my $input_dir = ".";
        $$nn++;
        $$npassed += test_case $dir, $input_dir, $test_name;

        # Only nbnxn test cases can execute GPU-based non-bonded
        # kernels. If GPU kernels were used to run this test case
        # (whether natively or in emulation), run it again in CPU-only
        # mode. This relies on test_case() not clearing up md.log.
        if ($test_name =~ /nbnxn/ && 0 < find_in_file("^Using .* 8x8 non-bonded kernels", "$dir/md.log")) {
            if ($mdparams =~ /-nb/) {
                print("Test case $test_name has been run using GPU non-bonded kernels,\n" .
                      "and would normally be run again using only CPU-based non-bonded\n" .
                      "kernels, but because you have set -nb in the -mdparam string,\n" .
                      "the test harness will stay out of your way.\n");
            } else {
                print("Re-running $test_name using only CPU-based non-bonded kernels\n");

                # Set up new variables that control where the test
                # case does its I/O
                my $new_dir = "$dir/cpu-only";
                my $new_input_dir = "..";
                my $new_test_name = "${test_name}-cpu-only";
                mkdir $new_dir;

                # Copy the files that limit parallelism
                foreach my $limiter_file ($max_openmp_threads_filename, $max_mpi_ranks_filename) {
                    if (-f "$dir/$limiter_file") {
                        copy("$dir/$limiter_file", $new_dir);
                    }
                }

                $$nn++;
                $$npassed += test_case $new_dir, $new_input_dir, $new_test_name;
            }
        }
    }
}

sub cleandirs {
    my $mydir = shift;
    chdir($mydir);
    foreach my $dir ( <*> ) {
	if ( -d $dir ) {
	    chdir($dir);
	    print "Cleaning $dir\n"; 
	    my @args = glob("#*# *~ *.out core.* field.xvg dgdl.xvg topol.tpr confout*.gro ener*.edr md.log traj*.trr *.tmp mdout.mdp step*.pdb *~ grompp[A-z]* state*.cpt *.xtc *.err confout*.gro cpu-only/*" );
	    unlink (@args);
	    chdir("..");
	}
    }
    chdir("..");
}

sub refcleandir {
    my $sdir = shift;
    
    if (-d $sdir ) {
	cleandirs($sdir);
	chdir($sdir);
	foreach my $dir ( <*> ) {
	    if ( -d $dir ) {
		chdir($dir);
		print "Removing reference files in $dir\n"; 
                my @extensions = ('log', 'warn', 'edr', 'tpr', 'trr');
                map { unlink("${ref}.$_") } @extensions;
		chdir("..");
	    }
	}
	chdir("..");
    }
}

sub test_dirs {
    my $dirs = shift;
    chdir($dirs);
    # glob all files, but retain only directories that match the regular expression
    my @subdirs = map { (-d $_ && $_ =~ $only_subdir) ? $_ : () } <*>;
    print XML "<testsuite name=\"$dirs\">\n" if ($xml);
    my ($nn, $npassed);
    test_systems(\$npassed, \$nn, @subdirs);
    print XML "</testsuite>\n" if ($xml);
    my $nerror=0;
    if ($npassed < $nn) {
	printf("%d out of $nn $dirs tests FAILED\n",$nn-$npassed);
	$nerror=1;
    }
    else {
	print("All $nn $dirs tests PASSED\n");
    }
    chdir("..");
    return $nerror;
}

#format for cfg files is:
#- command line 
#- list of output files which should be compared (one per line)
#- emtpy line
sub test_tools {
    chdir("tools");
    my $ncfg = 0;
    my $nerror_cfg = 0;

    foreach my $cfg ( glob("*.cfg") ) { #loop over config files
	$ncfg++;
	open(FIN,$cfg);
	my @cfg_name = split(".cfg",$cfg);
	my $cfg_name = $cfg_name[0];
	mkdir($cfg_name);
	chdir($cfg_name);
	my $ncmd = 0;
	my $nerror_cmd = 0;
	if ($verbose > 1) {
	    print "Testing $cfg_name . . . \n";
	}
	while(my $line=<FIN>) {  #loop over commands (seperated by empty line)
	    $ncmd++;
	    chomp($line);
	    my $cmdline = $line;
	    my @ofiles;
	    my $error = 0;
	    while(my $line=<FIN>) {  #add output fiels
		chomp($line);
		if ($line eq "") { last; }
		push(@ofiles, $line);
	    }
	    mkdir($ncmd);
	    chdir($ncmd);

	    do_system("$cmdline >$cfg_name.out 2>&1", 0,
		      sub { print "\n'".$cmdline."' failed"; $error = 1; });
	    
	    if ($error==0) {
		foreach my $of (@ofiles) {
		    if (! -f "$of.ref") {
			print "No file $of.ref. You are not really testing $cfg_name\n";
			copy($of, "$of.ref");
		    }
		    else {
			my $nerror = check_xvg("$of.ref",$of,1);  #TODO: check all columns
			if ( $nerror != 0 ) {
			    print "There were $nerror differences in $of output file\n";
			    $error += 1;
			}
		    }
		}
	    }
	    if ($error > 0) {
		print "FAILED $cfg_name test $ncmd\n";
		$nerror_cfg++;
		$nerror_cmd++;
	    }
	    chdir("..");
	}
	if ($nerror_cmd>0) {
	    print "$nerror_cmd out of $ncmd $cfg_name tests FAILED\n";
	} elsif ($verbose>1) {
	    print "All $ncmd $cfg_name tests PASSED\n";
	}
	chdir("..");
    }
    if ($nerror_cfg>0) {
	print "$nerror_cfg out of $ncfg tools tests FAILED\n";
    } else {
	print "All $ncfg tools tests PASSED\n";
    }
    return $nerror_cfg;
}

sub compare_xvg {
  my $xvg_ref = shift;
  my $xvg_test = shift;
  my $rel_toler = shift;
  
  open(REF, "$xvg_ref") || die("FAILED: Can not read $xvg_ref");
  my @ref = <REF>;
  close REF;
  open(TEST, "$xvg_test") || die("FAILED: Can not read $xvg_test");
  my @test = <TEST>;
  close TEST;
  die("FAILED: test $xvg_test has different size than reference $xvg_ref") if ($#ref != $#test);
  for(my $i = 0; ($i<=$#ref); $i++) {
    my @r = split(' ', $ref[$i]);
    my @t = split(' ', $test[$i]);
    die("FAILED: test $i = $t[1]  ref $i = $r[1]") if (abs($r[1]-$t[1])/$r[1] > $rel_toler);
  }
}

sub test_normal_modes {
  my $logfn = "normal_modes.log";
  
# Only run this test in double precision
  return if ($double == 0);
  
  chdir("normal_modes");
  open(LOG,">$logfn")  || die("FAILED: Opening $logfn for writing");
  foreach my $dir ( glob("*") ) {
    if ( -d $dir) {
      chdir $dir;
      open(PIPE,"$progs{'grompp'} -c conf.g96 2>&1 |");
      print LOG while <PIPE>;
      close PIPE;
      my $mtx = "nm.mtx";
      open(PIPE,"$progs{'mdrun'} -mtx $mtx 2>&1 |");
      print LOG while <PIPE>;
      close PIPE;
      open(PIPE,"$progs{'nmeig'} -f $mtx -first 7 -xvg no 2>&1 |");
      print LOG while <PIPE>;
      close PIPE;
      my $rel_tol = 1e-3;
      compare_xvg("eigenval_reference.xvg", "eigenval.xvg", $rel_tol);
      compare_xvg("eigenfreq_reference.xvg", "eigenfreq.xvg", $rel_tol);
      chdir "..";
    }
  }
  close LOG;
  chdir "..";
}

# Compares the data entries of two .xvg files at column 'index'
sub compare_xvg_column {
    my ( $refFile, $cmpFile, $index, $absTol ) = @_;

    if ($verbose > 1) {
        print("Comparing column #$index of .xvg files $refFile and $cmpFile");
    }
    open(REF,"$refFile") || die "Could not open file '$refFile'\n";
    open(CMP,"$cmpFile") || die "Could not open file '$cmpFile'\n";

    my $nLineRef = 0;
    my $nLineCmp = 0;

    my $nerr    = 0;
    my $bFirst  = 1;
    my $bFirst2 = 1;

    while (my $refLine = <REF>) {
        $nLineRef++;
        next if $refLine =~ /^[@#]/;  # Skip non-data entries in REF

        my $cmpLine;
        while ($cmpLine = <CMP> ) {
            $nLineCmp++;
            next if $cmpLine =~ /^[@#]/;  # Skip non-data entries in CMP
            last;
        }

        if (defined($refLine) && defined($cmpLine)) {
            chomp($refLine);
            chomp($cmpLine);

            my @refArray = split(' ', $refLine);
            my @cmpArray = split(' ',$cmpLine);
            my $xR       = $refArray[$index];
            my $xC       = $cmpArray[$index];
            my $absError = abs($xR - $xC);
 
            if ($absError > $absTol) {
                $nerr++;
                if ($bFirst) {
                    print("\nHere follows a list of the lines in $refFile and $cmpFile which did not\n");
                    print("pass the comparison test within a tolerance of $absTol\n");
                    print("Lines ref cmp  Reference   This test  Error\n");
                    $bFirst = 0;
                }
                printf("%5d  %5d  %10g  %10g  %10e\n", $nLineRef, $nLineCmp, $xR, $xC, $absError);
            }
        } else {
            if ($bFirst2) {
                printf("\nNumber of lines with valid data entries in $cmpFile does not match $refFile!\n");
                $bFirst2 = 0;
            }
            $nerr++;
        }
    }

    if (!$nerr && $verbose > 1) {
        print(" ... all within tolerance of $absTol\n");
    }

    close REF;
    close CMP;

    return $nerr;
}


# Special routine for the essential dynamics tests. It calculates the radius 
# from two eigenvectors and compares it to the expected radius.
sub check_radius {
    my ( $File, $ind1, $ind2, $indR, $absTol, $bRadcon ) = @_;

    if ($verbose > 1) {
        print("Checking radius for file $File ... ");
    }

    open(FILE,"$File") || die "Could not open file '$File'\n";

    my $nLine = 0;

    my $nerr   = 0;
    my $bFirst = 1;
    my $bFirst1 = 1;
    my $ev1start = 0;
    my $ev2start = 0;
    my $radstart = 0;
    my $rad;
    my $absError;

    while (my $Line = <FILE>) {
        $nLine++;
        next if $Line =~ /^[@#]/;  # Skip non-data entries

        if (defined($Line)) {
            chomp($Line);

            my @Array   = split(' ', $Line);
            my $ev1     = $Array[$ind1];
            my $ev2     = $Array[$ind2];
            my $rad_ref = $Array[$indR];
            if ($bFirst1) {
                $ev1start = $ev1;
                $ev2start = $ev2;
                $radstart = $rad_ref;
                $bFirst1  = 0;
            }

            if ($bRadcon) {
                # The explicit values given are the projections of the target
                # structure onto the eigenvectors one and two for the RADCON test case.
                $rad      = sqrt((-1.30289e-02 - $ev1)**2 + (1.81182e-02 - $ev2)**2);
            } else {
                $rad      = sqrt(($ev1 - $ev1start)**2 + ($ev2 - $ev2start)**2) + $radstart;
            }
           $absError = abs($rad - $rad_ref);

            if ($absError > $absTol) {
                $nerr++;
                if ($bFirst) {
                    print("\nHere follows a list of the lines in $File which did not\n");
                    print("pass the comparison test within a tolerance of $absTol\n");
                    print("Line     Expected      Actual  Error in radius (nm)\n");
                    $bFirst = 0;
                }
                printf("%5d  %10g  %10g  %6.2e\n", $nLine, $rad_ref, $rad, $absError);
            }
        } else {
            $nerr++;
        }
    }

    if (!$nerr && $verbose > 1) {
        print("all values within tolerance of $absTol\n");
    }

    close FILE;

    return $nerr;
}


# Checks whether the data entries in a .xvg column consistently increase or decrease with time
# Use $absTol > 0 to allow for slight deviations from monotonicity
sub check_monotonicity {
    my ( $File, $iCol, $dir, $absTol ) = @_;
    $absTol = abs($absTol);  # Ignore sign

    if ($verbose > 1) {
        print("Checking monotonicity for column $iCol of file $File ... ");
    }

    open(FILE,"$File") || die "Could not open file '$File'\n";

    my $nLine = 0;

    my $nerr   = 0;
    my $bFirst = 1;
    my $extremumSoFar = 0;

    while (my $Line = <FILE>) {
        $nLine++;
        next if $Line =~ /^[@#]/;  # Skip non-data entries

        if (defined($Line)) {
            chomp($Line);

            my @Array   = split(' ', $Line);
            my $value   = $Array[$iCol];
            if ($bFirst) {
                $extremumSoFar = $value;
                $bFirst = 0;
            } else {
                if ($dir >= 0) { # expect monotonically increasing values
                    $extremumSoFar = max($extremumSoFar, $value);
                    if ($extremumSoFar - $absTol > $value) {
                        $nerr++;
                        printf("ERROR: value in line %d is not increasing monotonically! ", $nLine);
                        printf("value=%10g  (but already reached %10g before)\n", $value, $extremumSoFar);
                    }
                } else { # expect monotonically decreasing values
                    $extremumSoFar = min($extremumSoFar, $value);
                    if ($extremumSoFar + $absTol < $value) {
                        $nerr++;
                        printf("ERROR: value in line %d is not decreasing monotonically! ", $nLine);
                        printf("value=%10g  (but already reached %10g)\n", $value, $extremumSoFar);
                    }
                }
            }
        }
    }

    close FILE;

    if (!$nerr && $verbose > 1) {
        my $buf = "increasing";
        if ($dir < 0) {
            $buf = "decreasing";
        }
        print("values of column $iCol are $buf monotonically.\n");
    }

    return $nerr;
}


sub compare_textfiles {
  my $txtref = shift;
  my $txtcmp = shift;

  my $nerror = 0;

  open(REF, "$txtref") || die("FAILED: Can not read $txtref");
  my @ref = <REF>;
  close REF;
  open(CMP, "$txtcmp") || die("FAILED: Can not read $txtcmp");
  my @cmp = <CMP>;
  close CMP;
  die("FAILED: test @cmp has different size than reference @ref") if ($#ref != $#cmp);
  for(my $i = 0; ($i<=$#ref); $i++) {
    my $r = $ref[$i];
    my $c = $cmp[$i];
    
    if ($r ne $c) {
      $nerror++;
      print("Diff line $i of ref ($txtref) and cmp ($txtcmp):\n");
      print("ref: $r"); 
      print("cmp: $c");
    }
  }
  return $nerror;
}

# Runs a short MD simulation with essential dynamics (ED) constraints.
# If $ediArgs is empty, an ED input file must already be present, otherwise
# gmx make_edi is used to produce one with the provided arguments.
sub run_single_ed_system {
  my $dir = shift;
  my $ediArgs = shift;
  my $inArgs = shift;


  my $nerror = 0;
  my $edifn  = "sam.edi";
  my $edofn  = "$dir.xvg";
  my $grompp_out = "grompp.out";
  my $grompp_err = "grompp.err";
  my $makeEdi_out = "make_edi.out";
  my $makeEdi_err = "make_edi.err";

  if ($verbose > 1) {
      print "Testing $dir . . .\n";
  }

  # Make the .tpr file
  $nerror += do_system("$progs{'grompp'} -maxwarn 1 >$grompp_out 2>$grompp_err");

  if ($ediArgs) {
      # Make the essential dynamics .edi input file
      unlink $edifn;  # delete old .edi file (if any)
      $nerror += do_system("echo $inArgs | $progs{'make_edi'} -f eigenvec.trr $ediArgs -o $edifn 1>$makeEdi_out 2>$makeEdi_err");
      # Check whether we get the expected .edi file
      my $nerr = compare_textfiles($edifn, "sam_reference.edi");
      if ( $nerr > 0 ) {
          $nerror += $nerr;
          printf("Essential dynamics '$dir' FAILED: Produced .edi file does not match the reference!\n");
      }
  } 

  # Make a short simulation with essential dynamics constraints:
  if ( $nerror == 0 ) { # If we already had errors, there is no use in going on ...
    unlink $edofn;  # delete old essential dynamics .xvg output file (if any)
    $nerror = run_mdrun($tmpi_ranks, $omp_threads, $mpi_ranks, $npme_ranks, $gpu_id, $mdprefix, "-ei $edifn -eo $edofn", "grompp.mdp");
  }

  return $nerror;
}


sub test_essentialdynamics {
  my $nerror = 0;
  my $retval = 0;
  my $ntest  = 0;
  my $edref  = "edsam_reference.xvg";
  my $dir = "";

  chdir("essentialdynamics");

  # ----------------------------------------------------------------------------
  # Test fixed-step linear expansion along the first eigenvector. In this test,
  # the projection on EV 1 should increase by 0.0013 nm at every step.
  $dir = "linfix";
  $ntest++;
  chdir $dir || die;
  $retval = run_single_ed_system($dir, "-outfrq 1 -linfix 1 -linstep 0.0013", "0");
  if ( 0 == $retval ) {
      # Compare the essential dynamics output to the reference. This
      # makes sense only if the system was successfully run, otherwise we will
      # get a huge bunch of errors.
      $nerror += compare_xvg_column($edref, "$dir.xvg", 0, 0.0   );  # 0 = Time (ps), must match exacly!
      $nerror += compare_xvg_column($edref, "$dir.xvg", 1, 0.05  );  # 1 = RMSD (nm), just making sure it does not diverge completely
      $nerror += compare_xvg_column($edref, "$dir.xvg", 2, 1.0e-6);  # 2 = EV1 projection (nm) LINFIX, last value is the tolerance
  } else {
  	  $nerror++;
  }
  chdir "..";

  # ----------------------------------------------------------------------------
  # Test acceptance linear expansion along the first eigenvector. In this test,
  # steps towards a larger value of the projection of EV 1 should be accepted,
  # others rejected.
  $dir = "linacc";
  $ntest++;
  chdir $dir || die;
  $retval = run_single_ed_system($dir, "-outfrq 2 -linacc 1 -accdir +1", "0");
  if ( 0 == $retval ) {
      # Compare the essential dynamics output to the reference:
      $nerror += compare_xvg_column($edref, "$dir.xvg", 0, 0.0   );  # 0 = Time (ps), must match exacly!
      $nerror += compare_xvg_column($edref, "$dir.xvg", 1, 0.005 );  # 1 = RMSD (nm)
      $nerror += compare_xvg_column($edref, "$dir.xvg", 2, 0.001 );  # 2 = EV projection (nm) LINACC, last value is the tolerance
      $nerror += check_monotonicity("$dir.xvg",     2, +1, 5e-7  );  # EV 1 projection must increase monotonically
  } else {
      $nerror++;
  }
  chdir "..";

  # ----------------------------------------------------------------------------
  # Test fixed-step radius expansion along eigenvectors 1-2. The reported radius
  # should increase by 0.002 nm per step in this test.
  $dir = "radfix";
  $ntest++;
  chdir $dir || die;
  $retval = run_single_ed_system($dir, "-outfrq 1 -radfix 1-2 -radstep 0.002", "0");
  if ( 0 == $retval ) {
      # Compare the essential dynamics output to the reference:
      $nerror += compare_xvg_column($edref, "$dir.xvg", 0, 0.0   );  # 0 = Time (ps), must match exacly!
      $nerror += compare_xvg_column($edref, "$dir.xvg", 1, 0.05  );  # 1 = RMSD (nm)
      $nerror += compare_xvg_column($edref, "$dir.xvg", 4, 0.0   );  # 4 = RADFIX radius (nm), must match exactly
      $nerror += check_radius("$dir.xvg", 2, 3, 4, 1e-5, 0);         # RADFIX radius is calculated from columns 2-3, must match column 4
  } else {
      $nerror++;
  }
  chdir "..";

  # ----------------------------------------------------------------------------
  # Test acceptance radius expansion along eigenvectors 1-2. Steps in the desired 
  # direction should be accepted, others should be rejected. In effect, the reported
  # radius should increase monotonically.
  $dir = "radacc";
  $ntest++;
  chdir $dir || die;
  $retval = run_single_ed_system($dir, "-outfrq 3 -radacc 1-2", "0");
  if ( 0 == $retval ) {
      # Compare the essential dynamics output to the reference:
      $nerror += compare_xvg_column($edref, "$dir.xvg", 0, 0.0    );  # 0 = Time (ps), must match exacly!
      $nerror += compare_xvg_column($edref, "$dir.xvg", 1, 0.05   );  # 1 = RMSD (nm)
      $nerror += check_radius("$dir.xvg", 2, 3, 4, 1e-5, 0);          # RADACC radius is calculated from columns 2-3, must match column 4
      $nerror += check_monotonicity("$dir.xvg", 4, +1, 0.0);          # RADACC radius must increase monotonically
  } else {
      $nerror++;
  }
  chdir "..";

  # ----------------------------------------------------------------------------
  # Test acceptance radius contraction along eigenvectors 1-2 towards target structure.
  $dir = "radcon";
  $ntest++;
  chdir $dir || die;
  $retval = run_single_ed_system($dir, "-outfrq 1 -radcon 1-2 -tar target.pdb", "0 0");
  if ( 0 == $retval ) {
      # Compare the essential dynamics output to the reference:
      $nerror += compare_xvg_column($edref, "$dir.xvg", 0, 0.0    );  # 0 = Time (ps), must match exacly!
      $nerror += compare_xvg_column($edref, "$dir.xvg", 1, 0.05   );  # 1 = RMSD (nm)
      $nerror += check_radius("$dir.xvg", 2, 3, 4, 1e-5, 1);          # RADCON radius is calculated from columns 2-3, must match column 4
      $nerror += check_monotonicity("$dir.xvg", 4, -1, 0.0);          # RADCON radius must decrease monotonically towards zero.
  } else {
      $nerror++;
  }
  chdir "..";

  # ----------------------------------------------------------------------------
  # Test flooding in a complex setting, using harmonic restraints on two eigenvectors.
  # This first of two flooding tests checks the flooding output values of a single step
  # using a small tolerance.
  $dir = "flooding1";
  $ntest++;
  chdir $dir || die;
  $retval = run_single_ed_system($dir, "", "");
  if ( 0 == $retval ) {
      # Compare the essential dynamics output to the reference:
      $nerror += compare_xvg_column($edref, "$dir.xvg", 0, 0.0    );  # 0 = Time (ps), must match exacly!
      $nerror += compare_xvg_column($edref, "$dir.xvg", 1, 1e-5   );  # 1 = RMSD (nm)
      $nerror += compare_xvg_column($edref, "$dir.xvg", 2, 1e-5   );  # 2 = EV 1 projection (nm)
      $nerror += compare_xvg_column($edref, "$dir.xvg", 3, 1      );  # 3 = EV 1 V_fl (kJ/mol) 
      $nerror += compare_xvg_column($edref, "$dir.xvg", 4, 1      );  # 4 = EV 1 fl_forces (kJ/mol/nm) 
      $nerror += compare_xvg_column($edref, "$dir.xvg", 5, 1e-5   );  # 5 = EV 2 projection (nm)
      $nerror += compare_xvg_column($edref, "$dir.xvg", 6, 1      );  # 6 = EV 2 V_fl (kJ/mol) 
      $nerror += compare_xvg_column($edref, "$dir.xvg", 7, 1      );  # 7 = EV 2 fl_forces (kJ/mol/nm) 
  } else {
      $nerror++;
  }
  chdir "..";

  # ----------------------------------------------------------------------------
  # This second test checks the evolution of the flooding characteristic output values
  # over 50 time steps, therefore using larger tolerances.
  $dir = "flooding2";
  $ntest++;
  chdir $dir || die;
  $retval = run_single_ed_system($dir, "", "");
  if ( 0 == $retval ) {
      # Compare the essential dynamics output to the reference:
      $nerror += compare_xvg_column($edref, "$dir.xvg", 0, 0.0    );  # 0 = Time (ps), must match exacly!
      $nerror += compare_xvg_column($edref, "$dir.xvg", 1, 0.05   );  # 1 = RMSD (nm)
      $nerror += compare_xvg_column($edref, "$dir.xvg", 2, 0.025  );  # 2 = EV 1 projection (nm)
      $nerror += compare_xvg_column($edref, "$dir.xvg", 3, 250    );  # 3 = EV 1 V_fl (kJ/mol) 
      $nerror += compare_xvg_column($edref, "$dir.xvg", 4, 1000   );  # 4 = EV 1 fl_forces (kJ/mol/nm) 
      $nerror += compare_xvg_column($edref, "$dir.xvg", 5, 0.025  );  # 5 = EV 2 projection (nm)
      $nerror += compare_xvg_column($edref, "$dir.xvg", 6, 250    );  # 6 = EV 2 V_fl (kJ/mol) 
      $nerror += compare_xvg_column($edref, "$dir.xvg", 7, 1000   );  # 7 = EV 2 fl_forces (kJ/mol/nm) 
  } else {
      $nerror++;
  }
  chdir "..";

  # ----------------------------------------------------------------------------

  if (0 == $nerror) {
    print "All $ntest essential dynamics tests PASSED\n";
  }
  else {
    print "Essential dynamics tests FAILED with $nerror errors!\n";
  }

  chdir "..";
}

sub test_pdb2gmx {
    my $logfn = "pdb2gmx.log";

    chdir("pdb2gmx");
    open (LOG,">$logfn") || die("FAILED: Opening $logfn for writing");
    my $npdb_dir = 0;
    my @pdb_dirs = ();
    my $ntest    = 0;
    my $nerror   = 0;
    my @pdb2gmx_test_names;
    foreach my $pdb ( glob("*.pdb") ) {
	my $pdir = "pdb-$pdb";
	my @kkk  = split('\.',$pdir);
	my $dir  = $kkk[0];
	$pdb_dirs[$npdb_dir++] = $dir;
	mkdir($dir);
	chdir($dir);
	foreach my $ff ( "gromos43a1", "oplsaa", "gromos53a6" ) {
	    mkdir("ff$ff");
	    chdir("ff$ff");
	    my @water = ();
	    my @vsite = ( "none", "h" );
	    if ( $ff eq "oplsaa"  ) {
		@water = ( "tip3p", "tip4p", "tip5p" );
	    }
	    elsif ( $ff eq "encads" ) {
		@vsite = ( "none" );
		@water = ( "spc" );
	    }
	    else {
		@water = ( "spc", "spce" );
	    }
	    foreach my $dd ( @vsite ) {
		mkdir("$dd");
		chdir("$dd");
		foreach my $ww ( @water ) {
		    $ntest++;
		    push @pdb2gmx_test_names, "$pdb with $ff using vsite=$dd and water=$ww";
		    my $line = "";
		    print(LOG "****************************************************\n");
		    print(LOG "** PDB = $pdb FF = $ff VSITE = $dd WATER = $ww\n");
		    printf(LOG "** Working directory = %s\n", getcwd());
		    print(LOG "****************************************************\n");
		    mkdir("$ww");
		    chdir("$ww");
		    print(LOG "****************************************************\n");
		    print(LOG "**  Running pdb2gmx\n");
		    print(LOG "****************************************************\n");
		    open(PIPE,"$progs{'pdb2gmx'} -f ../../../../$pdb -ff $ff -ignh -vsite $dd -water $ww -o conf.g96 2>&1 |");
		    print LOG while <PIPE>;
		    close PIPE;
		    print(LOG "****************************************************\n");
		    print(LOG "**  Running editconf\n");
		    print(LOG "****************************************************\n");
		    open(PIPE,"$progs{'editconf'} -o b4em.g96 -box 5 5 5 -c -f conf.g96 2>&1 |");
		    print LOG while <PIPE>;
		    close PIPE;
		    print(LOG "****************************************************\n");
		    print(LOG "**  Running grompp\n");
		    print(LOG "****************************************************\n");
		    open(PIPE,"$progs{'grompp'} -maxwarn 3 -f ../../../../em -c b4em.g96 2>&1 |");
		    print LOG while <PIPE>;
		    close PIPE;
		    print(LOG "****************************************************\n");
		    print(LOG "**  Running mdrun\n");
		    print(LOG "****************************************************\n");
		    open(PIPE,$mdprefix->($mpi_ranks)." $progs{'mdrun'} $mdparams 2>&1 |");
		    print LOG while <PIPE>;
		    close PIPE;
		    chdir("..");
		}
		chdir("..");
	    }
	    chdir("..");
	}
	chdir("..");
    }
    close LOG;
    
    my $only_energies_filename = 'ener.log';
    my $nsuccess = find_in_file('Potential Energy',$logfn,$only_energies_filename);

    if ( $nsuccess != $ntest ) {
	print "Error not all $ntest pdb2gmx tests have been done successfully\n";
	print "Only $nsuccess energies in the log file\n";
	$nerror = 1;
    }
    else {
	my $reflog = "${ref}.log";
	if (! -f $reflog) {
	    print "No file $reflog. You are not really testing pdb2gmx\n";
	    copy($only_energies_filename, $reflog);
	}
	else {
	    print XML "<testsuite name=\"pdb2gmx\">\n" if ($xml);
	    $nerror = check_xvg($reflog,$only_energies_filename,3,\@pdb2gmx_test_names);
	    print XML "</testsuite>\n" if ($xml);
	    if ( $nerror != 0 ) {
		print "There were $nerror/$ntest differences in final energy with the reference file\n";
	    }
	}
    }
    if (0 == $nerror) {
	print "All $ntest pdb2gmx tests PASSED\n";
	remove_tree(@pdb_dirs);
	unlink($logfn,$only_energies_filename);
    }
    else {
	print "pdb2gmx tests FAILED\n";
    }
    chdir("..");
    return $nerror;
}

sub clean_all {
    map { cleandirs("$_") } @all_dirs;
    chdir("pdb2gmx");
    unlink("pdb2gmx.log");
    remove_tree(glob "pdb-*");
    chdir("..");
    cleandirs('essentialdynamics')
}

sub usage {
    my $dirs = join(' | ', @all_dirs);
    print <<EOP;
Usage: ./gmxtest.pl [ -np N ] [ -nt 1 ] [ -npme n ]
                    [ -verbose ] [ -double ] [ -bluegene ]
                    [ -suffix xxx ] [ -reprod ] [ -mpirun mpirun_command ]
                    [ -mdrun mdrun_command ]
                    [ -crosscompile ] [ -relaxed ] [ -tight ] [ -mdparam xxx ]
                    [ $dirs | pdb2gmx | all ]
or:    ./gmxtest.pl clean | refclean | dist
EOP
    exit 1;
}

# Since Perl's File::Which is not yet a standard module, there's no portable way
# to see whether a GROMACS tool can be found in the path, so we only attempt to
# check for the tool with this routine when not cross compiling
sub test_gmx {
  foreach my $p ( values %progs ) {
    if (system("$p -h > help 2>&1") != 0) {
      print("ERROR: Can not find executable $p in your path.\nPlease source GMXRC and try again.\n");
      exit 1;
    }
    unlink("help");
  }
}

my $kk = 0;
my @work = ();
my $did_clean = 0;

for ($kk=0; ($kk <= $#ARGV); $kk++) {
    my $arg = $ARGV[$kk];
    # If $arg is a subdirectory of tests, test that subdirectory
    if (grep(/^$arg$/, @all_dirs)) {
        push @work, "test_dirs('$arg')";
    }
    elsif ($arg eq 'pdb2gmx' ) {
	push @work, "test_pdb2gmx()";
    }
    elsif ($arg eq 'nm' || $arg eq 'normal-modes') {
        push @work, "test_normal_modes()";
    }
    elsif ($arg eq 'ed' || $arg eq 'essentialdynamics') {
        push @work, "test_essentialdynamics()";
    }
#    elsif ($arg eq 'tools' ) {
#	push @work, "test_tools()";
#    }
    elsif ($arg eq 'all' ) {
        map { push @work, "test_dirs('$_')" } @all_dirs;
	push @work, "test_pdb2gmx()";
	push @work, "test_normal_modes()";
	push @work, "test_essentialdynamics()";
	#push @work, "test_tools()";
    }
    elsif ($arg eq 'clean' ) {
        clean_all();
        $did_clean = 1;
    }
    elsif ($arg eq 'refclean' ) {
        map { push @work, "refcleandir('$_')" } @all_dirs;
	push @work, "unlink('pdb2gmx/reference_s.log','pdb2gmx/reference_d.log')";
    }
    elsif ($arg eq 'dist' ) {
	push @work, "clean_all()";
	push @work, "chdir('..')";
	push @work, "system('tar --exclude .git --exclude .gitattributes --exclude foreach.sh -czvf regressiontests.tgz regressiontests')";
	push @work, "chdir('regressiontests')";
    }
    elsif ($arg eq 'help' ) {
	usage();
    }
    elsif ($arg eq '-verbose') {
	$verbose++;
    }
    elsif ($arg eq '-noverbose') {
	$verbose = 0;
    }
    elsif ($arg eq '-xml') {
	$verbose = 0;
        $xml = 1;
    }
    elsif ($arg eq '-double') {
	$double = 1;
    }
    elsif ($arg eq '-crosscompile') {
        $crosscompiling = 1;
    }
    elsif ($arg eq '-bluegene') {
        $crosscompiling = 1;
	$bluegene = 1;
	print "Will test BlueGene. You should probably set '-mpirun runjob' if you have not already.\n";
    }
    elsif ($arg eq '-np' ) {
	if ($kk <$#ARGV) {
	    $kk++;
	    $mpi_ranks = $ARGV[$kk];
	    if ($mpi_ranks <= 0) {
		$mpi_ranks = 0;
	    } else {
                print "Will test on $mpi_ranks MPI ranks (if possible)\n";
            }
	}
    }
    elsif ($arg eq '-nt' ) {
	if ($kk <$#ARGV) {
	    $kk++;
	    $tmpi_ranks = $ARGV[$kk];
	    if ($tmpi_ranks <= 1) {
		$tmpi_ranks = 1;
                # most of the tests don't scale at all well
	    } else {
                print "Will test on $tmpi_ranks thread-MPI ranks (if possible)\n";
	    }
	}
    }
    elsif ($arg eq '-ntomp' ) {
	if ($kk <$#ARGV) {
	    $kk++;
	    $omp_threads = $ARGV[$kk];
	    if ($omp_threads <= 1) {
		$omp_threads = 1;
                # most of the tests don't scale at all well
	    } else {
                print "Will test on $omp_threads OpenMP threads (if possible)\n";
	    }
	}
    }
    elsif ($arg eq '-npme' ) {
        if ($kk < $#ARGV) {
            $kk++;
            $npme_ranks = $ARGV[$kk];
            print "Will run PME tests using $npme_ranks separate PME ranks (if possible)\n";
        }
    }
    elsif ($arg eq '-suffix' ) {
	if ($kk <$#ARGV) {
	    $kk++;
	    $suffix = $ARGV[$kk];
	    print "Will test using executable suffix $suffix\n";
	}
    }
    elsif ($arg eq '-nosuffix' ) {
	$autosuffix = 0;
    }
    elsif ($arg eq '-prefix' ) {
        print "Option -prefix is deprecated. Please update your script\n";
    }
    elsif ($arg eq '-reprod' ) {
      $mdparams.=" -reprod"
    }
    elsif ($arg eq '-mpirun' ) {
	if ($kk <$#ARGV) {
	    $kk++;
	    $mpirun = $ARGV[$kk];
	}
    }
    elsif ($arg eq '-mdrun' ) {
	if ($kk <$#ARGV) {
	    $kk++;
	    $mdrun_cmd = $ARGV[$kk];
	}
    }
    elsif ($arg eq '-mdparam' ) {
	# The user is responsible for providing sensible values for
	# this flag when they want them
	if ($kk <$#ARGV) {
	    $kk++;
	    $mdparams .= $ARGV[$kk];
	    print "Will test using 'mdrun $mdparams'\n";
	}
    }
    elsif ($arg eq '-gpu_id' ) {
	# The user is responsible for providing sensible values for
	# this flag when they want them (ie. as for mdrun)
	if ($kk <$#ARGV) {
	    $kk++;
	    $gpu_id .= $ARGV[$kk];
	    print "Will test using 'mdrun -gpu_id ${gpu_id}' (if possible)\n";
	}
    }
    elsif ($arg eq '-only' ) {
	# only test a subdirectory if it matches the following
	# regular expression
	if ($kk <$#ARGV) {
	    $kk++;
	    print "Will only test subdirectories matching regular expression '$ARGV[$kk]'\n";
	    $only_subdir = qr/$ARGV[$kk]/;
	}
    }
    elsif ($arg eq '-relaxed' ) {
	$tolerance_factor *= 10.0;
	print "Will run tests with 10x larger variations allowed\n";
    }
    elsif ($arg eq '-tight' ) {
        $tolerance_factor *= 0.1;
        print "Will run tests with 10x smaller variations allowed\n";
    }
    elsif ($arg eq '-parse' ) {
	if ($kk <$#ARGV) {
	    $kk++;
	    $parse_cmd = $ARGV[$kk];
	}
	print "Will parse mdrun output with $parse_cmd\n";
    }
    else {
	usage();
    }
}

die "FAILED: Cannot specify both -nt and -np greater than zero!" if ($mpi_ranks > 0 && $tmpi_ranks > 0);

# Set up to write XML output if wanted

if ($xml) {
    my $xmlfn = "gmxtest.xml";
    open (XML,">$xmlfn") || die("FAILED: Opening $xmlfn for writing");
    print XML '<?xml version="1.1" encoding="UTF-8"?>';
    print XML "\n<testsuites>\n";
}

# Set up the initial state for doing actual testing when we are, and
# print the usage statement if the command line was unsuitable

if (scalar(@work)) {
  # Prepend standard things to do to the list of real work,
  # which would be empty in the case of 'gmxtest.pl clean', etc.
  unshift(@work, "test_gmx()") unless ($crosscompiling);
  unshift(@work, "setup_vars()");
  # setup_vars() is always the first work to do, so now
  # command line options will work correctly regardless of
  # order on the command line
}
elsif (!$did_clean) {
    # There's multiple reasons why @work might be empty; only print
    # the usage if it's warranted.
    usage();
}

# Now do the work, if there is any
my $nerror = sum 0, map
{
    eval $_;
    if ('' ne $@) {
        die "$@"
    }
} @work;

print XML "</testsuites>\n" if ($xml);
exit ($nerror>0);
