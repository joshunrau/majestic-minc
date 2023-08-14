#! /usr/bin/env perl

#----------------------------- MNI Header -----------------------------------
#@NAME       : print_all_executables.pl
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Lists all the executables of ILT, both as they are invoked,
#              and as full path that the shell finds, by using 'which'.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : May. 22, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

    use  strict;
    use  warnings "all";
    use  ILT::Executables;

    my( $rcsid ) = '$Header: /private-cvsroot/libraries/ILT/print_all_executables.pl,v 1.4 2006-05-05 00:48:25 claude Exp $';

    my( @execs, $exec, $actual_path );

    @execs = get_executable_list();

    foreach $exec ( @execs )
    {
        $actual_path = `which $exec`;

        print( "$exec :\n",
               "\n",
               "        $actual_path\n" );
    }
