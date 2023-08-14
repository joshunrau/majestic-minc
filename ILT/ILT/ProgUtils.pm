#!/usr/local/bin/perl5 -w

    package  ILT::ProgUtils;
    use  strict;
    use  Carp;
	use File::Temp qw/ tempdir tempfile /;
    use UNIVERSAL;


    BEGIN {
        use Exporter   ();
        use vars       qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

        @ISA         = qw(Exporter);
        %EXPORT_TAGS = ( );     # eg: TAG => [ qw!name1 name2! ],

        # your exported package globals go here,
        # as well as any optionally exported functions

        @EXPORT_OK      = ();
        @EXPORT         = qw(
                          &end_args
                          &arg_any &arg_int &arg_real &arg_enum &arg_object
                          &arg_array_of_ints
                          &arg_array_of_reals
                          &arg_array &arg_string
                          &arg_array_ref
                          &opt_arg_int &opt_arg_real &opt_arg_enum
                          &opt_arg_object &opt_arg_string
                          &opt_arg_array_of_reals
                          &get_tmp_file &system_call
                          &False
                          &True
                          &N_boolean_enums
                          &min
                          &max
                          &error
                          &warning
                          &fatal_error
                        );
    }

    my( $rcsid ) = '$Header: /private-cvsroot/libraries/ILT/ILT/ProgUtils.pm,v 1.3 2011-02-04 16:48:13 alex Exp $';

#---------------------------  error function

sub  error( @ )
{
    my( @args ) = @_;

    print( "Error: ", @args, "\n" );
}

sub  warning( @ )
{
    my( @args ) = @_;

    print( "Warning: ", @args, "\n" );
}

sub  fatal_error( @ )
{
    my( @args ) = @_;

    error( @args );
    confess( @_ , "\n" );
}

#---------------------------  Perl argument checking

sub  is_scalar( $ )
{
    my( $s ) = @_;

    return( !ref($s) );
}

sub  is_numeric( $ )
{
    my( $number ) = @_;

    return( 1 );       #--- insert numeric test here
}

sub  is_int( $ )
{
    my( $i ) = @_;

    return( is_numeric($i) && int($i) == $i );
}

sub  arg_object( $$ )
{
    my( $obj, $expected_class ) = @_;

    if( !UNIVERSAL::isa( $obj, $expected_class ) )
    {
        fatal_error(
           "Expected an object argument of class: $expected_class\n" );
    }

    return( $obj );
}

sub  arg_enum( $$ )
{
    my( $value, $n_enum ) = @_;

    if( !defined($value) )
        { fatal_error( "Expected an argument, but got <undef>\n" ); }

    if( !is_scalar( $value ) || $value < 0 || $value >= $n_enum ||
        !is_int( $value ) )
    {
        fatal_error( "Expected an enum, but got <$value>\n" );
    }

    return( $value );
}

sub  arg_int( $@ )
{
    my( $value, $min_value, $max_value ) = @_;

    if( !defined($value) )
        { fatal_error( "Expected an argument, but got <undef>\n" ); }

    if( !is_scalar( $value ) || !is_int($value) )
        { fatal_error( "Expected an integer, but got <$value>\n" ); }

    if( defined($min_value) && !defined($max_value) && $value < $min_value )
    {
        fatal_error( "Expected an int in the range > $min_value, " .
               "but got a value of $value\n" );
    }

    if( defined($min_value) && defined($max_value) &&
        ($value < $min_value || $value > $max_value) )
    {
        fatal_error(
               "Expected an integer in the range $min_value to $max_value, " .
               "but got a value of $value\n" );
    }

    return( $value );
}

sub  arg_real( $@ )
{
    my( $value, $min_value, $max_value ) = @_;

    if( !defined($value) )
        { fatal_error( "Expected an argument, but got <undef>\n" ); }

    if( !is_scalar($value) || !is_numeric($value) )
        { fatal_error( "Expected a real, but got a non-scalar <$value>\n" ); }

    if( defined($min_value) && !defined($max_value) && $value < $min_value )
    {
        fatal_error( "Expected a real in the range > $min_value, " .
               "but got a value of $value\n" );
    }

    if( defined($min_value) && defined($max_value) &&
        ($value < $min_value || $value > $max_value) )
    {
        fatal_error( "Expected a real in the range $min_value to $max_value, " .
               "but got a value of $value\n" );
    }

    return( $value );
}

sub  arg_any( $ )
{
    my( $value ) = @_;

    if( !defined($value) )
        { fatal_error( "Expected an argument, but got <undef>\n" ); }

    return( $value );
}

sub  arg_string( $ )
{
    my( $value ) = @_;

    if( !defined($value) )
        { fatal_error( "Expected an argument, but got <undef>\n" ); }

    if( !is_scalar($value) )
        { fatal_error( "Expected a string, but got a non-scalar <$value>\n" ); }

    return( $value );
}

sub  opt_arg_int( $@ )
{
    my( $value ) = @_;

    if( !defined($value) )
        { return( $value ); }
    else
        { return( arg_int( @_ ) ); }
}

sub  opt_arg_real( $@ )
{
    my( $value, @remaining ) = @_;

    if( !defined($value) )
        { return( $value ); }
    else
        { return( arg_real( $value, @remaining ) ); }
}

sub  opt_arg_string( $@ )
{
    my( $value ) = @_;

    if( !defined($value) )
        { return( $value ); }
    else
        { return( arg_string( $value ) ); }
}

sub  opt_arg_enum( $@ )
{
    my( $value, $max_enum ) = @_;

    if( !defined($value) )
        { return( $value ); }
    else
        { return( arg_enum( $value, $max_enum ) ); }
}

sub  opt_arg_object( $$ )
{
    my( $value, $class ) = @_;

    if( !defined($value) )
        { return( $value ); }
    else
        { return( arg_object( $value, $class ) ); }
}

sub  arg_array_of_reals( $@ )
{
    my( $arg_ref, $n ) = @_;

    my( @list, $count, $arg );

    if( ref($arg_ref) ne "ARRAY" )
        { fatal_error( "arg_array_of_reals() first argument should be \@_" ); }

    @list = ();

    $count = 0;
    while( defined($arg = shift(@$arg_ref)) && (!defined($n) || $count < $n) )
    {
        push( @list, arg_real($arg) );
        ++$count;
    }

    if( defined($n) && $count != $n )
    {
        { fatal_error( "arg_array_of_reals() incorrect number of args.\n" ); }
    }

    return( @list );
}

sub  opt_arg_array_of_reals( $@ )
{
    my( $arg_ref, $n ) = @_;

    if( @$arg_ref == 0 )
        { return(); }
    else
        { return( arg_array_of_reals( $arg_ref, $n ) ); }
}

sub  arg_array_of_ints( $@ )
{
    my( $arg_ref, $n ) = @_;

    my( @list, $count, $arg );

    if( ref($arg_ref) ne "ARRAY" )
        { fatal_error( "arg_list_of_ints() first argument should be \@_" ); }

    $count = 0;
    while( defined($arg = shift(@$arg_ref)) && (!defined($n) || $count < $n) )
    {
        push( @list, arg_int($arg) );
        ++$count;
    }

    if( defined($n) && $count != $n )
    {
        { fatal_error( "arg_array_of_ints() incorrect number of args.\n" ); }
    }

    return( @list );
}

sub  arg_array( $@ )
{
    my( $arg_ref, $n ) = @_;

    my( $arg, @list, $count );

    if( ref($arg_ref) ne "ARRAY" )
        { fatal_error( "arg_list_of_ints() first argument should be \@_" ); }

    $count = 0;

    while( defined($n) && $count < $n || !defined($n) && @{$arg_ref} )
    {
        $arg = shift( @$arg_ref );

        push( @list, $arg );
    }

    if( defined( $n ) && $count != $n )
        { fatal_error( "Expected $n args for arg_array()" ); }

    return( @list );
}

sub  arg_array_ref( $@ )
{
    my( $value, $n_elements ) = @_;

    if( ref($value) ne "ARRAY" )
        { fatal_error( "Expected an array ref, but got <$value>\n" ); }

    if( defined($n_elements) && @$value != $n_elements )
    {
        my( $size ) = @$value;
        fatal_error( "Expected an array of size $n_elements, but got $size\n" );
    }

    return( $value );
}

sub  end_args( @ )
{
    my( @args ) = @_;

    if( @args )
    {
        fatal_error(
             "Found extra arguments to a function call: <", qw(@args), ">\n");
    }
}

#--------------------------

sub  get_filename_base( $ )
{
    my( $file ) = arg_string( shift );
    end_args( @_ );

    $file =~ s/.*\///;
    $file =~ s/\..*//;

    return( $file );
}

sub  get_prefix( $ )
{
    my( $file ) = arg_string( shift );
    end_args( @_ );

    $file =~ s/\.[^\.]*$//;   

    return( $file );
}

sub  copy_file( $$ )
{
    my( $src ) = arg_string( shift );
    my( $dest ) = arg_string( shift );
    end_args( @_ );

    run_executable( "cp", "$src $dest" );
}

sub get_directory( $ )
{
    my( $filename ) = arg_string( shift );
    end_args( @_ );

    my( $dir );

    $dir = $filename;

    $dir =~ s/[^\/]*$//;

    if( ! $dir )  { $dir = "."; }

    if( substr( $dir, -1, 1 ) ne "/" )
        { $dir = $dir . "/"; }

    return( $dir );
}

sub get_filename_no_dirs( $ )
{
    my( $filename ) = arg_string( shift );
    end_args( @_ );

    my( $no_dirs );

    $no_dirs = $filename;

    $no_dirs =~ s/.*\///;

    if( ! $no_dirs )  { $no_dirs = "."; }

    return( $no_dirs );
}

#--------------- tmp files -------------------

my $tmp_dir=tempdir( "ILT_${$}_XXXXX", TMPDIR => 1, CLEANUP => 1 );

sub get_tmp_file( $ )
{
    my( $suffix ) = arg_string( shift );
    end_args( @_ );

    if( defined($suffix) )
        { $suffix = "." . $suffix; }
    else
        { $suffix = ""; }

    my($fh, $tmp_file) = tempfile( "ILT_XXXXX", DIR=>$tmp_dir, SUFFIX=>$suffix );

    return( $tmp_file );
}

#-------------------------

sub system_call( $ )
{
    my( $command ) = arg_string( shift );
    end_args( @_ );

    my( $ret, @separate, $com );

    print( "$command\n" );
    $ret = system( $command );
    if( $ret != 0 )
    {
        @separate = split( /\s+/, $command );
        $com = $separate[0];
        if( $ret == 2 )
            { confess( "System command <$com> was interrupted.\n" ); }
        elsif( $ret == 65280 )
            { confess( "System command <$com> was not found.\n" ); }
        else
            { confess( "System command <$com> failed with return value
 <$ret>.\n" ); }
    }
    return( $ret / 256 );
}

#----- Boolean enums

sub  False { return( 0 ); }
sub  True  { return( 1 ); }
sub  N_boolean_enums  { return( 2 ); }

#----- min and max

sub  min( @ )
{
    my( @args ) = arg_array_of_reals( \@_ );
    end_args( @_ );

    my( $arg, $min );

    $min = undef;

    foreach $arg ( @args )
    {
        if( !defined($min) || $arg < $min )
            { $min = $arg; }
    }

    return( $min );
}

sub  max( @ )
{
    my( @args ) = arg_array_of_reals( \@_ );
    end_args( @_ );

    my( $arg, $max );

    $max = undef;

    foreach $arg ( @args )
    {
        if( !defined($max) || $arg > $max )
            { $max = $arg; }
    }

    return( $max );
}

#---------------- return value for 'use'

1;

