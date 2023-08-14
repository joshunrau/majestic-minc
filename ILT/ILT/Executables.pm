#!/usr/local/bin/perl5 -w

    package  ILT::Executables;
    use  strict;
    use  ILT::ProgUtils;

    BEGIN {
        use Exporter   ();
        use vars       qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

        @ISA         = qw(Exporter);
        %EXPORT_TAGS = ( );     # eg: TAG => [ qw!name1 name2! ],

        # your exported package globals go here,
        # as well as any optionally exported functions

        @EXPORT_OK      = ();
        @EXPORT         = qw( &set_executable_path &find_executable 
                              &run_executable &get_output_of_command
                              &get_executable_list );
    }

    my( $rcsid ) = '$Header: /private-cvsroot/libraries/ILT/ILT/Executables.pm,v 1.5 1998-09-18 13:30:03 david Exp $';

#--------------------------------------------------------------------------
#  Define all the executables or shell commands used by the library
#--------------------------------------------------------------------------

my( %executables ) =
(
    ray_trace               => "ray_trace",          #--- ray_trace command

    compute_bounding_view   => "compute_bounding_view", #--- conglomerate's
    make_slice              => "make_slice",
    place_images            => "place_images",
    plane_polygon_intersect => "plane_polygon_intersect",
    transform_objects       => "transform_objects",
    find_image_bounding_box => "find_image_bounding_box",

    xfmconcat               => "xfmconcat",          #--- /usr/local/mni/bin

    convert                 => "convert",            #--- in path
    mogrify                 => "mogrify",
    imgcopy                 => "imgcopy",
    imginfo                 => "imginfo",

    cp                      => "/bin/cp",            #--- shell
    touch                   => "touch",
    cat                     => "cat",
);

#--------------------------------------------------------------------------
# Initially the path is empty
#--------------------------------------------------------------------------

my( $path ) = "";

#----------------------------- MNI Header -----------------------------------
#@NAME       : set_executable_path
#@INPUT      : path string
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Sets the path to be prefixed to any executables defined above
#              which are relative paths (probably all of them are).
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : May.  6, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  set_executable_path( $ )
{
    my( $value ) = arg_string( shift );
    end_args( @_ );

    $path = $value;
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : find_executable
#@INPUT      : executable_base_name
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Converts an executable base name into a command to pass to the
#              shell.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : May.  6, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  find_executable( $ )
{
    my( $name ) = arg_string( shift );
    end_args( @_ );

    my( $full_path );

    $full_path = $executables{$name};

    if( !defined($full_path) )
    {
        fatal_error( "Could not find executable for $name\n" );
    }

    if( substr( $full_path, 0, 1 ) ne "/" && $path ne "" )
    {
        $full_path = "${path}/$full_path";
    }

    return( $full_path );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : run_executable
#@INPUT      : name
#              argument_string
#@OUTPUT     : 
#@RETURNS    : status
#@DESCRIPTION: Finds the executable corresponding to name, and runs it
#              with the given arguments.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : May.  6, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  run_executable( $$ )
{
    my( $name ) = arg_string( shift );
    my( $arguments ) = arg_string( shift );
    end_args( @_ );

    my( $status, $full_path );

    $full_path = find_executable( $name );

    $status = system_call( "$full_path $arguments" );

    return( $status );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : get_output_of_command
#@INPUT      : name
#              argument_string
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Finds the executable corresponding to name, and runs it
#              with the given arguments, returning the output of the command.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : May.  6, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  get_output_of_command( $$ )
{
    my( $name ) = arg_string( shift );
    my( $arguments ) = arg_string( shift );
    end_args( @_ );

    my( $output, $full_path );

    $full_path = find_executable( $name );

    $output = `$full_path $arguments`;

    return( $output );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : get_executable_list
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Lists all the executables used by the ILT package.  Not
#              for general use.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : May. 22, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  get_executable_list()
{
    return( values(%executables) );
}

#---------------- return value for 'use'

1;

