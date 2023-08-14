use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'	=> 'ILT',
    'VERSION'   => '1.2.4',
    'EXE_FILES' => [qw(script/create_verify_image) ],
    'dist'      => { COMPRESS=>"gzip", SUFFIX=>"gz" }
);
