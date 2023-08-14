dnl Copyright (C) 1988 Eleftherios Gkioulekas <lf@amath.washington.edu>
dnl  
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2 of the License, or
dnl (at your option) any later version.
dnl 
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl 
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software 
dnl Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
dnl 
dnl As a special exception to the GNU General Public License, if you 
dnl distribute this file as part of a program that contains a configuration 
dnl script generated by Autoconf, you may include it under the same 
dnl distribution terms that you use for the rest of that program.
 
# ------------------------------------------------------------------------
# The following macro is useful for deep packages. It allows you to
# link all the header files *.h under a certain set of directories
# to be linked under an include directory from the toplevel.
# To use this feature in your configure.in call:
#   LF_LINK_HEADERS(dir1 dir2 dir3 .... , [directory] )
# where directory -> put links under include/directory
#       dir1 ...  -> the directories with header files we want to link
# WARNING: This macro will do  --> rm -rf include
# ------------------------------------------------------------------------

AC_DEFUN([LF_LINK_HEADERS],
[
  AC_REQUIRE([AC_PROG_LN_S])

  # Remove the include directory exactly once.
  if test -z "$lf_link_headers"
  then
    lf_link_headers="we are all Kosh"
    rm -rf "include"
  fi

  # Get the directory from the second argument which is optional
  ifelse([$2], ,  
         [lf_directory="include"] , 
         [lf_directory="include/$2"])
  AS_MKDIR_P( "$lf_directory" )

  # Link them
  lf_subdirs="`echo $1`"
  for lf_dir in $lf_subdirs
  do
    echo "linking headers from $srcdir/$lf_dir"

    if test -f "$srcdir/$lf_dir/Headers"
    then
      lf_headers=`(cd $srcdir/$lf_dir; cat Headers)`
    else
      lf_headers=`(cd $srcdir/$lf_dir; echo *.h)`
    fi

    for lf_file in $lf_headers
    do
      rm -f $lf_directory/$lf_file
      $LN_S "`pwd`/$srcdir/$lf_dir/$lf_file" "$lf_directory/$lf_file"
    done
  done
])
