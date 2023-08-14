#!/usr/local/bin/perl5 -w

# ----------------------------------------------------------------------------
#@COPYRIGHT  :
#              Copyright 1993,1994,1995,1996,1997,1998 David MacDonald,
#              McConnell Brain Imaging Centre,
#              Montreal Neurological Institute, McGill University.
#              Permission to use, copy, modify, and distribute this
#              software and its documentation for any purpose and without
#              fee is hereby granted, provided that the above copyright
#              notice appear in all copies.  The author and McGill University
#              make no representations about the suitability of this
#              software for any purpose.  It is provided "as is" without
#              express or implied warranty.
#-----------------------------------------------------------------------------


#----------------------------- MNI Header -----------------------------------
#@NAME       : LayoutInclude.pm
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: An include file which holds all object types for the 
#              Image Layout Tool
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------


    use  strict;

    use ILT::ImageLayout;
    use ILT::ImageInfo;
    use ILT::View;
    use ILT::SceneObject;
    use ILT::SceneObject::ColourObject;
    use ILT::SceneObject::GeometricObject;
    use ILT::SceneObject::IntersectionObject;
    use ILT::SceneObject::OneSubObject;
    use ILT::SceneObject::PlaneObject;
    use ILT::SceneObject::RenderObject;
    use ILT::SceneObject::TextObject;
    use ILT::SceneObject::TransformObject;
    use ILT::SceneObject::UnionObject;
    use ILT::SceneObject::VolumeObject;


1;
