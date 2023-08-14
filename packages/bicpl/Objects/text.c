/* ----------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1993,1994,1995 David MacDonald,
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

#include "bicpl_internal.h"

/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_text
@INPUT      : position
              colour
              text_font
              size
@OUTPUT     : text
@RETURNS    : 
@DESCRIPTION: Initializes the text structure to empty.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void   initialize_text(
    text_struct     *text,
    VIO_Point           *position,
    VIO_Colour          colour,
    Font_types      text_font,
    VIO_Real            size )
{
    if( position != NULL )
        text->origin = *position;
    else
        fill_Point( text->origin, 0.0, 0.0, 0.0 );

    text->colour = colour;
    text->font = text_font;
    text->size = size;

    text->string = create_string( NULL );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_text
@INPUT      : text
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the text object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sep. 25, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void   delete_text(
    text_struct     *text )
{
    delete_string( text->string );
}
