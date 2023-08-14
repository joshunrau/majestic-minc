#include  <volume_io.h>

  VIO_Colour  make_rgba_Colour(
    int    r,
    int    g,
    int    b,
    int    a )
{
    VIO_Colour          c;
    unsigned  char  *byte_ptr;

    c = 0;     /* to avoid used-before-set compiler messages */

    byte_ptr = (void *) &c;

    byte_ptr[0] = (unsigned char) a;
    byte_ptr[1] = (unsigned char) b;
    byte_ptr[2] = (unsigned char) g;
    byte_ptr[3] = (unsigned char) r;

    return( c );
}

  int  get_Colour_r(
    VIO_Colour   colour )
{
    unsigned  char  *b;

    b = (void *) &colour;

    return( (int) b[3] );
}

  int  get_Colour_g(
    VIO_Colour   colour )
{
    unsigned  char  *b;

    b = (void *) &colour;

    return( (int) b[2] );
}

  int  get_Colour_b(
    VIO_Colour   colour )
{
    unsigned  char  *b;

    b = (void *) &colour;

    return( (int) b[1] );
}

  int  get_Colour_a(
    VIO_Colour   colour )
{
    unsigned  char  *b;

    b = (void *) &colour;

    return( (int) b[0] );
}
