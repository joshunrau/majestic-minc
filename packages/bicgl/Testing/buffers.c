#include  <gl.h>

main()
{
    printf( "Double RGBA: %d %d %d %d\n",
            getgdesc( GD_BITS_NORM_DBL_RED ),
            getgdesc( GD_BITS_NORM_DBL_GREEN ),
            getgdesc( GD_BITS_NORM_DBL_BLUE ),
            getgdesc( GD_BITS_NORM_DBL_ALPHA ) );
    printf( "Single RGBA: %d %d %d %d\n",
            getgdesc( GD_BITS_NORM_SNG_RED ),
            getgdesc( GD_BITS_NORM_SNG_GREEN ),
            getgdesc( GD_BITS_NORM_SNG_BLUE ),
            getgdesc( GD_BITS_NORM_SNG_ALPHA ) );

    printf( "Double Colour map: %d\n",
            getgdesc( GD_BITS_NORM_DBL_CMODE ) );

    printf( "Single Colour map: %d\n",
            getgdesc( GD_BITS_NORM_SNG_CMODE ) );

    printf( "Depth Buffer: %d\n",
            getgdesc( GD_BITS_NORM_ZBUFFER ) );

    printf( "Blending: %d\n",
            getgdesc( GD_BLEND ) );

    printf( "A function: %d\n",
            getgdesc( GD_AFUNCTION ) );

    printf( "Texture mapping: %d\n",
            getgdesc( GD_TEXTURE ) );
}
