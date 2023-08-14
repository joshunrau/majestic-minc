#ifdef HAVE_CONFIG_H
#include "config.h"
#endif //HAVE_CONFIG_H  

#include  <volume_io.h>
#include  <graphics.h>


/* Modified: October 29, 2010  JGS

These numbers appear to be used to seed a linear feedback shift register type random number generator.  The 22bit case originally had poor run length.  The new value of 0x420000 is attributed to Wayne Stahnke while he was at Fairchild Semiconductor in 1970.  See http://www.xilinx.com/support/documentation/application_notes/xapp052.pdf */

static  Random_mask_type  random_order_masks[] = {
                                                 0x00000000,
                                                 0x00000000,
                                                 0x00000006,
                                                 0x0000000c,
                                                 0x00000014,
                                                 0x00000030,
                                                 0x00000060,
                                                 0x000000b8,
                                                 0x00000110,
                                                 0x00000240,
                                                 0x00000500,
                                                 0x00000ca0,
                                                 0x00001b00,
                                                 0x00003500,
                                                 0x00006000,
                                                 0x0000b400,
                                                 0x00012000,
                                                 0x00020400,
                                                 0x00072000,
                                                 0x00090000,
                                                 0x00140000,
                                                 0x00300000,
                                                 0x00420000,
                                                 0x00d80000,
                                                 0x01200000,
                                                 0x03880000,
                                                 0x07200000,
                                                 0x09000000,
                                                 0x14000000,
                                                 0x32800000,
                                                 0x48000000,
                                                 0xa3000000
                                               };

  Random_mask_type  get_random_order_mask(
    int  size )
{
    int  ind;

    ind = 0;

    while( size != 0 )
    {
        size >>= 1;
        ++ind;
    }

    return( random_order_masks[ind] );
}
