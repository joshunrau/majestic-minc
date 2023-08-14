/* ----------------------------- MNI Header -----------------------------------
@NAME       : ext_element_defs.h
@DESCRIPTION: Element definitions for extra elements needed for mosaics, etc.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : December 2001 (Rick Hoge)
@MODIFIED   : 
@COPYRIGHT  :
              Copyright 1993 Peter Neelin, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

/* Element id's for EXT */
/* bert- These appear to be completely nonstandard - I believe they were 
 * created by Peter and Rick to facilitate communication among the various
 * pieces of Siemens Mosaic handling code.  They do not represent any 
 * externally defined DICOM standard and therefore they may conflict with 
 * other manufacturer's proprietary fields.
 */
GLOBAL_ELEMENT(EXT_Mosaic_rows                        , 0x0023, 0x0001, LO);
GLOBAL_ELEMENT(EXT_Mosaic_columns                     , 0x0023, 0x0002, LO);
GLOBAL_ELEMENT(EXT_Slices_in_file                     , 0x0023, 0x0003, LO);
GLOBAL_ELEMENT(EXT_Sub_image_rows                     , 0x0023, 0x0004, US);
GLOBAL_ELEMENT(EXT_Sub_image_columns                  , 0x0023, 0x0005, US);
GLOBAL_ELEMENT(EXT_MrProt_dump                        , 0x0023, 0x0006, LO);
GLOBAL_ELEMENT(EXT_Diffusion_b_value                  , 0x0023, 0x0007, LO);
GLOBAL_ELEMENT(EXT_Delay_in_TR                        , 0x0023, 0x0008, LO); /*add for fMRI scans*/
/*Can't find slice acquisition in any standard fields; using ASCONV header sSliceArray.ucMode  ilana*/

/*0x1 ascending 0x2 descending 0x3 interleaved*/
GLOBAL_ELEMENT(EXT_Slice_order		      , 0x0023, 0x0009, LO);

/* Non-zero if slice order is inverted. */ 
GLOBAL_ELEMENT(EXT_Slice_inverted             , 0x0023, 0x000A, LO);

/* Slice orientation: 1 -> SAGITTAL, 2 -> CORONAL, 0 -> TRANSVERSE */
GLOBAL_ELEMENT(EXT_Slice_orientation          , 0x0023, 0x000B, LO);
