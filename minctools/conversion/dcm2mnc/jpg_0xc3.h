//This decodes DICOM images with Transfer Syntax 1.2.840.10008.1.2.4.70
// "JPEG Lossless, Nonhierarchical, First- Order Prediction"
// This format is described in http://www.w3.org/Graphics/JPEG/itu-t81.pdf
// It is identified with the 'Start of Frame' (SOF) code 0xC3
// It appears unique to medical imaging, and is not supported by most JPEG libraries
// http://www.dicomlibrary.com/dicom/transfer-syntax/
#ifndef _JPEG_SOF_0XC3_
#define _JPEG_SOF_0XC3_

uint8_t *decode_jpeg_sof_0xc3(uint8_t *lRawRA, int lRawSz, int verbose, int *dimX, int *dimY, int *bits, int *frames);

#endif
