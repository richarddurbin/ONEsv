/*  File: alncode.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description: IO for ONEcode .1aln files for Myers FASTGA package
 * Exported functions:
 * HISTORY:
 * Last edited: Aug 10 00:16 2024 (rd109)
 * Created: Sat Feb 24 12:19:16 2024 (rd109)
 *-------------------------------------------------------------------
 */

#include "ONElib.h"
#include "align.h"

OneSchema *alnMakeSchema();

// open the .1aln file for reading and read the header

OneFile *alnOpenRead (char *filename,  int nThreads,
		      I64  *nOverlaps, int *tspace,
		      char **db1_name, char **db2_name, char **cpath) ;

// next two routines read the records from the file

void alnReadOverlap (OneFile *of, Overlap *ovl);
int  alnReadTrace   (OneFile *of, U8 *trace);
void alnSkipTrace   (OneFile *of);

// and equivalents for writ1, ing

OneFile *alnOpenWrite (char *filename, int nThreads,
		       char *progname, char *version, char *commandLine,
		       int tspace, char *db1_name, char *db2_name, char *cpath);

void alnWriteOverlap (OneFile *of, Overlap *ovl);
void alnWriteTrace   (OneFile *of, U8 *trace, int tlen);

// end of file
