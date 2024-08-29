/*  File: alnseq.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Aug 18 11:14 2024 (rd109)
 * Created: Tue Aug 13 14:35:58 2024 (rd109)
 *-------------------------------------------------------------------
 */

#include "seqio.h"

typedef struct {
  SeqIO *si ;
  char  *seq ;
  int    inSeq ;
} AlnSeq ;

AlnSeq *alnSeqOpen (char *name, char *cpath, bool isIndexRequired) ; // open for read
char* alnSeqNext (AlnSeq *as, U64 *len) ; // DNA text (acgt) for next (contig) sequence
void alnSeqClose (AlnSeq *as) ;

/******** end of file *********/
