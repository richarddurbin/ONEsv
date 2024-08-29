/*  File: alnseq.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Aug 18 11:52 2024 (rd109)
 * Created: Tue Aug 13 14:34:13 2024 (rd109)
 *-------------------------------------------------------------------
 */

#include "alnseq.h"
#include "ONElib.h"

static char *gdbSchemaText =
  "1 3 def 1 0                 schema for genome skeleton\n"
  ".\n"
  "P 3 gdb                     GDB\n"
  "D f 4 4 REAL 4 REAL 4 REAL 4 REAL   global: base frequency vector\n"
  "O S 1 6 STRING              id for a scaffold\n"
  "D G 1 3 INT                 gap of given length\n"
  "D C 1 3 INT                 contig of given length\n"
;

static bool isACGT[256] ;

AlnSeq *alnSeqOpen (char *name, char *cpath, bool isIndexRequired) // open for read
{
  { bzero (isACGT, 256*sizeof(bool)) ;
    isACGT['a'] = true ; isACGT['c'] = true ; isACGT['g'] = true ; isACGT['t'] = true ;
    isACGT['A'] = true ; isACGT['C'] = true ; isACGT['G'] = true ; isACGT['T'] = true ;
  }    
  
  AlnSeq *as = new0 (1, AlnSeq) ;

  char   *fullPath = new(strlen(name) + strlen(cpath) + 2, char) ;
  strcpy (fullPath, cpath) ; strcat (fullPath, "/") ; strcat (fullPath, name) ;

  // first check whether this is a 1gdb file - if so find the parental DNA file
  OneFile *of = oneFileOpenRead (name, 0, "gdb", 1) ;
  if (!of) of = oneFileOpenRead (fullPath, 0, "gdb", 1) ;
  if (of)
    { int n = of->info['<']->accum.count ; // number of reference lines
      name = 0 ;
      OneReference *r = of->reference ;
      while (n--) if (r->count == 1) { name = r->filename ; break ; }
      if (!name) die ("failed to find reference name in GDB file %s", fullPath) ;
      free (fullPath) ; fullPath = new(strlen(name) + strlen(cpath) + 2, char) ;
      strcpy (fullPath, cpath) ; strcat (fullPath, "/") ; strcat (fullPath, name) ;
    }
  
  as->si = seqIOopenRead (name, dna2textConv, false) ;
  if (!as->si) as->si = seqIOopenRead (fullPath, dna2textConv, false) ;
  if (!as->si) die ("failed to open sequence file %s or %s", name, fullPath) ;

  if (seqIOread (as->si))
    as->seq = sqioSeq(as->si) ;
  else
    { alnSeqClose (as) ;
      as = 0 ;
    }

  return as ;
}

char* alnSeqNext (AlnSeq *as, U64 *len) // DNA text (acgt) for next (contig) sequence
{
  while (as->inSeq == as->si->seqLen)
    { if (!seqIOread (as->si)) return 0 ;
      as->inSeq = 0 ;
    }
  char *s = sqioSeq (as->si) + as->inSeq ;
  char *t = s, *tMax = sqioSeq(as->si) + as->si->seqLen ;
  while (t < tMax && isACGT[(int)*t]) ++t ;
  *len = t - s ;
  while (t < tMax && !isACGT[(int)*t]) ++t ;
  as->inSeq = t - sqioSeq (as->si) ;
  return s ;
}

void alnSeqClose (AlnSeq *as)
{
  seqIOclose (as->si) ;
  free (as) ;
}

/* the next two are only supported if there is an index - not implemented yet

char* alnSeq (AlnSeq *as, int i, int *len) // DNA text (acgt) for i'th (contig) sequence
{
  if (!as->nmax) die ("alnSeq requires index - must build with isIndexRequired true") ;
  if (iseq != as->parent[i])
    
}

bool alnSeqLoc (AlnSeq *as, int i, int x, int *s, int *sx) // source (scaffold) coords for 1aln i,x
{
  if (!as->nmax) die ("alnSeqLoc requires index - must build with isIndexRequired true") ;
  if (i < 0 || i >= as->nseq)
    { warn ("alnSeqLoc i %d is out of bounds [0,%d)", i, as->nseq) ; return false ; }
  if (s) *s = as->parent[i] ;
  if (x < 0 || x > as->len[i])
    { warn ("alnSeqLoc pos %d is out of bounds [0,%d]", i, as->len) ; return false ; }
  if (sx) *sx = as->offset[i] + x ;
  return true ;
}

*/

// end of file
