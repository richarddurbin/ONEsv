/*  File: svfind.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Aug 29 14:52 2024 (rd109)
 * Created: Fri Aug  9 22:41:44 2024 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "array.h"
#include "alncode.h" // includes ONElib.h and align.h
#include "alnseq.h"  // includes ONElib.h and align.h

#define PROG_NAME "svfind"
#define VERSION "0.1"

static char* schemaText =
  "1 3 def 2 1               schema for structural variants\n"
  ".                         expects the following reference lines\n"
  ". < a_file 1              source DNA file for insertions or duplications\n"
  ". < b_file 2              source DNA file for corresponding deletions/single copies (if not a)\n"
  ". < c_path 3              directory for a_file, b_file if their names are not absolute paths\n"
  ".                         \n"
  "P 3 seq                   SEQUENCE\n"
  "S 2 sv                    SEQUENCE VARIANT\n"
  "D o 1 3 INT               maximum overhang (global)\n"
  "D i 1 3 INT               maximum insert size (global)\n"
  ".                         \n"
  "O V 3 3 INT 3 INT 3 INT   variant: seqid, start, end (0-indexed, [start,end))\n"
  "D O 1 3 INT               overlap\n"
  "D B 3 3 INT 3 INT 3 INT   source, start-match, end-match\n"
  "D D 1 3 DNA               target site duplication (TSD): sequence\n"
  "D R 2 3 INT 3 INT         Terminal Inverted Repeat (TIR): length, number of mismatches\n"
  "D T 2 3 INT 3 INT         Long Terminal direct Repeat (LTR): length, number of mismatches\n"
  "G S                       insertions group sequences\n"
  ".\n"
  "O S 1 3 DNA               sequence of the insertion\n"
  "D I 1 6 STRING            identifier of the insertion\n"
  ;

static int MAX_OVERHANG = 50 ;
static int MAX_SIZE = 50000 ;

void usage (void)
{
  fprintf (stderr, "Usage: svfind [opts] <1alnFileName>\n") ;
  fprintf (stderr, "opts:     -w <int>         maximum overhang\n") ;
  fprintf (stderr, "          -m <int>         maximum length\n") ;
  fprintf (stderr, "          -a <filename>    outfile for insertions/duplications in a\n") ;
  fprintf (stderr, "          -b <filename>    outfile for insertions/duplications in b\n") ;
  
  exit (1) ;
}

static inline void flip (Overlap *o1, Overlap *o2) // must be safe for o2 == o1
{
  int t ;
  o2->flags = o1->flags ; // OK
  t = o1->aread ; o2->aread = o1->bread ; o2->bread = t ;
  t = o1->path.abpos ; o2->path.abpos = o1->path.bbpos ; o2->path.bbpos = t ; 
  t = o1->path.aepos ; o2->path.aepos = o1->path.bepos ; o2->path.bepos = t ; 
}

static int overlapOrder (const void *x, const void *y) // sort on b, a, bbpos
{ Overlap *ox = (Overlap*)x, *oy = (Overlap*)y ;
  int c = ox->bread - oy->bread ; if (c) return c ;
  c = ox->aread - oy->aread ; if (c) return c ;
  c = ox->path.bbpos - oy->path.bbpos ; return c ;
}

void insertionReport (OneFile *of, AlnSeq *as, Overlap *olap, int n) ;

int main (int argc, char *argv[])
{
  storeCommandLine (argc--, argv++) ;
  timeUpdate (0) ;

  OneSchema *schema = oneSchemaCreateFromText (schemaText) ;
  OneFile   *ofa = 0, *ofb = 0 ;
  char      *ofaName, *ofbName ;

  if (!argc) usage () ;
  
  while (argc > 1)
    if (!strcmp (*argv, "-w") && argc > 2)
      { if ((MAX_OVERHANG = atoi(argv[1])) <= 0)
	  die ("max_overhang %s must be a positive integer", argv[1]) ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-m") && argc > 2)
      { if ((MAX_SIZE = atoi(argv[1])) <= 0)
	  die ("max_size %s must be a positive integer", argv[1]) ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-a") && argc > 2)
      { if (!(ofa = oneFileOpenWriteNew (argv[1], schema, "sv", true, 1)))
	  die ("failed to open .1insert file %s to write", argv[1]) ;
	oneAddProvenance (ofa, PROG_NAME, VERSION, getCommandLine()) ;
	ofaName = argv[1] ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-b") && argc > 2)
      { if (!(ofb = oneFileOpenWriteNew (argv[1], schema, "sv", true, 1)))
	  die ("failed to open .1insert file %s to write", argv[1]) ;
	oneAddProvenance (ofb, PROG_NAME, VERSION, getCommandLine()) ;
	ofbName = argv[1] ;
	argc -= 2 ; argv += 2 ;
      }
    else
      { warn ("unknown option %s", *argv) ;
	usage () ;
      }

  I64      nOverlaps ;
  char    *db1Name = 0, *db2Name = 0, *cpath = 0 ;
  OneFile *ofIn = alnOpenRead (*argv, 1, &nOverlaps, 0, &db1Name, &db2Name, &cpath) ;
    
  if (!ofIn) die ("failed to open .1aln file %s", *argv) ;

  AlnSeq *as = 0, *bs = 0 ;
  
  if (ofa)
    { oneAddReference (ofa, db1Name, 1) ;
      if (db2Name) oneAddReference (ofa, db2Name, 2) ;
      if (cpath) oneAddReference (ofa, cpath, 3) ;
      if (!(as = alnSeqOpen (db1Name, cpath, false))) die ("failed to open %s", db1Name) ;
    }

  if (ofb)
    { if (!db2Name)
	die ("-b not possible: input %s has no b source (it has self-a alignments only)", *argv) ;
      oneAddReference (ofb, db2Name, 1) ; // NB change of order here
      oneAddReference (ofb, db1Name, 2) ;
      oneAddReference (ofb, cpath, 3) ;
      if (!(bs = alnSeqOpen (db2Name, cpath, false))) die ("failed to open %s", db2Name) ;
    }

  Overlap *olaps = new (nOverlaps, Overlap) ;
  int i ;
  for (i = 0 ; i < nOverlaps ; ++i)
    { alnReadOverlap (ofIn, olaps+i) ;
      alnSkipTrace (ofIn) ;
    }
  printf ("read %d overlaps\n", (int) nOverlaps) ;
  oneFileClose (ofIn) ;

  if (!db2Name) // add the reverse matches
    { resize (olaps, nOverlaps, 2*nOverlaps, Overlap) ;
      Overlap *o1 = olaps, *o2 = olaps + nOverlaps ;
      for (i = 0 ; i < nOverlaps ; ++i, ++o1, ++o2) flip (o1, o2) ;
      nOverlaps *= 2 ;
      printf ("self-alignment: doubled overlaps to %d\n", (int) nOverlaps) ;
    }
  timeUpdate (stdout) ;

  if (ofa)
    { qsort (olaps, nOverlaps, sizeof(Overlap), overlapOrder) ;
      insertionReport (ofa, as, olaps, nOverlaps) ;
      printf ("wrote %d insertions in %s to %s\n",
	      (int)ofa->info['V']->accum.count, db1Name, ofaName) ;
      oneFileClose (ofa) ;
      timeUpdate (stdout) ;
    }

  if (ofb)
    { Overlap *o1 = olaps ;
      for (i = 0 ; i < nOverlaps ; ++i, ++o1) flip (o1, o1) ;
      qsort (olaps, nOverlaps, sizeof(Overlap), overlapOrder) ;
      insertionReport (ofb, bs, olaps, nOverlaps) ;
      printf ("wrote %d insertions in %s to %s\n",
	      (int)ofb->info['V']->accum.count, db2Name, ofbName) ;
      oneFileClose (ofb) ;
      timeUpdate (stdout) ;
    }

  printf ("Total resources used: ") ; timeTotal (stdout) ;
}

typedef struct {
  int a, a_begin, a_end ;
  int b, b_match_begin, b_match_end ;
} Insertion ;

static int insertionOrder (const void *x, const void *y) // need complete sort because will compress
{
  Insertion *ix = (Insertion*)x, *iy = (Insertion*)y ;
  int c = ix->a - iy->a ; if (c) return c ;
  c = ix->a_begin - iy->a_begin ; if (c) return c ;
  c = ix->a_end - iy->a_end ; return c ;
}

void insertionReport (OneFile *of, AlnSeq *as, Overlap *olap, int n)
// look for insertions in a with respect to b, so olap is sorted on b, then a, then b_begin
// need to accumulate in an array and deduplicate
{
  int i,j ;
  Overlap *oi, *oj ;
  Array a = arrayCreate (4096, Insertion) ;

  for (i = 0, oi = olap ; i < n ; ++i, ++oi)
    for (j = i+1, oj = oi + 1 ; j < n ; ++j, ++oj)
      if (oj->aread != oi->aread || oj->bread != oi->bread) break ;
      else if (COMP(oj->flags) != COMP(oi->flags)) continue ;
      else if (oj->path.bbpos < oi->path.bepos - MAX_OVERHANG) continue ;
      else if (oj->path.bbpos > oi->path.bepos + MAX_OVERHANG) break ;
      else if (COMP(oj->flags) &&
	       oi->path.abpos > oj->path.aepos &&
	       oi->path.abpos < oj->path.aepos + MAX_SIZE)
	{ Insertion *ins = arrayp (a, arrayMax(a), Insertion) ;
	  ins->a = oj->aread ; ins->a_begin = oj->path.aepos ; ins->a_end = oi->path.abpos ;
	  ins->b = oj->bread ; ins->b_match_begin = oi->path.bepos ; ins->b_match_end = oj->path.bbpos ;
	}
      else if (!COMP(oj->flags) &&
	       oj->path.abpos > oi->path.aepos &&
	       oj->path.abpos < oi->path.aepos + MAX_SIZE)
	{ Insertion *ins = arrayp (a, arrayMax(a), Insertion) ;
	  ins->a = oj->aread ; ins->a_begin = oi->path.aepos ; ins->a_end = oj->path.abpos ;
	  ins->b = oj->bread ; ins->b_match_begin = oi->path.bepos ; ins->b_match_end = oj->path.bbpos ;
	}

  arraySort (a, insertionOrder) ;
  arrayCompress (a, insertionOrder) ;

  oneInt(of,0) = MAX_OVERHANG ; oneWriteLine (of, 'o', 0, 0) ;
  oneInt(of,0) = MAX_SIZE ; oneWriteLine (of, 'i', 0, 0) ;
  U64   ias = 0, sLen = 0 ;
  char *s = alnSeqNext (as, &sLen) ; // get 0'th sequence
  char *idBuf = new(256,char) ;
  for (i = 0 ; i < arrayMax (a) ; ++i)
    { Insertion *ins = arrp(a,i,Insertion) ;
      oneInt(of,0) = ins->a ; oneInt(of,1) = ins->a_begin ; oneInt(of,2) = ins->a_end ;
      oneWriteLine (of, 'V', 0, 0) ;
      oneInt(of,0) = ins->b ; oneInt(of,1) = ins->b_match_begin ; oneInt(of,2) = ins->b_match_end ;
      oneWriteLine (of, 'B', 0, 0) ;
      while (ias < ins->a)
	{ s = alnSeqNext (as, &sLen) ;
	  if (!s) die ("run out of contig sequences at %lld < %d", ias, ins->a) ;
	  ++ias ;
	}
      oneWriteLine (of, 'S', ins->a_end - ins->a_begin, s + ins->a_begin) ;
      sprintf (idBuf,"%d:%d-%d_%d:%d-%d",
	       ins->a, ins->a_begin, ins->a_end, ins->b, ins->b_match_begin, ins->b_match_end) ;
      oneWriteLine (of, 'I', strlen(idBuf), idBuf) ;
    }
  free (idBuf) ;
  arrayDestroy (a) ;
}
