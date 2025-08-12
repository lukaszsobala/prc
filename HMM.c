/* 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		   PRC, the profile comparer, version 1.5.6

	HMM.c: mundane struct HMM stuff (reading, allocation, reversal)

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Copyright (C) 2002-5 Martin Madera and MRC LMB, Cambridge, UK
   All Rights Reserved

   This source code is distributed under the terms of the GNU General Public 
   License. See the files COPYING and LICENSE for details.

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "prc.h"


// allocate an HMM
//
// hmm->i is *not* allocated, but set to i
//
// (to create an HMMinfo structure use alloc_HMMinfo() below)
//
HMM* alloc_HMM( HMMinfo *i )
{
  HMM* hmm;

  calloc_1D_array( hmm, HMM, 1);

  calloc_2D_array( hmm->Pmat,   float, i->M+4, 20 );
  calloc_2D_array( hmm->Pins,   float, i->M+4, 20 );
  calloc_2D_array( hmm->Smat,   float, i->M+4, 20 );
  calloc_2D_array( hmm->Sins,   float, i->M+4, 20 );
  calloc_2D_array( hmm->Pt,     float, i->M+4, N_PROF_HMM_TRANS );
  calloc_2D_array( hmm->St_lin, float, i->M+4, N_PROF_HMM_TRANS );
  calloc_2D_array( hmm->St_log, float, i->M+4, N_PROF_HMM_TRANS );

  hmm->i = i;

  return hmm;
}

// allocate HMMinfo
//
// the filename string is allocated and copied, the name is *not*
//
HMMinfo* alloc_HMMinfo( char *filename, char *name, int M )
{
  HMMinfo *info;
  int     len = strlen(filename)+1;

  calloc_1D_array( info, HMMinfo, 1 );
  calloc_1D_array( info->filename, char, len );
  strncpy(info->filename, filename, len );

  info->name = name;
  info->M = M;

  return info;
} 

// free an HMM
//
// hmm->i is *not* freed, and is returned
//
// (to free it, use free_HMMinfo below)
//
HMMinfo* free_HMM( HMM* hmm )
{
  HMMinfo *i = NULL;

  if( hmm != NULL )
    {
      i = hmm->i;

      free_2D_array( hmm->Pmat   );
      free_2D_array( hmm->Pins   );
      free_2D_array( hmm->Smat   );
      free_2D_array( hmm->Sins   );
      free_2D_array( hmm->Pt     );
      free_2D_array( hmm->St_lin );
      free_2D_array( hmm->St_log );
      
      free(hmm);
    };

  return i;
}

// fully free an HMMinfo structure
//
void free_HMMinfo( HMMinfo* info )
{
  if( info != NULL )
    {
      free_unless_null( info->filename );
      free_unless_null( info->name     );
      free_unless_null( info           );
    };
}

// get model name from a file name
//
char* get_HMM_name( char *filename )
{
  char *HMM_name;
  int  len = strlen(filename);
  int  start = 0, end = len;
  int  i;

  for(i=len-1; i>=0; i--)
    {
      if(filename[i]=='/') { start=i+1; break; };
    };

  for(i=len-1; i>=start; i--)
    {
      if(filename[i]=='.') { end = i-1; break; };
    };

  malloc_1D_array( HMM_name, char, end-start+2);
  strncpy( HMM_name, &filename[start], end-start+1 );
  HMM_name[end-start+1] = 0;

  return HMM_name;
}

// set St_xxx[tXX], where tXX should be either tBM or tME, which determine the
// LOCAL/GLOBAL alignment modes
//
// NB the B->M and M->E modes can differ! 
//
void set_HMM_align_mode( HMM* hmm, int tXX, int mode, int start, int end )
{
  float **St_lin = hmm->St_lin, **St_log = hmm->St_log;
  int   m;

  if( (tXX != tBM) && (tXX != tME) )
    die2("Unknown tXX=%d!", tXX);

  // no transitions up to start-1
  for(m=0; m<start; m++)
    {
      St_lin[m][tXX] = 0.0;
      St_log[m][tXX] = LOG_ZERO; 
    };

  // start ... end
  if( mode == LOCAL )
    {
      for(m=start; m<=end; m++)
	{
	  St_lin[m][tXX] = 1.0;
	  St_log[m][tXX] = 0.0;
	};

      if( tXX == tBM )
	{
	  St_lin[end][tXX] = 0.0;
	  St_log[end][tXX] = LOG_ZERO;
	}
      else if( tXX == tME )
	{
	  St_lin[start][tXX] = 0.0;
	  St_log[start][tXX] = LOG_ZERO;
	};
    }
  else if( mode == GLOBAL )
    {
      for(m=start; m<=end; m++)
	{
	  St_lin[m][tXX] = 0.0;
	  St_log[m][tXX] = LOG_ZERO;
	};

      if( tXX == tBM )
	{
	  St_lin[start][tXX] = 1.0;
	  St_log[start][tXX] = 0.0;
	}
      else if( tXX == tME )
	{
	  St_lin[end][tXX] = 1.0;
	  St_log[end][tXX] = 0.0;
	};
    }
  else
    {
      die2("Unknown alignment mode #%d!", mode);
    };
  
  // no transitions after end
  for(m=end+1; m<hmm->i->M+4; m++)
    {
      St_lin[m][tXX] = 0.0;
      St_log[m][tXX] = LOG_ZERO; 
    };

  /*
  printf( "set %s in HMM %s:%d-%d to %s\n",
	  (tXX==tBM)?"tBM":"tME",
	  hmm->name,
	  start,
	  end,
	  (mode==LOCAL)?"LOCAL":"GLOBAL" );
  */
}

// make sure probabilities sum up to one
//
// for PLAN7 this function redistributes any D->I or I->D probabilities among
// the remaining transitions
//
// NB see the comments in prc.h about struct HMM
//
void normalize_HMM( HMM* hmm )
{
#if PROF_HMM_TRANS == PLAN7
  float delta_DI, delta_ID;   // excess P assumed to be due to D->I & I->D
#endif
  float sum_M, sum_I, sum_D;  // sums of probalities; should be 1.0
  float **t = hmm->Pt;        // shortcut to save typing
  int   M = hmm->i->M;        // shortcut to save typing
  int   m, i;

  /* transition probabilities */

  // transitions out of M0, I0, D0, M_(M+2), I_(M+2) and D_(M+2) are handled by
  // the for-loop below; here we need to "manually" sort out the rest, namely:
  //
  // ... transitions out of M-1, I-1 and D-1
  t[0][tDM] = t[0][tMM] = t[0][tIM] = 0.0;
  t[0][tDD] = t[0][tMD] = 0.0;
#if PROF_HMM_TRANS == PLAN9
  t[0][tID] = 0.0;
#endif
  
  // ... transitions out of M_(M+3), D_(M+3) and I_(M+3)
  t[M+3][tMI] = t[M+3][tII] = 0.0;
#if PROF_HMM_TRANS == PLAN9
  t[M+3][tDI] = 0.0;
#endif

  // ... and transitions to D_(M+2)
  //
  // NB this part of the model will also be renormalized below
  //
  t[M+2][tDD] = t[M+2][tMD] = 0.0;
#if PROF_HMM_TRANS == PLAN9
  t[M+2][tID] = 0.0;
#endif

  // transitions out of Mm, Im and Dm
  for(m=0; m<=M+2; m++)
    {
#if PROF_HMM_TRANS == PLAN7
      // redistribute I->D and D->I between the rest 

      // get the excess Ps for Dm->Im and Im->Dm+1
      delta_DI = 1.0 - t[m+1][tDM] - t[m+1][tDD];
      delta_ID = 1.0 - t[m+1][tIM] - t[ m ][tII];

      // redistribute them

      // the "real" delete states are D2 ... D_(M+1)
      if( (2<=m) && (m<=M+1) )
	{
	  t[m+1][tDM] +=  delta_DI * t[m+1][tIM] / (1.0 - t[m][tII]);
	  t[m+1][tDD] +=  delta_DI *  delta_ID   / (1.0 - t[m][tII]);
	};

      // compared to the deletes, there is an extra insert state at the start
      // of the HMM, namely I1
      if( (1<=m) && (m<=M+1) )
	{
	  /*
	    this bit was removed in PRC 1.5.6 -- on reflection I think that the
	    M->I transition is too important to modify it 

	  t[m+1][tMD] += t[m][tMI] *  delta_ID   / (1.0 - t[m][tII]);
	  t[ m ][tMI] -= t[m][tMI] *  delta_ID   / (1.0 - t[m][tII]);
	  */

	  t[m+1][tIM] +=  delta_ID;
	};
#endif 

      /* now proper renormalization */

      // transitions out of Mm
      //
      // M1 = Begin; M2 ... M_(M+1) are the proper matches
      // no transitions out of the End state, M_(M+2)
      //
      if( (1<=m) && (m<=M+1) )
	{
	  sum_M = t[m+1][tMM] + t[m+1][tMD] + t[m][tMI];
	  t[m+1][tMM] /= sum_M;
	  t[m+1][tMD] /= sum_M;
	  t[ m ][tMI] /= sum_M;
	}
      else
	{
	  t[m+1][tMM] = t[m+1][tMD] = t[ m ][tMI] = 0.0;
	};

      // transitions out of Dm
      //
      // D2 ... D_(M+2) are the proper deletes
      //
      if( (2<=m) && (m<=M+1) )
	{
#if   PROF_HMM_TRANS == PLAN9
	  sum_D = t[m+1][tDD] + t[m+1][tDM] + t[m][tDI];
	  t[ m ][tDI] /= sum_D;
#elif PROF_HMM_TRANS == PLAN7
	  sum_D = t[m+1][tDD] + t[m+1][tDM];
#endif
	  t[m+1][tDM] /= sum_D;
	  t[m+1][tDD] /= sum_D;
	}
      else
	{
	  t[m+1][tDM] = t[m+1][tDD] = 0.0;
#if   PROF_HMM_TRANS == PLAN9
	  t[ m ][tDI] = 0.0;
#endif
	};

      // transitions out of Im
      //
      // I1 .. I_(M+1) are the proper inserts
      //
      if( (1<=m) && (m<=M+1) )
	{
#if   PROF_HMM_TRANS == PLAN9
	  sum_I = t[m+1][tID] + t[m+1][tIM] + t[m][tII];
	  t[m+1][tID] /= sum_I;
#elif PROF_HMM_TRANS == PLAN7
	  sum_I = t[m+1][tIM] + t[m][tII];
#endif
	  t[m+1][tIM] /= sum_I;
	  t[ m ][tII] /= sum_I;
	}
      else
	{
	  t[m+1][tIM] = t[m][tII] = 0.0;
#if   PROF_HMM_TRANS == PLAN9
	  t[m+1][tID] = 0.0;
#endif	  
	};
    };


  /* emission probabilities */

  for(m=0; m<=M+3; m++)
    {
      // insertions
      if( (1<=m) && (m<=M+1) ) 
	{
	  sum_I = 0;
	  for(i=0; i<20; i++) sum_I += hmm->Pins[m][i];

	  if(sum_I == 0)
	    {
	      // this happens for the initial & final insertions for HMMER

	      for(i=0; i<20; i++) hmm->Pins[m][i] = 0.0;

	      if( (m!=1) && (m!=M+1) )
		fprintf( stderr, 
			 "Warning: Model %s has zero insert emissions "
			 "in segment %d!\n\n",
			 hmm->i->name, m-1 );
	    }
	  else
	    {
	      for(i=0; i<20; i++) hmm->Pins[m][i] /= sum_I;
	    };
	} 
      else 
	{
	  for(i=0; i<20; i++) hmm->Pins[m][i] = 0.0;
	};

      // matches
      if( (2<=m) && (m<=M+1) ) 
	{
	  sum_M = 0;
	  for(i=0; i<20; i++) sum_M += hmm->Pmat[m][i];

	  if(sum_M == 0)
	    {
	      // this shouldn't happen, but just in case (see insertions above)

	      for(i=0; i<20; i++) hmm->Pmat[m][i] = 0.0;

	      fprintf( stderr, 
		       "Warning: model %s has zero match emissions "
		       "in segment %d!\n\n",
		       hmm->i->name, m-1 );
	    }
	  else
	    {
	      for(i=0; i<20; i++) hmm->Pmat[m][i] /= sum_M;
	    };
	}
      else
	{
	  for(i=0; i<20; i++) hmm->Pmat[m][i] = 0.0;
	};
    };
}

// set up hmm->Smat and hmm->Sins for dot2
//
void scorify_HMM_M_I( HMM *hmm, float *sqrt_null )
{
  int m, i;

  for(m=0; m<=hmm->i->M+3; m++) {
    for(i=0; i<20; i++)
      {
	hmm->Smat[m][i] = hmm->Pmat[m][i] / sqrt_null[i];
	hmm->Sins[m][i] = hmm->Pins[m][i] / sqrt_null[i];
      };
  };
};

// calculate the simple null, and transition scores
//
// match and insert scores are set independently by scorify_HMM_M_I() above
//
void scorify_HMM( HMM *hmm )
{
  int   i, m, M=hmm->i->M;
  float sum = 0.0;

  // PinsJ is just the arithmetic mean of Pmat[2..M+1]
  for(i=0; i<20; i++)
    hmm->PinsJ[i] = 0.0;

  for(m=2; m<=M+1; m++) 
    for(i=0; i<20; i++) 
      hmm->PinsJ[i] += hmm->Pmat[m][i];

  for(i=0; i<20; i++)
    {
      hmm->PinsJ[i] += ((float)M) / 1000.0;
      sum += hmm->PinsJ[i];
    };

  for(i=0; i<20; i++)
    hmm->PinsJ[i] /= sum;

  // set PtJJ so that mean emission length = M+1
  hmm->PtJJ = 1.0 - 1.0 / ((float)M);

  // sanity check: don't want PtJJ to be less than 0.25 (e.g. if M=1)
  if( hmm->PtJJ < 0.25 )
    hmm->PtJJ = 0.25;

  // got the null, now calculate the transition scores
  for(m=0; m<=M+3; m++) 
    {
      hmm->St_lin[m][tMM] = hmm->Pt[m][tMM] / hmm->PtJJ;
      hmm->St_lin[m][tMI] = hmm->Pt[m][tMI] / hmm->PtJJ;
      hmm->St_lin[m][tII] = hmm->Pt[m][tII] / hmm->PtJJ;
      hmm->St_lin[m][tIM] = hmm->Pt[m][tIM] / hmm->PtJJ;
      hmm->St_lin[m][tDM] = hmm->Pt[m][tDM] / hmm->PtJJ;
      hmm->St_lin[m][tDD] = hmm->Pt[m][tDD];
      hmm->St_lin[m][tMD] = hmm->Pt[m][tMD];
#if PROF_HMM_TRANS == PLAN9
      hmm->St_lin[m][tDI] = hmm->Pt[m][tDI] / hmm->PtJJ;
      hmm->St_lin[m][tID] = hmm->Pt[m][tID];
#endif

      for(i=2; i<N_PROF_HMM_TRANS; i++)
	hmm->St_log[m][i] = LOG( hmm->St_lin[m][i] );
    };
}

// Calculation of the reverse HMM
// ==============================
//
// This is the only vaguely interesting function in this file. The problem is
// as follows: for the purpose of the profile-profile reverse null it is
// necessary to calculate the "reverse HMM", defined as that for which
//
//     P(reverse seq|reverse HMM) = P(seq|HMM)     for all seqs
//
// The mapping of the match/insert states is obvious: we just reverse the
// order, like this:
//
//     reverse M_m = M_(M+3-m)
//     reverse D_m = D_(M+3-m)
//     reverse I_m = I_(M+2-m)
//
// The transition probabilities are more complicated; the correct equation
// appears to be:
//
//     reverse P(N->M) = P(M->N) * P(M) / P(N)
//
// where M,N are two adjacent nodes (in the original or "forward" HMM, M is
// before N), P(M->N) is the original "forward" transition probability, and
// P(M) and P(N) are the probabilities that a forward path goes through M and
// N, respectively.
//
// The formula expresses the fact that the probability of the reverse
// transition is just the probability that (in the forward model) we got to N
// from *M*, given that we got to N, i.e. this is basically Bayes' Theorem.
//
// Insert states are a pain, because they involve cyclical paths, but luckily
// there's a simple way around this. Firstly, one notes that P(I->I) must
// remain the same, because otherwise insertion lengths would differ. Then, one
// pretends that the I->I transition doesn't actually exist and normalizes it
// out:
//        P(I->D) /= 1 - P(I->I)
//        P(I->M) /= 1 - P(I->I)
//
// Now there are no cyclical paths, so the HMM can be reversed easily. Finally,
// one re-introduces the I->I transition and sets:
//
//       reverse P(I->D) *= 1 - P(I->I)
//       reverse P(I->M) *= 1 - P(I->I) .
//
// It can be shown that the reverse HMM constructed in this way does satisfy
// the obvious requirement that
//
//     reverse reverse HMM = HMM .
//
// Finally, we're assuming that the HMM we're passed is a nice one -- i.e. that
// it has been through normalize_HMM() above.
//
// See struct HMM in prc.h for PRC conventions.
//
HMM* reverse_HMM( HMM* hmm )
{
  HMM*  rev = alloc_HMM(hmm->i);
  float Pm_curr = 1.0, Pd_curr = 0.0, Pi_curr;
  float Pm_next, Pd_next, sum;
  float **t = hmm->Pt, f;
  int   m, i, M = hmm->i->M;

  // the emissions: simple
  for(m=2; m<=M+1; m++) 
    for(i=0; i<20; i++)
      rev->Pmat[M+3-m][i] = hmm->Pmat[m][i];

  for(m=1; m<=M+1; m++) 
    for(i=0; i<20; i++)
      rev->Pins[M+2-m][i] = hmm->Pins[m][i]; 

  // now the transitions
  for(m=1; m<=M+1; m++)
    {
      // NB _curr = _m, _next = _(m+1), in the forward HMM

      // f is the normalization factor for transitions out of I_m
      f = 1.0 - t[m][tII];

      // get Pi_curr
#if   PROF_HMM_TRANS == PLAN7
      Pi_curr = Pm_curr * t[m][tMI];
#elif PROF_HMM_TRANS == PLAN9
      Pi_curr = Pm_curr * t[m][tMI] + Pd_curr * t[m][tDI];
#endif

      // get Pm_next & Pd_next
      Pm_next = Pm_curr * t[m+1][tMM] + Pd_curr * t[m+1][tDM]
	+ Pi_curr * (t[m+1][tIM] / f);
#if   PROF_HMM_TRANS == PLAN7
      Pd_next = Pd_curr * t[m+1][tDD] + Pm_curr * t[m+1][tMD];
#elif PROF_HMM_TRANS == PLAN9
      Pd_next = Pd_curr * t[m+1][tDD] + Pm_curr * t[m+1][tMD]
	+ Pi_curr * (t[m+1][tID] / f);
#endif
      sum = Pm_next + Pd_next;
      Pm_next /= sum;
      Pd_next /= sum;

      // now get the reverse Ps
      if( Pm_next == 0.0 )
	{
	  // this shouldn't ever happen, but ...
	  rev->Pt[M+3-m][tMM] = 1.0;
	  rev->Pt[M+3-m][tMD] = 0.0;
	  rev->Pt[M+2-m][tMI] = 0.0;
	}
      else
	{
	  rev->Pt[M+3-m][tMM] = Pm_curr * t[m+1][tMM] / Pm_next;
	  rev->Pt[M+3-m][tMD] = Pd_curr * t[m+1][tDM] / Pm_next;
	  rev->Pt[M+2-m][tMI] = Pi_curr * (t[m+1][tIM] / f) / Pm_next; 
	};
      
      if( Pi_curr == 0.0 )
	{
	  // this does actually happen, e.g. for read_HMM_fasta() output
	  rev->Pt[M+3-m][tIM] = 1.0;
	  rev->Pt[M+2-m][tII] = 0.0;
#if PROF_HMM_TRANS == PLAN9
	  rev->Pt[M+3-m][tID] = 0.0;
#endif
	}
      else
	{
	  rev->Pt[M+2-m][tII] = t[m][tII];
	  rev->Pt[M+3-m][tIM] = f * Pm_curr * t[m][tMI] / Pi_curr;
#if PROF_HMM_TRANS == PLAN9
	  rev->Pt[M+3-m][tID] = f * Pd_curr * t[m][tDI] / Pi_curr;
#endif
	};
      
      if( Pd_next == 0.0 )
	{
	  // this does actually happen, e.g. for read_HMM_fasta() output
	  rev->Pt[M+3-m][tDM] = 1.0;
	  rev->Pt[M+3-m][tDD] = 0.0;
#if PROF_HMM_TRANS == PLAN9
	  rev->Pt[M+2-m][tDI] = 0.0;
#endif
	}
      else
	{
	  rev->Pt[M+3-m][tDM] = Pm_curr * t[m+1][tMD] / Pd_next;
	  rev->Pt[M+3-m][tDD] = Pd_curr * t[m+1][tDD] / Pd_next;
#if PROF_HMM_TRANS == PLAN9
	  rev->Pt[M+2-m][tDI] = Pi_curr * (t[m+1][tID] / f) / Pd_next; 
#endif
	};

      // prepare Px_curr for the next round
      Pm_curr = Pm_next;
      Pd_curr = Pd_next;
    };
  
  normalize_HMM(rev);
  scorify_HMM(rev);

  //print_HMM(hmm);
  //print_HMM(rev);

  return rev;
}

// print out the HMM in a vaguely SAM-ish, nice human-readable format
//
void print_HMM( FILE *stream, HMM* hmm )
{
  int  m, i, j;
  char order[]="ACDEFGHIKLMNPQRSTVWY";

  fprintf(stream, "Model name :  %s\n", hmm->i->name);
  fprintf(stream, "Length     :  %d\n\n", hmm->i->M);

  for(m=0; m<=hmm->i->M+3; m++)
    {
      fprintf( stream, "NODE %d\n", m );

      fprintf( stream, "Transitions:\n" );

      fprintf( stream, 
	       "D%d->D%d = %.3f\tM%d->D%d = %.3f\tI%d->D%d = %.3f\n", 
	       m-1, m, hmm->Pt[m][tDD], 
	       m-1, m, hmm->Pt[m][tMD],
	       m-1,m,
#if   PROF_HMM_TRANS == PLAN7
	       0.0
#elif PROF_HMM_TRANS == PLAN9
	       hmm->Pt[m][tID]
#endif
	     );
	  
      fprintf( stream, 
	       "D%d->M%d = %.3f\tM%d->M%d = %.3f\tI%d->M%d = %.3f\n", 
	       m-1, m, hmm->Pt[m][tDM], 
	       m-1, m, hmm->Pt[m][tMM], 
	       m-1, m, hmm->Pt[m][tIM] );
      
      fprintf( stream, 
	       "D%d->I%d = %.3f\tM%d->I%d = %.3f\tI%d->I%d = %.3f\n",
	       m, m,
#if   PROF_HMM_TRANS == PLAN7
	       0.0,
#elif PROF_HMM_TRANS == PLAN9
	       hmm->Pt[m][tDI],
#endif
	       m, m, hmm->Pt[m][tMI], 
	       m, m, hmm->Pt[m][tII] );
      

      fprintf( stream, "S(B->M%d) = %.3f\tS(M%d->E) = %.3f\n",
	       m, hmm->St_lin[m][tBM], m, hmm->St_lin[m][tME] );

      fprintf( stream, "\n" );


      fprintf( stream, "Match emissions:\n" );
      for(i=0; i<20; i+=5) {
	for(j=0; j<5; j++)
	  fprintf( stream, "%c = %.3f  ", order[i+j], hmm->Pmat[m][i+j] );
	fprintf( stream, "\n" );
      };
      fprintf( stream, "\n" );

      fprintf( stream, "Insert emissions:\n" );
      for(i=0; i<20; i+=5) {
	for(j=0; j<5; j++)
	  fprintf( stream, "%c = %.3f  ", order[i+j], hmm->Pins[m][i+j] );
	fprintf( stream, "\n" );
      };
      fprintf( stream, "\n\n" );
    };
}

// write the HMM to a binary file
//
// the file should be fopen()d before calling the function, and fclose()d
// afterwards 
//
void write_HMM_PRC_binary(FILE *file, HMM *hmm)
{
  char  version[7] = "bin005";
  int   len = strlen(hmm->i->name) + 1, n_trans = N_PROF_HMM_TRANS;
  int   M = hmm->i->M;

  fwrite( version,      sizeof(char),  7,             file );
  fwrite( &n_trans,     sizeof(int),   1,             file );
  fwrite( &M,           sizeof(int),   1,             file );
  fwrite( &len,         sizeof(int),   1,             file );
  fwrite( hmm->i->name, sizeof(char),  len,           file );
  fwrite( hmm->PinsJ,   sizeof(float), 20,            file );
  fwrite( &hmm->PtJJ,   sizeof(float), 1,             file );
  fwrite( hmm->Pmat[0], sizeof(float), (M+4)*20,      file );
  fwrite( hmm->Pins[0], sizeof(float), (M+4)*20,      file );
  fwrite( hmm->Pt[0],   sizeof(float), (M+4)*n_trans, file );
}

// read the HMM from a binary PRC file written by the routine above
//
// supports bin004 and bin005
//
HMM* read_HMM_PRC_binary(char *filename)
{
  FILE*    file;
  char     version[7];
  HMMinfo* info;
  HMM*     hmm;
  int      len, M, n_trans, i;
  float    dummy;

  open_file_or_die(file, filename, "r");

  fread(version, sizeof(char), 7, file);

  if(strncmp(version, "bin004", 6) == 0) 
    {
      fread(&n_trans, sizeof(int), 1, file);

      // silly bug in old versions of PRC -- N_PROF_HMM_TRANS was one bigger
      // than necessary
      //
      // this is the only difference between bin004 and bin005
      //
      if(n_trans != N_PROF_HMM_TRANS+1)
	{
	  fclose(file);
	  die2( "Sorry, the file '%s' was generated by a version of PRC\n"
		"compiled with a different PROF_HMM_TRANS setting!\n",
		filename );
	};
      
      fread(&M,   sizeof(int), 1, file);
      fread(&len, sizeof(int), 1, file);
  
      info = alloc_HMMinfo(filename, NULL, M);
      hmm  = alloc_HMM(info);
      malloc_1D_array(info->name, char, len);
      
      fread( info->name,   sizeof(char),  len,      file );
      fread( hmm->PinsJ,   sizeof(float), 20,       file );
      fread( &hmm->PtJJ,   sizeof(float), 1,        file );
      fread( hmm->Pmat[0], sizeof(float), (M+4)*20, file );
      fread( hmm->Pins[0], sizeof(float), (M+4)*20, file );

      for(i=0; i<M+4; i++)
	{
	  fread( hmm->Pt[i], sizeof(float), n_trans-1, file );
	  fread( &dummy,     sizeof(float), 1,         file );
	};
    }
  else if(strncmp(version, "bin005", 6) == 0)
    {
      fread(&n_trans, sizeof(int), 1, file);

      if(n_trans != N_PROF_HMM_TRANS)
	{
	  fclose(file);
	  die2( "Sorry, the file '%s' was generated by a version of PRC\n"
		"compiled with a different PROF_HMM_TRANS setting!\n",
		filename );
	};
      
      fread(&M,   sizeof(int), 1, file);
      fread(&len, sizeof(int), 1, file);
  
      info = alloc_HMMinfo(filename, NULL, M);
      hmm  = alloc_HMM(info);
      malloc_1D_array(info->name, char, len);
      
      fread( info->name,   sizeof(char),  len,           file );
      fread( hmm->PinsJ,   sizeof(float), 20,            file );
      fread( &hmm->PtJJ,   sizeof(float), 1,             file );
      fread( hmm->Pmat[0], sizeof(float), (M+4)*20,      file );
      fread( hmm->Pins[0], sizeof(float), (M+4)*20,      file );
      fread( hmm->Pt[0],   sizeof(float), (M+4)*n_trans, file );
    }
  else
    {      
      fclose(file);
      die2( "Sorry, the file '%s' is not in the bin004 or bin005 formats!",
	    filename );
    };

  fclose(file);
  scorify_HMM(hmm);

  return hmm;
}


// the SAM order of transition Ps 
static int SAM_map[] = {
#if   PROF_HMM_TRANS == PLAN9
  tDD, tMD, tID,   
  tDM, tMM, tIM,
  tDI, tMI, tII 
#elif PROF_HMM_TRANS == PLAN7
  tDD, tMD, -1,   
  tDM, tMM, tIM,
  -1,  tMI, tII 
#endif
};

// read 4 bytes (= sizeof(float), sizeof(int)) from a file; this
// function is needed to read the binary SAM format, which does this
// to make the model files architecture-independent
//
// blatantly based on Richard Hughey's SAM code
//
// returns 1 if success, 0 otherwise
//
int SAM_binary_read_4B (FILE* file, void* ptr)
{
  unsigned char b[4];
  unsigned int  y = 0x0;
  int i, dummy;

  for ( i=0; i < 4; i++)
    {
      if( (dummy = getc(file)) == EOF )
	return 0;
      b[i] = (unsigned char) dummy;
    }

  for (i = 3; i >=1; i-- ) 
    {
      y |= b[i];
      y <<=8;
    }
  y |=b[i];
 
  *(unsigned int  *)ptr = y;

  return 1;
}

// converting SAM logs to probabilities
// blatantly based on Richard Hughey's SAM code
//
#define SAM_BASE -0.001
#define SAM_NLOG0 80000
#define SAM_ZERO 2.0612e-9
static inline float SAM_NEXP(int log)
{
  return ( log>=SAM_NLOG0 ? SAM_ZERO : exp(SAM_BASE*(double)log) );
}

#undef DEBUG_SAM_BINARY

// need float data[49]
// return 1 if success, 0 otherwise
//
int SAM_binary_read_node(FILE* file, float* data, int n)
{
  int i, dummy;

#ifdef DEBUG_SAM_BINARY
  printf("read_node 1\n");
#endif

  if( (!SAM_binary_read_4B(file, &dummy)) || (dummy==0) )
    return 0;

#ifdef DEBUG_SAM_BINARY
  printf("read_node 2\n");
#endif

  // test in binary_r_NodeEnv()
  SAM_binary_read_4B(file, &dummy);
  if( dummy != 0 )
    {
      SAM_binary_read_4B(file, &dummy);
      SAM_binary_read_4B(file, &dummy);
      SAM_binary_read_4B(file, &dummy);
    }

  // node type
  getc(file);

  for(i=0; i<n; i++)
    {
      if( ! SAM_binary_read_4B(file, &data[i]) )
	return 0;

#ifdef DEBUG_SAM_BINARY
      printf("%2i: %08X %f\n", i, *(unsigned int*)&data[i], data[i]);
#endif
    }

  return 1;
}

// read in a profile HMM in the SAM 3.x binary format
//
// many details were worked out by looking at Richard Hughey's SAM code
//
HMM* read_HMM_SAM3x_binary( char* filename ) 
{
  FILE    *file;
  HMM     *hmm;
  HMMinfo *info;
  char    line[200] = "";
  float   data[49], itable[20], ftable[20];
  int     i, dummy, iindex, findex, m=0;

  open_file_or_die(file, filename, "r");

  // get to the start of the file
  while( toupper(line[0]) != 'A' )
    fgets_null_check(fgets(line, 200, file), filename);

  // ntot, npos and nneg, whatever that means
  // we're using npos
  SAM_binary_read_4B(file, &dummy);
  if( ! SAM_binary_read_4B(file, &m) ) return NULL;
  SAM_binary_read_4B(file, &dummy);

  // findex, iindex from binary_read_tables()
  SAM_binary_read_4B(file, &findex);
  if( findex != -1 )
    for(i=0; i<20; i++)
      SAM_binary_read_4B(file, &ftable[i]);

  SAM_binary_read_4B(file, &iindex);
  if( iindex != -1 )
    for(i=0; i<20; i++)
      SAM_binary_read_4B(file, &itable[i]);

  // generic, lettcount, freqave
  // ignore the first two; freqave is needed sometimes
  SAM_binary_read_node(file, data, 49);
  SAM_binary_read_node(file, data, 49);
  SAM_binary_read_node(file, data, 49);

  // m is the number of "positive nodes", whatever that means!
  info = alloc_HMMinfo( filename, get_HMM_name(filename), m );
  hmm  = alloc_HMM( info );

  for( m=1; m<=hmm->i->M+2; m++ )
    {
      if( (iindex !=-1) || (findex !=-1) )
	{
	  if( ! SAM_binary_read_node(file, data, 29) )
	    die3("Error parsing node %d of '%s'", m-2, filename); 
	}
      else
	{
	  if( ! SAM_binary_read_node(file, data, 49) )
	    die3("Error parsing node %d of '%s'", m-2, filename); 
	}

      for(i=0; i<9; i++ )
        {
          if(SAM_map[i]==-1) continue;
          hmm->Pt[m][SAM_map[i]] = data[i];
        };

      for(i=0; i<20; i++ )
	{
	  hmm->Pmat[m][i] = data[i+9];

	  if ( findex != -1 )      hmm->Pins[m][i] = ftable[i];
	  else if ( iindex != -1 ) hmm->Pins[m][i] = itable[i];
	  else                     hmm->Pins[m][i] = data[i+29];
        };
    };

  fclose( file );  

  normalize_HMM(hmm);
  scorify_HMM(hmm);  

  return hmm;
}

// read in a profile HMM in the SAM 3.x ASCII format
//
HMM* read_HMM_SAM3x_ascii( char* filename ) 
{
  FILE    *file;
  HMM     *hmm;
  HMMinfo *info;
  float   **data;
  char    line[200] = "";
  int     i, m=0;

  open_file_or_die(file, filename, "r");

  // get to the start of the file
  while( strncmp(line, "Begin", 5) )
    {
      fgets_null_check(fgets(line, 200, file), filename);
      
      if( strncmp(line, "BINARY", 6) == 0)
	{
	  fclose(file);
	  return read_HMM_SAM3x_binary(filename);
	}
    }

  // bug spotted by Alejandro Ochoa -- alloc data only once we know it's ASCII
  malloc_2D_array(data, float, MAX_PROF_HMM_LENGTH, 49);

  // raw read into **data -- don't know M
  do 
    {
      // beware, FIMs = extra line, eg TYPE 174 F
      // (should handle this better but this will do for now)
      if(strncmp(line, "TYPE", 4)==0)
	fgets_null_check(fgets(line, 200, file), filename);

      for(i=0; i<49; i++)
	if( fscanf(file, "%f", &data[m][i]) != 1 )
	    die4( "Error parsing record %d, node %d in file '%s'!",
		  i, m, filename );

      // \n after the last float
      fgets_null_check(fgets(line, 200, file), filename); 

      m++;

      if( m == MAX_PROF_HMM_LENGTH )
	die3( "The model in file '%s' is too long!\n"
	      "Please reset MAX_PROF_HMM_LENGTH in prc.h (currently %d) "
	      "and recompile ...", filename, MAX_PROF_HMM_LENGTH );
     
      // node number
      if( fscanf(file, "%s", line) != 1 )
	die3( "Error parsing node number, node %d in file '%s'!",
	      m, filename );
    }
  while( strncmp(line, "ENDMODEL", 8) );
  fclose( file );  

  // now store the data properly
  info = alloc_HMMinfo( filename, get_HMM_name(filename), m-2 );
  hmm  = alloc_HMM( info );
    
  for( m=1; m<=hmm->i->M+2; m++ )
    {
      for(i=0; i<9; i++ )
	{
	  if(SAM_map[i]==-1) continue;
	  hmm->Pt[m][SAM_map[i]] = data[m-1][i];
	};
      
      for(i=0; i<20; i++ )
	{
	  hmm->Pmat[m][i] = data[m-1][i+9];
	  hmm->Pins[m][i] = data[m-1][i+29];
	};
    };
  
  free_2D_array(data);

  normalize_HMM(hmm);
  scorify_HMM(hmm);  

  return hmm;
}

// score -> probability converter for HMMER parser
//
// returns -1.0 when there are problems with the conversion
//
float score2prob( char* score, float null_P )
{
  int int_score;

  if(!score)
    return -1.0;
  else if(score[0]=='*')
    return 0.0;
  else
    {
      if(sscanf(score, "%d", &int_score) != 1)
	return -1.0;

      return null_P * exp(log(2.0)*((float)int_score)/1000.0);
    };
}

// read in a profile HMM in the HMMER 2.x ASCII format
//
HMM* read_HMM_HMMER2x_ascii( char* filename )
{
#define LINELEN 512
  FILE     *file;
  HMM*     hmm;
  HMMinfo* info;
  char     line[LINELEN]="", *name=NULL, *s;
  int      M=-1, got_nule=0, len, m, i;
  float    nule[20], tBD;
  
  // the order of transition Ps in HMMER files
  static int map[]  = { tMM, tMI, tMD, tIM, tII, tDM, tDD };

  // HMMER stores P(M1->M2) in [1], SAM & PRC in [2]
  static int diff[] = {   1,   0,   1,   1,   0,   1,   1 };

  open_file_or_die(file, filename, "r");

  // the HMMER2.0 [2.2g] line 
  fgets_null_check(fgets(line, LINELEN, file), filename);
  if( strncmp(line, "HMMER2.0", 8) != 0 )
    die2( "The file '%s' is not in the HMMER2.0 format!", filename );

  // parse the header
  while( fgets(line, LINELEN, file) )
    {
      if(  (strncmp(line, "ACC  ", 5)==0) ||
	  ((strncmp(line, "NAME ", 5)==0) && (name==NULL)) )
	{
	  free_unless_null(name);
	  len = strcspn(&line[6], " \n\t\r");

	  // sanity check
	  if((len==0) || (len>500)) {
	    if( strncmp(line, "ACC  ", 5)==0 )
	      {
		die2( "Strange ACC line in file '%s', please correct it!", 
		      filename );
	      }
	    else
	      {
		die2( "Strange NAME line in file '%s', please correct it!", 
		      filename );
	      };
	  };

	  malloc_1D_array(name, char, len+1);
	  strncpy(name, &line[6], len);
	  name[len] = 0;
	}
      else if( strncmp(line, "LENG ", 5)==0 )
	{
	  sscanf(&line[6], "%d", &M);
	  if(M<0)
	    die2( "LENG record < 0 in file '%s', please correct it!",
		  filename );
	}
      else if( strncmp(line, "NULE ", 5)==0 )
	{
	  got_nule = 1;
	  for(i=0; i<20; i++)
	    {
	      nule[i] = score2prob(&line[5+7*i], 1.0/20.0);
	      if( nule[i] == -1 )
		die3( "Error parsing NULE score number %d in file '%s'!",
		      i+1, filename );
	    };
	}
      else if( strncmp(line, "HMM  ", 5)==0 )
	{
	  break;
	}
    };

  if( M<0 )
    die2( "No LENG record in file '%s', please add it!", filename );

  if( !got_nule )
    die2( "No NULE record in file '%s', please add it!", filename );

  if( name==NULL )
    name = get_HMM_name(filename);

  info = alloc_HMMinfo(filename, name, M);
  hmm  = alloc_HMM(info);


  // the m->m m->i etc. line
  fgets_null_check(fgets(line, LINELEN, file), filename);

  // the B->D1 line
  fgets_null_check(fgets(line, LINELEN, file), filename);

  strtok(line, " \t\n");      // 1 - tBD
  strtok(NULL, " \t\n");      // *
  s = strtok(NULL, " \t\n");  // tBD
  tBD = score2prob(s, 1.0);

  if( tBD == -1 )
    die2( "Error parsing B->D score in file %s!", filename );

  // the main loop
  for(m=2; m<=M+1; m++)
    {
      // match emissions
      fgets_null_check(fgets(line, LINELEN, file), filename);

      // throw away the first token (line number)
      strtok(line, " \t\n");

      for(i=0; i<20; i++)
	{
	  s = strtok(NULL, " \t\n");
	  hmm->Pmat[m][i] = score2prob(s, nule[i]);
	  if( hmm->Pmat[m][i] == -1 )
	    die4( "Error parsing match emmision score %d, segment %i "
		  "in file '%s'!", i+1, m-1, filename );
	};
		  
      // insert emissions
      fgets_null_check(fgets(line, LINELEN, file), filename);

      // throw away the first token ("-")
      strtok(line, " \t\n");

      for(i=0; i<20; i++)
	{
	  s = strtok(NULL, " \t\n");
	  hmm->Pins[m][i] = score2prob(s, nule[i]);
	  if( hmm->Pins[m][i] == -1 )
	    die4( "Error parsing insert emmision score %d, segment %i "
		  "in file '%s'!", i+1, m-1, filename );
	};

      // transitions
      fgets_null_check(fgets(line, LINELEN, file), filename);

      // throw away the first token ("-")
      strtok(line, " \t\n");

#if PROF_HMM_TRANS == PLAN9
      hmm->Pt[m][tDI] = 0.0;
      hmm->Pt[m][tID] = 0.0;
#endif

      for(i=0; i<7; i++)
	{
	  s = strtok(NULL, " \t\n");
	  hmm->Pt[m+diff[i]][map[i]] = score2prob(s, 1.0);
	  
	  if( hmm->Pt[m+diff[i]][map[i]] == -1 )
	    die4( "Error parsing transition score %d, segment %i "
		  "in file '%s'!", i+1, m-1, filename );
	};
    };
  
  // Pins[1] needs setting to 0 (Pins[M+1]=0 is in the HMMER file)
  for(i=0; i<20; i++)
    hmm->Pins[1][i] = 0.0;

  // route traffic around I1 and I_(M+1)
  hmm->Pt[1][tMI] = hmm->Pt[M+1][tMI] = 0.0;
  hmm->Pt[1][tII] = hmm->Pt[M+1][tII] = 0.0;
  hmm->Pt[2][tIM] = hmm->Pt[M+2][tIM] = 1.0;
#if PROF_HMM_TRANS == PLAN9
  hmm->Pt[1][tDI] = hmm->Pt[M+1][tDI] = 0.0;
  hmm->Pt[2][tID] = hmm->Pt[M+2][tID] = 0.0;
#endif

  // the end is all set to 0.0 by HMMER, but we need M->E and D->E to be 1.0
  hmm->Pt[M+2][tDM] = 1.0;
  hmm->Pt[M+2][tMM] = 1.0;

  // set P(B1->D2) and P(B1->M2)
  hmm->Pt[2][tMM] = 1.0 - tBD;
  hmm->Pt[2][tMD] = tBD;

      
  // finally, check that the '//' line is there
  fgets_null_check(fgets(line, LINELEN, file), filename);
  if( strncmp(line, "//", 2) !=0 )
    die2("Error parsing file '%s': '//' line seems missing!", filename);

  fclose(file);  

  normalize_HMM(hmm);
  scorify_HMM(hmm);  

  return hmm;
}

// read in a PSIBLAST checkpoint file
//
HMM* read_PSIBLAST_binary( char* filename ) 
{
  FILE    *file;
  HMM     *hmm;
  HMMinfo *info;

  // the order in which PSI-BLAST saves its 20 emission probabilities
  // appears to be ARNDC QEGHI LKMFP STWYV (Ala Arg Asn Asp Cys ...)
  //               
  // our order is: ACDEF GHIKL MNPQR STVWY
  //               01234 56789 01234 56789
  //
  int     map[] = { 0,  14, 11, 2,  1,
		    13, 3,  5,  6,  7,
		    9,  8,  10, 4,  12,
		    15, 16, 18, 19, 17  };
  
  double  Pmat[20];
  float   Pins[20], Pt[N_PROF_HMM_TRANS], sum=0.0;

  int     M, m, i;


  // these are all wild guesses, take them with a pinch of salt
  // (I should probably read more about PSI-BLAST ...)
  //
  Pt[tBM] = 0.000;
  Pt[tME] = 0.000;
  Pt[tMM] = 0.965;
  Pt[tMI] = 0.015;
  Pt[tMD] = 0.020;
  Pt[tIM] = 0.350;
  Pt[tII] = 0.650;
  Pt[tDM] = 0.150;
  Pt[tDD] = 0.850;
#if      PROF_HMM_TRANS==PLAN9
  Pt[tID] = 0.000;
  Pt[tDI] = 0.000;
#endif


  open_file_or_die(file, filename, "r");

  // model length
  fread(&M, sizeof(int), 1, file);

  // skip the seed sequence
  fseek(file, M, SEEK_CUR);

  // allocate the HMM
  info = alloc_HMMinfo(filename, get_HMM_name(filename), M);
  hmm  = alloc_HMM(info);

  // read in the emission probabilities
  for(m=1; m<=M; m++)
    {
      fread(Pmat, sizeof(double), 20, file);

      for(i=0; i<20; i++)
	hmm->Pmat[m+1][map[i]] = (float) Pmat[i];
    };

  // calculate inserts: arithmetic average of match state emissions
  // (another wild guess)
  for(i=0; i<20; i++) 
    {
      Pins[i] = 0.0;
      for(m=1; m<=M; m++)
	Pins[i] += hmm->Pmat[m+1][i];

      Pins[i] /= (float)M;
      sum += Pins[i];
    };

  if( sum==0.0 )
    {
      // assuming none of the match probabilities are negative, the only way
      // to get 0 here is if ALL the match probabilities are zero
      // 
      // this happened once in an example Julian Gough sent me (probably a
      // strange PSI-BLAST bug) ...

      fprintf( stderr, 
	       "Warning: All probabilities in '%s' are zero!\n"
	       "(I can handle it, but it's dodgy -- double-check the file.)\n"
	       "\n",
	       filename );

      for(i=0; i<20; i++)
	Pins[i] = 0.05;

      // better set up that matches as well
      for(m=1; m<=M; m++)
	hmm->Pmat[m+1][i] = 0.05;
    }
  else
    {
      for(i=0; i<20; i++)
	Pins[i] /= sum;
    };

  // sort out the rest of the model
  for(m=1; m<=M+2; m++)
    {
      memcpy(hmm->Pt[m],   Pt,   sizeof(float)*N_PROF_HMM_TRANS);
      memcpy(hmm->Pins[m], Pins, sizeof(float)*20);
    };

  normalize_HMM(hmm);
  scorify_HMM(hmm);  
  //print_HMM(hmm);

  return hmm;
}


// construct a simple profile HMM from a single sequence read from a FASTA file
//
// this is useful if you want to compare the performance of PRC against that of
// standard HMM programs like SAM or HMMER
//
// WARNING: QUICK HACK. I should do a Sean and write a dynamically allocated
// fgets, or just grab his sre_fgets ...
//
HMM* read_HMM_FASTA( char* filename ) 
{
  FILE    *file;
  HMM     *hmm;
  HMMinfo *info;
  char    *name, *line, *string, aa[] = "ACDEFGHIKLMNPQRSTVWY";
  int     n=0, m, i, map[26];
  float   null[20], sum = 0.0;

#define NAMELEN  500
#define SEQLEN   10000

  calloc_1D_array(name,   char, NAMELEN + 1);
  calloc_1D_array(line,   char, SEQLEN  + 1);
  calloc_1D_array(string, char, SEQLEN  + 1);

  open_file_or_die(file, filename, "r");

  // get the name from the header
  while( fgets(line, SEQLEN, file ) )
    {
      if( line[0] == '>' )
	{
	  line[SEQLEN-1] = 0;
	  strtok(&line[1], " ;,\t\n");
	  strncpy(name, &line[1], NAMELEN);
	  name[NAMELEN-1] = 0;
	  break;
	};
    };

  // now parse the sequence
  while( fgets(line, SEQLEN, file ) )
    {
      if( line[0] == '>' )
	break;

      strtok(line, " \t\n");
      strncpy(&string[n], line, SEQLEN-n);
      n+=strlen(line);

      if( n > SEQLEN - 100 )
	{
	  die4( "Error parsing file '%s':\n"
		"The sequence %s is too long! The limit is %d amino-acids.",
		filename, name, SEQLEN - 100 );
	};
    };

  fclose(file);


  // create an amino acid -> index map
  for(i=0; i<26; i++)
    map[i]=-1;

  for(i=0; i<20; i++)
    map[aa[i]-'A']=i;

  // check the sequence & calculate the null model
  for(i=0; i<20; i++)
    null[i] = 0.0;

  for(m=0; m<n; m++)
    {
      i = toupper((int)string[m]) - ((int)'A');

      if( (i<0) || (i>= 26) )
	{
	  die5( "Error parsing file '%s':\n"
		"Strange character '%c' at position %d in sequence %s!",
		filename, string[m], m+1, name );
	};

      i = map[i];

      if( i == -1 )
	fprintf( stderr,
		 "Warning parsing file '%s':\n"
		 "Strange amino-acid '%c' at position %d in sequence %s!\n"
		 "Setting emissions to null probabilities ...\n"
		 "\n",
		 filename, string[m], m+1, name );
      else
	  null[i] += 1.0;
    };

  for(i=0; i<20; i++)
    sum += null[i];

  for(i=0; i<20; i++)
    null[i] /= sum;

  
  // construct the HMM
  info = alloc_HMMinfo(filename, NULL, n);
  hmm  = alloc_HMM(info);
  malloc_1D_array(info->name, char, strlen(name)+1);
  strcpy(info->name, name);

  // match and insert emissions
  for(m=0; m<n; m++)
    {
      i = map[ toupper((int)string[m]) - ((int)'A') ];

      memset( hmm->Pmat[m+2], 0, sizeof(float)/sizeof(int)*20 );
      memset( hmm->Pins[m+2], 0, sizeof(float)/sizeof(int)*20 );

      if( i == -1 )
	{
	  // non-standard amino-acid, X, or something
	  // just set it to the null 
	  for(i=0; i<20; i++)
	    {
	      hmm->Pmat[m+2][i] = null[i];
	      hmm->Pins[m+2][i] = null[i];
	    };
	}
      else
	{
	  // wooof, a normal amino-acid
	  //
	  // the inserts aren't necessary, because the DP routine never gets
	  // there, but that's how the null is calculated in scorify_HMM()
	  // above 
	  //
	  hmm->Pmat[m+2][i] = 1.0;
	  hmm->Pins[m+2][i] = 1.0;
	};
    };

  // transitions
  for(m=2; m<n+3; m++)
    {
      memset( hmm->Pt[m], 0, sizeof(float)/sizeof(int)*N_PROF_HMM_TRANS  );

      hmm->Pt[m][tMM] = 1.0;
      hmm->Pt[m][tIM] = 1.0;
      hmm->Pt[m][tDM] = 1.0;
    };

  normalize_HMM(hmm);
  scorify_HMM(hmm);  

  free_unless_null(name);
  free_unless_null(line);
  free_unless_null(string);


  //print_HMM(hmm);

  return hmm;
}

// Function read_HMM: read in an HMM stored in a file
//
// Depending on the extension, the file is assumed to be of the following 
// format:
//
//   .prc     :  PRC binary
//   .mod     :  SAM 3.x (ASCII or binary)
//   .hmm     :  HMMER 2.x ASCII
//   .profile :  PSI-BLAST (alternatives: .prof .pbla .psi .psiblast)
//   .fa      :  a simple FASTA sequence (also .seq)
//
HMM* read_HMM(char *filename)
{
  HMM  *hmm;
  int  len = strlen(filename);

  if( strncmp(&filename[len-4],".prc",4)==0 )
    {
      hmm = read_HMM_PRC_binary(filename);
    }
  else if( strncmp(&filename[len-4],".mod",4)==0 ) 
    {
      hmm = read_HMM_SAM3x_ascii(filename);
    }
  else if( strncmp(&filename[len-4],".hmm",4)==0 ) 
    {
      hmm = read_HMM_HMMER2x_ascii(filename);
    }
  else if( (strncmp(&filename[len-3],".fa" ,3)==0) || 
	   (strncmp(&filename[len-4],".seq",4)==0)    )
    {
      hmm = read_HMM_FASTA(filename);
    }
  else if( (strncmp(&filename[len-5],".pbla"    ,5)==0) ||
	   (strncmp(&filename[len-4],".psi"     ,4)==0) ||
	   (strncmp(&filename[len-9],".psiblast",9)==0) ||
	   (strncmp(&filename[len-5],".prof"    ,5)==0) ||
	   (strncmp(&filename[len-8],".profile" ,8)==0)    )
    {
      hmm = read_PSIBLAST_binary(filename);
    }
  else 
    {
      die2( "File '%s' has an unknown extension!\n\n"
	    "Known extensions:\n\n" KNOWN_EXTENSIONS,
	    filename );
    };

  if(hmm == NULL)
    die2("Weird, failed to read '%s', not sure why!", filename);

  return hmm;
}

// check that a file has a known extension
//
// known extension => return 1, unknown extension => return 0
//
// see read_HMM() above for more details
//
int known_extension(char *filename)
{
  int  len = strlen(filename);

  if( (strncmp(&filename[len-4],".prc"     ,4)==0) ||
      (strncmp(&filename[len-4],".mod"     ,4)==0) ||
      (strncmp(&filename[len-4],".hmm"     ,4)==0) ||
      (strncmp(&filename[len-3],".fa"      ,3)==0) ||
      (strncmp(&filename[len-4],".seq"     ,4)==0) ||
      (strncmp(&filename[len-5],".pbla"    ,5)==0) ||
      (strncmp(&filename[len-4],".psi"     ,4)==0) ||
      (strncmp(&filename[len-9],".psiblast",9)==0) ||
      (strncmp(&filename[len-5],".prof"    ,5)==0) ||
      (strncmp(&filename[len-8],".profile" ,8)==0)    )
    return 1;

  return 0;
}
