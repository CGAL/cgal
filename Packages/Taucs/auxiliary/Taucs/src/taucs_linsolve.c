/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*#include <stdarg.h>*/
#include <assert.h>
#include <math.h>

#define TAUCS_CORE_CILK

#include "taucs.h"

#define TRUE 1
#define FALSE 0

#ifdef TAUCS_CORE_GENERAL

/*********************************************************/
/* utility routines                                      */
/*********************************************************/

static int element_size(int flags)
{
  if (flags & TAUCS_SINGLE)   return sizeof(taucs_single);
  if (flags & TAUCS_DOUBLE)   return sizeof(taucs_double);
  if (flags & TAUCS_SCOMPLEX) return sizeof(taucs_scomplex);
  if (flags & TAUCS_DCOMPLEX) return sizeof(taucs_dcomplex);
  if (flags & TAUCS_INT)      return sizeof(int);
  assert(0);
  return -1;
}

/*********************************************************/
/* argument parsing                                      */
/*********************************************************/

int taucs_getopt_boolean(char* cmd, void* args[], char* name, int* x) {
  int lc = strlen(cmd);
  int ln = strlen(name);
  if (!strncmp(cmd,name,ln)) {
    if (lc > ln && cmd[ln] == '.') return 0;
    if (lc > ln && cmd[ln] == '=') {
      if (cmd[ln+1] == '#') {
	unsigned int p;
	if (sscanf(cmd+ln+2,"%u",&p) == 1) {
	  unsigned int i;
	  for (i=0; args[i]; i++) {
	    if (i==p) { *x = *( (int*) args[i] ); return 1; }
	  }
	  taucs_printf("taucs: WARNING, pointer argument out of range in [%s]\n",cmd);
	}
	taucs_printf("taucs: WARNING, illegal pointer argument in [%s]\n",cmd);
	return 0;
      }
      if (!strcmp(cmd+ln+1,"true")) {
	*x = TRUE;
	return 1;
      }
      if (!strcmp(cmd+ln+1,"false")) {
	*x = FALSE;
	return 1;
      }
    }
    taucs_printf("taucs: WARNING, illegal argument in [%s]\n",cmd);
  }

  return 0;
}

int taucs_getopt_double(char* cmd, void* args[], char* name, double* x) {
  int lc = strlen(cmd);
  int ln = strlen(name);
  if (!strncmp(cmd,name,ln)) {
    if (lc > ln && cmd[ln] == '.') return 0;
    if (lc > ln && cmd[ln] == '=') {
      if (cmd[ln+1] == '#') {
	unsigned int p;
	if (sscanf(cmd+ln+2,"%u",&p) == 1) {
	  unsigned int i;
	  for (i=0; args[i]; i++) {
	    if (i==p) { *x = *( (double*) args[i] ); return 1; }
	  }
	  taucs_printf("taucs: WARNING, pointer argument out of range in [%s]\n",cmd);
	}
	taucs_printf("taucs: WARNING, illegal pointer argument in [%s]\n",cmd);
	return 0;
      }
      if (sscanf(cmd+ln+1,"%le",x) == 1) {
	return 1;
      }
    }
    taucs_printf("taucs: WARNING, illegal argument in [%s]\n",cmd);
  }

  return 0;
}


int taucs_getopt_pointer(char* cmd, void* args[], char* name, void** x) {
  int lc = strlen(cmd);
  int ln = strlen(name);
  if (!strncmp(cmd,name,ln)) {
    if (lc > ln && cmd[ln] == '.') return 0;
    if (lc > ln && cmd[ln] == '=') {
      if (cmd[ln+1] == '#') {
	unsigned int p;
	if (sscanf(cmd+ln+2,"%u",&p) == 1) {
	  unsigned int i;
	  for (i=0; args[i]; i++)
	    if (i==p) { *x = *( (void**) args[i] ); return 1; }
	  taucs_printf("taucs: WARNING, pointer argument out of range in [%s]\n",cmd);
	}
	taucs_printf("taucs: WARNING, illegal pointer argument in [%s]\n",cmd);
	return 0;
      }
    }
    taucs_printf("taucs: WARNING, illegal argument in [%s]\n",cmd);
  }

  return 0;
}

int taucs_getopt_string(char* cmd, void* args[], char* name, char** x) {
  int lc = strlen(cmd);
  int ln = strlen(name);
  if (!strncmp(cmd,name,ln)) {
    if (lc > ln && cmd[ln] == '.') return 0;
    if (lc > ln && cmd[ln] == '=') {
      if (cmd[ln+1] == '#') {
	unsigned int p;
	if (sscanf(cmd+ln+2,"%u",&p) == 1) {
	  unsigned int i;
	  for (i=0; args[i]; i++)
	    if (i==p) { *x = *( (char**) args[i] ); return 1; }
	  taucs_printf("taucs: WARNING, pointer argument out of range in [%s]\n",cmd);
	}
	taucs_printf("taucs: WARNING, illegal pointer argument in [%s]\n",cmd);
	return 0;
      }
      *x = cmd+ln+1;
      return 1;
    }
    taucs_printf("taucs: WARNING, illegal argument in [%s]\n",cmd);
  }

  return 0;
}

/*********************************************************/
/* Generic Factor routines                               */
/* (Experimental, unstable interface)                    */
/*********************************************************/

#define TAUCS_FACTORTYPE_NONE           0
#define TAUCS_FACTORTYPE_LLT_SUPERNODAL 1
#define TAUCS_FACTORTYPE_LLT_CCS        2
#define TAUCS_FACTORTYPE_LDLT_CCS       3
#define TAUCS_FACTORTYPE_LLT_OOC        4
#define TAUCS_FACTORTYPE_LU_OOC         5

typedef struct {
  int   n;
  int   flags;
  int   type;
  int*  rowperm;
  int*  colperm;
  void* L;
} taucs_factorization;

static void taucs_linsolve_free(void* vF)
{
  taucs_factorization* F = (taucs_factorization*) vF;

  if (!F) return;

  if (F->type == TAUCS_FACTORTYPE_LLT_SUPERNODAL)
    taucs_supernodal_factor_free(F->L);
  if (F->type == TAUCS_FACTORTYPE_LLT_CCS)
    taucs_ccs_free(F->L);
  taucs_free(F->rowperm);
  taucs_free(F->colperm);
  taucs_free(F);
}

int taucs_linsolve(taucs_ccs_matrix* A, 
		   void**            F,
		   int               nrhs,
		   void*             X,
		   void*             B,
		   char*             options[],
		   void*             opt_arg[])
{
  int retcode = TAUCS_SUCCESS;
  double tw,tc;

  int i;
  taucs_ccs_matrix*    PAPT    = NULL;
  int*                 rowperm = NULL;
  int*                 colperm = NULL;
  taucs_factorization* f       = NULL;

  void* PX = NULL;
  void* PB = NULL;

  void*  opt_context   = NULL;
  double opt_cilk_nproc= -1.0;
#ifdef TAUCS_CILK
  int    local_context = FALSE;
#endif

  int    opt_factor    =  1;
  int    opt_symbolic  =  1;
  int    opt_numeric   =  1;

  int    opt_llt       =  0;
  int    opt_lu        =  0;

  int    opt_mf        =  0;
  int    opt_ll        =  0;

  double opt_maxdepth  = 0.0; /* default meaning no limit */

  int    opt_ooc       =  0;
  char*            opt_ooc_name   = NULL;
  void*            opt_ooc_handle = NULL;
  int              local_handle_open   = FALSE;
  int              local_handle_create = FALSE;
  double           opt_ooc_memory = -1.0;

  char*            opt_ordering   = NULL;

  int    opt_cg          = 0;
  int    opt_minres      = 0;
  double opt_maxits      = 300.0;
  double opt_convergetol = 1e-6;

  int    opt_sg          = 0;
  int    opt_amwb        = 0;
  double opt_amwb_sg     = 1;
  double opt_amwb_rnd    = 170566;
  taucs_ccs_matrix* M    = NULL;
  taucs_ccs_matrix* PMPT = NULL;

  if (!A && nrhs==0) {
    if (F) taucs_linsolve_free(*F);
    *F = NULL;
    return TAUCS_SUCCESS;
  }

  if (options) {
    for (i=0; options[i]; i++) {
      int understood = FALSE;
      
      understood |= taucs_getopt_pointer(options[i],opt_arg,"taucs.cilk.context",&opt_context);
      understood |= taucs_getopt_double(options[i],opt_arg,"taucs.cilk.nproc",&opt_cilk_nproc);

      understood |= taucs_getopt_boolean(options[i],opt_arg,"taucs.approximate.sg",&opt_sg); 
      understood |= taucs_getopt_boolean(options[i],opt_arg,"taucs.approximate.amwb",&opt_amwb); 
      understood |= taucs_getopt_double(options[i],opt_arg,"taucs.approximate.amwb.randomseed",&opt_amwb_rnd); 
      understood |= taucs_getopt_double(options[i],opt_arg,"taucs.approximate.amwb.subgraphs",&opt_amwb_sg); 
      understood |= taucs_getopt_boolean(options[i],opt_arg,"taucs.factor",&opt_factor); 
      understood |= taucs_getopt_boolean(options[i],opt_arg,"taucs.factor.symbolic",&opt_symbolic); 
      understood |= taucs_getopt_boolean(options[i],opt_arg,"taucs.factor.numeric",&opt_numeric); 
      understood |= taucs_getopt_boolean(options[i],opt_arg,"taucs.factor.LLT",&opt_llt); 
      understood |= taucs_getopt_boolean(options[i],opt_arg,"taucs.factor.LU",&opt_lu); 
      understood |= taucs_getopt_boolean(options[i],opt_arg,"taucs.factor.mf",&opt_mf); 
      understood |= taucs_getopt_boolean(options[i],opt_arg,"taucs.factor.ll",&opt_ll); 
      understood |= taucs_getopt_string(options[i],opt_arg,"taucs.factor.ordering",&opt_ordering); 
      understood |= taucs_getopt_double(options[i],opt_arg,"taucs.maxdepth",&opt_maxdepth); 

      understood |= taucs_getopt_boolean(options[i],opt_arg,"taucs.ooc",&opt_ooc); 
      understood |= taucs_getopt_string (options[i],opt_arg,"taucs.ooc.basename",&opt_ooc_name); 
      understood |= taucs_getopt_pointer(options[i],opt_arg,"taucs.ooc.iohandle",&opt_ooc_handle); 
      understood |= taucs_getopt_double (options[i],opt_arg,"taucs.ooc.memory",  &opt_ooc_memory); 

      understood |= taucs_getopt_boolean(options[i],opt_arg,"taucs.solve.cg",&opt_cg); 
      understood |= taucs_getopt_boolean(options[i],opt_arg,"taucs.solve.minres",&opt_minres); 
      understood |= taucs_getopt_double(options[i],opt_arg,"taucs.solve.maxits",&opt_maxits); 
      understood |= taucs_getopt_double(options[i],opt_arg,"taucs.solve.convergetol",&opt_convergetol); 

      if (!understood) taucs_printf("taucs_linsolve: illegal option [[%s]]\n",
				    options[i]);
    }
  }

  /* First, construct a preconditioner if one is needed */

  if (opt_amwb) {
    M = taucs_amwb_preconditioner_create(A,(int) opt_amwb_rnd,opt_amwb_sg,0 /* stretch flag */);
    if (!M)
      taucs_printf("taucs_linsolve: AMWB preconditioner construction failed, using A\n");
  }

  /* First, decide on the kind of factorization */

  if (opt_factor) {
    taucs_printf("taucs_linsolve: preparing to factor\n");
    f = (taucs_factorization*) taucs_malloc(sizeof(taucs_factorization));
    if (!f) {
      taucs_printf("taucs_factor: memory allocation\n");
      retcode = TAUCS_ERROR_NOMEM;
      goto release_and_return;
    }
    f->n       = A->n;
    f->type    = TAUCS_FACTORTYPE_NONE;
    f->flags   = A->flags; /* remember data type */

    if (!opt_numeric && (nrhs > 0)) {
      taucs_printf("taucs_linsolve: WARNING, you can't solve without a numeric factorization\n");
      opt_numeric = 1;
    }

    /* decide on ordering and order */  

    if (!opt_ordering)
      opt_ordering = opt_lu ? 
	"colamd" : 
#if defined(TAUCS_CONFIG_METIS)
	"metis"
#elif defined(TAUCS_CONFIG_GENMMD)
	"genmmd"
#elif defined(TAUCS_CONFIG_AMD)
	"amd"
#endif
	;
  
    taucs_printf("taucs_linsolve: ordering (llt=%d, lu=%d, ordering=%s)\n",
		 opt_llt,opt_lu,opt_ordering);
    tw = taucs_wtime();
    tc = taucs_ctime();
    taucs_ccs_order(M ? M : A,&rowperm,&colperm,opt_ordering);
    if (!rowperm) {
      taucs_printf("taucs_factor: ordering failed\n");
      retcode = TAUCS_ERROR_NOMEM;
      goto release_and_return;
    } else
      taucs_printf("taucs_linsolve: ordering time %.02e seconds (%.02e seconds CPU time)\n",taucs_wtime()-tw,taucs_ctime()-tc);

    f->rowperm = rowperm;
    f->colperm = colperm;

    if (opt_llt) {
      taucs_printf("taucs_linsolve: starting LLT factorization\n");
      if (M) {
	taucs_printf("taucs_linsolve: pre-factorization permuting of M\n");
	PMPT = taucs_ccs_permute_symmetrically(M,rowperm,colperm);
	if (!PMPT) {
	  taucs_printf("taucs_factor: permute rows and columns failed\n");
	  retcode = TAUCS_ERROR_NOMEM;
	  goto release_and_return;
	}
      } else {
	taucs_printf("taucs_linsolve: pre-factorization permuting of A\n");
	PAPT = taucs_ccs_permute_symmetrically(A,rowperm,colperm);
	if (!PAPT) {
	  taucs_printf("taucs_factor: permute rows and columns failed\n");
	  retcode = TAUCS_ERROR_NOMEM;
	  goto release_and_return;
	}
      }

      if (opt_ooc) {
	taucs_printf("taucs_linsolve: starting OOC LLT factorization\n");
	if ((!opt_ooc_name && !opt_ooc_handle)
	    || (opt_ooc_name && opt_ooc_handle)) {
	  taucs_printf("taucs_linsolve: ERROR, you must specify either a basename or an iohandle for an out-of-core factorization\n");
	  retcode = TAUCS_ERROR_BADARGS;
	  goto release_and_return;
	}

	if (opt_ooc_name) {
	  opt_ooc_handle = taucs_io_open_multifile(opt_ooc_name);
	  if (opt_ooc_handle) {
	    local_handle_open = TRUE;
	  } else {
	    opt_ooc_handle = taucs_io_create_multifile(opt_ooc_name);
	    if (opt_ooc_handle) {
	      local_handle_create = TRUE;
	    } else {
	      taucs_printf("taucs_linsolve: ERROR, could neither open nor create file [%s]\n",
			   opt_ooc_name);
	      retcode = TAUCS_ERROR;
	      goto release_and_return;
	    }
	  }
	}
	taucs_printf("taucs_linsolve: ooc file created?=%d opened?=%d\n",
		     local_handle_create,local_handle_open);
	if (opt_ooc_memory < 0.0) opt_ooc_memory = taucs_available_memory_size();
	if (taucs_ooc_factor_llt(PMPT ? PMPT : PAPT, 
				 opt_ooc_handle, opt_ooc_memory) == TAUCS_SUCCESS)
	  f->type = TAUCS_FACTORTYPE_LLT_OOC;
	else {
	  retcode = TAUCS_ERROR;
	  goto release_and_return;
	}
      } else { /* in-core */
	taucs_printf("taucs_linsolve: starting IC LLT factorization\n");
	if (opt_mf) {
	  taucs_printf("taucs_linsolve: starting IC LLT MF factorization\n");

#ifdef TAUCS_CILK
	  if (!opt_context) {
	    char* argv[16]  = {"program_name" };
	    char  bufs[16][16];
	    int   p = 0;
	    int   argc;

	    for (argc=1; argc<16; argc++) argv[argc] = 0;
	    argc = 1;
	    
	    argv[argc++] = "--pthread-stacksize";
	    argv[argc++] = "2000000";
	    argv[argc++] = "--stack";
	    argv[argc++] = "2000000";

	    if (opt_cilk_nproc > 0) {
	      argv[argc++] = "--nproc";
	      sprintf(bufs[p],"%d",(int) opt_cilk_nproc);
	      argv[argc++] = bufs[p++];
	    }

	    taucs_printf("taucs_ccs_linsolve:_cilk_init\n");
	    opt_context = Cilk_init(&argc,argv);
	    local_context = TRUE;
	  }
#endif

	  if (!opt_numeric && opt_symbolic)
	    f->L = taucs_ccs_factor_llt_symbolic_maxdepth(PMPT ? PMPT : PAPT,(int) opt_maxdepth);

	  if (opt_numeric && !opt_symbolic) {
	    int rc;
	    if (!F 
		|| !(*F) 
		|| ((taucs_factorization*)*F)->type != TAUCS_FACTORTYPE_LLT_SUPERNODAL) {
	      taucs_printf("taucs_linsolve: ERROR, you need to provide a symbolic factorization for a numeric factorization\n");
	      retcode = TAUCS_ERROR_BADARGS;
	      goto release_and_return;
	    }
	    f->L = ((taucs_factorization*)*F)->L;
	    taucs_supernodal_factor_free_numeric(f->L);

#ifdef TAUCS_CILK	  
	    rc = EXPORT(taucs_ccs_factor_llt_numeric)(opt_context, PMPT ? PMPT : PAPT, f->L);
#else
	    rc = taucs_ccs_factor_llt_numeric(PMPT ? PMPT : PAPT, f->L);
#endif
	  }

	  if (opt_numeric && opt_symbolic) {
#ifdef TAUCS_CILK	  
	    f->L = EXPORT(taucs_ccs_factor_llt_mf_maxdepth)(opt_context,
							    PMPT ? PMPT : PAPT,
							    (int) opt_maxdepth);
#else
	    f->L = taucs_ccs_factor_llt_mf_maxdepth(PMPT ? PMPT : PAPT,(int) opt_maxdepth);
#endif
	  }

	  if (! (f->L) ) {
	    taucs_printf("taucs_factor: factorization failed\n");
	    retcode = TAUCS_ERROR;
	    goto release_and_return;
	  } else {
	    f->type = TAUCS_FACTORTYPE_LLT_SUPERNODAL;
	  }

#ifdef TAUCS_CILK
	  if (local_context) Cilk_terminate((CilkContext*) opt_context);
#endif
	} else if (opt_ll || TRUE) { /* this will be the default LLT */
	  taucs_printf("taucs_linsolve: starting IC LLT LL factorization\n");
	  f->L = taucs_ccs_factor_llt_ll_maxdepth(PMPT ? PMPT : PAPT,(int) opt_maxdepth);
	  if (! (f->L) ) {
	    taucs_printf("taucs_factor: factorization failed\n");
	    retcode = TAUCS_ERROR;
	    goto release_and_return;
	  } else {
	    f->type = TAUCS_FACTORTYPE_LLT_SUPERNODAL;
	  }
	} /* left-looking */
      } /* in-core */
    } /* llt */
  }

  if (nrhs > 0) {

    int             (*precond_fn)(void*,void* x,void* b) = NULL;
    void*           precond_arg = NULL;
    int    j;

    if (!f) {
      if (!F || !(*F)) {
	taucs_printf("taucs_linsolve: can't solve, no factorization\n");
	retcode = TAUCS_ERROR;
	goto release_and_return;
      } else {
	if (F && *F)
	  f = (taucs_factorization*) *F;
	else {
	  taucs_printf("taucs_linsolve: can't solve, no factorization\n");
	  retcode = TAUCS_ERROR;
	  goto release_and_return;
	} 
      }
    }

    taucs_printf("taucs_linsolve: preparing to solve\n");
    PX = (void*) taucs_malloc(element_size(A->flags)*nrhs*(A->n));
    PB = (void*) taucs_malloc(element_size(A->flags)*nrhs*(A->n));
    if (!PB || !PX) {
      taucs_printf("taucs_linsolve: memory allocation\n");
      retcode = TAUCS_ERROR_NOMEM;
      goto release_and_return;
    }      

    switch (f->type) {
    case TAUCS_FACTORTYPE_NONE:
      taucs_printf("taucs_linsolve: WARNING, no preconditioner\n");
      precond_fn  = NULL;
      precond_arg = NULL;
      break;
    case TAUCS_FACTORTYPE_LLT_SUPERNODAL:
      precond_fn  = taucs_supernodal_solve_llt;
      precond_arg = f->L;
      break;
    case TAUCS_FACTORTYPE_LLT_CCS:
      precond_fn  = taucs_ccs_solve_llt;
      precond_arg = f->L;
      break;
    case TAUCS_FACTORTYPE_LDLT_CCS:
      precond_fn  = taucs_ccs_solve_ldlt;
      precond_arg = f->L;
      break;
    case TAUCS_FACTORTYPE_LLT_OOC:
      precond_fn  = taucs_ooc_solve_llt;
      precond_arg = opt_ooc_handle;
      break;
    default:
      assert(0);
    }
    
    taucs_printf("taucs_linsolve: pre-solve permuting of A\n");
    if (!PAPT) PAPT = taucs_ccs_permute_symmetrically(A,f->rowperm,f->colperm);
    if (!PAPT) {
      taucs_printf("taucs_factor: permute rows and columns failed\n");
      retcode = TAUCS_ERROR_NOMEM;
      goto release_and_return;
    }

    for (j=0; j<nrhs; j++) {
      int ld = (A->n) * element_size(A->flags);
      
      taucs_vec_permute (A->n,A->flags,(char*)B+j*ld,(char*)PB+j*ld,f->rowperm);

      if (opt_cg) {
	taucs_conjugate_gradients (PAPT,
				   precond_fn, precond_arg,
				   (char*)PX+j*ld, (char*)PB+j*ld,
				   (int) opt_maxits,
				   opt_convergetol);
	
      } else if (opt_minres) {
	taucs_minres              (PAPT,
				   precond_fn, precond_arg,
				   (char*)PX+j*ld, (char*)PB+j*ld,
				   (int) opt_maxits,
				   opt_convergetol);
      } else if (precond_fn) {
	(*precond_fn)(precond_arg,(char*)PX+j*ld,(char*)PB+j*ld);
      } else {
	taucs_printf("taucs_linsolve: I don't know how to solve!\n");
	retcode = TAUCS_ERROR;
	goto release_and_return;
      }
   
      taucs_vec_ipermute(A->n,A->flags,(char*)PX+j*ld,(char*)X+j*ld,f->rowperm);
    }

    taucs_free(PB);
    taucs_free(PX);
  }

  if (F) {
    if (local_handle_open)   taucs_io_close(opt_ooc_handle);
    if (local_handle_create) taucs_io_close(opt_ooc_handle);
    
    *F = f;
  } else {
    if (f->type == TAUCS_FACTORTYPE_LLT_OOC) {
      if (local_handle_open)   taucs_io_close(opt_ooc_handle);
      if (local_handle_create) taucs_io_delete(opt_ooc_handle);
    }
    taucs_linsolve_free(f);
  }

  taucs_ccs_free(PMPT);
  taucs_ccs_free(PAPT);
  taucs_ccs_free(M);

  return retcode;

release_and_return:
  taucs_printf("taucs_linsolve: an error occured, releasing resources and bailing out\n");
#ifdef TAUCS_CILK
  if (local_context)  Cilk_terminate((CilkContext*) opt_context);
#endif
  taucs_free(rowperm);
  taucs_free(colperm);
  taucs_ccs_free(PMPT);
  taucs_ccs_free(PAPT);
  taucs_ccs_free(M);
  taucs_free(PB);
  taucs_free(PX);
  taucs_free(f);
  return retcode;
}

#endif /* TAUCS_CORE_GENERAL */

/*********************************************************/
/*                                                       */
/*********************************************************/

