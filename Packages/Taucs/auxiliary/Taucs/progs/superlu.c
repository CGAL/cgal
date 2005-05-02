/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <taucs.h>

/*********************************************************/
/*                                                       */
/*********************************************************/

#define max(x,y) ( ((x) > (y)) ? (x) : (y) )

int        M,N;					/* Size of matrix */

/*********************************************************/
/*                                                       */
/*********************************************************/

static double twonorm(int n, double* v)
{
  /*
  double norm;
  int i;

  for (i=0, norm=0.0; i<n; i++) norm += v[i]*v[i];

  norm = sqrt(norm);
  return norm;
  */

  double norm, ssq, scale, absvi;
  int i;

  if (n==1) return fabs(v[0]);

  scale = 0.0;
  ssq   = 1.0;

  for (i=0; i<n; i++) {
    if ( v[i] != 0 ) {
      absvi = fabs(v[i]);
      if (scale < absvi) {
	ssq   = 1.0 + ssq * (scale/absvi)*(scale/absvi);
	scale = absvi;
      } else
	ssq   = ssq + (absvi/scale)*(absvi/scale);
    }
  }
  return scale * sqrt( ssq );
}
/*********************************************************/
/*                                                       */
/*********************************************************/

void usage(int argc, char* argv[])
{
  printf("%s OPTIONS\n",argv[0]);
  printf("  -help or -h or no arguments at all\n");
  printf("        this help\n");
  printf("  MATRIX OPTIONS:\n");
  printf("        -hb       Harwll-Boeing-filename (matrix from file)\n");
  printf("        -ijv      ijvfilename (matrix from file)\n");
  printf("        -mtx      mtxfilename (matrix from file)\n");
  printf("        -ccs      ccsfilename (matrix from file)\n");
  printf("        -mesh3d      X Y Z (X-by-Y-by-Z poisson problem)\n");
  printf("        -mesh2d          n (n-by-n poisson problem)\n");
  printf("        -mat_type [anisotropic_y/anisotropic_x/dirichlet/neumann]\n");
  printf("                  mat_type only applicable to -mesh2d\n");
  printf("                  (3D meshes are always Neumann BC)\n");
  printf("                  anisotropic problems are Neumann\n");
  printf("  RIGHT-HAND SIDE OPTIONS:\n");
  printf("        -n+rhs filename (order n followed by rhs, ascii)\n");
  printf("  FACTOR/SOLVE OPTIONS:\n");
  printf("  -snmf\n");
  printf("        supernodal-multifrontal\n");
  printf("  -snll\n");
  printf("        supernodal-left-looking\n");
  printf("  -symb\n");
  printf("        supernodal-multifrontal, separate symbolic & numeric phases\n");
  printf("  -ooc\n");
  printf("        out-of-core supernodal left-looking\n");
  printf("  OOC OPTIONS:\n");
  printf("  -matrixfile basename\n");
  printf("        base name for files; default is /tmp/taucs\n");
  printf("  -memory M\n");
  printf("        amount of memory in megabytes to use;\n");
  printf("        if not given, TAUCS tries to figure an optimal amount;\n");
  printf("  OTHER OPTIONS:\n");
  printf("  -log filename\n");
  printf("        write log into filename; can be stdout and stderr\n");
  exit(0);
}

/*********************************************************/
/*                                                       */
/*********************************************************/

main(int argc, char* argv[])
{
  double wtime_order;
  double wtime_permute;
  double wtime_factor;
  double wtime_solve;
  double wtime_precond_create;
  double ctime_factor;

  int i,j;
  double NormErr;
  taucs_ccs_matrix*  A;
  taucs_ccs_matrix*  PAPT;
  taucs_ccs_matrix*  L;

  double*      X = NULL;
  double*      B = NULL;
  double*      PX;
  double*      PB;
  double*      NX;

  char*        ordering = "metis";
  char*        mat_type = "neumann";
  int*         perm;
  int*         invperm;
  int          superlu_flag = 0;
  int          snmf_flag = 0;
  int          snll_flag = 0;
  int          symb_flag = 0;
  int          mesh2d_flag = 0,mesh2d_size = 0;

  int          ooc_flag = 0;
  int          panelize = 0;
  double       memory_mb = -1.0;
  char*        matrixfile = "/tmp/taucs.L";
  taucs_io_handle* oocL;
  
  int             (*precond_fn)(void*,double* x,double* b);
  void*           precond_args;


  /***********************************************************/
  /* Read arguments: log file, matrix, memory size, ooc name */
  /***********************************************************/

  if ((argc == 1) ||((argc == 2) && !strncmp(argv[1],"-h",2)))
    usage(argc,argv);

  A = NULL;
  for (i=1; i<argc; i++) {
    if (!strcmp(argv[i],"-superlu")) superlu_flag = 1;
    if (!strcmp(argv[i],"-snmf")) snmf_flag = 1;
    if (!strcmp(argv[i],"-snll")) snll_flag = 1;
    if (!strcmp(argv[i],"-symb")) symb_flag = 1;
    if (!strcmp(argv[i],"-ooc"))  ooc_flag = 1;
    if (!strcmp(argv[i],"-log") && i <= argc-1 ) {
      i++;
      taucs_logfile(argv[i]);
    }

    if (!strcmp(argv[i],"-panelize") && i <= argc-1) {
      i++;
      if (sscanf(argv[i],"%d",&panelize) != 1) {
	taucs_printf("0 (smart), 1 (in-core), or 2 (single supernode) follow -panelize argument\n");
	exit(1);
      }
    }

    if (!strcmp(argv[i],"-memory") && i <= argc-1) {
      i++;
      if (sscanf(argv[i],"%lf",&memory_mb) != 1) {
	taucs_printf("memory size in MB must follow -memory argument\n");
	exit(1);
      }
    }

    if (!strcmp(argv[i],"-matrixfile") && i <= argc-1 ) {
      i++;
      matrixfile = argv[i];
    }
    if (!strcmp(argv[i],"-ordering") && i <= argc-1 ) {
      i++;
      ordering = argv[i];
    }
    if (!strcmp(argv[i],"-mat_type") && i <= argc-1 ) {
      i++;
      mat_type = argv[i];
    }

    if (!strcmp(argv[i],"-hb") && i <= argc-1 ) {
      int  nrows,ncols,nnz,j;
      char fname[256];
      char type[3];

      i++;
      for (j=0; j<256; j++) fname[j] = ' ';
      strcpy(fname,argv[i]);
      taucs_printf("main: reading HB matrix %s\n",argv[i]);
      ireadhb_(fname,type,&nrows,&ncols,&nnz);
      A = taucs_ccs_create(nrows,ncols,nnz);
      if (type[1] == 's' || type[1] == 'S')
	A->flags = TAUCS_SYMMETRIC | TAUCS_LOWER;
      dreadhb_(fname,&nrows,&ncols,&nnz,
	       A->colptr,A->rowind,A->values);
      /* make indices 0-based */
      for (j=0; j<=ncols; j++) ((A->colptr)[j])--;
      for (j=0; j<nnz;    j++) ((A->rowind)[j])--;
      taucs_printf("main: done reading\n");
    }

    if (!strcmp(argv[i],"-mtx") && i <= argc-1) {
      i++;
      taucs_printf("main: reading mtx matrix %s\n",argv[i]);
      A = taucs_ccs_read_mtx (argv[i],TAUCS_SYMMETRIC);
      taucs_printf("main: done reading\n");
    }

    if (!strcmp(argv[i],"-ijv") && i <= argc-1) {
      i++;
      taucs_printf("main: reading ijv matrix %s\n",argv[i]);
      A = taucs_ccs_read_ijv (argv[i],TAUCS_SYMMETRIC);
      taucs_printf("main: done reading\n");
    }

    if (!strcmp(argv[i],"-ccs") && i <= argc-1) {
      i++;
      taucs_printf("main: reading ccs matrix %s\n",argv[i]);
      A = taucs_ccs_read_ccs (argv[i],TAUCS_SYMMETRIC);
      taucs_printf("main: done reading\n");
    }

    if (!strcmp(argv[i],"-mesh2d") && i <= argc-1) {
      mesh2d_flag = 1;
      taucs_printf("A is a mesh2d\n");
      i++;
      if (sscanf(argv[i],"%d",&mesh2d_size) != 1) {
	taucs_printf("mesh size (n, where the mesh is n-by-n) must follow -mesh2d argument\n");
	exit(1);
      }
    }
    if (!strcmp(argv[i],"-mesh3d") && i <= argc-3) {
      int X,Y,Z;
      taucs_printf("A is a mesh3d\n");
      if (sscanf(argv[i+1],"%d",&X) != 1 
	  || sscanf(argv[i+2],"%d",&Y) != 1 
	  || sscanf(argv[i+3],"%d",&Z) != 1) {
	taucs_printf("mesh size (X Y Z must follow -mesh3d argument\n");
	exit(1);
      }
      i += 3;
      taucs_printf("main: creating mesh\n");
      A = taucs_ccs_generate_mesh3d(X,Y,Z);
    }

    if (!strcmp(argv[i],"-n+rhs") && i <= argc-1) {
      FILE* f;
      int n,j,nnz;

      i++;

      taucs_printf("main: reading right-hand side %s\n",argv[i]);
      f=fopen(argv[i],"r");
      assert(f);
      fscanf(f,"%d",&n);
      B=(double*) malloc(n*sizeof(double));
      nnz = 0;
      for (j=0; j<n; j++) {
	fscanf(f,"%lg",B+j);
	if (B[j]) nnz++;
      }
      fclose(f);
      taucs_printf("main: done reading rhs, %d nonzeros\n",nnz);
    }

  }

  taucs_printf("Chosen Ordering: %s\n",ordering);

  if (mesh2d_flag)
    {
      taucs_printf("Matrix type is %s\n",mat_type);
      taucs_printf("Grid Size is %d\n",mesh2d_size);
      A = taucs_ccs_generate_mesh2d(mesh2d_size,mat_type);
    }

  if (!A) {
    taucs_printf("matrix argument not given or matrix file not found\n");
    usage(argc,argv);
  }
  N = M = A->n;
  
  /***********************************************************/
  /* Create exact solution, compute right-hand-side          */
  /***********************************************************/

  if (!X) {
    X=(double*)malloc(N*sizeof(double));
    for(i=0; i<N; i++)  X[i]=(double)random()/RAND_MAX;
  } else
    taucs_printf("iter: not using a random X, already allocated\n");

  if (!B) {
    B=(double*)malloc(N*sizeof(double));
    taucs_ccs_times_vec(A,X,B);
  } else {
    for(i=0; i<N; i++)  X[i]= 0.0 / 0.0;
  }

  NX=(double*)malloc(N*sizeof(double));

  PX=(double*)malloc(N*sizeof(double));
  PB=(double*)malloc(N*sizeof(double));

  /***********************************************************/
  /* Compute column ordering                                 */
  /***********************************************************/

  /***********************************************************/
  /* factor                                                  */
  /***********************************************************/

  {
    int n;
    double unit,curr;

    n = A->n;
    unit = (n-1.)+n;

    wtime_order = taucs_wtime();
    taucs_ccs_order(A,&perm,&invperm,ordering);
    wtime_order = taucs_wtime() - wtime_order;
    taucs_printf("\tOrdering time = % 10.3f seconds\n",wtime_order);

    if (!perm) {
      taucs_printf("\tOrdering Failed\n");
      exit(1);
    }

    if (0) {
      int i;
      FILE* f;
      f=fopen("p.ijv","w");
      for (i=0; i<n; i++) fprintf(f,"%d\n",perm[i]+1);
      fclose(f);
    }

    wtime_permute = taucs_wtime();
    PAPT = taucs_ccs_permute_symmetrically(A,perm,invperm);
    wtime_permute = taucs_wtime() - wtime_permute;
    taucs_printf("\tPermute time  = % 10.3f seconds\n",wtime_permute);

    wtime_factor = taucs_wtime();
    ctime_factor = taucs_ctime();

    if (superlu_flag) {
      void* taucs_ccs_factor_superlu(void*);
      int taucs_ccs_solve_superlu(void*,double*,double*);
      L = taucs_ccs_factor_superlu(PAPT);
      precond_args = L;
      precond_fn   = taucs_ccs_solve_superlu;
    }
    else if (snmf_flag) {

      taucs_ccs_matrix* C;
      L = taucs_ccs_factor_llt_mf(PAPT);
      precond_args = L;
      precond_fn   = taucs_supernodal_solve_llt;
      
      /*
      C = taucs_supernodal_factor_to_ccs(L);
      taucs_ccs_write_ijv(C,"C.ijv");
      precond_args = C;
      precond_fn   = taucs_ccs_solve_llt;
      */
    } else if (ooc_flag) {
#define TESTING
#ifdef TESTING
      int taucs_ooc_factor_llt_panelchoice(taucs_ccs_matrix* A, 
					   taucs_io_handle* handle,
					   double memory,
					   int panelization_method);

      int c;
      oocL = taucs_io_create_multifile(matrixfile);
      assert(oocL);
      if (memory_mb == -1.0) memory_mb = taucs_available_memory_size()/1048576.0;
      taucs_ooc_factor_llt_panelchoice(PAPT, oocL, memory_mb*1048576.0,panelize);
      precond_args = oocL;
      precond_fn   = taucs_ooc_solve_llt;
#else
      int c;
      oocL = taucs_io_create_multifile(matrixfile);
      assert(oocL);
      if (memory_mb == -1.0) memory_mb = taucs_available_memory_size()/1048576.0;
      taucs_ooc_factor_llt(PAPT, oocL, memory_mb*1048576.0);
      precond_args = oocL;
      precond_fn   = taucs_ooc_solve_llt;
#endif
    } else if (snll_flag) {
      L = taucs_ccs_factor_llt_ll(PAPT);
      precond_args = L;
      precond_fn   = taucs_supernodal_solve_llt;
    } else if (symb_flag) {
      L = taucs_ccs_factor_llt_symbolic(PAPT);
      taucs_ccs_factor_llt_numeric(PAPT,L); /* should check error code */
      precond_args = L;
      precond_fn   = taucs_supernodal_solve_llt;
    } else {
      L = taucs_ccs_factor_llt(PAPT,0.0,0);
      precond_args = L;
      precond_fn   = taucs_ccs_solve_llt;
    }

    wtime_factor = taucs_wtime() - wtime_factor;
    ctime_factor = taucs_ctime() - ctime_factor;
    taucs_printf("\tFactor time   = % 10.3f seconds  ",wtime_factor);
    taucs_printf("(%.3f cpu time)\n",ctime_factor);
  }

  if (!L && !ooc_flag /* no L in ooc */) {
    taucs_printf("\tFactorization Failed\n");
    exit(1);
  }

  /*taucs_ccs_write_ijv(PAPT,"A.ijv",1);*/ /* 1 = complete the upper part */
  /*taucs_ccs_write_ijv(L,"L.ijv",0);*/

  /***********************************************************/
  /* solve                                                   */
  /***********************************************************/

  taucs_vec_permute(A->n,B,PB,perm);

  wtime_solve = taucs_wtime();

  precond_fn(precond_args,PX,PB); /* direct solver */

  wtime_solve = taucs_wtime() - wtime_solve;
  taucs_printf("\tSolve time    = % 10.3f seconds\n",wtime_solve);

  taucs_vec_ipermute(A->n,PX,NX,perm);

  /***********************************************************/
  /* delete out-of-core matrices                             */
  /***********************************************************/

  if (ooc_flag) {
    taucs_io_delete(oocL);
    /*taucs_io_close(oocL);*/
  }

  /***********************************************************/
  /* Compute norm of forward error                           */
  /***********************************************************/

  NormErr = 0.0;
  for(i=0; i<N; i++) NormErr = max(NormErr,fabs((NX[i]-X[i])/X[i]));
  taucs_printf("main: max relative error = %1.6e \n",NormErr); 

  /***********************************************************/
  /* Exit                                                    */
  /***********************************************************/

  taucs_printf("main: done\n");
} 

