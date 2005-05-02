/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo and Doron Chen                   */
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/*
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/timeb.h>
*/

#ifndef WIN32
#include <unistd.h>
#endif

#include <sys/types.h>                                                 
#include <assert.h>
#include <taucs.h>


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

  double ssq, scale, absvi;
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

int        M,N;					/* Size of matrix */

/*********************************************************/
/*                                                       */
/*********************************************************/

void usage(int argc, char* argv[])
{
  printf("%s OPTIONS\n",argv[0]);
  printf("  -help or -h or no arguments at all\n");
  printf("        this help\n");
  printf("  MATRIX OPTIONS:\n");
  printf("        -hb  filename    (matrix from Harwell-Boeing file)\n");
  printf("        -ijv ijvfilename (matrix from file)\n");
  printf("        -mtx mtxfilename (matrix from file)\n");
  printf("        -ccs ccsfilename (matrix from file)\n");
  printf("        -discont     X Y Z jump (X-by-Y-by-Z poisson with discontinuous coeff.)\n");
  printf("        -rrn         X Y Z dp rmin (X-by-Y-by-Z random resistor network)\n");
  printf("        -mesh3d      X Y Z (X-by-Y-by-Z poisson problem)\n");
  printf("        -mesh2d          n (n-by-n poisson problem)\n");
  printf("        -mat_type [anisotropic_y/anisotropic_x/dirichlet/neumann]\n");
  printf("                  mat_type only applicable to -mesh2d\n");
  printf("                  (3D meshes are always Neumann BC)\n");
  printf("                  anisotropic problems are Neumann\n");
  printf("  RIGHT-HAND SIDE OPTIONS:\n");
  printf("        -n+rhs filename (order n followed by rhs, ascii)\n");
  printf("  GLOBAL PRECOND OPTIONS:\n");
  printf("  -icc\n");
  printf("        uses incomplete Cholesky\n");
  printf("  -vaidya\n");
  printf("        uses Vaidya\n");
  printf("  -stretch\n");
  printf("        uses Vaidya's preconditioner - starting with a low stretch tree\n");
  printf("  -recvaidya\n");
  printf("        uses recursive Vaidya\n");
  printf("  -sg argument\n");
  printf("        uses support-graph preconditioning\n");
  printf("  ICC PRECOND OPTIONS:\n");
  printf("  -modified\n");
  printf("        modified ICC; maintains rowsums; only for -icc\n");
  printf("  -droptol droptolerance\n");
  printf("        drop tolerance for icc and micc\n");
  printf("  -ordering [tree/metis/natural/genmmd/md/mmd/amd]\n");
  printf("        ordering to use\n");
  printf("  VAIDYA PRECOND OPTIONS:\n");
  printf("  -nomf\n");
  printf("        use column Cholesky instead of multifrontal\n");
  printf("  -subgraphs t\n");
  printf("        desired number of subgraphs for Vaidya\n");
  printf("  -seed s\n");
  printf("        seed for random-number generator; affects decomposition to subtrees\n");
  printf("  -rnd x\n");
  printf("        uses x INSTEAD of a random number to decompose to subtrees\n");
  printf("        we use this to repeat experiments deterministically\n");
  printf("  RECURSIVE VAIDYA PRECOND OPTIONS:\n");
  printf("  -c c -epsilon e -nsmall nsmall -maxlevels l -convratio r -innerits m\n");
  printf("        splits tree into about k=c*n^(1/(1+e)) subgraphs,\n");
  printf("        matrices smaller than nsmall are factored directly\n");
  printf("        preconditioner has at most l levels\n");
  printf("        inner solves reduce their residual by a factor of r\n");
  printf("        using at most m iterations\n");
  printf("  SUPPORT GRAPH PRECOND OPTIONS:\n");
  printf("  argument\n");
  printf("        use regular:XX:N\n");
  printf("        where XX is GM (Gremban-Miller), CT (Chen-Toledo) or VA (Vaidya)\n");
  printf("        and N is an integer, the number of subdomains in each partition\n");
  printf("  -subgraphs t\n");
  printf("        desired number of subgraphs for Vaidya\n");
  printf("  -seed s\n");
  printf("        seed for random-number generator; affects decomposition to subtrees\n");
  printf("  -rnd x\n");
  printf("        uses x INSTEAD of a random number to decompose to subtrees\n");
  printf("        we use this to repeat experiments deterministically\n");
  printf("  OTHER OPTIONS:\n");
  printf("  -log filename\n");
  printf("        write log into filename; can be stdout and stderr\n");
  exit(0);
}

/*********************************************************/
/*                                                       */
/*********************************************************/

int main(int argc, char* argv[])
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
  taucs_ccs_matrix*  A_orig = NULL;
  taucs_ccs_matrix*  PAPT;
  taucs_ccs_matrix*  L;
  taucs_ccs_matrix*  V;
  taucs_ccs_matrix*  PVPT;

  double*      X = NULL;
  double*      B = NULL;
  double*      PX;
  double*      PB;
  double*      NX;
  double*      X_orig;
  double*      B_orig;

  char*        ordering = "metis";
  char*        sg_command = "qqq";
  char*        mat_type = "neumann";
  int*         perm;
  int*         invperm;
  int          multifrontal_flag = 1;
  int          modified_flag = 0;
  int          icc_flag = 0;
  int          sg_flag = 0;
  int          vaidya_flag = 0;
  int          stretch_flag = 0;
  int          vaidya_icc_flag = 0;
  int          contx_flag = 0;
  int          recvaidya_flag = 0;
  int          mesh2d_flag = 0,mesh2d_negative_flag=0,mesh2d_size = 0;
  int          trick_flag = 0;
  double       subgraphs = 1.0;
  
  double droptol = 0.0;
  int rnd,force_rnd=0;
  int seed = 123;

  double C = 0.25;    /* constant factor reduction in rec vaidya */
  double epsilon=0.2; 
  int    nsmall=10000;
  int    maxlevels=2;
  int    innerits=2;
  double innerconv=0.01;

  double restol = 1e-8;
  int    maxits = 10000;

  int             (*precond_fn)(void*,void* x,void* b);
  void*           precond_args;


  /***********************************************************/
  /* Read arguments: log file, matrix, memory size, ooc name */
  /***********************************************************/

  if ((argc == 1) ||((argc == 2) && !strncmp(argv[1],"-h",2)))
    usage(argc,argv);

  A = NULL;
  for (i=1; i<argc; i++) {
    if (!strcmp(argv[i],"-nomf"))      multifrontal_flag = 0;
    if (!strcmp(argv[i],"-trick"))     trick_flag = 1;
    if (!strcmp(argv[i],"-icc"))       icc_flag = 1;
    if (!strcmp(argv[i],"-vaidya"))    vaidya_flag = 1;
    if (!strcmp(argv[i],"-stretch"))    stretch_flag = 1;
    if (!strcmp(argv[i],"-cont-x"))    contx_flag = 1;
    if (!strcmp(argv[i],"-vaidya-icc"))    vaidya_icc_flag = 1;
    if (!strcmp(argv[i],"-recvaidya")) recvaidya_flag = 1;
    if (!strcmp(argv[i],"-modified")) {modified_flag = 1;taucs_printf("Factoring will be Modified ICC\n");}
    if (!strcmp(argv[i],"-log") && i <= argc-1 ) {
      i++;
      taucs_logfile(argv[i]);
    }
    if (!strcmp(argv[i],"-ordering") && i <= argc-1 ) {
      i++;
      ordering = argv[i];
    }
    if (!strcmp(argv[i],"-sg") && i <= argc-2 ) {
      sg_flag = 1;
      i++;
      sg_command = argv[i];
    }
    if (!strcmp(argv[i],"-mat_type") && i <= argc-1 ) {
      i++;
      mat_type = argv[i];
    }
    if (!strcmp(argv[i],"-droptol") && i <= argc-1 ) {
      i++;
      if (sscanf(argv[i],"%lg",&droptol) != 1) {
	taucs_printf("main: -droptol must be followed by a real value\n");
	exit(1);
      }
    }

    if (!strcmp(argv[i],"-restol") && i <= argc-1 ) {
      i++;
      if (sscanf(argv[i],"%lg",&restol) != 1) {
	taucs_printf("main: -restol must be followed by a real value\n");
	exit(1);
      }
    }
    if (!strcmp(argv[i],"-maxits") && i <= argc-1 ) {
      i++;
      if (sscanf(argv[i],"%d",&maxits) != 1) {
	taucs_printf("main: -maxits must be followed by an integer value\n");
	exit(1);
      }
    }

    if (!strcmp(argv[i],"-c") && i <= argc-1 ) {
      i++;
      if (sscanf(argv[i],"%lg",&C) != 1) {
	taucs_printf("main: -c must be followed by a real value\n");
	exit(1);
      }
    }
    if (!strcmp(argv[i],"-epsilon") && i <= argc-1 ) {
      i++;
      if (sscanf(argv[i],"%lg",&epsilon) != 1) {
	taucs_printf("main: -epsilon must be followed by a real value\n");
	exit(1);
      }
    }
    if (!strcmp(argv[i],"-nsmall") && i <= argc-1 ) {
      i++;
      if (sscanf(argv[i],"%d",&nsmall) != 1) {
	taucs_printf("main: -nsmall must be followed by an integer\n");
	exit(1);
      }
    }

    if (!strcmp(argv[i],"-maxlevels") && i <= argc-1 ) {
      i++;
      if (sscanf(argv[i],"%d",&maxlevels) != 1) {
	taucs_printf("main: -maxlevels must be followed by an integer\n");
	exit(1);
      }
    }
    if (!strcmp(argv[i],"-innerits") && i <= argc-1 ) {
      i++;
      if (sscanf(argv[i],"%d",&innerits) != 1) {
	taucs_printf("main: -innerits must be followed by an integer\n");
	exit(1);
      }
    }
    if (!strcmp(argv[i],"-innerconv") && i <= argc-1 ) {
      i++;
      if (sscanf(argv[i],"%lg",&innerconv) != 1) {
	taucs_printf("main: -innerconv must be followed by a real value\n");
	exit(1);
      }
    }

    if (!strcmp(argv[i],"-seed") && i <= argc-1 ) {
      i++;
      if (sscanf(argv[i],"%d",&seed) != 1) {
	taucs_printf("main: -seed must be followed by an integer value\n");
	exit(1);
      }
    }
    if (!strcmp(argv[i],"-rnd") && i <= argc-1 ) {
      i++;
      if (sscanf(argv[i],"%d",&force_rnd) != 1) {
	taucs_printf("main: -rnd must be followed by an integer value\n");
	exit(1);
      }
    }
    if (!strcmp(argv[i],"-subgraphs") && i <= argc-1 ) {
      i++;
      if (sscanf(argv[i],"%lg",&subgraphs) != 1) {
	taucs_printf("main: -subgraphs must be followed by a real value\n");
	exit(1);
      }
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
      A = taucs_dccs_create(nrows,ncols,nnz);
      if (type[1] == 's' || type[1] == 'S')
	A->flags = TAUCS_SYMMETRIC | TAUCS_LOWER | TAUCS_DOUBLE;
      dreadhb_(fname,&nrows,&ncols,&nnz,
	       A->colptr,A->rowind,((double*)A->values.d));
      /* make indices 0-based */
      for (j=0; j<=ncols; j++) ((A->colptr)[j])--;
      for (j=0; j<nnz;    j++) ((A->rowind)[j])--;
      taucs_printf("main: done reading\n");
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

    if (!strcmp(argv[i],"-mesh2d_negative") && i <= argc-1) {
      mesh2d_negative_flag = 1;
      taucs_printf("A is a mesh2d_negative\n");
      i++;
      if (sscanf(argv[i],"%d",&mesh2d_size) != 1) {
	taucs_printf("mesh size (n, where the mesh is n-by-n) must follow -mesh2d_negative argument\n");
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
    if (!strcmp(argv[i],"-rrn") && i <= argc-5) {
      int X,Y,Z;
      double dp, rmin;

      if (sscanf(argv[i+1],"%d",&X) != 1 
	  || sscanf(argv[i+2],"%d",&Y) != 1 
	  || sscanf(argv[i+3],"%d",&Z) != 1
	  || sscanf(argv[i+4],"%lg",&dp) != 1
	  || sscanf(argv[i+5],"%lg",&rmin) != 1) {
	taucs_printf("mesh size (X Y Z dp rmin must follow -rrn argument\n");
	exit(1);
      }
      i += 5;
      taucs_printf("main: creating mesh %d\n",Y);
      A = taucs_ccs_generate_rrn(X,Y,Z,dp,rmin);

      taucs_printf("A is a random resistor network %d-by-%d-by-%d, drop prob=%.2e, rmin=%.2e\n",
		   X,Y,Z,dp,rmin);
    }
    if (!strcmp(argv[i],"-discont") && i <= argc-4) {
      int x,y,z;
      double jump;

      if (sscanf(argv[i+1],"%d",&x) != 1 
	  || sscanf(argv[i+2],"%d",&y) != 1 
	  || sscanf(argv[i+3],"%d",&z) != 1
	  || sscanf(argv[i+4],"%lg",&jump) != 1) {
	taucs_printf("mesh size (X Y Z jump must follow -discont argument\n");
	exit(1);
      }
      i += 4;
      taucs_printf("main: creating mesh\n",y);
      A = taucs_ccs_generate_discontinuous(x,y,z,jump);

      taucs_printf("A is a %d-by-%d-by-%d discontinuous Poisson problem, jump=%.2e\n",
		   x,y,z,jump);

      if (contx_flag) {
	X = taucs_vec_generate_continuous(x,y,z,"default");
	taucs_printf("iter: using continuous solution\n");
      }
    }
  }

  srand(seed);
  rnd = rand();
  if (force_rnd)
    rnd = force_rnd;

  taucs_printf("Chosen Ordering: %s\n",ordering);

  if (mesh2d_flag)
    {
      taucs_printf("Matrix type is %s\n",mat_type);
      taucs_printf("Grid Size is %d\n",mesh2d_size);
      A = taucs_ccs_generate_mesh2d(mesh2d_size,mat_type);
    }

  if (mesh2d_negative_flag)
    {
      taucs_printf("Grid Size is %d\n",mesh2d_size);
      A = taucs_ccs_generate_mesh2d_negative(mesh2d_size);
    }

  /* A = symccs_read_ijv("/tmp/A.ijv"); */

  if (!A) {
    taucs_printf("matrix argument not given or matrix file not found\n");
    usage(argc,argv);
  }
  
  N = A->n;
  if (!X) {
    X=(double*)malloc(N*sizeof(double));
    for(i=0; i<N; i++)  X[i]=(double)rand()/RAND_MAX;
  } else
    taucs_printf("iter: not using a random X, already allocated\n");

  if (!B) {
    B=(double*)malloc(N*sizeof(double));
    taucs_ccs_times_vec(A,X,B);
  } else {
    for(i=0; i<N; i++)  X[i]= taucs_get_nan();/*omer 0.0 / 0.0*/;
  }

  if ((sg_flag)&&(trick_flag==0))
    {
      for(i=0;i<A->n;i++)
	{
	  for (j=A->colptr[i];j<A->colptr[i+1];j++)
	    if ((A->rowind[j] != i)&&(((double*)A->values.d)[j]>0))
	      {
		trick_flag = 1;
		break;
	      }
	  if (trick_flag) break;
	}
    }

  if (trick_flag) {
    int n;
    int *tmp;

    taucs_printf("Warning: augmenting to ensure nonpositive offdiagonals\n");

    B_orig = B;
    X_orig = X;
    A_orig = A;
    
    X = (double *)malloc(2*N*sizeof(double));
    for(i=0;i<N;i++) {
      X[i] = X_orig[i];
      X[i+N] = -X_orig[i];
    }

    B = (double *)malloc(2*N*sizeof(double));
    for(i=0;i<N;i++) {
      B[i] = B_orig[i];
      B[i+N] = -B_orig[i];
    }

#if 1
    {
    taucs_ccs_matrix* A_tmp;
    n=A->n;
    A_tmp = taucs_dccs_create(2*n,2*n,2*(A->colptr[n]));
    A_tmp->flags |= TAUCS_SYMMETRIC | TAUCS_LOWER;

    tmp = (int *)calloc((2*n+1),sizeof(int));

    for(i=0;i<n;i++)
	{
	  for(j=A->colptr[i];j<A->colptr[i+1];j++)
	    {
	      if ((i == A->rowind[j])||(((double*)A->values.d)[j] < 0))
		{
		  tmp[i]++;
		  tmp[i+n]++;
		}
	      else
		{
		  /* printf("WWW\n");exit(345); */
		  tmp[i]++;
		  tmp[A->rowind[j]]++;
		}
	    }
	}
      A_tmp->colptr[0]=0;
      for(i=0;i<2*n;i++)
	A_tmp->colptr[i+1] = A_tmp->colptr[i] + tmp[i];
      for(i=0;i<2*n;i++)
	tmp[i] = A_tmp->colptr[i];

      /* for(i=0;i<n;i++) */
	/* printf("ZZZ %d %d\n",tmp[i],tmp[i+n]); */
      /* exit(345); */

      for(i=0;i<n;i++)
	{
	  for(j=A->colptr[i];j<A->colptr[i+1];j++)
	    {
	      if ((i == A->rowind[j])||(((double*)A->values.d)[j] < 0))
		{
 		  A_tmp->rowind[tmp[i]]=A->rowind[j];
		  ((double*)A_tmp->values.d)[tmp[i]++]=((double*)A->values.d)[j];
		  A_tmp->rowind[tmp[i+n]]=A->rowind[j]+n;
		  ((double*)A_tmp->values.d)[tmp[i+n]++]=((double*)A->values.d)[j];
		}
	      else
		{
		  /* printf("WWW\n");exit(345); */
 		  A_tmp->rowind[tmp[i]]=A->rowind[j]+n;
		  ((double*)A_tmp->values.d)[tmp[i]++]=-((double*)A->values.d)[j];
		  A_tmp->rowind[tmp[A->rowind[j]]]=i+n;
		  ((double*)A_tmp->values.d)[tmp[A->rowind[j]]++]=-((double*)A->values.d)[j];
		}
	    }
	}
      /*
      taucs_ccs_write_ijv(A,"A.ijv");
      taucs_ccs_write_ijv(A_tmp,"Aaug.ijv");
      */

      free(tmp);
      taucs_ccs_free(A);
      A = A_tmp;

    }
#else
    A = taucs_ccs_augment_nonpositive_offdiagonals(A_orig);
#endif
    }
  
  N = M = A->n;
  
  /***********************************************************/
  /* Create exact solution, compute right-hand-side          */
  /***********************************************************/

  NX=(double*)malloc(N*sizeof(double));
  PX=(double*)malloc(N*sizeof(double));
  PB=(double*)malloc(N*sizeof(double));
  
  /***********************************************************/
  /* Compute column ordering                                 */
  /***********************************************************/

  /***********************************************************/
  /* factor                                                  */
  /***********************************************************/

  if (icc_flag) {
    int n;
    double unit,curr;

    n = A->n;
    unit = (n-1.)+n;

    taucs_printf("main: using incomplete Cholesky droptol=%lg modified=%d\n",
	       droptol,modified_flag);

    wtime_order = taucs_wtime();
    taucs_ccs_order(A,&perm,&invperm,ordering);
    wtime_order = taucs_wtime() - wtime_order;
    taucs_printf("\tOrdering time = % 10.3f seconds\n",wtime_order);

    wtime_permute = taucs_wtime();
    PAPT = taucs_ccs_permute_symmetrically(A,perm,invperm);
    wtime_permute = taucs_wtime() - wtime_permute;
    taucs_printf("\tPermute time  = % 10.3f seconds\n",wtime_permute);

    wtime_factor = taucs_wtime();
    ctime_factor = taucs_ctime();
    L = taucs_ccs_factor_llt(PAPT,droptol,modified_flag);
    /* taucs_ccs_write_ijv(PAPT,"PAPT.ijv"); */
    /* taucs_ccs_write_ijv(L,"L.ijv"); */
    wtime_factor = taucs_wtime() - wtime_factor;
    ctime_factor = taucs_ctime() - ctime_factor;
    taucs_printf("\tFactor time   = % 10.3f seconds  ",wtime_factor);
    taucs_printf("(%.3f cpu time)\n",ctime_factor);

    curr = L->colptr[n];
    taucs_printf("droptol = %1.20lf curr = %lf curr_ratio = %lf\n",droptol,curr,curr/unit);

    precond_args = L;
    precond_fn   = taucs_ccs_solve_llt;

  } else if (vaidya_flag) {
    int n;
    double unit,curr;
    /*
    int taucs_supernodal_solve_llt(void*,
				   double*, double*);
    */


    n = A->n;
    unit = (n-1.)+n;

    taucs_printf("main: using Vaidya # subgraphs requested=%lf\n",subgraphs);
    taucs_printf("Chosen rnd: %d\n",rnd);

    wtime_precond_create = taucs_wtime();
#if 1
    V = taucs_amwb_preconditioner_create(A,
					 rnd /* random seed */,
					 subgraphs,
					 stretch_flag);
#else
    V = taucs_amst_preconditioner_create(A,
					 rnd /* random seed */,
					 subgraphs);
#endif
    /* taucs_ccs_write_ijv(A,"A.ijv"); */
    /* taucs_ccs_write_ijv(V,"V.ijv"); */
    wtime_precond_create = taucs_wtime() - wtime_precond_create;
    taucs_printf("\tCreation time = % 10.3f seconds\n",wtime_precond_create);

    wtime_order = taucs_wtime();
    taucs_ccs_order(V,&perm,&invperm,ordering);
    wtime_order = taucs_wtime() - wtime_order;
    taucs_printf("\tOrdering time = % 10.3f seconds\n",wtime_order);

    wtime_permute = taucs_wtime();
    PAPT = taucs_ccs_permute_symmetrically(A,perm,invperm);
    PVPT = taucs_ccs_permute_symmetrically(V,perm,invperm);
    wtime_permute = taucs_wtime() - wtime_permute;
    taucs_printf("\tPermute time  = % 10.3f seconds\n",wtime_permute);

    taucs_ccs_free(V);

    wtime_factor = taucs_wtime();

    if (multifrontal_flag) {
      void* snL;
      snL = taucs_ccs_factor_llt_mf(PVPT);
#if 1
      L = taucs_supernodal_factor_to_ccs(snL);
      taucs_supernodal_factor_free(snL);
      precond_args = L;
      precond_fn   = taucs_ccs_solve_llt;
#else
      L = snL;
      precond_args = L;
      precond_fn   = taucs_supernodal_solve_llt;
#endif
    } else {
      L = taucs_ccs_factor_llt(PVPT,0.0,0);
      precond_args = L;
      precond_fn   = taucs_ccs_solve_llt;
      curr = L->colptr[n];
      taucs_printf("subgraphs = %6.15lf curr = %lf curr_ratio = %lf\n",subgraphs,curr,curr/unit);fflush(stdout);
    }

    wtime_factor = taucs_wtime() - wtime_factor;
    taucs_printf("\tFactor time   = % 10.3f seconds\n",wtime_factor);

    taucs_ccs_free(PVPT);
    ((taucs_ccs_matrix *)precond_args)->flags |= TAUCS_DOUBLE;

  } else if (recvaidya_flag) {

    wtime_precond_create = taucs_wtime();
    precond_args = taucs_recursive_amwb_preconditioner_create(A,
							      C,
							      epsilon,
							      nsmall,
							      maxlevels,
							      innerits,
							      innerconv,
							      &perm,&invperm);
    precond_fn = taucs_recursive_amwb_preconditioner_solve;
    wtime_precond_create = taucs_wtime() - wtime_precond_create;
    taucs_printf("\tCreation time = % 10.3f seconds\n",wtime_precond_create);

    wtime_permute = taucs_wtime();
    PAPT = taucs_ccs_permute_symmetrically(A,perm,invperm);
    wtime_permute = taucs_wtime() - wtime_permute;
    taucs_printf("\tOrdering time = % 10.3f seconds\n",wtime_permute);
  } 
#ifdef SIVAN_5_DEC_2002_XXX
  else if (sg_flag) {
    int
      taucs_sg_preconditioner_solve(void*  vP,
				    double* Z, 
				    double* R);


    wtime_precond_create = taucs_wtime();
    /* precond_args = taucs_multilevel_preconditioner_create(A, */
							  /* subgraphs, */
							  /* &perm,&invperm); */

    precond_args = taucs_sg_preconditioner_create(A,&perm,&invperm,
						  ordering,sg_command);
    if (!precond_args)
      {
	printf("Creating the preconditioner failed. Exiting.\n");
	exit(345);
      }

    precond_fn = taucs_sg_preconditioner_solve;
    wtime_precond_create = taucs_wtime() - wtime_precond_create;
    taucs_printf("\tCreation time = % 10.3f seconds\n",wtime_precond_create);

    wtime_permute = taucs_wtime();
    PAPT = taucs_ccs_permute_symmetrically(A,perm,invperm);
    wtime_permute = taucs_wtime() - wtime_permute;
    taucs_printf("\tOrdering time = % 10.3f seconds\n",wtime_permute);
  }
#endif
  else {
    precond_args = NULL;
    precond_fn   = NULL;
    taucs_ccs_order(A,&perm,&invperm,"identity");
    PAPT = A;
  }
  /*taucs_ccs_write_ijv(PAPT,"A.ijv",1);*/ /* 1 = complete the upper part */
  /*taucs_ccs_write_ijv(L,"L.ijv",0);*/

  /***********************************************************/
  /* solve                                                   */
  /***********************************************************/



  for (i=0; i<A->n; i++) PB[i] = B[perm[i]];
  for (i=0; i<A->n; i++) PX[i] = 0.0; /* initial guess */

  wtime_solve = taucs_wtime();
#if 1

  if (vaidya_icc_flag && vaidya_flag) {
    taucs_ccs_matrix* L;
    double* R;
    double init_residual_reduction;
    
    wtime_factor = taucs_wtime();
    ctime_factor = taucs_ctime();

    L = taucs_ccs_factor_llt(PAPT,1.0,0); /* unmodified ICC(0) */
    wtime_factor = taucs_wtime() - wtime_factor;
    ctime_factor = taucs_ctime() - ctime_factor;
    taucs_printf("\tFactor time   = % 10.3f seconds  ",wtime_factor);
    taucs_printf("(%.3f cpu time)\n",ctime_factor);

    taucs_conjugate_gradients (PAPT,
			       taucs_ccs_solve_llt, L,
			       PX, PB,
			       25, 1e-15);

    R = (double*) malloc(A->n * sizeof(double));
    taucs_ccs_times_vec(PAPT,PX,R);
    for (i=0; i<A->n; i++) R[i] = PB[i] - R[i];
    init_residual_reduction = twonorm(A->n,R)/twonorm(A->n,B);
    taucs_printf("reduction by a factor of %.2e so far\n",init_residual_reduction);
    free(R);

    taucs_conjugate_gradients (PAPT,
			       precond_fn, precond_args,
			       PX, PB,
			       10000, 1e-15 / init_residual_reduction);
  } else
    {
      taucs_conjugate_gradients (PAPT,
				 precond_fn, precond_args,
				 PX, PB,
				 maxits, restol);
    }
#else
  precond_fn(precond_args,PX,PB); /* direct solver */
#endif

  wtime_solve = taucs_wtime() - wtime_solve;
  taucs_printf("\tSolve time    = % 10.3f seconds\n",wtime_solve);

  for (i=0; i<A->n; i++) NX[i] = PX[invperm[i]];

  /***********************************************************/
  /* delete out-of-core matrices                             */
  /***********************************************************/

  /***********************************************************/
  /* extract solution from augmented system                  */
  /***********************************************************/

  if (0 && A_orig) {
    /* 
       the first half of the elements of NX contains the
       solution of the original system
    */

    taucs_ccs_free(A);
    A = A_orig;

    free(X);
    X = X_orig;

    free(B);
    B = B_orig;

    N = M = A->n;
  }

  /***********************************************************/
  /* Compute norm of forward error                           */
  /***********************************************************/

  NormErr = 0.0;
  for(i=0; i<N; i++) NormErr = max(NormErr,fabs((NX[i]-X[i])/X[i]));
  taucs_printf("main: max relative error = %1.6e \n",NormErr); 

  for(i=0; i<N; i++) B[i] = NX[i]-X[i];
  taucs_printf("main: ||computed X - X||/||X||=%.2e\n",twonorm(N,B)/twonorm(N,X)); 
  
  /* for(i=0; i<N; i++) */
    /* printf("QQQ %i %lf\n",i,NX[i]); */

  /***********************************************************/
  /* Exit                                                    */
  /***********************************************************/
  
  taucs_printf("main: done\n");
	return 0;
} 

