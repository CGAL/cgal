/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/

#if defined(TAUCS_CORE_CILK) && defined(TAUCS_CILK)
#pragma lang -C
#endif

/*** taucs_ccs_factor.c ***/
int taucs_getopt_boolean(char* cmd, void* args[], char* name, int*    x);
int taucs_getopt_double (char* cmd, void* args[], char* name, double* x);
int taucs_getopt_pointer(char* cmd, void* args[], char* name, void**  x);
int taucs_getopt_string (char* cmd, void* args[], char* name, char**  x);

int taucs_linsolve(taucs_ccs_matrix* A, 
		   void**            F,
		   int               nrhs,
		   void*             X,
		   void*             B,
		   char*             options[],
		   void*             opt_arg[]);

/*** taucs_ccs_base.c ***/

extern taucs_datatype taucs_dtl(zero_const);
extern taucs_datatype taucs_dtl(one_const);

#ifndef TAUCS_C99_COMPLEX
taucs_datatype taucs_dtl(complex_create_fn)(taucs_real_datatype re, 
					    taucs_real_datatype im);
#endif
taucs_datatype taucs_dtl(add_fn)(taucs_datatype a, taucs_datatype b);
taucs_datatype taucs_dtl(sub_fn)(taucs_datatype a, taucs_datatype b);
taucs_datatype taucs_dtl(mul_fn)(taucs_datatype a, taucs_datatype b);
taucs_datatype taucs_dtl(div_fn)(taucs_datatype a, taucs_datatype b);
taucs_datatype taucs_dtl(neg_fn)(taucs_datatype a);
taucs_datatype taucs_dtl(sqrt_fn)(taucs_datatype a);
taucs_datatype taucs_dtl(conj_fn)(taucs_datatype a);
double         taucs_dtl(abs_fn)(taucs_datatype a);

/*** taucs_ccs_base.c ***/

taucs_ccs_matrix* taucs_dtl(ccs_create)          (int m, int n, int nnz);
taucs_ccs_matrix* taucs_ccs_create               (int m, int n, int nnz, int flags);
void              taucs_dtl(ccs_free)            (taucs_ccs_matrix* matrix);
void              taucs_ccs_free                 (taucs_ccs_matrix* matrix);

/*** taucs_ccs_ops.c ***/

void              taucs_dtl(ccs_split)           (taucs_ccs_matrix* A, 
						  taucs_ccs_matrix** L, 
						  taucs_ccs_matrix** R, 
						  int p);
void              taucs_ccs_split                (taucs_ccs_matrix* A, 
						  taucs_ccs_matrix** L, 
						  taucs_ccs_matrix** R, 
						  int p);

taucs_ccs_matrix* taucs_dtl(ccs_permute_symmetrically)(taucs_ccs_matrix* A, 
						       int* perm, int* invperm);
taucs_ccs_matrix*     taucs_ccs_permute_symmetrically (taucs_ccs_matrix* A, 
						       int* perm, int* invperm);

void              taucs_dtl(ccs_times_vec)       (taucs_ccs_matrix* m, 
						  taucs_datatype* X,
						  taucs_datatype* B);
void                        taucs_ccs_times_vec  (taucs_ccs_matrix* m, 
						  void* X,
						  void* B);

/* matrix-vector with double-precision accumulator for iterative refinement */
void              taucs_sccs_times_vec_dacc      (taucs_ccs_matrix* m, 
						  taucs_single* X,
						  taucs_single* B);

taucs_ccs_matrix* taucs_dtl(ccs_augment_nonpositive_offdiagonals)(taucs_ccs_matrix* A);
taucs_ccs_matrix*     taucs_ccs_augment_nonpositive_offdiagonals (taucs_ccs_matrix* A);

/*** taucs_ccs_io.c ***/

int               taucs_dtl(ccs_write_ijv)       (taucs_ccs_matrix* matrix, 
						  char* filename);
int                    taucs_ccs_write_ijv       (taucs_ccs_matrix* matrix, 
						  char* filename);
taucs_ccs_matrix* taucs_dtl(ccs_read_ijv)        (char* filename,int flags);
taucs_ccs_matrix*     taucs_ccs_read_ijv         (char* filename,int flags);
taucs_ccs_matrix* taucs_dtl(ccs_read_mtx)        (char* filename,int flags);
taucs_ccs_matrix*     taucs_ccs_read_mtx         (char* filename,int flags);
taucs_ccs_matrix* taucs_dtl(ccs_read_ccs)        (char* filename,int flags);
taucs_ccs_matrix*     taucs_ccs_read_ccs         (char* filename,int flags);
taucs_ccs_matrix* taucs_ccs_read_binary          (char* filename);
void*             taucs_vec_read_binary (int n, int flags,          char* filename);
int               taucs_vec_write_binary(int n, int flags, void* v, char* filename);
taucs_ccs_matrix* taucs_ccs_read_hb              (char* filename,int flags);

/*** taucs_ccs_order.c ***/

void              taucs_ccs_order                (taucs_ccs_matrix* matrix, 
						  int** perm, int** invperm,
						  char* which);

/*** taucs_ccs_factor_llt.c ***/

taucs_ccs_matrix* taucs_dtl(ccs_factor_llt)      (taucs_ccs_matrix* A,
						  double droptol, int modified);
taucs_ccs_matrix*     taucs_ccs_factor_llt       (taucs_ccs_matrix* A,
						  double droptol, int modified);
taucs_ccs_matrix* taucs_ccs_factor_llt_partial   (taucs_ccs_matrix* A, 
						  int p);
taucs_ccs_matrix* taucs_dtl(ccs_factor_llt_partial)(taucs_ccs_matrix* A, 
						    int p);

taucs_ccs_matrix* taucs_dtl(ccs_factor_ldlt)     (taucs_ccs_matrix* A);
taucs_ccs_matrix*     taucs_ccs_factor_ldlt      (taucs_ccs_matrix* A);

taucs_ccs_matrix* taucs_dtl(ccs_factor_xxt)      (taucs_ccs_matrix* A);

/*** taucs_ccs_solve_llt.c ***/

int               taucs_ccs_solve_llt            (void* L, void* x, void* b);
int               taucs_dtl(ccs_solve_llt)       (void* L, taucs_datatype* x, taucs_datatype* b);
int               taucs_ccs_solve_ldlt           (void* L, void* x, void* b);
int               taucs_dtl(ccs_solve_ldlt)      (void* L, taucs_datatype* x, taucs_datatype* b);
int               taucs_ccs_solve_schur          (taucs_ccs_matrix* L,
						  taucs_ccs_matrix* schur_comp,
						  int    (*schur_precond_fn)(void*,void* x,void* b),
						  void*  schur_precond_args,
						  int    maxits,
						  double convratio,
						  void* x, void* b);
int               taucs_dtl(ccs_solve_schur)     (taucs_ccs_matrix* L,
						  taucs_ccs_matrix* schur_comp,
						  int    (*schur_precond_fn)(void*,void* x,void* b),
						  void*  schur_precond_args,
						  int    maxits,
						  double convratio,
						  taucs_datatype* x, taucs_datatype* b);

/***  ***/

taucs_ccs_matrix* taucs_ccs_factor_xxt           (taucs_ccs_matrix* A);
int               taucs_ccs_solve_xxt            (void* X, double* x, double* b);

taucs_ccs_matrix* taucs_ccs_generate_mesh2d      (int n,char *which);
taucs_ccs_matrix* taucs_ccs_generate_mesh2d_negative(int n);
taucs_ccs_matrix* taucs_ccs_generate_mesh3d      (int X, int Y, int Z);
taucs_ccs_matrix* taucs_ccs_generate_dense       (int m,int n, int flags);
taucs_ccs_matrix* taucs_ccs_generate_rrn         (int X, int Y, int Z, 
						  double drop_probability, 
						  double rmin);
taucs_ccs_matrix* taucs_ccs_generate_discontinuous(int X, int Y, int Z, 
						   double jump);
double* taucs_vec_generate_continuous            (int X, int Y, int Z, char* which);

int taucs_conjugate_gradients                    (taucs_ccs_matrix*  A,
						  int               (*precond_fn)(void*,void* x,void* b),
						  void*             precond_args,
						  void*             X,
						  void*             B,
						  int               itermax,
						  double            convergetol);

int taucs_minres                                 (taucs_ccs_matrix*  A,
						  int               (*precond_fn)(void*,void* x,void* b),
						  void*             precond_args,
						  void*             X,
						  void*             B,
						  int               itermax,
						  double            convergetol);

int taucs_sg_preconditioner_solve                (void*   P,
						  double* z, 
						  double* r);

void *taucs_sg_preconditioner_create             (taucs_ccs_matrix *A,
						  int **perm,
						  int **invperm,
						  char* ordering,
						  char *specification);
void taucs_sg_preconditioner_free                (void* P);

taucs_ccs_matrix*
taucs_amwb_preconditioner_create                 (taucs_ccs_matrix *symccs_mtxA, 
						  int rnd,
						  double subgraphs,
						  int stretch_flag);

void* 
taucs_recursive_amwb_preconditioner_create       (taucs_ccs_matrix* A, 
						  double c, 
						  double epsilon, 
						  int nsmall,
						  int maxlevels,
						  int innerits,
						  double convratio,
						  int** perm, 
						  int** invperm);

int
taucs_recursive_amwb_preconditioner_solve        (void* P, 
						  void* Z, 
						  void* R);

int taucs_dtl(ccs_etree)                         (taucs_ccs_matrix* A,
						  int* parent,
						  int* l_colcount,
						  int* l_rowcount,
						  int* l_nnz);

int      
taucs_dtl(ccs_symbolic_elimination)              (taucs_ccs_matrix* A,
						  void* L,
						  int do_order,
						  int max_depth
						  );

void* taucs_dtl(ccs_factor_llt_symbolic)         (taucs_ccs_matrix* A);
void* taucs_dtl(ccs_factor_llt_symbolic_maxdepth)(taucs_ccs_matrix* A,int max_depth);
taucs_cilk int   taucs_dtl(ccs_factor_llt_numeric)          (taucs_ccs_matrix* A,void* L);

taucs_cilk void* taucs_dtl(ccs_factor_llt_mf)               (taucs_ccs_matrix* A);
taucs_cilk void* taucs_dtl(ccs_factor_llt_mf_maxdepth)      (taucs_ccs_matrix* A,int max_depth);
void* taucs_dtl(ccs_factor_llt_ll)               (taucs_ccs_matrix* A);
void* taucs_dtl(ccs_factor_llt_ll_maxdepth)      (taucs_ccs_matrix* A,int max_depth);
int   taucs_dtl(supernodal_solve_llt)            (void* vL, void* x, void* b);
void taucs_dtl(supernodal_factor_free)                (void* L);
void taucs_dtl(supernodal_factor_free_numeric)        (void* L);
taucs_ccs_matrix* taucs_dtl(supernodal_factor_to_ccs) (void* L);
taucs_datatype* taucs_dtl(supernodal_factor_get_diag) (void* L);

int taucs_ccs_etree                              (taucs_ccs_matrix* A,
						  int* parent,
						  int* l_colcount,
						  int* l_rowcount,
						  int* l_nnz);

int      
taucs_ccs_symbolic_elimination                   (taucs_ccs_matrix* A,
						  void* L,
						  int do_order,
						  int max_depth
						  );

void* taucs_ccs_factor_llt_symbolic              (taucs_ccs_matrix* A);
void* taucs_ccs_factor_llt_symbolic_maxdepth     (taucs_ccs_matrix* A,int max_depth);
taucs_cilk int   taucs_ccs_factor_llt_numeric               (taucs_ccs_matrix* A,void* L);

taucs_cilk void* taucs_ccs_factor_llt_mf                    (taucs_ccs_matrix* A);
taucs_cilk void* taucs_ccs_factor_llt_mf_maxdepth           (taucs_ccs_matrix* A,int max_depth);
void* taucs_ccs_factor_llt_ll                    (taucs_ccs_matrix* A);
void* taucs_ccs_factor_llt_ll_maxdepth           (taucs_ccs_matrix* A,int max_depth);
int   taucs_supernodal_solve_llt                 (void* vL, void* x, void* b);
void taucs_supernodal_factor_free                (void* L);
void taucs_supernodal_factor_free_numeric        (void* L);
taucs_ccs_matrix* taucs_supernodal_factor_to_ccs (void* L);
void* taucs_supernodal_factor_get_diag           (void* L);

taucs_double taucs_vec_norm2(int n, int flags, void* x);
void* taucs_dtl(vec_create)  (int n);
void* taucs_vec_create       (int n, int flags);
void  taucs_dtl(vec_axpby)   (int n,
			      taucs_real_datatype a,taucs_datatype* x,
			      taucs_real_datatype b,taucs_datatype* y,
			      taucs_datatype* axpby);
void  taucs_vec_axpby        (int n,int flags,
			      taucs_double a,void* x,
			      taucs_double b,void* y,
			      void* axpby);
void taucs_dtl(vec_permute)  (int n, taucs_datatype v[],  taucs_datatype pv[], int p[]);
void taucs_dtl(vec_ipermute) (int n, taucs_datatype pv[], taucs_datatype v[],  int invp[]);
void taucs_vec_permute(int n, int flags, void* v, void* pv, int p[]);
void taucs_vec_ipermute(int n, int flags, void* v, void* pv, int p[]);


/*********************************************************/
/* Utilities                                             */
/*********************************************************/

void   taucs_logfile(char* file_prefix);
int    taucs_printf(char *fmt, ...);
int    taucs_maximize_stacksize(void);
double taucs_system_memory_size(void);
double taucs_available_memory_size(void);
double taucs_wtime(void);
double taucs_ctime(void);

/*********************************************************/
/* Out-of-core IO routines                               */
/*********************************************************/

taucs_io_handle* taucs_io_create_singlefile(char* filename);
taucs_io_handle* taucs_io_open_singlefile(char* filename);

taucs_io_handle* taucs_io_create_multifile(char* filename);
taucs_io_handle* taucs_io_open_multifile(char* filename);

int              taucs_io_close (taucs_io_handle* f);
int              taucs_io_delete(taucs_io_handle* f);

int              taucs_io_append(taucs_io_handle* f,
				 int   index,
				 int   m,int n,
				 int   flags,
				 void* data
				 );
int              taucs_io_write(taucs_io_handle* f,
				int   index,
				int   m,int n,
				int   flags,
				void* data
				);
int              taucs_io_read(taucs_io_handle* f,
			       int   index,
			       int   m,int n,
			       int   flags,
			       void* data
			       );

char*            taucs_io_get_basename(taucs_io_handle* f);

/*********************************************************/
/* Out-of-core Sparse Choleksy routines                  */
/*********************************************************/

int taucs_dtl(ooc_factor_llt)(taucs_ccs_matrix* A, 
			      taucs_io_handle*  L,
			      double memory);
/*added omer*/
int taucs_dtl(ooc_factor_llt_panelchoice)(taucs_ccs_matrix* A, 
					  taucs_io_handle* handle,
					  double memory,
					  int panelization_method);
/* end omer*/
int taucs_dtl(ooc_solve_llt) (void* L /* actual type: taucs_io_handle* */,
			      void* x, void* b);

int taucs_ooc_factor_llt(taucs_ccs_matrix* A, 
			 taucs_io_handle*  L,
			 double memory);
int taucs_ooc_solve_llt (void* L /* actual type: taucs_io_handle* */,
			 void* x, void* b);

/*********************************************************/
/* Out-of-core Sparse LU                                 */
/*********************************************************/

void taucs_dtl(ooc_factor_lu)(taucs_ccs_matrix* A_in,
		              int    colperm[],
                              taucs_io_handle* LU,
  	                      double memory);

int  taucs_ooc_factor_lu     (taucs_ccs_matrix* A_in,
		              int*   colperm,
                              taucs_io_handle* LU,
  	                      double memory);

int taucs_dtl(ooc_solve_lu)(taucs_io_handle*   LU,
			    taucs_datatype* x, 
                            taucs_datatype* b);

int taucs_ooc_solve_lu     (taucs_io_handle*   LU,
			    void* x, 
                            void* b);


/*********************************************************/
/* Utilities                                             */
/*********************************************************/

void   taucs_logfile(char* file_prefix);
int    taucs_printf(char *fmt, ...);
double taucs_system_memory_size(void);
double taucs_available_memory_size(void);
double taucs_wtime(void);
double taucs_ctime(void);

/*********************************************************/
/*                                                       */
/*********************************************************/



