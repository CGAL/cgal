/*********************************************************/
/* TAUCS                                                 */
/* Author: Doron Chen and Sivan Toledo                   */
/* File  : taucs_vaidya.c                                */
/* Description: constructs Vaidya's preconditioners      */
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "taucs.h"

/*long int random(void); omer*/

/*********************************************************/
/*                                                       */
/*********************************************************/

#ifdef TAUCS_CORE_DOUBLE

typedef unsigned char byte;
/*typedef int byte;*/

#define Do(i,n) for ((i)=0;(i)<(n);(i)++)

typedef struct {
  int n;
  int nent;
  int max_size;
  int *ivec1;
  int *ivec2;
  double *dvec;
} sym_matrix;

typedef struct {
  int    i;
  int    j;
  double v;
} wedge; /* weighted edge */

typedef struct {
  int n;
  int nent;
  int max_size;
  wedge* edges;
} graph;

/************** UNION FIND ***********/

static char *label = NULL;
static int *p      = NULL;
static int *rank   = NULL;

static
int unionfind_init(int size)
{
  int i;
  p = (int *)taucs_malloc(size * sizeof(int));
  rank = (int *)taucs_malloc(size * sizeof(int));
  label = (char *)taucs_malloc(size * sizeof(char));
  if (!p || !rank || !label) {
    taucs_free(p);
    taucs_free(rank);
    taucs_free(label);
    return -1;
  }

  Do(i,size)
    {
      p[i] = i;
      rank[i] = 0;
      label[i] = 0;
    }

  return 0;
}

static void 
unionfind_free (void)
{
  taucs_free(p);
  taucs_free(rank);
  taucs_free(label);
}

static
int Union(int a,int b,int x,int y,int l) /* unite a's and b's trees, whose roots are x and y. returns the root
					  of the united tree */
{
  if (rank[x] > rank[y])
    {
      p[y] = x;
      label[y] = label[a] ^ label[b] ^ l;
      return(x);
    }
  else
    {
      p[x] = y;
      label[x] = label[a] ^ label[b] ^ l;
      if (rank[x] == rank[y])
	rank[y]++;
      return(y);
    }
}

static
int find_set(int x)
{
  int tmp;

  if (x != p[x])
    {
      tmp = find_set(p[x]);
      label[x] ^= label[p[x]];
      p[x] = tmp;
    }
  return(p[x]);
}

/*********** HEAP OPERATIONS ************/

typedef struct hea {
  int heap_size;
  int alloc_size;
  int *edges;
  double *key;
} heap;

#define INF 100000000.0

#define Parent(i) ((((i)+1)/2) - 1)
#define Left(i) ((((i)+1)*2) - 1)
#define Right(i) ((((i)+1)*2 + 1) - 1)

static
void exchange(heap A,int a,int b)
{
  int tmp1;
  double tmp2;

  tmp1 = A.edges[a];
  A.edges[a] = A.edges[b];
  A.edges[b] = tmp1;

  tmp2 = A.key[a];
  A.key[a] = A.key[b];
  A.key[b] = tmp2;

}

static
void Heapify(heap A,int i)
{
  int l,r,smallest;
  
  l = Left(i);
  r = Right(i);
  
  if ((l < A.heap_size) && (A.key[l] < A.key[i]))
    smallest = l;
  else
    smallest = i;

  if ((r < A.heap_size) && (A.key[r] < A.key[smallest]))
    smallest = r;

  if (smallest != i)
    {
      exchange(A,i,smallest);
      Heapify(A,smallest);
    }
}

#if 0
static
int build_heap(int size,heap *h,graph *A)
{
  int i,k=0;

  h->heap_size = 0;
  h->edges = (int *)taucs_malloc(size * sizeof(int));
  h->key = (double *)taucs_malloc(size * sizeof(double));
  if (!(h->edges) || !(h->key)) {
    taucs_free(h->edges);
    taucs_free(h->key);
    return -1;
  }

  Do(i,size)
    if ((A->edges)[i].i != (A->edges)[i].j)
      {
	h->edges[k] = i;
	h->key[k] = fabs((A->edges)[i].v);
	k++;
      }
  
  h->heap_size = k;
  size = k;

  /*
  for(i=(size/2)-1;i>=0;i--)
    Heapify(*h,i);
  */

  return 0;
}
#endif

static
void free_heap(heap h)
{
  taucs_free(h.edges);
  taucs_free(h.key);
}


static int partition(heap h, int p, int r)
{
  int pivot;
  double x;
  int i,j;

  if (r-p < 16) pivot = 0;
  /* Sivan: chnaged random() to rand() to remove warning (random is not ansi C) */
  else if ((r - p) < 8) pivot = p + (rand() % (r-p+1)); 
  else {
    int c[3]; /* candidates */
    int t;
    c[0] = p + (rand() % (r-p+1));
    c[1] = p + (rand() % (r-p+1));
    c[2] = p + (rand() % (r-p+1));
    
    if (h.key[c[1]] < h.key[c[0]]) {t=c[0]; c[0]=c[1]; c[1]=t;}
    if (h.key[c[2]] < h.key[c[1]]) {t=c[1]; c[1]=c[2]; c[2]=t;}
    if (h.key[c[1]] < h.key[c[0]]) {t=c[0]; c[0]=c[1]; c[1]=t;}

    pivot = c[1];
  }

  x = h.key[pivot];
  /*
  i = p-1;
  j = r+1;

  while (1) {
    do {
      j--;
    } while ( h.key[j] > x );
    
    do {
      i++;
    } while ( h.key[i] < x );
    
    if (i < j)
      exchange(h,i,j);
    else 
      return j;
  }
  */
  i = p-1;
  j = r+1;

  while (1) {
    for (j--; h.key[j] > x; j--);
    for (i++; h.key[i] < x; i++);

    if (i < j)
      exchange(h,i,j);
    else 
      return j;
  }
}

#if 0
static void insertion_sort(heap h, int p, int r)
{
  int i,j;

  for (j=p+1; j<r; j++) {
    double key  = h.key[j];
    int    edge = h.edges[j];

    for (i=j-1; i>=p && h.key[i] > key; i--) {
      h.key[i+1] = h.key[i];
      h.edges[i+1] = h.edges[i];
    }

    h.key[i+1] = key;
    h.edges[i+1] = edge;
  }
}
#endif /* 0, we don't need insertion sort, heap sort */

static
void heapify_offset(heap A,int p,int r,int i)
{
  int L,R,smallest;
  int size = r-p+1;
  
  L = Left(i);
  R = Right(i);
  
  if ((L < size) && (A.key[p+L] < A.key[p+i]))
    smallest = L;
  else
    smallest = i;

  if ((R < size) && (A.key[p+R] < A.key[p+smallest]))
    smallest = R;

  if (smallest != i)
    {
      exchange(A,p+i,p+smallest);
      heapify_offset(A,p,r,smallest);
    }
}


static void heapsort_sort(heap h, int p, int r)
{
  int size = r-p+1;
  int i;

  for(i=(size/2)-1;i>=0;i--)
    heapify_offset(h,p,r,i);

  for(i=size-1;i>=1;i--)
    {
      exchange(h,0,i);
      heapify_offset(h,0,i,i);
    }
}

static void quick_sort(heap h, int p, int r)
{
  int q;
  if (p >= r) return;
  if (r - p < 100) {
    /*insertion_sort(h,p,r);*/
    heapsort_sort(h,p,r);
    return;
  }
  q = partition(h,p,r);
  quick_sort(h,p,q);
  quick_sort(h,q+1,r);
}

#if 0
static
int heap_sort(int size,heap *h,graph *A)
{
  int i;
  /*double *dvec;*/

  if (build_heap(size,h,A) == -1) 
    return -1;
  size = h->heap_size;

#define noQSORT
#ifdef QSORT
  quick_sort(*h,0,size-1);
#else

  for(i=(size/2)-1;i>=0;i--)
    Heapify(*h,i);


  for(i=size-1;i>=1;i--)
    {
      exchange(*h,0,i);
      h->heap_size--;
      Heapify(*h,0);
    }
#endif

  return(size); /* cannot be -1, so -1 is an error */
}
#endif /* 0, no heap_sort */

static int pqueue_fill(heap* h, graph* G)
{
  int i,size;

  size=0;
  Do(i,G->nent) {
    if ((G->edges)[i].i != (G->edges)[i].j) {
      assert(size <= h->alloc_size);
      h->edges[size] = i;
      h->key[size] = fabs((G->edges)[i].v);
      size++;
    }
  }
  
  h->heap_size = size;

#define noQSORT
#ifdef QSORT
  quick_sort(*h,0,size-1);
#else

  for(i=(size/2)-1;i>=0;i--)
    Heapify(*h,i);

  for(i=size-1;i>=1;i--) {
    exchange(*h,0,i);
    h->heap_size--;
    Heapify(*h,0);
  }
#endif

  h->heap_size = size;

  return size;
}

static int pqueue_create(heap* h, int size)
{
  h->heap_size  = 0;
  h->alloc_size = size;
  h->edges = (int *)taucs_malloc(size * sizeof(int));
  h->key = (double *)taucs_malloc(size * sizeof(double));
  if (!(h->edges) || !(h->key)) {
    taucs_free(h->edges);
    taucs_free(h->key);
    return -1;
  }
  return 0;
}

/***************************************************/
#ifdef GRAPHSORT

#define Parent(i) ((((i)+1)/2) - 1)
#define Left(i) ((((i)+1)*2) - 1)
#define Right(i) ((((i)+1)*2 + 1) - 1)

static
void new_exchange(wedge* e,int i,int j)
{
  wedge t;
  t = e[i];
  e[i] = e[j];
  e[j] = t;
  /*
  int tmp1;
  double tmp2;

  tmp1 = A.edges[a];
  A.edges[a] = A.edges[b];
  A.edges[b] = tmp1;

  tmp2 = A.key[a];
  A.key[a] = A.key[b];
  A.key[b] = tmp2;
  */
}

static
void heapify(wedge* e,int n,int i)
{
  int l,r,smallest;
  
  l = Left(i);
  r = Right(i);
  
  if ((l < n) && (fabs(e[l].v) < fabs(e[i].v)))
    smallest = l;
  else
    smallest = i;

  if ((r < n) && (fabs(e[r].v) < fabs(e[smallest].v)))
    smallest = r;

  if (smallest != i)
    {
      new_exchange(e,i,smallest);
      heapify(e,n,smallest);
    }
}

static
int new_heap_sort(wedge* e, int n)
{
  int i;

  for (i=(n/2)-1; i>=0; i--)
    heapify(e,n,i);

  for(i=n-1; i>=1; i--) {
    new_exchange(e,0,i);
    n--;
    heapify(e,n,0);
  }
}


int wedge_compare(const void* ve1, const void* ve2)
{
  wedge* e1 = (wedge*) ve1;
  wedge* e2 = (wedge*) ve2;

  double k1, k2;

  /*k1 = fabs(e1->v);*/
  /*k2 = fabs(e2->v);*/

  k1 = fabs(e1->v);
  k2 = fabs(e2->v);

  if (k1 < k2) return -1;
  if (k1 > k2) return  1;
  return 0;
}



static insertion_sort(wedge* e, int n)
{
  int i,j;

  for (j=1; j<n; j++) {
    double key = fabs(e[j].v);
    wedge  ej = e[j];
    for (i=j-1; i>=0 && fabs(e[i].j) > key; i--) {
      /*e[i+1] = e[i];*/

      e[i+1].i = e[i].i;
      e[i+1].j = e[i].j;
      e[i+1].v = e[i].v;
    }
    e[i+1].i = ej.i;
    e[i+1].j = ej.j;
    e[i+1].v = ej.v;
    /*    e[i+1] = ej;*/
  }
}

static int partition(wedge* e, int n)
{
  int pivot = (rand() % n);
  double x = fabs(e[pivot].v);
  int i,j;

  i = -1;
  j = n;

  while (1) {
    do {
      j--;
    } while ( fabs(e[j].v) > x );
    
    do {
      i++;
    } while ( fabs(e[i].v) < x );
    
    if (i < j) {
      /*
      wedge t;
      t = e[i];
      e[i] = e[j];
      e[j] = t;
      */

      int ti,tj; double tv;
      ti = e[i].i;
      tj = e[i].j;
      tv = e[i].v;
      e[i].i = e[j].i;
      e[i].j = e[j].j;
      e[i].v = e[j].v;
      e[j].i = ti;
      e[j].j = tj;
      e[j].v = tv;
    } else return j;
  }
}

static quick_sort(wedge* e, int n)
{
  int q;
  if (n <= 1) return;
  if (n < 32) {
    insertion_sort(e,n);
    return;
  }
  q = partition(e,n);
  quick_sort(e    ,q+1);
  quick_sort(e+q+1,n-q-1);
}

static int
graph_sort(graph* G) 
{
  /*
  qsort(G->edges,
	G->nent,
	sizeof(wedge),
	wedge_compare);
  */
  quick_sort(G->edges,G->nent);
}

#endif /* GRAPHSORT */
/************ VAIDYA'S PRECONDITIONERS *************/

#define swap(a,b) {int TMP; TMP = a; a = b; b = TMP;}
#define EPSILON 0.00000001

/************ GRAPHS *************/

static
graph* construct_graph(int size)
{
  graph *out;
  
  out = (graph *)taucs_malloc(sizeof(graph));
  if (!out) return NULL;

  out->edges = (wedge*) taucs_malloc(size*sizeof(wedge));
  if (!(out->edges)) {
    taucs_free(out);
    return NULL;
  }
  
  out->max_size = size;

  return out;
}

static
void free_graph(graph *a)
{
  if(a)
    {
      taucs_free(a->edges);
      taucs_free(a);
    }
}

void free_ccs_matrix(taucs_ccs_matrix *a)
{
  if (a)
    {
      taucs_free(a->rowind);
      taucs_free(a->colptr);
      taucs_free(a->values.d/*taucs_values*/);
      taucs_free(a);
    }
}

static
taucs_ccs_matrix* construct_ccs_matrix(int nent,int n)
{
  taucs_ccs_matrix *out;
  
  out = (taucs_ccs_matrix *)taucs_malloc(sizeof(taucs_ccs_matrix));
  if (!out) return NULL;
  out->colptr = (int *)taucs_malloc((n+1)*sizeof(int));
  out->rowind = (int *)taucs_malloc(nent*sizeof(int));
  out->values.d/*taucs_values*/ = (double *)taucs_malloc(nent*sizeof(double));
  if (!(out->colptr) || !(out->rowind) || !(out->values.d/*taucs_values*/)) {
    taucs_free(out->colptr);
    taucs_free(out->rowind);
    taucs_free(out->values.d/*taucs_values*/);
    taucs_free(out);
    return NULL;
  }
  
  out->n = n;
  out->m = n;
  out->flags = TAUCS_SYMMETRIC | TAUCS_LOWER | TAUCS_DOUBLE;
  
  return out;
}

#if 0
static
graph *ccs_matrix_to_graph(taucs_ccs_matrix *in)
{
  graph *out;
  int nent,n;
  int j,ip;

  nent = in->colptr[in->n];
  out = construct_graph(nent);
  if (!out) return NULL;

  n = in->n;

  out->n = n;
  out->nent = nent;
  
  for(j=0;j<n;j++) {
    for(ip=in->colptr[j];ip<in->colptr[j+1];ip++) {
      (out->edges)[ip].i = (in->rowind)[ip];
      (out->edges)[ip].j = j;
      (out->edges)[ip].v = (in->values.d/*taucs_values*/)[ip];
    }
  }
  
  return(out);
}
#endif /* 0, we don't need this routine */

static
int graph_resize(graph *a,int new_size)
{
  wedge* edges;
  
  assert(new_size > a->max_size);
  
  edges = (wedge*) taucs_malloc(new_size*sizeof(wedge));
  if (!edges) {
    return -1;
  }
  
  memcpy(edges,a->edges,a->max_size*sizeof(wedge));
  
  taucs_free(a->edges);

  a->edges=edges;

  a->max_size = new_size;

  return 0;
}

/************ LINKED LISTS *************/

typedef struct edg {
  int entry_no;
  struct edg *next;
} edge;

typedef struct linke {
  edge **point;
  edge *array;
} linked;

typedef struct thre {
  int group_1;
  int group_2;
  int a;
  int b;
  double c;
  byte already_connected;
  byte completed_to_basis;
  struct thre *next;
} three;

typedef struct si {
  int group_1;
  int group_2;
  int a[2];
  int b[2];
  double c[2];
  byte cross[2];
  byte no_edges; /* number of edges connecting group_1 and group_2 (0,1 or 2) */
  struct si *next;
} six;

int taucs_check_diag_dominant_matrix(graph *A, int force_diagonal_dominance)
{
  int i;
  double *sum;
  int n;
  int diagonally_dominant, all_nonpositive;

  n = A->n;

  sum = (double *)taucs_calloc(n,sizeof(double));
  if (!sum) return -1;

  Do(i,A->nent)
    {
      if ((A->edges)[i].i != (A->edges)[i].j)
	{
	  sum[(A->edges)[i].i]-=fabs((A->edges)[i].v);
	  sum[(A->edges)[i].j]-=fabs((A->edges)[i].v);
	}
      else
	{
	  sum[(A->edges)[i].i]+=(A->edges)[i].v;
	  if ((A->edges)[i].v < 0)
	    {
	      taucs_printf("ERROR! This matrix is not diagonally dominant. It has negative diagonals.\n");
	      /* taucs_free(sum); */
	      /* return -2; */
	    }
	}
    }
  
  diagonally_dominant = 1; /* until proven otherwise */
  all_nonpositive = 1;
  Do(i,n)
    {
      if (sum[i] < -EPSILON) diagonally_dominant = 0;
      if (sum[i] > EPSILON)  all_nonpositive     = 0;
    }

  if ((force_diagonal_dominance)&&(diagonally_dominant == 0)) {
    int first_time = 1;
    for(i=0;i<A->nent;i++)
      {
	if ((A->edges)[i].i == (A->edges)[i].j && sum[(A->edges)[i].i] <= EPSILON)
	  {
	    if (first_time) {
	      first_time=0; 
	      taucs_printf("\t\tAMWB warning: perturbing to force diagonal dominance\n");
	    }
	    (A->edges)[i].v -= sum[ (A->edges)[i].i ];
	    if (all_nonpositive && (A->edges)[i].i == 0) {
	      taucs_printf("taucs warning: perturbing to ensure strict diagonal dominance\n");
	      (A->edges)[i].v += 0.1; /* arbitrary perturbation */
	    }
	  }
      }
  } else
    if (diagonally_dominant == 0)
      {
	taucs_printf("ERROR! This matrix is not diagonally dominant. sum[%d] = %lf\n",i,sum[i]);
	taucs_free(sum);
	return -2;
      }
  
  taucs_free(sum);
  return 0;
}

#if 0
static
double *analyze_graph(graph *A)
{
  int i;
  int t1,t2;
  int n;
  double *diag,t3;
  
  n = A->n;
  diag = (double *)taucs_calloc(n,sizeof(double));
  if (!diag) return NULL;

  Do(i,A->nent)
    {
      t1=(A->edges)[i].i;
      t2=(A->edges)[i].j;
      t3=(A->edges)[i].v;

      if (t1 == t2)
	diag[t1] += fabs(t3);
      else
	{
	  diag[t1] -= fabs(t3);
	  diag[t2] -= fabs(t3);
	}
    }
  return(diag);
}
#endif /* 0, we don't need this routine */

/*********************************************************/
/* ccs diagnostics, row sums, and conversion to a graph  */
/*********************************************************/

#define TAUCS_SYM_NOT_SYMLOWER     1
#define TAUCS_SYM_POS_OFFDIAGONALS 2
#define TAUCS_SYM_NEG_DIAGONALS    4
#define TAUCS_SYM_NOT_DIAGDOMINANT 8

static
graph *ccs_matrix_to_graph_plus(taucs_ccs_matrix *in,
				int*   diagnostics,
				double diag[],
				int    force_diagonal_dominance)
{
  graph *out;
  int nent,n;
  int i,j,k,ip;
  double v;
  int negative_on_diagonal;
  int positive_off_diagonal;
  int not_diagonally_dominant;

  *diagnostics = 0;

  if (!(in->flags & TAUCS_SYMMETRIC) || !(in->flags & TAUCS_LOWER)) {
    *diagnostics = TAUCS_SYM_NOT_SYMLOWER;
    return NULL;
  }

  nent = in->colptr[in->n];
  out = construct_graph(nent);
  if (!out) return NULL;

  n = in->n;

  out->n = n;
  out->nent = nent;
  
  for (i=0; i<n; i++) diag[i] = 0.0;

  negative_on_diagonal  = 0;
  positive_off_diagonal = 0;

  for(j=0;j<n;j++) {
    for(ip=in->colptr[j];ip<in->colptr[j+1];ip++) {
      i = (in->rowind)[ip];
      v = (in->values.d/*taucs_values*/)[ip];
      (out->edges)[ip].i = i;
      (out->edges)[ip].j = j;
      (out->edges)[ip].v = v;

      if (i == j) {
	negative_on_diagonal |= (v < 0.0);

	diag[i] += fabs(v);
      } else {
	positive_off_diagonal |= (v > 0.0);

	diag[i] -= fabs(v);
	diag[j] -= fabs(v);
      }
    }
  }
  
  if (force_diagonal_dominance) {
    int strict_diagdominance = 0;
    int first_time = 1;

    for (i=0; i<n; i++) 
      strict_diagdominance |= (diag[i] > 0.0);

    for(k=0;k<out->nent;k++) {
      i = (out->edges)[k].i;
      j = (out->edges)[k].j;
      v = (out->edges)[k].v;

      if (i == j && diag[i] < 0.0) {
	if (first_time) {
	  first_time=0; 
	  taucs_printf("taucs warning: perturbing to force diagonal dominance\n");
	}
	(out->edges)[k].v -= diag[i];
	diag[i] = 0.0;
	if (strict_diagdominance == 0 && i == 0) {
	  taucs_printf("taucs warning: perturbing to ensure strict diagonal dominance\n");
	  (out->edges)[k].v += 1e-8; /* arbitrary perturbation */
	}
      }
    }

    not_diagonally_dominant = 0;

  } else {

    not_diagonally_dominant = 0;
    for (i=i; i<n; i++) 
      not_diagonally_dominant |= (diag[i] < -1e-12); /* arbitrary threashold */

  }

  *diagnostics = 0;
  
  if (not_diagonally_dominant) *diagnostics |= TAUCS_SYM_NOT_DIAGDOMINANT;
  if (negative_on_diagonal   ) *diagnostics |= TAUCS_SYM_NEG_DIAGONALS;
  if (positive_off_diagonal  ) *diagnostics |= TAUCS_SYM_POS_OFFDIAGONALS;

  return(out);
}


/*********************************************************/
/*                                                       */
/*********************************************************/


static
void free_linked_list(linked* a)
{
  taucs_free(a->point);
  taucs_free(a->array);
  taucs_free(a);
}

static
void free_linked_list_2(six *a)
{
  if (a!=NULL)
    {
      free_linked_list_2(a->next);
      taucs_free(a);
    }

}

static
linked* create_linked_list(graph *A,int n,int Anent,double *min,double *max)
{
  /* creates linked list which holds the off-diagonal entries of the sparse graph. */
  int i;
  edge *tmp;
  /*linked out; */
  linked* out;
  int free_place = 0;

  *min = INF;
  *max = -INF;

  out = (linked*) taucs_malloc(sizeof(linked));
  if (!out) {
    return NULL;
  }
  out->point = (edge **)taucs_calloc(n,sizeof(edge *));
  out->array = (edge *)taucs_calloc(2*Anent,sizeof(edge));
  if (!(out->point) || !(out->array)) {
    taucs_free(out->point);
    taucs_free(out->array);
    taucs_free(out);
    return NULL;
  }
  
  Do(i,Anent)
    {
      if ((A->edges)[i].i != (A->edges)[i].j)
	{
	  if (-(A->edges)[i].v < *min)
	    *min = -(A->edges)[i].v;
	  if (-(A->edges)[i].v > *max)
	    *max = -(A->edges)[i].v;
	  
	  tmp = &(out->array[free_place++]);
	  tmp->entry_no = i;
	  tmp->next = out->point[(A->edges)[i].i];
	  out->point[(A->edges)[i].i] = tmp;

	  tmp = &(out->array[free_place++]);
	  tmp->entry_no = i;
	  tmp->next = out->point[(A->edges)[i].j];
	  out->point[(A->edges)[i].j] = tmp;
	}
    }

  return(out);
}

static
linked* create_linked_list_cluster(graph *A,int n,int Anent,double *min,double *max,int *partition,int *new_partition)
{
  /* creates linked list which holds off-diagonal entries of the sparse graph.
   This linked list contains all the edges which connect vertices whose endpoints are
   in different sections in partition, but in the same section in new_partition.
   This will help us build trees within each section in new_partition. Each vertex
   in these trees will be a contracted section in partition */
  int i;
  edge *tmp;

  linked* out = NULL;
  int free_place = 0;

  *min = INF;
  *max = -INF;

  out = (linked*) taucs_malloc(sizeof(linked));
  if (!out) {
    return NULL;
  }
  out->point = (edge **)taucs_calloc(n,sizeof(edge *));
  out->array = (edge *)taucs_calloc(2*Anent,sizeof(edge));
  if (!(out->point) || !(out->array)) {
    taucs_free(out->point);
    taucs_free(out->array);
    taucs_free(out);
    return NULL;
  }
  
  Do(i,Anent)
    {
      if ((partition[(A->edges)[i].i] != partition[(A->edges)[i].j]) &&
	  (new_partition[(A->edges)[i].i] == new_partition[(A->edges)[i].j]))
	{
	  if (-(A->edges)[i].v < *min)
	    *min = -(A->edges)[i].v;
	  if (-(A->edges)[i].v > *max)
	    *max = -(A->edges)[i].v;
	  
	  tmp = &(out->array[free_place++]);
	  tmp->entry_no = i;
	  tmp->next = out->point[partition[(A->edges)[i].i]];
	  out->point[partition[(A->edges)[i].i]] = tmp;

	  tmp = &(out->array[free_place++]);
	  tmp->entry_no = i;
	  tmp->next = out->point[partition[(A->edges)[i].j]];
	  out->point[partition[(A->edges)[i].j]] = tmp;
	}
    }

  return(out);
}


static
taucs_ccs_matrix *graph_to_ccs_matrix(graph *A)
{
  taucs_ccs_matrix *out;
  int n,nent,i,j1,j2;
  /*int count=0;*/
  int *tmp;

  n = A->n;
  nent = A->nent;

  tmp = (int *)taucs_malloc(n*sizeof(int));
  if (!tmp) return NULL;

  out=construct_ccs_matrix(nent,n);
  if (!out) {
    taucs_free(tmp);
    return NULL;
  }
  out->flags = TAUCS_SYMMETRIC | TAUCS_LOWER | TAUCS_DOUBLE;

  Do(i,n)
    tmp[i] = 0;
  Do(i,nent)
    tmp[min((A->edges)[i].i,(A->edges)[i].j)]++;
  out->colptr[0] = 0;
  Do(i,n)
    out->colptr[i+1] = out->colptr[i] + tmp[i];

  Do(i,n)
    tmp[i] = out->colptr[i];

  Do(i,nent)
    {
      j1 = min((A->edges)[i].i , (A->edges)[i].j);
      j2 = max((A->edges)[i].i , (A->edges)[i].j);
      out->rowind[tmp[j1]]=j2;
      out->values.d/*taucs_values*/[tmp[j1]]=(A->edges)[i].v;
      tmp[j1]++;
    }

  taucs_free(tmp);
  return(out);
}

static
int compute_sub_tree_sizes(int ver,int *first_child,int *next_child,int *sizes)
{
  int sum = 1,v;

  if (first_child[ver] == -1)
    {
      sizes[ver] = 1;
      return(1);
    }
  else
    {
      v=first_child[ver];
      while(v != -1)
	{
	  sum += compute_sub_tree_sizes(v,first_child,next_child,sizes);
	  v = next_child[v];
	}
    }
  sizes[ver] = sum;
  return(sum);
}

static
void assign_group(int ver,int gr,int *first_child,int *next_child,int *groups)
{
  int v;

  groups[ver] = gr;
  if (first_child[ver] != -1)
    {
      v=first_child[ver];
      while(v != -1)
	{
	  assign_group(v,gr,first_child,next_child,groups);
	  v = next_child[v];
	}
    }
}

static
int create_children_arrays(int *pi,int n,int **fc,int **nc)
{
  int *first_child,*next_child;
  int father,child,ch;
  int i;
  
  first_child = (int *)taucs_malloc(n*sizeof(int));
  next_child = (int *)taucs_malloc(n*sizeof(int));
  if (!first_child || !next_child) {
    taucs_free(first_child);
    taucs_free(next_child);
    return -1;
  }

  Do(i,n)
    first_child[i] = next_child[i] = -1;


  Do(i,n)
    {
      child = i;
      father = pi[i];
      
      if (father != -1)
	{
	  if (first_child[father] == -1)
	    first_child[father] = child;
	  else
	    {
	      ch = first_child[father];
	      while(next_child[ch] != -1)
		ch = next_child[ch];
	      next_child[ch] = child;
	      
	    }
	}
    }

  
  *fc = first_child;
  *nc = next_child;

  return 0; /* success */
}

static
void disconnect(int father,int child,
		int *first_child,int *next_child,int *pi)
     /* disconnect subtree whose root is 'child', from the tree */
{

  int oldfather;
  int v;
  /* int tmp;*/
  
  oldfather = father;

  assert(first_child[father] != -1);
  
  if (first_child[father] == child)
    first_child[father] = next_child[child];
  else
    {
      v = first_child[father];
      while(next_child[v] != child)
	v = next_child[v];

      next_child[v] = next_child[next_child[v]];

    }

}

static
void divide_to_groups(int r,int *first_child,int *next_child,int *pi,int *curr_group,
		      int *groups,int *sub_tree_sizes,int root,double subgraphs,
		      int n)
     /* divides the vertices into different groups (divides the tree is subtrees) */
{
  int v;
  double low;

  low = max(1,((double)n/subgraphs));

  if(first_child[r] != -1)
    {
      v=first_child[r];
      sub_tree_sizes[r] = 1;
      while(v != -1)
	{
	  if (sub_tree_sizes[v] > low) 
	    divide_to_groups(v,first_child,next_child,pi,curr_group,groups,
			     sub_tree_sizes,root,subgraphs,n);
	  
	  if (sub_tree_sizes[v] >= low)
	    {
	      assign_group(v,*curr_group,first_child,next_child,groups);
	      disconnect(r,v,first_child,next_child,pi);
	      (*curr_group)++;
	    }
	  else
	    sub_tree_sizes[r] += sub_tree_sizes[v];
	  v = next_child[v];
	}
      
    }
  
}


static
void DFS_visit(graph *precond,int r,byte *color,linked l,int *pi,int *visited)
{
  edge *p;
  int r1;
  
  color[r] = 1;
  (*visited)++;

  p = l.point[r];
  while (p != NULL)
    {
      /* this looks strange. Sivan */
      /*r1 = ivec1[ p->entry_no ]+ivec2[ p->entry_no ]-r;*/
      r1 = (precond->edges)[ p->entry_no ].i + (precond->edges)[ p->entry_no ].j - r;
      if (color[r1]==0)
	{
	  DFS_visit(precond,r1,color,l,pi,visited);
	  pi[r1] = r;
	}
      p = p->next;
    }
}

static
void make_perm(int perm[],int k)
{
  int i,tmp,tmp1;
 
 for(i=0;i<k;i++)
    perm[i]=i;

  for(i=0;i<k;i++)
    {
      tmp = rand()%(k-i);
      tmp1 = perm[i+tmp];
      perm[i+tmp]=perm[i];
      perm[i]=tmp1;
    }
}

#define USE_HEAPSORT

static
taucs_ccs_matrix*
amwb_preconditioner_create(graph *mtxA, double* diag,
			   int rnd,
			   double subgraphs)
     /*
amwb_preconditioner_create(taucs_ccs_matrix *taucs_ccs_mtxA, 
			   int rnd,
			   double subgraphs)
     */
{
  taucs_ccs_matrix *out;

  /*  graph *mtxA,*mtxA_tmp;*/
  graph *mtxA_tmp;
  graph *precond;

  int i;
  /*int j,tmp,entry;*/
  /*
  int *ivec1, *ivec2 ;
  int *ivec1_p, *ivec2_p ;
  double *dvec; 
  double *dvec_p; 
  */

  int size;
  int *pi;
  /*edge **array;*/
  /*edge *p;*/
#ifdef USE_HEAPSORT
  heap Ah;
  heap Bh;
  heap h;
#endif
  int Bent;
  /*int row,col;*/
  double weight;
  int *first_child,*next_child;
  int *groups,*sub_tree_sizes;
  int curr_group;
  int n,Anent;
  int precond_maxnent,chunk;
  /*double *diag;*/
  byte *closed_cycle,closed_cycle_x,closed_cycle_y;
  byte *color;
  byte *already_added;
  int count=0;
  int edge_sign; /*  byte edge_sign;*/
  int u,v,x,y,un_root,r;
  linked* l;
  int *perm;
  int visited,*roots,rcount;
  int basis_Bent=0,step2_Bent=0;
  three *complete_subgraph;
  six **pairs;
  /*char bool;*/
  /* FILE *graph_file0,*graph_file1,*graph_file2,*graph_file3,*group_file; */

  double wtime;
  double wtime_sort  = 0.0;
  double wtime_global_basis;
  double wtime_treepartition;
  double wtime_component_bases;
  double wtime_pair_bases;
  double wtime_total;
  double dummy;

  wtime_total = taucs_wtime();

  /* graph_file0 = fopen("graphfile0.txt","w"); */
  /* graph_file1 = fopen("graphfile1.txt","w"); */
  /* graph_file2 = fopen("graphfile2.txt","w"); */
  /* graph_file3 = fopen("graphfile3.txt","w"); */
  /* group_file  = fopen("groupfile.txt" ,"w"); */


  /********************************************************/
  /*                                                      */
  /********************************************************/

  /********************************************************/
  /* convert matrix to a graph                            */
  /********************************************************/

  /*** ALLOCATED: NONE ***/

  /*
  wtime = taucs_wtime();
  mtxA = ccs_matrix_to_graph(taucs_ccs_mtxA);
  if (!mtxA) {
    return NULL;
  }
  wtime = taucs_wtime() - wtime;
  taucs_printf("\t\tAMWB matrix-to-graph = %.3f seconds\n",wtime);
  */

  /********************************************************/
  /* check that the matrix is diagonally dominant         */
  /********************************************************/

  /*** ALLOCATED: mtxA ***/
  
#if 0
  wtime = taucs_wtime();
  i = taucs_check_diag_dominant_matrix(mtxA,1 /* force diagonal dominance */);
  if (i == -1) {
    free_graph(mtxA);
    return NULL;
  }
  if (i == -2) {
    free_graph(mtxA);
    return taucs_ccs_mtxA; /* not diagonally dominant */
  }
  wtime = taucs_wtime() - wtime;
  taucs_printf("\t\tAMWB check-diag-dominance = %.3f seconds\n",wtime);
#endif

  n = mtxA->n;

  /********************************************************/
  /* generate random permutation and permute vertices     */
  /********************************************************/

  wtime = taucs_wtime();
  perm = (int *)taucs_malloc(mtxA->nent*sizeof(int));
  if (!perm) {
    free_graph(mtxA);
    return NULL;
  }

  /*** ALLOCATED: mtxA,perm ***/

  make_perm(perm,mtxA->nent);

  mtxA_tmp = construct_graph(mtxA->max_size);
  if (!mtxA_tmp) {
    free_graph(mtxA);
    taucs_free(perm);
    return NULL;
  }

  /*** ALLOCATED: mtxA,perm,mtxA_tmp ***/

  mtxA_tmp->nent = mtxA->nent;
  mtxA_tmp->n = mtxA->n;
  Do(i,mtxA->nent)
    {
      /*
      mtxA_tmp->ivec1[i] = mtxA->ivec1[perm[i]];
      mtxA_tmp->ivec2[i] = mtxA->ivec2[perm[i]];
      mtxA_tmp->dvec[i]  = mtxA->dvec[perm[i]];
      */

      mtxA_tmp->edges[i].i = mtxA->edges[perm[i]].i;
      mtxA_tmp->edges[i].j = mtxA->edges[perm[i]].j;
      mtxA_tmp->edges[i].v = mtxA->edges[perm[i]].v;
    }

  taucs_free(perm);
  free_graph(mtxA);

  wtime = taucs_wtime() - wtime;
  taucs_printf("\t\tAMWB random permute = %.3f seconds\n",wtime);

  /********************************************************/
  /* compute and remember row weights                     */
  /********************************************************/

  wtime = taucs_wtime();

  /*** ALLOCATED: mtxA_tmp ***/

  /*
  diag = analyze_graph(mtxA_tmp);
  if (!diag) {
    free_graph(mtxA_tmp);
    return NULL;
  }
  */

  wtime = taucs_wtime() - wtime;
  taucs_printf("\t\tAMWB row weights = %.3f seconds\n",wtime);

  /********************************************************/
  /* allocate vectors                                     */
  /********************************************************/

  /*** ALLOCATED: mtxA_tmp,diag ***/

  Anent = mtxA_tmp->nent;

  already_added = (byte *)taucs_calloc(Anent,sizeof(byte));
  pi            = (int *) taucs_malloc(n*sizeof(int));
  if (!already_added || !pi) {
    free_graph(mtxA_tmp);
    taucs_free(diag);
    taucs_free(already_added);
    taucs_free(pi);
    return NULL;
  }

  /*** ALLOCATED: mtxA_tmp,diag,already_added,pi ***/

  Do(i,n)
    pi[i] = -1;

  /********************************************************/
  /* construct empty preconditioner                       */
  /********************************************************/

  precond_maxnent = 3*n;

  taucs_printf("allocating space for %d entries in precond\n",precond_maxnent);fflush(stdout);
  
  precond = construct_graph(precond_maxnent);
  if (!precond) {
    free_graph(mtxA_tmp);
    taucs_free(diag);
    taucs_free(already_added);
    taucs_free(pi);
    return NULL;
  }

  /*** ALLOCATED: mtxA_tmp,diag,already_added,pi,precond ***/

  precond->n=mtxA_tmp->n;

  /*
  ivec1_p = precond->ivec1 ;
  ivec2_p = precond->ivec2 ;
  dvec_p = precond->dvec ;
  */
	 
  Bent = 0;

  /********************************************************/
  /* allocate vectors                                     */
  /********************************************************/

  wtime_global_basis = taucs_wtime();

  closed_cycle = (byte *)taucs_calloc(n,sizeof(byte));
  if (!closed_cycle) {
    free_graph(mtxA_tmp);
    free_graph(precond);
    taucs_free(diag);
    taucs_free(already_added);
    taucs_free(pi);
    return NULL;
  }

  /*** ALLOCATED: mtxA_tmp,diag,already_added,pi,precond,closed_cycle ***/

  /* Variation on Kruskal - Introduction to Algorithms page 505 */

  /********************************************************/
  /* initialize union-find                                */
  /********************************************************/

  if (unionfind_init(n) == -1) {
    free_graph(mtxA_tmp);
    free_graph(precond);
    taucs_free(diag);
    taucs_free(already_added);
    taucs_free(pi);
    taucs_free(closed_cycle);
    return NULL;
  }

  /********************************************************/
  /* sort edges of matrix                                 */
  /********************************************************/

  /*** ALLOCATED: mtxA_tmp,diag,already_added,pi,precond,closed_cycle,UF ***/

  wtime = taucs_wtime();
#ifdef USE_HEAPSORT


  /*
  size = heap_sort(Anent,&h,mtxA_tmp);
  if (size == -1) {
  */
  if (pqueue_create(&Ah,Anent) == -1) {
    free_graph(mtxA_tmp);
    free_graph(precond);
    taucs_free(diag);
    taucs_free(already_added);
    taucs_free(pi);
    taucs_free(closed_cycle);
    unionfind_free();
    return NULL;
  }
  size = pqueue_fill(&Ah,mtxA_tmp);
#else
  assert(Anent == mtxA_tmp->nent);
  size = mtxA_tmp->nent;
  graph_sort(mtxA_tmp);
#endif
  wtime = taucs_wtime() - wtime;
  taucs_printf("\t\tAMWB sort(%d) = %.3f seconds\n",Anent,wtime);
  wtime_sort += wtime;

  /*** ALLOCATED: mtxA_tmp,diag,already_added,pi,precond,closed_cycle,UF,heap ***/

  /********************************************************/
  /* build a basis for the entire graph                   */
  /********************************************************/

  Do(i,size)
    {
      if (count == n)
	break;

#ifdef USE_HEAPSORT
      u = mtxA_tmp->edges[Ah.edges[i]].i;
      v = mtxA_tmp->edges[Ah.edges[i]].j;
      weight = mtxA_tmp->edges[Ah.edges[i]].v;
#else
      u      = mtxA_tmp->edges[i].i;
      v      = mtxA_tmp->edges[i].j;
      weight = mtxA_tmp->edges[i].v;
      if (u==v) continue;
#endif
      
      edge_sign = (weight>0);

      x = find_set(u);
      y = find_set(v);

      if (x!=y)
	{
	  /* printf("different trees\n"); */
	  if (!((closed_cycle[x])&&(closed_cycle[y])))
	    {
	      count++;
	      /* printf("(%d,%d) - %lf\n",u,v,weight); */
	      /*
	      ivec1_p[Bent] = u; 
	      ivec2_p[Bent] = v; 
	      dvec_p[Bent] = weight; 
	      */
	      precond->edges[Bent].i = u; 
	      precond->edges[Bent].j = v; 
	      precond->edges[Bent].v = weight; 

	      Bent++; 

	      diag[u] += fabs(weight);
	      diag[v]+=fabs(weight);

#ifdef USE_HEAPSORT
	      already_added[Ah.edges[i]] = 1;
#else
	      already_added[i] = 1;
#endif
	      un_root = Union(u,v,x,y,edge_sign);
	      closed_cycle[un_root] = closed_cycle[x] | closed_cycle[y];
	    }
	  /* else
	    {
	      printf("cannot add (%d,%d) - %lf - both trees already have a cycle\n",u,v,weight); 
	    }
	  */
	}
      else
	{
	  /* printf("same tree\n"); */
	  if ((edge_sign != (label[u]^label[v])) && (closed_cycle[x]==0))
	    {
	      count++;
	      /* printf("(%d,%d) - %lf\n",u,v,weight); */
	      /*
	      ivec1_p[Bent] = u; 
	      ivec2_p[Bent] = v; 
	      dvec_p[Bent] = weight; 
	      */
	      precond->edges[Bent].i = u; 
	      precond->edges[Bent].j = v; 
	      precond->edges[Bent].v = weight; 

	      Bent++; 

	      diag[u] += fabs(weight);
	      diag[v]+=fabs(weight);

#ifdef USE_HEAPSORT
	      already_added[Ah.edges[i]] = 1;
#else
	      already_added[i] = 1;
#endif
	      closed_cycle[x] = 1;
	    }
	  /* else
	    {
	      if (closed_cycle[x]==1)
		printf("cannot add (%d,%d) - %lf - tree already contains cycle\n",u,v,weight); 
	      else
		printf("cannot add (%d,%d) - %lf - it closes a positive cycle\n",u,v,weight); 
	    }
	  */
	}
    }

  /********************************************************/
  /* the preconditioner is now a max-weight-basis         */
  /********************************************************/

  wtime_global_basis = taucs_wtime() - wtime_global_basis;
  taucs_printf("\t\tAMWB global basis = %.3f seconds\n",wtime_global_basis);
  
  taucs_free(closed_cycle);
#ifdef USE_HEAPSORT
  /*free_heap(h);*/
#endif
  unionfind_free();

  precond->nent = Bent;
  basis_Bent = Bent;

  /********************************************************/
  /* break into subgraphs                                 */
  /********************************************************/

  /*** ALLOCATED: mtxA_tmp,diag,already_added,pi,precond ***/

  l = create_linked_list(precond,n,Bent,&dummy,&dummy);
  if (!l) {
    free_graph(mtxA_tmp);
    free_graph(precond);
    taucs_free(diag);
    taucs_free(already_added);
    taucs_free(pi);
    return NULL;
  }

  /*** ALLOCATED: mtxA_tmp,diag,already_added,pi,precond,linked ***/

  color = (byte *)taucs_calloc(n,sizeof(byte));
  roots = (int *)taucs_malloc(n*sizeof(int));
  if (!color || !roots) {
    free_graph(mtxA_tmp);
    free_graph(precond);
    taucs_free(diag);
    taucs_free(already_added);
    taucs_free(pi);
    free_linked_list(l);
    return NULL;
  }

  /*** ALLOCATED: mtxA_tmp,diag,already_added,pi,precond,linked,color,roots ***/

  wtime_treepartition = taucs_wtime();

  visited = 0;
  rcount = 0;
  /*
  while(visited<n)
    {
      r = rand()%n;
      while (color[r])
	r = rand()%n;
      roots[rcount++] = r;
      pi[r] = -1;
      DFS_visit(precond,r,color,*l,pi,&visited);
    }
  */

  for (r=0; r<n; r++) {
    if (color[r] != 0) continue;
    roots[rcount++] = r;
    pi[r] = -1;
    DFS_visit(precond,r,color,*l,pi,&visited);
  }

  taucs_free(color);

  /*** ALLOCATED: mtxA_tmp,diag,already_added,pi,precond,linked,roots ***/

  /* pi now contains the parent array of the tree */

  if (create_children_arrays(pi,n,&first_child,&next_child) == -1) {
    free_graph(mtxA_tmp);
    free_graph(precond);
    taucs_free(diag);
    taucs_free(already_added);
    taucs_free(pi);
    taucs_free(roots);
    free_linked_list(l);
    return NULL;
  }

  /*** ALLOCATED: mtxA_tmp,diag,already_added,pi,precond,linked,roots,FC,NC ***/

  /* the first_child/next_child arrays enable us to find the children
     of any vertex in the tree */
  
  groups = (int *)taucs_malloc(n*sizeof(int));
  sub_tree_sizes = (int *)taucs_malloc(n*sizeof(int));
  if(!groups || !sub_tree_sizes) {
    free_graph(mtxA_tmp);
    free_graph(precond);
    taucs_free(diag);
    taucs_free(already_added);
    taucs_free(pi);
    taucs_free(roots);
    taucs_free(groups);
    taucs_free(sub_tree_sizes);
    taucs_free(first_child);
    taucs_free(next_child);
    free_linked_list(l);
    return NULL;
  }

  /*** ALLOCATED: mtxA_tmp,diag,already_added,pi,precond,linked,roots,FC,NC ***/
  /*** ALLOCATED: groups,sub_tree_sizes ***/

  Do(i,n)
    groups[i] = -1;
  
  curr_group = 0;
  
  Do(i,rcount)
    {
      r = roots[i];
      compute_sub_tree_sizes(r,first_child,next_child,sub_tree_sizes);
      /* now for every vertex v in the tree, sub_tree_sizes[v] is the size
	 of the subtree whose root is v */
      
      divide_to_groups(r,first_child,next_child,pi,&curr_group,groups,
		       sub_tree_sizes,r,subgraphs,n);
      if ((sub_tree_sizes[r]<((double)n/subgraphs))&&(curr_group>0))
	curr_group--;
      assign_group(r,curr_group,first_child,next_child,groups);
      curr_group++;
    }

  taucs_printf("actual number of subgraphs = %ld\n",curr_group);fflush(stdout);
  /* now the tree is divided into linked groups */

  chunk = max(min((curr_group*(curr_group-1)/2)/10,5000),10000);

  taucs_free(roots);
  taucs_free(first_child);
  taucs_free(next_child);
  taucs_free(sub_tree_sizes);

  wtime_treepartition = taucs_wtime() - wtime_treepartition;
  taucs_printf("\t\tAMWB treepartition = %.3f seconds\n",wtime_treepartition);




  wtime_component_bases = taucs_wtime();

  /*** ALLOCATED: mtxA_tmp,diag,already_added,pi,precond,linked ***/
  /*** ALLOCATED: groups ***/

  precond->nent = Bent;

  /********************************************************/
  /* complete each subgraph into a basis                  */
  /********************************************************/

  complete_subgraph = (three *)taucs_calloc(curr_group,sizeof(three));
  /*** ALLOCATED: mtxA_tmp,diag,already_added,pi,precond,linked ***/
  /*** ALLOCATED: groups ***/
  if (!complete_subgraph) {
    free_graph(mtxA_tmp);
    free_graph(precond);
    taucs_free(diag);
    taucs_free(already_added);
    taucs_free(pi);
    taucs_free(groups);
    free_linked_list(l);
    return NULL;
  }
  
  /*** ALLOCATED: mtxA_tmp,diag,already_added,pi,precond,linked ***/
  /*** ALLOCATED: groups,complete_subgraph ***/
  
  /* For each subgraph, complete_subgraph will store the edge needed
     to complete the subgraph into a basis (if such an edge exists) */
  
  /* Variation on Kruskal - Introduction to Algorithms page 505 */
  if (unionfind_init(n) == -1) {
    free_graph(mtxA_tmp);
    free_graph(precond);
    free_linked_list(l);
    taucs_free(diag);
    taucs_free(already_added);
    taucs_free(pi);
    taucs_free(groups);
    taucs_free(complete_subgraph);
    return NULL;
  }
  
  /*** ALLOCATED: mtxA_tmp,diag,already_added,pi,precond,linked ***/
  /*** ALLOCATED: groups,complete_subgraph,UF ***/
  
  closed_cycle = (byte *)taucs_calloc(n,sizeof(byte));
  if (!closed_cycle) {
    free_graph(mtxA_tmp);
    free_graph(precond);
    free_linked_list(l);
    unionfind_free();
    taucs_free(diag);
    taucs_free(already_added);
    taucs_free(pi);
    taucs_free(groups);
    taucs_free(complete_subgraph);
    return NULL;
  }
  
  /*** ALLOCATED: mtxA_tmp,diag,already_added,pi,precond,linked ***/
  /*** ALLOCATED: groups,complete_subgraph,UF,closed_cycle ***/
  
  wtime = taucs_wtime();
#ifdef USE_HEAPSORT
  /*
    size = heap_sort(Bent,&h,precond);
    if (size == -1) {
    free_graph(mtxA_tmp);
    free_graph(precond);
    free_linked_list(l);
    unionfind_free();
    taucs_free(diag);
    taucs_free(already_added);
    taucs_free(pi);
    taucs_free(groups);
    taucs_free(complete_subgraph);
    taucs_free(closed_cycle);
    return NULL;
    }
  */
  if (pqueue_create(&Bh,Anent) == -1) {
    free_graph(mtxA_tmp);
    free_graph(precond);
    taucs_free(diag);
    taucs_free(already_added);
    taucs_free(pi);
    taucs_free(closed_cycle);
    unionfind_free();
    return NULL;
  }
  size = pqueue_fill(&Bh,precond);
#else
  assert(Bent == precond->nent);
  size = precond->nent;
  graph_sort(precond);
#endif
  wtime = taucs_wtime() - wtime;
  taucs_printf("\t\tAMWB sort 2(%d) = %.3f seconds\n",Bent,wtime);
  wtime_sort += wtime;
  
  /*** ALLOCATED: mtxA_tmp,diag,already_added,pi,precond,linked ***/
  /*** ALLOCATED: groups,complete_subgraph,UF,closed_cycle,heap ***/

  Do(i,size)
    {
#ifdef USE_HEAPSORT
      u = precond->edges[Bh.edges[i]].i;
      v = precond->edges[Bh.edges[i]].j;
#else
      u = precond->edges[i].i;
      v = precond->edges[i].j;
      if (u==v) continue;
#endif
      
      
      if (groups[u] == groups[v])
	{
#ifdef USE_HEAPSORT
	  weight = precond->edges[Bh.edges[i]].v;
#else
	  weight = precond->edges[i].v;
#endif
	  
	  edge_sign = (weight>0);
	  
	  x = find_set(u);
	  y = find_set(v);
	  
	  if (x!=y)
	    {
	      if (!((closed_cycle[x])&&(closed_cycle[y])))
		{
		  un_root = Union(u,v,x,y,edge_sign);
		  closed_cycle[un_root] = closed_cycle[x] | closed_cycle[y];
		}
	      else
		assert(0); /* this is a subgraph of a basis */
	    }
	  else
	    {
	      if ((edge_sign != (label[u]^label[v])) && (closed_cycle[x]==0))
		closed_cycle[x] = 1;
	    }
	  
	}
    }
  
  wtime = taucs_wtime();
#ifdef USE_HEAPSORT
  /*
    free_heap(h);
    size = heap_sort(Anent,&h,mtxA_tmp);
    if (size == -1) {
    free_graph(mtxA_tmp);
    free_graph(precond);
    free_linked_list(l);
    unionfind_free();
    taucs_free(diag);
    taucs_free(already_added);
    taucs_free(pi);
    taucs_free(groups);
    taucs_free(complete_subgraph);
    taucs_free(closed_cycle);
    return NULL;
    }
  */
  /*size = pqueue_fill(&Ah,mtxA_tmp);*/
  size = Ah.heap_size;
#else
  assert(Anent == mtxA_tmp->nent);
  size = mtxA_tmp->nent;
  graph_sort(mtxA_tmp);
#endif
  wtime = taucs_wtime() - wtime;
  taucs_printf("\t\tAMWB sort(%d) = %.3f seconds\n",Anent,wtime);
  wtime_sort += wtime;
  
  Do(i,size)
    {
#ifdef USE_HEAPSORT
      u = mtxA_tmp->edges[Ah.edges[i]].i;
      v = mtxA_tmp->edges[Ah.edges[i]].j;
#else
      u      = mtxA_tmp->edges[i].i;
      v      = mtxA_tmp->edges[i].j;
      if (u==v) continue;
#endif
      
      if (groups[u] == groups[v])
	{
#ifdef USE_HEAPSORT
	  weight = mtxA_tmp->edges[Ah.edges[i]].v;
#else
	  weight = mtxA_tmp->edges[i].v;
#endif
	  
	  edge_sign = (weight>0);
	  
	  x = find_set(u);
	  y = find_set(v);
	  
	  if (x!=y)
	    {
	      if (!((closed_cycle[x])&&(closed_cycle[y])))
		{
		  /*
		    ivec1_p[Bent] = u; 
		    ivec2_p[Bent] = v; 
		    dvec_p[Bent] = weight; 
		  */

		  precond->edges[Bent].i = u; 
		  precond->edges[Bent].j = v; 
		  precond->edges[Bent].v = weight; 
		  
		  Bent++; 

		  diag[u] += fabs(weight);
		  diag[v]+=fabs(weight);
		  
		  assert(complete_subgraph[groups[u]].completed_to_basis==0);
		  complete_subgraph[groups[u]].completed_to_basis = 1;
		  complete_subgraph[groups[u]].a = u;
		  complete_subgraph[groups[u]].b = v;
		  complete_subgraph[groups[u]].c = weight;

		  un_root = Union(u,v,x,y,edge_sign);
		  closed_cycle[un_root] = closed_cycle[x] | closed_cycle[y];
		}
	    }
	  else
	    {
	      if ((edge_sign != (label[u]^label[v])) && (closed_cycle[x]==0))
		{
		  /*
		    ivec1_p[Bent] = u; 
		    ivec2_p[Bent] = v; 
		    dvec_p[Bent] = weight; 
		  */
		  precond->edges[Bent].i = u; 
		  precond->edges[Bent].j = v; 
		  precond->edges[Bent].v = weight; 
		  
		  Bent++; 
		  
		  diag[u] += fabs(weight);
		  diag[v] += fabs(weight);
		  
		  assert(complete_subgraph[groups[u]].completed_to_basis==0);
		  complete_subgraph[groups[u]].completed_to_basis = 1;
		  complete_subgraph[groups[u]].a = u;
		  complete_subgraph[groups[u]].b = v;
		  complete_subgraph[groups[u]].c = weight;
		  
		  closed_cycle[x] = 1;
		}
	    }
	  
	}
    }
#ifdef USE_HEAPSORT
  /*    free_heap(h);*/
#endif
  taucs_free(closed_cycle);
  unionfind_free();
  
  /*** ALLOCATED: mtxA_tmp,diag,already_added,pi,precond,linked ***/
  /*** ALLOCATED: groups,complete_subgraph ***/
  
  precond->nent = Bent;
  
  wtime_component_bases = taucs_wtime() - wtime_component_bases;
  taucs_printf("\t\tAMWB component bases = %.3f seconds\n",wtime_component_bases);



  wtime_pair_bases = taucs_wtime();

  step2_Bent = Bent;
  
  /* COMPLETE EACH PAIR OF SUBGRAPHS INTO A BASIS */
  if (curr_group>1) {
    pairs = (six **)taucs_calloc(curr_group,sizeof(six *));
    closed_cycle = (byte *)taucs_calloc(n,sizeof(byte));
    if(!pairs || !closed_cycle) {
      free_graph(mtxA_tmp);
      free_graph(precond);
      free_linked_list(l);
      taucs_free(diag);
      taucs_free(already_added);
      taucs_free(pi);
      taucs_free(groups);
      taucs_free(complete_subgraph);
      taucs_free(pairs);
      taucs_free(closed_cycle);
    }

    /*** ALLOCATED: mtxA_tmp,diag,already_added,pi,precond,linked ***/
    /*** ALLOCATED: groups,complete_subgraph,pairs,closed_cycle ***/

    if (unionfind_init(n) == -1) {
      free_graph(mtxA_tmp);
      free_graph(precond);
      free_linked_list(l);
      taucs_free(diag);
      taucs_free(already_added);
      taucs_free(pi);
      taucs_free(groups);
      taucs_free(complete_subgraph);
      taucs_free(pairs);
      taucs_free(closed_cycle);
      return NULL;
    }

    /*** ALLOCATED: mtxA_tmp,diag,already_added,pi,precond,linked,UF ***/
    /*** ALLOCATED: groups,complete_subgraph,pairs,closed_cycle ***/

    wtime = taucs_wtime();
#ifdef USE_HEAPSORT
    /*
    size = heap_sort(basis_Bent,&h,precond);
    if (size == -1) {
      free_graph(mtxA_tmp);
      free_graph(precond);
      free_linked_list(l);
      unionfind_free();
      taucs_free(diag);
      taucs_free(already_added);
      taucs_free(pi);
      taucs_free(groups);
      taucs_free(complete_subgraph);
      taucs_free(pairs);
      taucs_free(closed_cycle);
      return NULL;
    }
    */
    size = pqueue_fill(&Bh,precond);

#else
    assert(basis_Bent == precond->nent);
    size = precond->nent;
    graph_sort(precond);
#endif
    wtime = taucs_wtime() - wtime;
    taucs_printf("\t\tAMWB sort(%d) = %.3f seconds\n",basis_Bent,wtime);
    wtime_sort += wtime;

    /*** ALLOCATED: mtxA_tmp,diag,already_added,pi,precond,linked,UF ***/
    /*** ALLOCATED: groups,complete_subgraph,pairs,closed_cycle,heap ***/

    Do(i,size)
      {
#ifdef USE_HEAPSORT
	u = precond->edges[Bh.edges[i]].i;
	v = precond->edges[Bh.edges[i]].j;
#else
	u = precond->edges[i].i;
	v = precond->edges[i].j;
	if (u==v) continue;
#endif

	if (groups[u] == groups[v])
	  {
#ifdef USE_HEAPSORT
	    weight = precond->edges[Bh.edges[i]].v;
#else
	    weight = precond->edges[i].v;
#endif
	    
	    edge_sign = (weight>0);
	    
	    x = find_set(u);
	    y = find_set(v);
	    
	    if (x!=y)
	      {
		if (!((closed_cycle[x])&&(closed_cycle[y])))
		  {
		    un_root = Union(u,v,x,y,edge_sign);
		    closed_cycle[un_root] = closed_cycle[x] | closed_cycle[y];
		  }
	      }
	    else
	      {
		if ((edge_sign != (label[u]^label[v])) && (closed_cycle[x]==0))
		  closed_cycle[x] = 1;
	      }
	    
	  }
      }

    wtime = taucs_wtime();
#ifdef USE_HEAPSORT
    /*
    free_heap(h);
    size = heap_sort(Anent,&h,mtxA_tmp);
    if (size == -1) {
      free_graph(mtxA_tmp);
      free_graph(precond);
      free_linked_list(l);
      unionfind_free();
      taucs_free(diag);
      taucs_free(already_added);
      taucs_free(pi);
      taucs_free(groups);
      taucs_free(complete_subgraph);
      taucs_free(pairs);
      taucs_free(closed_cycle);
      return NULL;
    }
    */
    /*size = pqueue_fill(&Ah,mtxA_tmp);*/
    size = Ah.heap_size;
#else
    assert(Anent == mtxA_tmp->nent);
    size = mtxA_tmp->nent;
    graph_sort(mtxA_tmp);
#endif
    wtime = taucs_wtime() - wtime;
    taucs_printf("\t\tAMWB sort(%d) = %.3f seconds\n",Anent,wtime);
    wtime_sort += wtime;

    Do(i,size)
      {
	int g1,g2;
	six *p;
	six *tmp;


#ifdef USE_HEAPSORT
	u = mtxA_tmp->edges[Ah.edges[i]].i;
	v = mtxA_tmp->edges[Ah.edges[i]].j;
#else
	u = mtxA_tmp->edges[i].i;
	v = mtxA_tmp->edges[i].j;
	if (u==v) continue;
#endif

	if (Bent == precond_maxnent)
	  {
	    precond_maxnent += chunk;
	    taucs_printf("adding space for %d entries in precond\n",chunk);fflush(stdout);
	    if (graph_resize(precond,precond_maxnent) == -1) {
	      int i_local;
	      Do(i_local,curr_group)
		free_linked_list_2(pairs[i_local]);
	      /* precond has not been freed */
	      free_graph(mtxA_tmp);
	      free_graph(precond);
	      free_linked_list(l);
	      unionfind_free();
#ifdef USE_HEAPSORT
	      free_heap(h);
#endif
	      taucs_free(diag);
	      taucs_free(already_added);
	      taucs_free(pi);
	      taucs_free(groups);
	      taucs_free(complete_subgraph);
	      taucs_free(pairs);
	      taucs_free(closed_cycle);
	      return NULL;
	    }
	    
	    /*
	    ivec1_p = precond->ivec1;
	    ivec2_p = precond->ivec2;
	    dvec_p = precond->dvec;
	    */
	  }

	g1 = min(groups[u],groups[v]); g2 = max(groups[u],groups[v]);
	
	if (g1 != g2)
	  {
	    p = pairs[(g1+g2)%curr_group];
	    while (p!=NULL)
	      {
		if ((p->group_1 == g1)&&(p->group_2 == g2))
		  goto after2;
		p = p->next;
	      }
	  after2:
	    
	    if (p == NULL)
	      {
		tmp = (six *)taucs_calloc(1,sizeof(six));
		if (tmp == NULL) {
		  int i_local;
		  Do(i_local,curr_group)
		    free_linked_list_2(pairs[i_local]);
		  free_graph(mtxA_tmp);
		  free_graph(precond);
		  free_linked_list(l);
		  unionfind_free();
#ifdef USE_HEAPSORT
		  free_heap(h);
#endif
		  taucs_free(diag);
		  taucs_free(already_added);
		  taucs_free(pi);
		  taucs_free(groups);
		  taucs_free(complete_subgraph);
		  taucs_free(pairs);
		  taucs_free(closed_cycle);
		  return NULL;
		}
		tmp->group_1 = g1;
		tmp->group_2 = g2;
		tmp->no_edges = 0;
		tmp->next = pairs[(g1+g2)%curr_group];
		pairs[(g1+g2)%curr_group] = tmp;
		p = tmp;
	      }

	    /* first check if the next edges to be added are inner to one of the subgraphs
	       (not cross edges) */
	    if (p->no_edges < 2)
	      {
#ifdef USE_HEAPSORT
		weight = mtxA_tmp->edges[Ah.edges[i]].v;
#else
		weight = mtxA_tmp->edges[i].v;
#endif
		
		edge_sign = (weight>0);
		
		x = find_set(u);
		y = find_set(v);
		closed_cycle_x = closed_cycle[x];
		closed_cycle_y = closed_cycle[y];

		if (complete_subgraph[g1].completed_to_basis)
		  if (fabs(weight)<=fabs(complete_subgraph[g1].c))
		    if ((p->no_edges==0)||
			((groups[p->a[0]]!=g1)||
			 (groups[p->b[0]]!=g1)))
		    {
		      p->a[p->no_edges] = complete_subgraph[g1].a;
		      p->b[p->no_edges] = complete_subgraph[g1].b;
		      p->c[p->no_edges] = complete_subgraph[g1].c;
		      p->cross[p->no_edges] = 0;
		      (p->no_edges)++;
		    }

		if (p->no_edges < 2)		
		  if (complete_subgraph[g2].completed_to_basis)
		    if (fabs(weight)<=fabs(complete_subgraph[g2].c))
		      if ((p->no_edges==0)||
			  ((groups[p->a[0]]!=g2)||
			   (groups[p->b[0]]!=g2)))
			{
			  p->a[p->no_edges] = complete_subgraph[g2].a;
			  p->b[p->no_edges] = complete_subgraph[g2].b;
			  p->c[p->no_edges] = complete_subgraph[g2].c;
			  p->cross[p->no_edges] = 0;
			  (p->no_edges)++;
			}
		
		/* if p->no_edges == 2 get out of if */

		if (p->no_edges == 1)
		  {
		    if (groups[p->a[0]] == groups[p->b[0]])
		      {
			if (groups[p->a[0]] == groups[u])
			  closed_cycle_x = 1;
			else
			  closed_cycle_y = 1;
		      }
		      
		  }

		
		if (p->no_edges == 0) {
		  /* sivan: added this brace below to avoid warning. */
		  /* I have no idea why there are two identical ifs. */
		  /* I hope I added the matching brace in the right place */
		  if (p->no_edges == 0) { 
		    if (x != y) {
		      if (!((closed_cycle[x])&&(closed_cycle[y]))) {
#ifdef USE_HEAPSORT
			if (already_added[Ah.edges[i]]==0)
#else
			  if (already_added[i]==0)
#endif
			    {
			      /*
				ivec1_p[Bent] = u; 
				ivec2_p[Bent] = v; 
				dvec_p[Bent] = weight; 
			      */
			      precond->edges[Bent].i = u; 
			      precond->edges[Bent].j = v; 
			      precond->edges[Bent].v = weight; 
			      
			      Bent++; 
			      
			      diag[u] += fabs(weight);
			      diag[v]+=fabs(weight);
			    }
			p->a[p->no_edges] = u;
			p->b[p->no_edges] = v;
			p->c[p->no_edges] = weight;
			p->cross[p->no_edges] = 1;
			(p->no_edges)++;
		      }
		    } else {
		      if ((edge_sign != (label[u]^label[v])) && (closed_cycle[x]==0)) {
#ifdef USE_HEAPSORT
			if (already_added[Ah.edges[i]]==0)
#else
			  if (already_added[i]==0)
#endif
			    {
			      /*
				ivec1_p[Bent] = u; 
				ivec2_p[Bent] = v; 
				dvec_p[Bent] = weight; 
			      */
			      precond->edges[Bent].i = u; 
			      precond->edges[Bent].j = v; 
			      precond->edges[Bent].v = weight; 
			      
			      Bent++; 
			      
			      diag[u] += fabs(weight);
			      diag[v] += fabs(weight);
			    }
			p->a[p->no_edges] = u;
			p->b[p->no_edges] = v;
			p->c[p->no_edges] = weight;
			p->cross[p->no_edges] = 1;
			(p->no_edges)++;
		      }
		    }
		  }
		}

		      if (p->no_edges == 1) {
			if (p->cross[0]==0) {
			  if (groups[p->a[0]]==groups[u])
			    closed_cycle_x = 1;
			  else
			    closed_cycle_y = 1;
			}
		    
		    if ((!((closed_cycle_x)&&(closed_cycle_y)))&&(p->cross[0]==0))
		      {
#ifdef USE_HEAPSORT
			if (already_added[Ah.edges[i]]==0)
#else
			if (already_added[i]==0)
#endif
			  {
			    /*
			    ivec1_p[Bent] = u; 
			    ivec2_p[Bent] = v; 
			    dvec_p[Bent] = weight; 
			    */
			    precond->edges[Bent].i = u;
			    precond->edges[Bent].j = v;
			    precond->edges[Bent].v = weight;

			    Bent++; 

			    diag[u] += fabs(weight);
			    diag[v]+=fabs(weight);
			  }
			p->a[p->no_edges] = u;
			p->b[p->no_edges] = v;
			p->c[p->no_edges] = weight;
			p->cross[p->no_edges] = 1;
			(p->no_edges)++;
		      }
		    
		    if ((!closed_cycle_x)&&(!closed_cycle_y)&&(p->cross[0]==1))
		      {
			int x1,y1,edge_sign1;
			x1 = find_set(p->a[0]);
			y1 = find_set(p->b[0]);
			edge_sign1 = (p->c[0]>0);			    
			
			if ((edge_sign^edge_sign1^label[u]^label[v]^label[p->a[0]]^label[p->b[0]]) ==1)
			  {
#ifdef USE_HEAPSORT
			    if (already_added[Ah.edges[i]]==0)
#else
			    if (already_added[i]==0)
#endif
			      {
				/*
				ivec1_p[Bent] = u; 
				ivec2_p[Bent] = v; 
				dvec_p[Bent] = weight; 
				*/

				precond->edges[Bent].i = u; 
				precond->edges[Bent].j = v; 
				precond->edges[Bent].v = weight; 

				Bent++; 

				diag[u] += fabs(weight);
				diag[v]+=fabs(weight);
			      }
			    p->a[p->no_edges] = u;
			    p->b[p->no_edges] = v;
			    p->c[p->no_edges] = weight;
			    p->cross[p->no_edges] = 1;
			    (p->no_edges)++;
			  }
		      }
		  }
	      }
	  }
      }

    Do(i,curr_group)
      free_linked_list_2(pairs[i]);

    taucs_free(pairs);
#ifdef USE_HEAPSORT
    /*    free_heap(h);*/
    free_heap(Ah);
    free_heap(Bh);
#endif
    taucs_free(closed_cycle);
    unionfind_free();
    
    precond->nent = Bent;

    /*** ALLOCATED: mtxA_tmp,diag,already_added,pi,precond,linked ***/
    /*** ALLOCATED: groups,complete_subgraph ***/
  }

  taucs_free(complete_subgraph);
  taucs_free(already_added);

  wtime_pair_bases = taucs_wtime() - wtime_pair_bases;
  taucs_printf("\t\tAMWB pair bases = %.3f seconds\n",wtime_pair_bases);


  /*** ALLOCATED: mtxA_tmp,diag,pi,precond,linked ***/
  /*** ALLOCATED: groups ***/

  /* allocate more memory to the preconditioner if needed */

  wtime = taucs_wtime();
  if (precond_maxnent < Bent + n)
    {
      taucs_printf("adding space for %d entries in precond for diagonal entries\n",Bent+n-precond_maxnent);
      precond_maxnent = Bent + n;


      if (graph_resize(precond,precond_maxnent) == -1) {
	/* precond has not been freed */
	free_graph(mtxA_tmp);
	free_graph(precond);
	free_linked_list(l);
	taucs_free(diag);
	taucs_free(pi);
	taucs_free(groups);
	return NULL;
      }

      precond->nent = Bent;
      /*
      ivec1_p = precond->ivec1 ;
      ivec2_p = precond->ivec2 ;
      dvec_p = precond->dvec ;
      */
    }

  Do(i,n)
    {
      /*
      ivec1_p[Bent] = i;
      ivec2_p[Bent] = i;
      dvec_p[Bent] = diag[i];
      */

      precond->edges[Bent].i = i; 
      precond->edges[Bent].j = i; 
      precond->edges[Bent].v = diag[i]; 

      Bent++;
    }
  precond->nent = Bent;

  wtime = taucs_wtime() - wtime;
  taucs_printf("\t\tAMWB precond resize = %.3f seconds\n",wtime);


  taucs_printf("actual number of entries in preconditioner = %d\n",Bent);fflush(stdout);

  taucs_free(diag);
  taucs_free(groups);
  taucs_free(pi);
  free_linked_list(l);
  free_graph(mtxA_tmp);

  /*** ALLOCATED: precond ***/

  wtime = taucs_wtime();
  out = graph_to_ccs_matrix(precond);
  if (!out) {
    free_graph(precond);
    return NULL;
  }
  wtime = taucs_wtime() - wtime;
  taucs_printf("\t\tAMWB graph-to-matrix = %.3f seconds\n",wtime);

  free_graph(precond);

  
  wtime_total = taucs_wtime() - wtime_total;
  taucs_printf("\t\tAMWB time = %.3f seconds (%.3f sort)\n",
	       wtime_total,
	       wtime_sort);
  
  return out;
}

/*********************************************************/
/* MST-specific routines                                 */
/*********************************************************/

typedef struct msthea {
  int     heap_size;
  int*    vertices;
  double* key;
} mstheap;

#define INF 100000000.0

#define Parent(i) ((((i)+1)/2) - 1)
#define Left(i) ((((i)+1)*2) - 1)
#define Right(i) ((((i)+1)*2 + 1) - 1)

static
void mstheap_exchange(mstheap A,int a,int b,int *point_to_heap)
{
  int tmp1;
  double tmp2;

  tmp1 = A.vertices[a];
  A.vertices[a] = A.vertices[b];
  A.vertices[b] = tmp1;

  tmp2 = A.key[a];
  A.key[a] = A.key[b];
  A.key[b] = tmp2;

  point_to_heap[A.vertices[a]] = a;
  point_to_heap[A.vertices[b]] = b;
}

static
void Mstheapify(mstheap A,int i,int *point_to_heap)
{
  int l,r,largest;
  
  l = Left(i);
  r = Right(i);
  
  if ((l < A.heap_size) && (A.key[l] > A.key[i]))
    largest = l;
  else
    largest = i;

  if ((r < A.heap_size) && (A.key[r] > A.key[largest]))
    largest = r;

  if (largest != i)
    {
      mstheap_exchange(A,i,largest,point_to_heap);
      Mstheapify(A,largest,point_to_heap);
    }
}

static
int build_mstheap(int r,int size,mstheap *h,int *point_to_heap,char alloc_flag)
{
  int i;

  h->heap_size = size;

  if (alloc_flag)
    {
      h->vertices = (int *)taucs_malloc(size * sizeof(int));
      h->key = (double *)taucs_malloc(size * sizeof(double));
      
      if (h->key == NULL || h->vertices == NULL) return -1;
    }

  Do(i,size)
    {
      h->vertices[i] = i;
      h->key[i] = (- INF);
      point_to_heap[i] = i;
    }
  
  h->key[r] = 0;

  mstheap_exchange((*h),r,0,point_to_heap);
  
  return 0;
}

static
int Mstheap_Extract_Max(mstheap *A,int *point_to_heap)
{
  int out;

  assert (A->heap_size >= 1);
  
  out = A->vertices[0];
  A->vertices[0] = A->vertices[A->heap_size - 1];
  A->key[0] = A->key[A->heap_size - 1];
  point_to_heap[A->vertices[0]] = 0;
  
  A->heap_size --;
  
  Mstheapify((*A),0,point_to_heap);
  
  return(out);
}

static
void mstheap_increase_key(mstheap h,int vv,double val,int *point_to_heap)
{
  int i;
  /*int count=0;*/
  double key;
  int ver, v;
  
  v = point_to_heap[vv];
  
  h.key[v] = val;

  key = h.key[v];
  ver = h.vertices[v];

  i = v;

  while((i>0) && (h.key[Parent(i)] < key))
    {
      h.key[i]      = h.key[Parent(i)];
      h.vertices[i] = h.vertices[Parent(i)];
      point_to_heap[h.vertices[i]] = i;
      i = Parent(i);
    }
  
  h.key[i] = key;
  h.vertices[i] = ver;
  point_to_heap[h.vertices[i]] = i;  

}

static
void free_mstheap(mstheap h)
{
  taucs_free(h.vertices);
  taucs_free(h.key);
}

static
int add_heavy_edges(graph* mtxA,
		    graph* precond,
		    int Bent,
		    edge **array,
		    int *groups,
		    int no_groups,
		    int *pi,
		    double *diag)
{
  three **pairs;
  int i,j,a,b,k;
  double w;
  /*
  int *ivec1_p, *ivec2_p ;
  double *dvec_p; 
  int *ivec1, *ivec2 ;
  double *dvec; 
  */
  int orig_Bent;
  three *p,*tmp;
  int n,nent;
  int precond_maxnent,chunk;

  three *pool;
  int   next_in_pool;

  n = mtxA->n;
  nent = mtxA->nent;

  precond_maxnent = precond->max_size;
  chunk = max(min((no_groups*(no_groups-1)/2)/10,5000),10000);
  orig_Bent = Bent;
  
  pairs = (three **) taucs_calloc(no_groups,sizeof(three *));
  pool  = (three*)   taucs_malloc((n+nent) * sizeof(three));
  if (!pairs || !pool) {
    taucs_free(pairs);
    taucs_free(pool);
    return -1;
  }
  next_in_pool = 0;

  /*
  ivec1 = mtxA->ivec1;
  ivec2 = mtxA->ivec2;
  dvec = mtxA->dvec;
  */

  Do(k,n) {
    i = k;
    j = pi[k];
    
    if (j != (-1)) {
      a = min(groups[i],groups[j]);
      b = max(groups[i],groups[j]);
      
      if (a != b) {
	p = pairs[(a+b)%no_groups];
	while (p!=NULL) {
	  if ((p->group_1 == a)&&(p->group_2 == b))
	    goto after;
	  p = p->next;
	}
      after:
	if (p == NULL) {
	  /*tmp = (three *)taucs_malloc(sizeof(three));*/
	  /*if (tmp == NULL)  {taucs_printf("ERROR! OUT OF MEMORY\n");exit(234);}*/
	  assert(next_in_pool < n+nent);
	  tmp = pool + next_in_pool; 
	  next_in_pool ++;
	  tmp->group_1 = a;
	  tmp->group_2 = b;
	  tmp->already_connected = 1;
	  tmp->next = pairs[(a+b)%no_groups];
	  pairs[(a+b)%no_groups] = tmp;
	}
      }
    }
  }

  Do(k,nent) {
    i = (mtxA->edges)[k].i;
    j = (mtxA->edges)[k].j;

    a = min(groups[i],groups[j]);
    b = max(groups[i],groups[j]);
    
    if (a != b)	{
      w = - (mtxA->edges)[k].v;
      if (w) {
	p = pairs[(a+b)%no_groups];
	while (p!=NULL) {
	  if ((p->group_1 == a)&&(p->group_2 == b)) {
	    if (p->already_connected==0)
	      if (w > p->c) {
		p->a = i;
		p->b = j;
		p->c = w;
	      }
	    goto after1;
	  }
	  p = p->next;
	}
      after1:
	if (p == NULL) {
 	  /*tmp = (three *)taucs_malloc(sizeof(three));*/
	  /*if (tmp == NULL)  {taucs_printf("ERROR! OUT OF MEMORY\n");exit(234);}*/
	  assert(next_in_pool < n+nent);
	  tmp = pool + next_in_pool; 
	  next_in_pool ++;
	  tmp->group_1 = a;
	  tmp->group_2 = b;
	  tmp->a = i;
	  tmp->b = j;
	  tmp->c = w;
	  tmp->already_connected = 0;
	  tmp->next = pairs[(a+b)%no_groups];
	  pairs[(a+b)%no_groups] = tmp;
	}
      }
    }
  }
  
  /*
  ivec1_p = precond->ivec1;
  ivec2_p = precond->ivec2;
  dvec_p = precond->dvec;
  */
  
  Do(i,no_groups) {
    p = pairs[i];
    while(p!=NULL) {
      if(p->already_connected == 0) {
	if (p->a > p->b)
	  swap(p->a,p->b);
	
	(precond->edges)[Bent].i = p->a;
	(precond->edges)[Bent].j = p->b;
	(precond->edges)[Bent].v = -(p->c);

	Bent++;
	
	if (Bent == precond_maxnent) {
	  precond_maxnent += chunk;
	  taucs_printf("adding space for %d entries in precond\n",chunk);
	  graph_resize(precond,precond_maxnent);
	  precond->nent = orig_Bent;
		
	  /*  
	  ivec1_p = precond->ivec1;
	  ivec2_p = precond->ivec2;
	  dvec_p = precond->dvec;
	  */
	  
	}
	diag[p->a] -= (-p->c);
	diag[p->b] -= (-p->c);
      }

      p = p->next;
    }
  }
  
  /*Do(i,no_groups) free_linked_list_2(pairs[i]);*/
  taucs_free(pool);
  taucs_free(pairs);

  return(Bent);
}

static int Dijkstra(graph *mtxA,int r,int *pi,linked *l,int *d,int *maxdist,int *partition,double min,double y,int j);


static int Av_Part_W(graph *mtxA,int *partition,int *new_partition,int *parts,graph *out,int nparts)
{
  /* Peleg, Noga et al.'s partition. From Peleg's book, page 217 */
   
  linked *l = NULL, *l_c = NULL;
  int *pi=NULL,*d=NULL,i,k,*findrho=NULL,minrho,maxdist,classes,n,nent,root,j,curr_partition=0;
  int row, col;
  int *pi1 = 0; /* warning */
  double x, y, min, max, not;
  byte bool=1;
  edge *p,*dummy, *pe ,*max_pe;
  int count = 0;
 
  n = mtxA->n;
  nent = mtxA->nent;
 
  x = exp(sqrt(log(n)*log(log(n))))/3;
  y = x * 9*log(n) * (floor(3*log(n)/log(x))+1);

  pi      = (int *)taucs_malloc(n*sizeof(int));
  d       = (int *)taucs_malloc(n*sizeof(int));
  l       = create_linked_list(mtxA,n,mtxA->nent,&min,&max);
  if (!pi || !d || !l)
    goto exit_Av_Part_W;
  
  classes = (int)(log(max/min)/log(y))+1;

  for(i=0;i<n;i++)
    new_partition[i] = -1;

  j = 1;

  while (count < n)
    {
      root = rand() % n;
      if (new_partition[root] == -1)
	{
	  for(i=0;i<n;i++)
	    d[i] = -1;
	  Dijkstra(mtxA,root,pi,l,d,&maxdist,partition,min,y,j);
	  
	  findrho     = (int *)taucs_calloc((maxdist+1)*classes,sizeof(int));
	  if (!findrho)
	    goto exit_Av_Part_W;
	  
	  for(i=0;i<n;i++)
	    {
	      if (d[i] != -1)
		{
		  p = (l->point)[i];
		  while (p != NULL)
		    {
		      if ((d[mtxA->edges[p->entry_no].i]!=-1) && (d[mtxA->edges[p->entry_no].j]!=-1))
			if ((d[mtxA->edges[p->entry_no].i] - d[mtxA->edges[p->entry_no].j] >= -1) &&
			    (d[mtxA->edges[p->entry_no].i] - d[mtxA->edges[p->entry_no].j] <= 1))
			  findrho[(max(d[mtxA->edges[p->entry_no].i],d[mtxA->edges[p->entry_no].j]))*classes+
				 /*(int)(log(abs(mtxA->edges[p->entry_no].v)/min)/log(y))]++; omer*/
				 (int)(log(abs((int)(mtxA->edges[p->entry_no].v))/min)/log(y))]++;
		      p = p->next;
		    }
		}
	    }

	  for(i=0;i<min(j,classes);i++)
	    findrho[i] = 0; /* ignore edges connecting two vertices whose distance is 0 */
	
	  /* At this point, findrho[k,i], or findrho[k*classes+i], contains the number of edges in E_i,
	     connecting two vertices whose distance is k, or a vertex of distance k with a vertex of distance k-1 */
	  
	  for(k=1;k<maxdist;k++)
	    for(i=0;i<min(j,classes);i++)
	      findrho[k*classes+i] += findrho[(k-1)*classes+i];
	  
	  /* At this point, findrho[k,i], or findrho[k*classes+i], contains the number of edges in E_i,
	     connecting two vertices whose distance is j, or a vertex of distance j with a vertex of distance j-1
	     for j=1,...,k */
	  
	  for(minrho=1;minrho<maxdist;minrho++)
	    {
	      bool = 1;
	      for(k=0;k<min(j,classes);k++)
		{
		  if ((double)(findrho[(minrho+1)*classes+k]-findrho[minrho*classes+k]) > (findrho[minrho*classes+k])/x)
		    bool = 0;
		}
	      if (bool)
		goto afterr;
	    }
	
	afterr:
	  if (bool)
	    {
	      for(i=0;i<n;i++)
		if ((d[i] <= minrho) && (d[i] != -1) )
		  if (new_partition[i] == -1)
		    {
		      count ++;
		      new_partition[i] = curr_partition;
		    }
	    }
	  else
	    {
	      for(i=0;i<n;i++)
		if ((new_partition[i] == -1) && (d[i] != -1))
		  {
		    count ++;
		    new_partition[i] = curr_partition;
		  }
	    }
	  
	  for(i=0;i<n;i++)
	    {
	      if (new_partition[i] == curr_partition)
		l->point[i] = NULL;
	      else
		{
		  p = l->point[i];
		  l->point[i] = NULL;
		  while (p != NULL)
		    {
		      if ((new_partition[mtxA->edges[p->entry_no].i] != curr_partition) && (new_partition[mtxA->edges[p->entry_no].j] != curr_partition))
			{
			  dummy = l->point[i];
			  l->point[i] = p;
			  p = p->next;
			  (l->point[i])->next = dummy;
			}
		      else
			p = p->next;
		    }
		}
	    }
	  
	  curr_partition ++;
	  j++;
	  taucs_free(findrho);
	}
    }

  *parts = curr_partition;

  l_c = create_linked_list_cluster(mtxA,nparts,mtxA->nent,&not,&not,partition,new_partition);
  pi1 = (int *)taucs_malloc(nparts*sizeof(int));

  if (!l_c || !pi1)
    goto exit_Av_Part_W;
  
  for(i=0;i<n;i++)
    if (pi[i] != -1)
      {
	if (partition[i] != partition[pi[i]])
	  pi1[partition[i]] = partition[pi[i]];
      }
    else
      pi1[partition[i]] = -1;

  Do(i,nparts) {
    row = i;
    col = pi1[i];
    if (col != (-1)) {
      pe = l_c->point[row];
      max_pe = NULL;
      while (pe != NULL) {
	if (   (partition[(mtxA->edges)[pe->entry_no].j] == col) 
	       || (partition[(mtxA->edges)[pe->entry_no].i] == col))
	  {
	    if (!max_pe)
	      max_pe = pe;
	    else
	      if (-mtxA->edges[pe->entry_no].v > -mtxA->edges[max_pe->entry_no].v)
		max_pe = pe;
	  }
	pe = pe->next;
      }
      
      assert(max_pe);
      out->edges[out->nent].i = (mtxA->edges)[max_pe->entry_no].i;
      out->edges[out->nent].j = (mtxA->edges)[max_pe->entry_no].j;
      out->edges[out->nent].v = (mtxA->edges)[max_pe->entry_no].v;
      out->nent++;
      
    }
  }


  taucs_free(pi);taucs_free(d);free_linked_list(l);free_linked_list(l_c);taucs_free(pi1);
  return 1;
  
 exit_Av_Part_W:
  taucs_free(pi);taucs_free(d);free_linked_list(l);taucs_free(l);taucs_free(findrho);free_linked_list(l_c);taucs_free(l_c);taucs_free(pi1);
  return 0;

  
}

static taucs_ccs_matrix *amst_preconditioner_create(graph *mtxA, double* diag,int rnd,double subgraphs,int stretch_flag);
static int Prim(graph *mtxA,int r,int *pi,int *d,linked *l);
/*static int Prim_cluster(graph *mtxA,int r,int *pi,linked *l,int *partition,int *new_partition,int nparts,int *point_to_heap,char *in_Q,mstheap h);*/

static int Dijkstra(graph *mtxA,int r,int *pi,linked *l,int *d,int *maxdist,int *partition,double min,double y,int j)
{
  /* Dijkstra is used in order to compute the distance of vertices from the root.
     The distance of two vertices within the same partition is 0 */

  int n,entry,u,v,i;
  mstheap h;
  int *point_to_heap = NULL;
  char *in_Q = NULL;
  int size_of_Q;
  edge *p;
  double weight;

  *maxdist=0;
  
  n = mtxA->n;

  point_to_heap = (int*)  taucs_malloc(n*sizeof(int));
  in_Q          = (char*) taucs_malloc(n*sizeof(char));

  Do(i,n) in_Q[i] = 1;
  size_of_Q = n;
  
  /* Prim's Algorithm - Introduction to Algorithms page 509 */
  
  if (build_mstheap(r,n,&h,point_to_heap,1) == -1) 
    goto exit_Dijkstra;

  pi[r] = -1;
  
  while(size_of_Q > 0) {
    if (h.key[0] == -INF)
      goto after_Dijkstra;
    assert(h.key[0] != -INF);
    weight = -h.key[0];
    u = Mstheap_Extract_Max(&h,point_to_heap);
    d[u] = (int)weight;
    if (weight>*maxdist)
      *maxdist = (int)weight;

    in_Q[u] = 0;
    size_of_Q --;

    p = (l->point)[u];
    while (p != NULL) {
      entry = p->entry_no;
      /*if ((int)(log(abs(mtxA->edges[entry].v/min))/log(y)) < j) omer*/
			if ((int)(log(abs((int)(mtxA->edges[entry].v/min)))/log(y)) < j)
	{      
	  /* v belongs to Adj[u] */

	  if ((mtxA->edges)[entry].j != u)
	    v = (mtxA->edges)[entry].j;
	  else
	    v = (mtxA->edges)[entry].i;
	  
	  if (in_Q[v] && (-h.key[point_to_heap[v]] > ((partition[v]==partition[u])?0:1) + weight))
	    {
	      pi[v] = u;
	      mstheap_increase_key(h,v,-(((partition[v]==partition[u])?0:1) + weight),point_to_heap);
	    }
	}
      p = p->next;
    }
  }

 after_Dijkstra:
  free_mstheap(h);
  taucs_free(point_to_heap);
  taucs_free(in_Q);
  return 1;

 exit_Dijkstra:
  free_mstheap(h);
  taucs_free(point_to_heap);
  taucs_free(in_Q);
  return 0;
  
}

void stupid_part(int *partition,int n,int j,int *nparts)
{
  int i,k,q;
  
  k = 1<<j;
  
  q = ((n%k == 0)?(n/k):(n/k+1));
  
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      {
	partition[i*n+j] = q*(i/k)+(j/k);
      }
  *nparts=partition[n*n-1]+1;

  /* for(i=0;i<n*n;i++) */
    /* printf("ASD %d %d - %d %d\n",i,partition[i],n,k); */
  /* exit(345); */
}

graph *low_stretch(graph *mtxA)
{
  int *partition = NULL,*new_partition = NULL,*choose_root = NULL;
  int i,n,nparts,new_nparts;/* r omer*/
  /*double dummy; omer*/
  int *pi = NULL;
  /*linked *l = NULL;*/
  /*int k=0,j=0;*/ /* row,col, omer*/
  /*edge *pe,*max_pe; omer*/
  graph *out;
  /*int *point_to_heap = NULL;*/
  /*char *in_Q = NULL;*/
  /*mstheap h; omer*/

  n = mtxA->n;
  
  out = construct_graph(n-1);
  partition = (int *)taucs_malloc(n*sizeof(int));
  new_partition = (int *)taucs_malloc(n*sizeof(int));
  choose_root = (int *)taucs_malloc(n*sizeof(int));
  pi = (int *)taucs_malloc(n*sizeof(int));

  if (!out || !partition || !new_partition || !choose_root || !pi)
    goto exit_low_stretch;

  out->n = n;
  out->nent = 0;
  
  for (i=0;i<n;i++)
    partition[i] = i;

  nparts = n;

  while (nparts > 1)
    {
      if (Av_Part_W(mtxA,partition,new_partition,&new_nparts,out,nparts) == 0)
	goto exit_low_stretch;

#if 0
      for(i=0;i<n;i++)
	pi[i] = -1;

      out->n = n;
      out->nent = n-1;

      j++;

      /* stupid_part(new_partition,(int)(sqrt(n)),j,&new_nparts); */

      for(i=0;i<new_nparts;i++)
	choose_root[i] = 0;
      
      l = create_linked_list_cluster(mtxA,nparts,mtxA->nent,&dummy,&dummy,partition,new_partition);
      /* This linked list contains all the edges which connect vertices whose endpoints are
	 in different sections in partition, but in the same section in new_partition.
	 This will help us build trees within each section in new_partition. */

      if (l == NULL)
	goto exit_low_stretch;


      point_to_heap = (int *)taucs_malloc(n*sizeof(int));
      in_Q = (char *)taucs_malloc(n*sizeof(char));
      h.vertices = NULL;
      h.key = NULL;
      h.vertices = (int *)taucs_malloc(n * sizeof(int));
      h.key = (double *)taucs_malloc(n * sizeof(double));

      if (!point_to_heap || !in_Q || !h.vertices || !h.key)
	{
	  taucs_free(point_to_heap);taucs_free(in_Q);taucs_free(h.vertices);taucs_free(h.key);
	  goto exit_low_stretch;
	}

      for(i=0;i<n;i++)
	{
	  if (choose_root[new_partition[i]] == 0)
	    {
	      /* new_partition[i] is a section for which a tree has not yet been found */
	      choose_root[new_partition[i]] = 1;
	      r = partition[i];
	      if (Prim_cluster(mtxA,r,pi,l,partition,new_partition,nparts,point_to_heap,in_Q,h) == 0)
		{
		  free_linked_list(l);
		  goto exit_low_stretch;
		}
	    }
	}
      
      taucs_free(point_to_heap);
      taucs_free(in_Q);
      free_mstheap(h);

      Do(i,nparts) {
	row = i;
	col = pi[i];
	if (col != (-1)) {
	  pe = l->point[row];
	  max_pe = NULL;
	  while (pe != NULL) {
	    if (   (partition[(mtxA->edges)[pe->entry_no].j] == col) 
		   || (partition[(mtxA->edges)[pe->entry_no].i] == col))
	      {
		if (!max_pe)
		  max_pe = pe;
		else
		  if (-mtxA->edges[pe->entry_no].v > -mtxA->edges[max_pe->entry_no].v)
		    max_pe = pe;
	      }
	    pe = pe->next;
	  }
	  
	  assert(max_pe);
	  out->edges[k].i = (mtxA->edges)[max_pe->entry_no].i;
	  out->edges[k].j = (mtxA->edges)[max_pe->entry_no].j;
	  out->edges[k].v = (mtxA->edges)[max_pe->entry_no].v;
	  k++;
	  
	}
      }

      out->nent = k;
      free_linked_list(l);
#endif

      for(i=0;i<n;i++)
	partition[i] = new_partition[i];
      
      nparts = new_nparts;

    }
  
  assert(out->nent==(n-1)); /* helps verify that out is a tree */

  taucs_free(partition);taucs_free(new_partition);taucs_free(choose_root);taucs_free(pi);
  return out;
  
 exit_low_stretch:
  free_graph(out);taucs_free(partition);taucs_free(new_partition);taucs_free(choose_root);taucs_free(pi);
  return 0;
}

static
int Prim(graph *mtxA,int r,int *pi,int *d,linked *l)
{
  /* Prim's Algorithm - Introduction to Algorithms page 509 */
  
  int n,entry,u,v,i;
  mstheap h;
  int *point_to_heap;
  char *in_Q;
  int size_of_Q;
  edge *p;

  n = mtxA->n;

  point_to_heap = (int*)  taucs_malloc(n*sizeof(int));
  in_Q          = (char*) taucs_malloc(n*sizeof(char));

  Do(i,n) in_Q[i] = 1;
  size_of_Q = n;
  
  /* Prim's Algorithm - Introduction to Algorithms page 509 */

  if ((build_mstheap(r,n,&h,point_to_heap,1) == -1)) {
    taucs_free(in_Q);
    taucs_free(pi);
    taucs_free(point_to_heap);
    /* free linked_list; */
    return 0;
  }

  pi[r] = -1;
  d[r] = 0;
  
  while(size_of_Q > 0) {
    u = Mstheap_Extract_Max(&h,point_to_heap);
    in_Q[u] = 0;
    size_of_Q --;

    p = (l->point)[u];
    while (p != NULL) {
      entry = p->entry_no;
      /* v belongs to Adj[u] */

      if ((mtxA->edges)[entry].j != u)
	v = (mtxA->edges)[entry].j;
      else
	v = (mtxA->edges)[entry].i;
      
      if (in_Q[v] && ((-((mtxA->edges)[entry].v)) > h.key[point_to_heap[v]])) {
	pi[v] = u;
	d[v] = d[u] + 1;
	mstheap_increase_key(h,v,-((mtxA->edges)[entry].v),point_to_heap);
      }
      p = p->next;
    }
  }

  free_mstheap(h);
  taucs_free(point_to_heap);
  taucs_free(in_Q);
  return 1;
}

#if 0
static
int Prim_cluster(graph *mtxA,int r,int *pi,linked *l,int *partition,int *new_partition,int nparts,int *point_to_heap,char *in_Q,mstheap h)
{
  /* Prim's Algorithm - Introduction to Algorithms page 509 */
  
  int n,entry,u,v,i;
  /* mstheap h; */
  /* int *point_to_heap = NULL; */
  /* char *in_Q = NULL; */
  int size_of_Q;
  edge *p;

  n = mtxA->n;

  /* point_to_heap = (int*)  taucs_malloc(nparts*sizeof(int)); */
  /* in_Q          = (char*) taucs_malloc(nparts*sizeof(char)); */
  /* if (!point_to_heap || !in_Q) */
    /* goto exit_Prim_cluster; */

  Do(i,nparts) in_Q[i] = 1;
  size_of_Q = nparts;
  
  /* Prim's Algorithm - Introduction to Algorithms page 509 */

  if (build_mstheap(r,nparts,&h,point_to_heap,0) == -1)
    goto exit_Prim_cluster;

  pi[r] = -1;
  
  while(size_of_Q > 0) {
    if (h.key[0] == -INF)
      goto after_Prim_cluster;
    u = Mstheap_Extract_Max(&h,point_to_heap);

    in_Q[u] = 0;
    size_of_Q --;

    p = (l->point)[u];
    while (p != NULL) {
      entry = p->entry_no;
      /* v belongs to Adj[u] */

      if (partition[(mtxA->edges)[entry].j] != u)
	v = partition[(mtxA->edges)[entry].j];
      else
	v = partition[(mtxA->edges)[entry].i];

      if (in_Q[v] && ((-((mtxA->edges)[entry].v)) > h.key[point_to_heap[v]])) {
	pi[v] = u;
	mstheap_increase_key(h,v,-((mtxA->edges)[entry].v),point_to_heap);
      }
      p = p->next;
    }
  }

 after_Prim_cluster:
  /* free_mstheap(h); */
  /* taucs_free(point_to_heap); */
  /* taucs_free(in_Q); */
  return 1;

 exit_Prim_cluster:
  /* free_mstheap(h); */
  /* taucs_free(point_to_heap); */
  /* taucs_free(in_Q); */
  return 0;
}
#endif /* 0, we don't need this routine */

static double dist(int i, int j, double w, graph *mtxA,int *pi,int *d,linked *l,double *dilation,double *congestion)
{
  double out=0;
  int tmp;
  edge *e;
  double q = 0;
  
  if (d[i] < d[j])
    {tmp = i;i=j;j=tmp;}

  /* now we know that d[i] >= d[j] */
  
  while (d[i] > d[j])
    {
      e = l->point[i];
      while ((mtxA->edges[e->entry_no].i != pi[i]) && (mtxA->edges[e->entry_no].j != pi[i]))
	e = e->next;
      out += mtxA->edges[e->entry_no].v;
      q += mtxA->edges[e->entry_no].v;
      congestion[i] += mtxA->edges[e->entry_no].v/w;
      i = pi[i];
    }
  
  /* now we know that d[i] == d[j] */
  
  while (i != j)
    {
      e = l->point[i];
      while ((mtxA->edges[e->entry_no].i != pi[i]) && (mtxA->edges[e->entry_no].j != pi[i]))
	e = e->next;
      out += mtxA->edges[e->entry_no].v;
      q += mtxA->edges[e->entry_no].v;
      congestion[i] += mtxA->edges[e->entry_no].v/w;
      i = pi[i];

      e = l->point[j];
      while ((mtxA->edges[e->entry_no].i != pi[j]) && (mtxA->edges[e->entry_no].j != pi[j]))
	e = e->next;
      out += mtxA->edges[e->entry_no].v;
      q += mtxA->edges[e->entry_no].v;
      congestion[j] += mtxA->edges[e->entry_no].v/w;
      j = pi[j];
    }

  *dilation = max(*dilation,q/w);
  return(out);
}

static double find_stretch(graph *mtxA,int *pi,int *d)
{
  int i,E=0;
  double stretch=0,dummy;
  linked *l;
  double dilation;
  double *congestion,cong=0;
  
  congestion = (double *)taucs_calloc(mtxA->n,sizeof(double));

  l = create_linked_list(mtxA,mtxA->n,mtxA->nent,&dummy,&dummy);
  assert(l && congestion);

  for (i=0;i<mtxA->nent;i++)
    {
      if (mtxA->edges[i].i != mtxA->edges[i].j)
	{
	  E++;
	  stretch += dist(mtxA->edges[i].i,mtxA->edges[i].j,mtxA->edges[i].v,mtxA,pi,d,l,&dilation,congestion)/mtxA->edges[i].v;
	}
    }
  
  Do(i,mtxA->n)
    cong = max(cong,congestion[i]);

  printf("Cong-Dil = %f\n",cong*dilation);

  free_linked_list(l);
  taucs_free(congestion);
  return(stretch/E);
}

static
taucs_ccs_matrix*
amst_preconditioner_create(graph *mtxA, double* diag,
			   int rnd,
			   double subgraphs,
			   int stretch_flag)
{
  taucs_ccs_matrix* out;
  /*
  sym_matrix *mtxA;
  sym_matrix *precond;
  */
  /*graph *mtxA;*/
  graph *precond;
  int i;
  /*int tmp;*/
  /*
  int *ivec1, *ivec2 ;
  int *ivec1_p, *ivec2_p ;
  double *dvec; 
  double *dvec_p; 
  */
  int *pi,*d;
  edge **array;
  edge *p;
  int r;
  int Bent,row,col;
  double weight = 0;
  int *first_child,*next_child;
  int *groups,*sub_tree_sizes;
  int curr_group;
  int n;
  int precond_maxnent,chunk;
  /*double *diag;*/
  /*linked l;*/
  linked* lp;
  double dummy;
  double wtime;

  n = mtxA->n;

  wtime = taucs_wtime();

  pi            = (int*)  taucs_malloc(n*sizeof(int));
  d            = (int*)  taucs_malloc(n*sizeof(int));

  /*l = create_linked_list_old(mtxA,n,mtxA->nent); */ /* THIS MAY RUN OUT OF MEMORY ! */
  lp = create_linked_list(mtxA,n,mtxA->nent,&dummy,&dummy); /* THIS MAY RUN OUT OF MEMORY ! */

  if (!pi || !d) {
    taucs_free(pi);
    taucs_free(d);
    taucs_free(diag);
    free_graph(mtxA);
  }

#if 0
  taucs_free(diag);
  diag = analyze_graph(mtxA); /* should change! */
#endif

  /* array is an array of linked lists, which hold the
     off-diagonal entries of mtxA */

  array = lp->point;

  Do(i,n) pi[i] = -1;

  /*
  ivec1 = mtxA->ivec1 ;
  ivec2 = mtxA->ivec2 ;
  dvec = mtxA->dvec ;
  */


  wtime = taucs_wtime() - wtime;
  taucs_printf("\t\tAMST prepare for mst = %.3f seconds\n",wtime);

  wtime = taucs_wtime();
  
  r = rnd % n;

  if (stretch_flag)
    {
      graph *low_stretch_tree;
      linked *lp_low;

      if ((low_stretch_tree = low_stretch(mtxA)) == NULL)
	{
	  taucs_free(pi);
	  taucs_free(diag);
	  free_graph(mtxA);
	  return 0;
	}

      lp_low = create_linked_list(low_stretch_tree,n,low_stretch_tree->nent,&dummy,&dummy); /* THIS MAY RUN OUT OF MEMORY ! */
      if (lp_low == NULL)
	{
	  taucs_free(pi);
	  taucs_free(diag);
	  free_graph(mtxA);
	  free_graph(low_stretch_tree);
	  return 0;
	}
      Do(i,n)
	pi[i] = -2;
      Prim(low_stretch_tree,r,pi,d,lp_low);
      
      free_graph(low_stretch_tree);
      free_linked_list(lp_low);
    }

  else
    {
      Prim(mtxA,r,pi,d,lp);
    }

  /* pi now contains the parent array of the tree, d the distance from the root */
  printf("Stretch = %f\n",find_stretch(mtxA,pi,d));
  

  groups         = (int *) taucs_malloc(n*sizeof(int));
  first_child    = (int *) taucs_malloc(n*sizeof(int));
  next_child     = (int *) taucs_malloc(n*sizeof(int));
  sub_tree_sizes = (int *) taucs_malloc(n*sizeof(int));

  if (!groups || !first_child || !next_child || !sub_tree_sizes) {
    taucs_free(groups);
    taucs_free(first_child);
    taucs_free(next_child);
    taucs_free(sub_tree_sizes);

    taucs_free(pi);
    taucs_free(diag);
    free_graph(mtxA);
    /* free linked_list; */
  }


  if (create_children_arrays(pi,n,&first_child,&next_child) == -1) {
    taucs_free(groups);
    taucs_free(first_child);
    taucs_free(next_child);
    taucs_free(sub_tree_sizes);

    taucs_free(pi);
    taucs_free(diag);
    free_graph(mtxA);
    /* free linked_list; */
  }

  wtime = taucs_wtime() - wtime;
  taucs_printf("\t\tAMST mst = %.3f seconds\n",wtime);

  wtime = taucs_wtime();

  Do(i,n) groups[i] = -1;
  
  curr_group = 0;
  
  compute_sub_tree_sizes(r,first_child,next_child,sub_tree_sizes);

  /* now for every vertex v in the tree, sub_tree_sizes[v] is the size
     of the subtree whose root is v */
  
  divide_to_groups(r,first_child,next_child,pi,&curr_group,groups,
		   sub_tree_sizes,r,subgraphs,n);
  assign_group(r,curr_group,first_child,next_child,groups);
  curr_group++;

  taucs_printf("actual number of subgraphs = %ld\n",curr_group);

  
  taucs_free(first_child);
  taucs_free(next_child);
  taucs_free(sub_tree_sizes);

  wtime = taucs_wtime() - wtime;
  taucs_printf("\t\tAMST partition = %.3f seconds\n",wtime);

  /* now the tree is devided into linked groups */

  wtime = taucs_wtime();

  chunk = max(min((curr_group*(curr_group-1)/2)/10,5000),100);
  precond_maxnent = (n-1) + n + chunk;

  taucs_printf("allocating space for %d entries in precond\n",precond_maxnent);
  
  precond = construct_graph(precond_maxnent);
  if (!precond) {
    taucs_free(pi);
    taucs_free(diag);
    free_graph(mtxA);
    /* free linked_list; */
    return NULL;
  }
  precond->n=mtxA->n;

  /*
  ivec1_p = precond->ivec1 ;
  ivec2_p = precond->ivec2 ;
  dvec_p = precond->dvec ;
  */
	 
  Bent = 0;

  wtime = taucs_wtime() - wtime;
  taucs_printf("\t\tAMST allocating mst precond = %.3f seconds\n",wtime);

  wtime = taucs_wtime();

  /* adds the tree edges to the preconditioner */
  Do(i,n) {
    row = i;
    col = pi[i];
    if (col != (-1)) {
      p = array[row];
      
      while (p != NULL) {
	if (   ((mtxA->edges)[p->entry_no].j == col) 
	    || ((mtxA->edges)[p->entry_no].i == col)) {
	  weight = (mtxA->edges)[p->entry_no].v;
	  break;
	}
	p = p->next;
      }
      
      (precond->edges)[Bent].i = row;
      (precond->edges)[Bent].j = col;
      (precond->edges)[Bent].v = weight;

      Bent++;

      diag[row] -= weight;
      diag[col] -= weight;
    }
  }
  
  precond->nent = Bent;


  wtime = taucs_wtime() - wtime;
  taucs_printf("\t\tAMST adding tree edges = %.3f seconds\n",wtime);

  /* 
     add the heavy edges between every two subgraphs, if such an edge
     exists, and if the subraphs were not already connected through the
     tree 
  */

  wtime = taucs_wtime();

  if (curr_group>1) /* more than 1 group */
    Bent = add_heavy_edges(mtxA,precond,Bent,array,groups,curr_group,pi,diag);

  if (Bent == -1) { /* memory allocation failure in add_heavy_edges */
    taucs_free(pi);
    taucs_free(diag);
    taucs_free(groups);
    free_linked_list(lp);
    /* taucs_free(point_to_heap); */
    free_graph(mtxA);
    free_graph(precond);
    return NULL;
  }

  wtime = taucs_wtime() - wtime;
  taucs_printf("\t\tAMST finding heavy edges = %.3f seconds\n",wtime);
  
  /* allocate more memory to the preconditioner if needed */

  wtime = taucs_wtime();

  if (precond_maxnent < Bent + n) {
    taucs_printf("adding space for %d entries in precond for diagonal entries\n",Bent+n-precond_maxnent);
    precond_maxnent = Bent + n;

    graph_resize(precond,precond_maxnent);
    precond->nent = Bent;
    /*
    ivec1_p = precond->ivec1 ;
    ivec2_p = precond->ivec2 ;
    dvec_p = precond->dvec ;
    */
  }
  
  Do(i,n) {

    (precond->edges)[Bent].i = i;
    (precond->edges)[Bent].j = i;
    (precond->edges)[Bent].v = diag[i];

    Bent++;
  }

  precond->nent = Bent;

  wtime = taucs_wtime() - wtime;
  taucs_printf("\t\tAMST resize and add heavy edges = %.3f seconds\n",wtime);

  wtime = taucs_wtime();

  taucs_free(pi);
  taucs_free(diag);
  taucs_free(groups);
  free_linked_list(lp);
  taucs_printf("actual number of entries in preconditioner = %d\n",Bent);
  free_graph(mtxA);

  /* out could be NULL, but we return out anyway */

  out = graph_to_ccs_matrix(precond); 

  assert(out);

  free_graph(precond);

  wtime = taucs_wtime() - wtime;
  taucs_printf("\t\tAMST free memory and convert to ccs = %.3f seconds\n",wtime);

#define TAUCS_VAIDYA_RELAX_NONONO
#ifdef TAUCS_VAIDYA_RELAX
  {
    taucs_ccs_matrix* A = ccs_mtxA;
    taucs_ccs_matrix* M = out;
    int j;
    
    for (j=0; j<A->n; j++) {
      double dA, dM, modification;
      int i,ip;
      for (ip=(A->colptr)[j]; ip<(A->colptr)[j+1]; ip++) {
	i = (A->rowind)[ ip ];
	if (i==j) {
	  dA = (A->values.d/*taucs_values*/)[ ip ];
	  break;
	}
      }

      for (ip=(M->colptr)[j]; ip<(M->colptr)[j+1]; ip++) {
	i = (M->rowind)[ ip ];
	if (i==j) {
	  dM = (M->values.d/*taucs_values*/)[ ip ];
	  break;
	}
      }
      
      assert(dA >= dM);
      
      modification = dA - dM;
      if (j < 30) printf(">>> %.4e %.4e\n",dA,dM);
      
      (M->values.d/*taucs_values*/)[ ip ] += 0.01 * modification;
    }
  }
#endif

  return out;
}


taucs_ccs_matrix*
taucs_amwb_preconditioner_create(taucs_ccs_matrix *A, 
				 int rnd,
				 double subgraphs,
				 int stretch_flag)
{
  double  wtime;
  double* diag;
  graph*  G_A;
  int     diagnostics;
  int     n;

  if (!(A->flags & TAUCS_DOUBLE)) {
    taucs_printf("taucs_amwb_preconditioner_create: matrix must be double-precision real\n");
    return NULL;
  }

  n = A->n;

  diag = (double*) taucs_malloc(n*sizeof(double));
  if (diag == NULL) return NULL;

  wtime = taucs_wtime();
  G_A = ccs_matrix_to_graph_plus(A,&diagnostics,diag,1 /* force diag dominance */);
  if (!G_A) {
    taucs_free(diag);
    return NULL;
  }
  wtime = taucs_wtime() - wtime;
  taucs_printf("\t\tAMWB matrix-to-graph + analysis = %.3f seconds\n",wtime);

  if (diagnostics & TAUCS_SYM_NOT_SYMLOWER) {
    taucs_printf("taucs_amwb_preconditioner_create: matrix must be symmetrix & lower\n");
    /* in this case, G_A == NULL, no need to free */
    taucs_free(diag);
    return A;
  }
  if (diagnostics & TAUCS_SYM_NOT_DIAGDOMINANT) {
    taucs_printf("taucs_amwb_preconditioner_create: matrix not diagonally dominant\n");
    taucs_free(diag);
    free_graph(G_A);
    return A;
  }
  if (diagnostics & TAUCS_SYM_NEG_DIAGONALS) {
    taucs_printf("taucs_amwb_preconditioner_create: negative diagonal elements\n");
    taucs_free(diag);
    free_graph(G_A);
    return A;
  }

  if (diagnostics & TAUCS_SYM_POS_OFFDIAGONALS)
    return amwb_preconditioner_create(G_A, diag, rnd, subgraphs);
  else
    return amst_preconditioner_create(G_A, diag, rnd, subgraphs,stretch_flag);
}

#endif /* TAUCS_CORE_DOUBLE */

/*********************************************************/
/*                                                       */
/*********************************************************/

