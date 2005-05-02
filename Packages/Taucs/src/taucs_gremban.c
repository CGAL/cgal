/*********************************************************/
/* TAUCS                                                 */
/* Author: Doron Chen                                    */
/* File  : taucs_gremban.c                               */
/* Description: constructs multilevel support            */
/*        reconditioners (including Gremban-Miller)      */
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "taucs.h"

/*#include <unistd.h>*/

/*long int random() omer*/

#ifdef TAUCS_CORE_DOUBLE

/* #include "../metis-4.0/Lib/defs.h" */
/* #include "../metis-4.0/Lib/struct.h" */
/* #include "../metis-4.0/Lib/proto.h" */

typedef int idxtype;

typedef struct {
  taucs_ccs_matrix* L;
  int n,k;
  double* Ztilde;
  double* Rtilde;
} multilevel_args;

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

#define Do(i,n) for ((i)=0;(i)<(n);(i)++)

static
taucs_ccs_matrix* construct_ccs_matrix(int nent,int n)
{
  taucs_ccs_matrix *out;
  
  out = (taucs_ccs_matrix *)taucs_malloc(sizeof(taucs_ccs_matrix));
  if (!out) return NULL;
  out->colptr = (int *)taucs_malloc((n+1)*sizeof(int));
  out->rowind = (int *)taucs_malloc(nent*sizeof(int));
  out->taucs_values = (double *)taucs_malloc(nent*sizeof(double));
  if (!(out->colptr) || !(out->rowind) || !(out->taucs_values)) {
    taucs_free(out->colptr);
    taucs_free(out->rowind);
    taucs_free(out->taucs_values);
    taucs_free(out);
    return NULL;
  }
  
  out->n = n;
  out->m = n;
  out->flags = TAUCS_SYMMETRIC | TAUCS_LOWER | TAUCS_DOUBLE;
  
  return out;
}

static
taucs_ccs_matrix *graph_to_ccs_matrix(graph *A)
{
  taucs_ccs_matrix *out;
  int n,nent,i,j1,j2;
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
      out->taucs_values[tmp[j1]]=(A->edges)[i].v;
      tmp[j1]++;
    }

  taucs_free(tmp);
  return(out);
}

#if 0
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
#endif /* 0, we don't need this function */

static
void free_graph(graph *a)
{
  if(a)
    {
      taucs_free(a->edges);
      taucs_free(a);
    }
}

/* we use the version in taucs_vaidya.c */
#if 1
extern int taucs_check_diag_dominant_matrix(graph *A, int force_diagonal_dominance);
#else
extern int check_diag_dominant_matrix(graph *A, int force_diagonal_dominance);
#define EPSILON 0.00000001
static
int check_diag_dominant_matrix(graph *A)
{
  int i;
  double *sum;
  int n;

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
	  sum[(A->edges)[i].i]+=fabs((A->edges)[i].v);
	  if ((A->edges)[i].v < 0)
	    {
	      taucs_printf("ERROR! This matrix is not diagonally dominant. It has negative diagonals.\n");
	      taucs_free(sum);
	      return -2;
	    }
	}
      
    }
  
  Do(i,n)
    {
      if (sum[i] < -EPSILON)
	{
	  taucs_printf("ERROR! This matrix is not diagonally dominant. sum[%d] = %lf\n",i,sum[i]);
	  taucs_free(sum);
	  return -2;
	}
    }
  
  taucs_free(sum);
  return 0;
}
#endif


int
taucs_sg_preconditioner_solve(void*  vP,
			      double* Z, 
			      double* R)
{
  multilevel_args* P = (multilevel_args*) vP;
  /*int nplusk = (P->L)->n;*/
  int i;
  int n = P->n;
  int k = P->k;

  for (i=0; i<n;   i++) (P->Rtilde)[i] = R[i];
  for (i=n; i<n+k; i++) (P->Rtilde)[i] = 0.0;

  taucs_ccs_solve_llt(P->L,
		      P->Ztilde,
		      P->Rtilde);

  for (i=0; i<n;   i++) Z[i] = (P->Ztilde)[i];

  return 0;
}

typedef struct {
  int n;
  idxtype *colptr;
  idxtype *rowind;
  idxtype *values;
} Metis_struct;

Metis_struct *Metis_struct_create(int n,int nent)
{
  Metis_struct *out;

  out=(Metis_struct *)taucs_malloc(sizeof(Metis_struct));
  
  if (!out)
    return NULL;
  
  out->n = n;
  out->colptr = (idxtype *)taucs_malloc((n+1)*sizeof(idxtype));
  out->rowind = (idxtype *)taucs_malloc(nent*sizeof(idxtype));
  out->values = (idxtype *)taucs_malloc(nent*sizeof(idxtype));
  
  if ((out->colptr==NULL)||(out->rowind==NULL)||(out->values==NULL))
    {
      taucs_free(out->colptr);
      taucs_free(out->rowind);
      taucs_free(out->values);
      return NULL;
    }

  return(out);
}

void Metis_struct_free(Metis_struct *A)
{
  if (A)
    {
      taucs_free(A->colptr);
      taucs_free(A->rowind);
      taucs_free(A->values);
      taucs_free(A);
    }
}

Metis_struct *taucs_ccs_matrix_to_Metis_struct(taucs_ccs_matrix *A)
{
  Metis_struct *out;
  int n,nent,i,j,j1,j2;
  int *tmp;

  n = A->n;
  nent = 0;

  tmp = (int *)taucs_malloc(n*sizeof(int));
  if (!tmp) return NULL;

  Do(i,n)
    tmp[i] = 0;
  Do(i,n)
    {

      for(j=A->colptr[i];j<A->colptr[i+1];j++)
	if (i!=(A->rowind)[j])
	  {
	    tmp[i]++;
	    tmp[(A->rowind)[j]]++;
	    nent+=2;
	  }
    }

  out = Metis_struct_create(n,nent);
  if (out == NULL)
    {
      taucs_free(tmp);
      return NULL;
    }
  
  out->colptr[0] = 0;
  Do(i,n)
    out->colptr[i+1] = out->colptr[i] + tmp[i];

  Do(i,n)
    tmp[i] = out->colptr[i];

  Do(i,n)
    for(j=A->colptr[i];j<A->colptr[i+1];j++)    
      if (i!=A->rowind[j])
	{
	  j1 = i;
	  j2 = A->rowind[j];
	  out->rowind[tmp[j1]]=j2;
	  out->rowind[tmp[j2]]=j1;
	  /* out->values[tmp[j1]]=(idxtype)min(10000,-A->values[j]); */
	  /* out->values[tmp[j2]]=(idxtype)min(10000,-A->values[j]); */
	  out->values[tmp[j1]]=1;
	  out->values[tmp[j2]]=1;
	  tmp[j1]++;
	  tmp[j2]++;
	}

  taucs_free(tmp);
  return(out);
}

taucs_ccs_matrix *taucs_ccs_matrix_to_taucs_ccs_matrix(taucs_ccs_matrix *A,double *diag)
{
  taucs_ccs_matrix *out;
  int n,nent,i,j,j1,j2;
  int *tmp;

  n = A->n;
  nent = 0;

  tmp = (int *)taucs_malloc(n*sizeof(int));
  if (!tmp) return NULL;

  Do(i,n)
    tmp[i] = 0;

  for(i=0;i<n;i++)
    {
      for(j=(A->colptr[i]);j<(A->colptr[i+1]);j++)
	{
	  if (i!=(A->rowind)[j])
	    {
	      tmp[i]++;
	      tmp[(A->rowind)[j]]++;
	      nent+=2;
	    }
	  else
	    diag[i] = (A->taucs_values)[j];
	}
    }

  out = taucs_dtl(ccs_create)(n,n,nent);
  if (out == NULL)
    {
      taucs_free(tmp);
      return NULL;
    }
  
  out->colptr[0] = 0;
  Do(i,n)
    out->colptr[i+1] = out->colptr[i] + tmp[i];

  Do(i,n)
    tmp[i] = out->colptr[i];

  Do(i,n)
    for(j=A->colptr[i];j<A->colptr[i+1];j++)    
      if (i!=A->rowind[j])
	{
	  j1 = i;
	  j2 = A->rowind[j];
	  out->rowind[tmp[j1]]=j2;
	  out->rowind[tmp[j2]]=j1;
	  out->taucs_values[tmp[j1]]=(idxtype)A->taucs_values[j];
	  out->taucs_values[tmp[j2]]=(idxtype)A->taucs_values[j];
	  tmp[j1]++;
	  tmp[j2]++;
	}

  taucs_free(tmp);
  return(out);
}

void Metis_struct_print(Metis_struct *A)
{
  int i;
  int j;
  int n,nent;

  n=A->n;
  nent = A->colptr[n];

  Do(i,n)
    for(j=A->colptr[i];j<A->colptr[i+1];j++)
      printf("%d %d %d\n",i,A->rowind[j],A->values[j]);
  exit(345);
}

graph *graph_create(int size)
{
  graph *out;

  out = (graph *)taucs_malloc(sizeof(graph));
  if (out == NULL)
    return NULL;
  
  out->edges = (wedge*) taucs_malloc(size*sizeof(wedge));
  if (!out->edges) {
    taucs_free (out);
    return NULL;
  }
  
  out->max_size = size;
  
  return out;
}

typedef unsigned char byte;

typedef struct {
  byte type; /* 0 - gremban, 1 - toledo, 2 - vaidya */
  int k;
} instruction;

int partition(int *quicksort_array_nodes_1,int *quicksort_array_nodes_2,double *quicksort_array_values,int p,int r)
{
  int x1,x2,i,j,tmpi1,tmpi2;
  double tmpd;
  
  x1 = quicksort_array_nodes_1[p];
  x2 = quicksort_array_nodes_2[p];
  i = p-1;
  j = r+1;
  while(1)
    {
      do
	j--;
      while ((quicksort_array_nodes_1[j]>x1)||((quicksort_array_nodes_1[j]==x1)&&(quicksort_array_nodes_2[j]>x2)));
      do
	i++;
      while ((quicksort_array_nodes_1[i]<x1)||((quicksort_array_nodes_1[i]==x1)&&(quicksort_array_nodes_2[i]<x2)));
      
      if (i<j)
	{
	  tmpi1 = quicksort_array_nodes_1[i];
	  tmpi2 = quicksort_array_nodes_2[i];
	  tmpd = quicksort_array_values[i];
	  quicksort_array_nodes_1[i] = quicksort_array_nodes_1[j];
	  quicksort_array_nodes_2[i] = quicksort_array_nodes_2[j];
	  quicksort_array_values[i] = quicksort_array_values[j];
	  quicksort_array_nodes_1[j] = tmpi1;
	  quicksort_array_nodes_2[j] = tmpi2;
	  quicksort_array_values[j] = tmpd;
	}
      else
	return(j);
    }
}

void quicksort(int *quicksort_array_nodes_1,int *quicksort_array_nodes_2,double *quicksort_array_values,int p,int r)
{
  int q;
  
  if (p<r)
    {
      q = partition(quicksort_array_nodes_1,quicksort_array_nodes_2,quicksort_array_values,p,r);
      quicksort(quicksort_array_nodes_1,quicksort_array_nodes_2,quicksort_array_values,p,q);
      quicksort(quicksort_array_nodes_1,quicksort_array_nodes_2,quicksort_array_values,q+1,r);
    }
  
}

/* from metis.h */
extern
void METIS_PartGraphRecursive(int *, 
			      idxtype *, 
			      idxtype *, 
			      idxtype *, 
			      idxtype *, 
			      int *, 
			      int *, 
			      int *, 
			      int *, 
			      int *, 
			      idxtype *); 
void METIS_PartGraphKway(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 

int quicksort_and_shrink(int *quicksort_array_nodes_1,int *quicksort_array_nodes_2,double *quicksort_array_values,int quicksort_index)
{
  int i,outindex=0,curr_pair1,curr_pair2,tmp,tmp1,tmp2;
  double acc=0,tmp3;

  for(i=0;i<quicksort_index;i++)
    {
      tmp = rand()%(quicksort_index-i);
      tmp1 = quicksort_array_nodes_1[i+tmp];
      tmp2 = quicksort_array_nodes_2[i+tmp];
      tmp3 = quicksort_array_values[i+tmp];
      quicksort_array_nodes_1[i+tmp]=quicksort_array_nodes_1[i];
      quicksort_array_nodes_2[i+tmp]=quicksort_array_nodes_2[i];
      quicksort_array_values[i+tmp]=quicksort_array_values[i];
      quicksort_array_nodes_1[i]=tmp1;
      quicksort_array_nodes_2[i]=tmp2;
      quicksort_array_values[i]=tmp3;
    }

  quicksort(quicksort_array_nodes_1,quicksort_array_nodes_2,quicksort_array_values,0,quicksort_index-1);

  curr_pair1 = quicksort_array_nodes_1[0];
  curr_pair2 = quicksort_array_nodes_2[0];
  
  for(i=0;i<quicksort_index;i++)
    {
      if ((quicksort_array_nodes_1[i]!=curr_pair1)||(quicksort_array_nodes_2[i]!=curr_pair2))
	{
	  quicksort_array_nodes_1[outindex]=curr_pair1;
	  quicksort_array_nodes_2[outindex]=curr_pair2;
	  quicksort_array_values[outindex++]=acc;
	  acc=quicksort_array_values[i];
	  curr_pair1=quicksort_array_nodes_1[i];
	  curr_pair2=quicksort_array_nodes_2[i];
	}
      else
	acc += quicksort_array_values[i];
    }
  
  quicksort_array_nodes_1[outindex]=curr_pair1;
  quicksort_array_nodes_2[outindex]=curr_pair2;
  quicksort_array_values[outindex++]=acc;
  
  return(outindex);
  
}

int create_recursive_preconditioner(graph *out, /* output: graph of the preconditioner. */
				    int curr_vertex, 
				    int *next_unused_vertex,
				    int *curr_entry,
				    Metis_struct *father,
				    int *vertex_perm, /* translates the vertex numbers in mtxA to the actual vertex numbers */
				    int *inv_perm,  /* inv_perm is an array of -1 when this function is entered */
				    double diag,
				    double *diagonal, /* diagonal values of mtxA */
				    instruction *inst, /* array of instructions : how the preconditioner is built at each level */
				    int max_inst, /* number of levels of instructions */
				    int curr_inst, /* current level of instructions */
				    taucs_ccs_matrix *taucs_ccs_mtxA,   /* NOTICE - taucs_ccs_mtxA should not contain the diagonal values,
									   and it should contains each entry twice !!! */
				    char *ordering,
				    int *p1
				    )
{
  int i,j,n,k,orig_n,p,t1,p2;
  /* int t2 */
  int options[5]={0};
  int wgtflag = 1;
  int numflag = 0;
  int nparts,edgecut;
  idxtype *part;
  Metis_struct **sons;
  int *tmp,*tmp1,*vertices_in_subgraphs;
  double *weights,*diags;
  int *perm_tmp=NULL,*inv_perm_tmp=NULL;
  int **vertex_perms;
  int local_next_unused_vertex;
  int success = 1;
  graph *order_graph;
  taucs_ccs_matrix *order_ccs=NULL;
  int is_root=0;
  static int ordering_counter,first=1;
  static int ordering_counter_leaves = 0;
  int *quicksort_array_nodes_1=NULL,*quicksort_array_nodes_2=NULL,quicksort_index=0;
  double *quicksort_array_values=NULL;
  taucs_ccs_matrix *vaidya_ccs=NULL;

  if (first)
    {
      first = 0;
      ordering_counter = curr_vertex;
      is_root = 1;
    }

  local_next_unused_vertex = *next_unused_vertex;
  
  n = father->n;
  orig_n = taucs_ccs_mtxA->n;

  part = (idxtype *)taucs_malloc(n*sizeof(idxtype));
  if (!part)
    return 0;

  nparts = min(inst[curr_inst].k,n);
  if (curr_inst == max_inst)
    nparts = n;
  
  if (nparts == n)
    for(i=0;i<n;i++)
      part[i]=i;
  else
    if (nparts == 1)
      for(i=0;i<n;i++)
	part[i]=0;
    else
      {
	int *visited; /* helps determine how many parts the graph was actually divided into */
	int actual_nparts=0;
	/* taucs_printf("calling METIS\n"); */
#ifdef NOMETIS
	/* omer - for warning*/
	edgecut = 0;
#else
	if (nparts < 8)
	  METIS_PartGraphRecursive(&n,father->colptr,father->rowind,
				   NULL,father->values,&wgtflag,&numflag,
				   &nparts,options,&edgecut,part);
	else
	  METIS_PartGraphKway(&n,father->colptr,father->rowind,
			      NULL,father->values,&wgtflag,&numflag,
			      &nparts,options,&edgecut,part);
#endif	
	/* taucs_printf("calling METIS: done\n"); */
	visited=taucs_calloc(nparts,sizeof(int));
	for(i=0;i<n;i++)
	  if(visited[part[i]]==0)
	    {
	      visited[part[i]]=1;
	      actual_nparts++;
	    }
	
	if (actual_nparts!=nparts)
	  {
	    actual_nparts=0;
	    for(i=0;i<nparts;i++)
	      if (visited[i])
		visited[i]=actual_nparts++;
	    
	    for(i=0;i<n;i++)
	      part[i] = visited[part[i]];
	    nparts = actual_nparts;
	  }
	taucs_free(visited);
	
      }
  
  for(i=0;i<n;i++)
    inv_perm[vertex_perm[i]] = i;

  weights = (double *)taucs_calloc(n,sizeof(double));
  diags = (double *)taucs_calloc(nparts,sizeof(double));
  if ((!weights)||(!diags))
    {
      taucs_free(part);
      taucs_free(weights);
      taucs_free(diags);
      return(0);
    }

  if (inst[curr_inst].type == 0) /* Gremban */
    {
      for(i=0;i<n;i++)
	{
	  p = part[i];
	  t1 = vertex_perm[i];
	  for(j=taucs_ccs_mtxA->colptr[t1];j<taucs_ccs_mtxA->colptr[t1+1];j++)
	    if ((inv_perm[taucs_ccs_mtxA->rowind[j]] == -1) || (part[inv_perm[taucs_ccs_mtxA->rowind[j]]] != p))
	      weights[p] += taucs_ccs_mtxA->taucs_values[j]; /* weights[p] is the sum of all weights of all edges between part p and the rest the graph */
	}
      if (curr_inst < max_inst)
	{
	  for(i=0;i<nparts;i++)
	    if (weights[i]!=0)
	      {
		out->edges[(*curr_entry)].i=curr_vertex;
		out->edges[(*curr_entry)].j=local_next_unused_vertex+i;
		out->edges[(*curr_entry)].v=weights[i];
		diags[i] -= weights[i];
		diag     -= weights[i];
		(*curr_entry)++;
	      }
	}
      else /* curr_inst == max_inst */
	{
	  for(i=0;i<nparts;i++)
	    if (weights[i]!=0)
	      {
		out->edges[(*curr_entry)].i=curr_vertex;
		out->edges[(*curr_entry)].j=vertex_perm[i];
		out->edges[(*curr_entry)].v=weights[i];
		diag -= weights[i];
		(*curr_entry)++;
	      }
	  for(i=0;i<nparts;i++)
	    {
	      out->edges[(*curr_entry)].i=vertex_perm[i];
	      out->edges[(*curr_entry)].j=vertex_perm[i];
	      if (diagonal[vertex_perm[i]])
		out->edges[(*curr_entry)].v=diagonal[vertex_perm[i]];
	      else
		out->edges[(*curr_entry)].v=1;
	      (*curr_entry)++;
	      p1[ordering_counter_leaves++] = vertex_perm[i];
	    }
	}
    }
  else /* Toledo or Vaidya */
    {
      int order_count;
      
      for(i=0;i<n;i++)
	{
	  p = part[i];
	  t1 = vertex_perm[i];
	
	  for(j=taucs_ccs_mtxA->colptr[t1];j<taucs_ccs_mtxA->colptr[t1+1];j++)
	    {
	      if (inv_perm[taucs_ccs_mtxA->rowind[j]] == -1)
		weights[p] += taucs_ccs_mtxA->taucs_values[j]; /* weights[p] is the sum of all weights of all edges between part p and vertices outside current subgraph */
	      else
		if (part[inv_perm[taucs_ccs_mtxA->rowind[j]]] > p)
		  quicksort_index++;
	    }
	  
	}

      if (quicksort_index)
	{
	  
	  quicksort_array_nodes_1 = (int *)taucs_malloc(quicksort_index*sizeof(int));
	  quicksort_array_nodes_2 = (int *)taucs_malloc(quicksort_index*sizeof(int));
	  quicksort_array_values = (double *)taucs_malloc(quicksort_index*sizeof(double));
	  if ((!quicksort_array_nodes_1)||(!quicksort_array_nodes_2)||(!quicksort_array_values))
	    {
	      taucs_free(quicksort_array_nodes_1);
	      taucs_free(quicksort_array_nodes_2);
	      taucs_free(quicksort_array_values);
	      taucs_free(part);
	      taucs_free(weights);
	      taucs_free(diags);
	      return(0);
	    }
	  
	  quicksort_index = 0;
	  
	  for(i=0;i<n;i++)
	    {
	      p = part[i];
	      t1 = vertex_perm[i];
	      
	      for(j=taucs_ccs_mtxA->colptr[t1];j<taucs_ccs_mtxA->colptr[t1+1];j++)
		{
		  if (inv_perm[taucs_ccs_mtxA->rowind[j]] != -1)
		    if (part[inv_perm[taucs_ccs_mtxA->rowind[j]]] > p)
		      /* since each entry appears twice in taucs_ccs_mtxA,
			 we need not update pairs when part[...]<p */
		      {
			p2=part[inv_perm[taucs_ccs_mtxA->rowind[j]]];
			quicksort_array_nodes_1[quicksort_index] = p;
			quicksort_array_nodes_2[quicksort_index] = p2;
			quicksort_array_values[quicksort_index++] = taucs_ccs_mtxA->taucs_values[j];
		      }
		}
	    }
	  
	  quicksort_index = quicksort_and_shrink(quicksort_array_nodes_1,quicksort_array_nodes_2,quicksort_array_values,quicksort_index);

	}
      
      order_count = quicksort_index;

      order_graph = graph_create(order_count+nparts);
      order_graph->nent = order_count+nparts;
      order_count = 0;
      order_graph->n = nparts;
      for(i=0;i<nparts;i++)
	{
	  order_graph->edges[order_count].i=i;
	  order_graph->edges[order_count].j=i;
	  order_graph->edges[order_count].v=0;
	  order_count++;
	}
      if (curr_inst < max_inst)
	{
	  for(i=0;i<nparts;i++)
	    if (weights[i]!=0)
	      {
		out->edges[(*curr_entry)].i=curr_vertex;
		out->edges[(*curr_entry)].j=local_next_unused_vertex+i;
		out->edges[(*curr_entry)].v=weights[i];
		diags[i] -= weights[i];
		diag     -= weights[i];
		(*curr_entry)++;
	      }

	  for(k=0;k<quicksort_index;k++)
	    {
	      i = quicksort_array_nodes_1[k];
	      j = quicksort_array_nodes_2[k];
	      order_graph->edges[order_count].i=i;
	      order_graph->edges[order_count].j=j;
	      order_graph->edges[order_count].v=quicksort_array_values[k];
	      order_graph->edges[i].v-=quicksort_array_values[k];
	      order_graph->edges[j].v-=quicksort_array_values[k];
	      order_count++;
	    }

	  order_ccs = graph_to_ccs_matrix(order_graph);
	  free_graph(order_graph);
	  if (!order_ccs)
	    {
	      taucs_free(part);
	      taucs_free(weights);
	      taucs_free(diags);
	      taucs_free(quicksort_array_values);
	      taucs_free(quicksort_array_nodes_1);
	      taucs_free(quicksort_array_nodes_2);
	      return(0);
	    }

	  if (inst[curr_inst].type == 1) /* Toledo */
	    {
	      for(k=0;k<quicksort_index;k++)
		{
		  i = quicksort_array_nodes_1[k];
		  j = quicksort_array_nodes_2[k];
		  out->edges[(*curr_entry)].i=local_next_unused_vertex+i;
		  out->edges[(*curr_entry)].j=local_next_unused_vertex+j;
		  out->edges[(*curr_entry)].v=quicksort_array_values[k];
		  diags[i] -= quicksort_array_values[k];
		  diags[j] -= quicksort_array_values[k];
		  weights[i] += quicksort_array_values[k];
		  weights[j] += quicksort_array_values[k];
		  (*curr_entry)++;
		}

	    }
	  else /* Vaidya */
	    {
	      vaidya_ccs = taucs_amwb_preconditioner_create(order_ccs,1,(order_ccs->n)/8,0);
	      for(i=0;i<vaidya_ccs->n;i++)
		for(j=vaidya_ccs->colptr[i];j<vaidya_ccs->colptr[i+1];j++)
		  if (i != (vaidya_ccs->rowind[j]))
		    {
		      out->edges[(*curr_entry)].i=local_next_unused_vertex+i;
		      out->edges[(*curr_entry)].j=local_next_unused_vertex+vaidya_ccs->rowind[j];
		      out->edges[(*curr_entry)].v=vaidya_ccs->taucs_values[j];
		      diags[i] -= vaidya_ccs->taucs_values[j];
		      diags[vaidya_ccs->rowind[j]] -= vaidya_ccs->taucs_values[j];
		      weights[i] += quicksort_array_values[k];
		      weights[vaidya_ccs->rowind[j]] += quicksort_array_values[k];
		      (*curr_entry)++;
		    }
	      
	      taucs_ccs_free(vaidya_ccs);
	    }

	  taucs_ccs_order(order_ccs,&perm_tmp,&inv_perm_tmp,ordering);
	  taucs_free(inv_perm_tmp);
	  taucs_ccs_free(order_ccs);
	}
      else /* curr_inst == max_inst */
	{
	  for(i=0;i<nparts;i++)
	    if (weights[i]!=0)
	      {
		out->edges[(*curr_entry)].i=curr_vertex;
		out->edges[(*curr_entry)].j=vertex_perm[i];
		out->edges[(*curr_entry)].v=weights[i];
		diag     -= weights[i];
		(*curr_entry)++;
	      }
	  for(k=0;k<quicksort_index;k++)
	    {
	      i = quicksort_array_nodes_1[k];
	      j = quicksort_array_nodes_2[k];
	      order_graph->edges[order_count].i=i;
	      order_graph->edges[order_count].j=j;
	      order_graph->edges[order_count].v=quicksort_array_values[k];
	      order_graph->edges[i].v-=quicksort_array_values[k];
	      order_graph->edges[j].v-=quicksort_array_values[k];
	      order_count++;
	    }
	  order_ccs = graph_to_ccs_matrix(order_graph);
	  free_graph(order_graph);
	  if (!order_ccs)
	    {
	      taucs_free(part);
	      taucs_free(weights);
	      taucs_free(diags);
	      taucs_free(quicksort_array_values);
	      taucs_free(quicksort_array_nodes_1);
	      taucs_free(quicksort_array_nodes_2);
	      return(0);
	    }
	  
	  if (inst[curr_inst].type == 1) /* Toledo */
	    {
	      for(k=0;k<quicksort_index;k++)
		{
		  i = quicksort_array_nodes_1[k];
		  j = quicksort_array_nodes_2[k];
		  out->edges[(*curr_entry)].i=vertex_perm[i];
		  out->edges[(*curr_entry)].j=vertex_perm[j];
		  out->edges[(*curr_entry)].v=quicksort_array_values[k];
		  weights[i] += quicksort_array_values[k];
		  weights[j] += quicksort_array_values[k];
		  (*curr_entry)++;
		}
	    }
	  else /* Vaidya */
	    {
	      vaidya_ccs = taucs_amwb_preconditioner_create(order_ccs,1,1,0);
	      for(i=0;i<vaidya_ccs->n;i++)
		for(j=vaidya_ccs->colptr[i];j<vaidya_ccs->colptr[i+1];j++)
		  if (i != (vaidya_ccs->rowind[j]))
		    {
		      out->edges[(*curr_entry)].i=vertex_perm[i];
		      out->edges[(*curr_entry)].j=vertex_perm[vaidya_ccs->rowind[j]];
		      out->edges[(*curr_entry)].v=vaidya_ccs->taucs_values[j];
		      weights[i] += vaidya_ccs->taucs_values[j];
		      weights[vaidya_ccs->rowind[j]] += vaidya_ccs->taucs_values[j];
		      (*curr_entry)++;
		    }
	      taucs_ccs_free(vaidya_ccs);
	    }

	  taucs_ccs_order(order_ccs,&perm_tmp,&inv_perm_tmp,ordering);
	  taucs_free(inv_perm_tmp);
	  taucs_ccs_free(order_ccs);
	  for(i=0;i<nparts;i++)
	    {
	      out->edges[(*curr_entry)].i=vertex_perm[i];
	      out->edges[(*curr_entry)].j=vertex_perm[i];
	      if (diagonal[vertex_perm[i]])
		out->edges[(*curr_entry)].v=diagonal[vertex_perm[i]];
	      else
		out->edges[(*curr_entry)].v=1;
	      (*curr_entry)++;
	      p1[ordering_counter_leaves++] = vertex_perm[perm_tmp[i]];	      
	    }
	}
    }
  

  taucs_free(quicksort_array_values);
  taucs_free(quicksort_array_nodes_1);
  taucs_free(quicksort_array_nodes_2);

  out->edges[(*curr_entry)].i=curr_vertex;
  out->edges[(*curr_entry)].j=curr_vertex;
  if (diag)
    out->edges[(*curr_entry)].v=diag;
  else
    out->edges[(*curr_entry)].v=1;
  (*curr_entry)++;
  
  /* return inv_perm to its original state - contains only -1-s */
  for(i=0;i<n;i++)
    inv_perm[vertex_perm[i]] = -1;
  
  vertices_in_subgraphs = (int *)taucs_malloc(n*sizeof(int)); /* for each vertex i, tmp[i] is the vertex's position in the subgraph*/
  tmp = (int *)taucs_calloc(nparts,sizeof(int));  /* tmp[i] is the next free vertex in part i */
  tmp1 = (int *)taucs_calloc(nparts,sizeof(int)); /* helps compute number of vertices in each part i */
  
  if ((!vertices_in_subgraphs)||(!tmp)||(!tmp1))
    {
      taucs_free(part);
      taucs_free(weights);
      taucs_free(diags);
      taucs_free(vertices_in_subgraphs);
      taucs_free(tmp);
      taucs_free(tmp1);
      return(0);
    }

  for(i=0;i<n;i++)
    {
      vertices_in_subgraphs[i]=tmp[part[i]]++;
      tmp1[part[i]]++;
    }

  for(i=0;i<nparts;i++)
    tmp[i] = 0;
  /* now tmp[i] will help compute number of edges in each part */
  for(i=0;i<n;i++)
    {
      p = part[i];
      for(j=father->colptr[i];j<father->colptr[i+1];j++)
	{
	  if (part[father->rowind[j]]==p)
	    tmp[p]++;
	}
    }

  sons = (Metis_struct **)taucs_malloc(nparts*sizeof(Metis_struct *));
  if (!sons)
    {
      taucs_free(part);
      taucs_free(weights);
      taucs_free(diags);
      taucs_free(vertices_in_subgraphs);
      taucs_free(tmp);
      taucs_free(tmp1);
      return(0);
    }

  for(i=0;i<nparts;i++)
    {
      sons[i] = Metis_struct_create(tmp1[i],tmp[i]);
      if (sons[i] == 0)
	{
	  taucs_free(part);
	  taucs_free(diags);
	  taucs_free(weights);
	  taucs_free(vertices_in_subgraphs);
	  taucs_free(tmp);
	  taucs_free(tmp1);
	  for(j=0;j<i;j++)
	    Metis_struct_free(sons[j]);
	  taucs_free(sons);
	  return(0);
	}
    }

  for(i=0;i<nparts;i++)
    {
      tmp[i] = 0;
      sons[i]->colptr[0]=0;
    }
  
  /* tmp[p] is the next free extry in sons[p]->colptr[p] */
  for(i=0;i<n;i++)
    {
      p = part[i];
      for(j=father->colptr[i];j<father->colptr[i+1];j++)
	if(part[father->rowind[j]]==p)
	  {
	    sons[p]->rowind[tmp[p]] = vertices_in_subgraphs[father->rowind[j]];
	    sons[p]->values[tmp[p]] = father->values[j];
	    tmp[p]++;
	  }
      sons[p]->colptr[vertices_in_subgraphs[i]+1]=tmp[p];
    }

  Metis_struct_free(father);
  taucs_free(tmp);

  if (curr_inst < max_inst)
    (*next_unused_vertex) += nparts;

  out->n = max(out->n,*next_unused_vertex);

  vertex_perms=(int **)taucs_malloc(nparts*sizeof(int *));
  if (!vertex_perms)
    {
      taucs_free(part);
      taucs_free(weights);
      taucs_free(diags);
      taucs_free(vertices_in_subgraphs);
      taucs_free(tmp1);
      for(j=0;j<nparts;j++)
	Metis_struct_free(sons[j]);
      taucs_free(sons);
      return(0);
    }
  
  for(i=0;i<nparts;i++)
    {
      vertex_perms[i] = (int *)taucs_malloc(tmp1[i]*sizeof(int));
      if (!vertex_perms[i])
	{
	  taucs_free(part);
	  taucs_free(diags);
	  taucs_free(weights);
	  taucs_free(vertices_in_subgraphs);
	  taucs_free(tmp1);
	  for(j=0;j<nparts;j++)
	    Metis_struct_free(sons[j]);
	  taucs_free(sons);
	  for(j=0;j<i;j++)
	    taucs_free(vertex_perms[i]);
	  taucs_free(vertex_perms);
	  return(0);
	}
    }
  for(i=0;i<n;i++)
    vertex_perms[part[i]][vertices_in_subgraphs[i]] = vertex_perm[i];
  
  taucs_free(vertices_in_subgraphs);

  taucs_free(part);

  if (curr_inst < max_inst)
    {
      for(i=0;i<nparts;i++)
	{
	  if (tmp1[i] > 1) /* if the subgraph contains more than one vertex */
	    {
	      success = create_recursive_preconditioner(out,local_next_unused_vertex+i,next_unused_vertex,curr_entry,sons[i],vertex_perms[i],inv_perm,
							diags[i],diagonal,inst, max_inst,curr_inst+1,taucs_ccs_mtxA,ordering,p1);
	      if (success == 0)
		{
		  taucs_free(diags);
		  taucs_free(weights);
		  taucs_free(tmp1);
		  for(j=0;j<nparts;j++)
		    Metis_struct_free(sons[j]);
		  taucs_free(sons);
		  for(j=0;j<i;j++)
		    taucs_free(vertex_perms[i]);
		  taucs_free(vertex_perms);
		  return(0);
		}
	      
	    }
	  else
	    {
	      
	      /* I sure hope there is no bug in here */
	      out->edges[(*curr_entry)].i=local_next_unused_vertex+i;
	      out->edges[(*curr_entry)].j=vertex_perms[i][0];
	      out->edges[(*curr_entry)].v=weights[i];
	      (*curr_entry)++;
	      
	      out->edges[(*curr_entry)].i=local_next_unused_vertex+i;
	      out->edges[(*curr_entry)].j=local_next_unused_vertex+i;
	      out->edges[(*curr_entry)].v=-2*weights[i];
	      (*curr_entry)++;
	      
	      out->edges[(*curr_entry)].i=vertex_perms[i][0];
	      out->edges[(*curr_entry)].j=vertex_perms[i][0];
	      out->edges[(*curr_entry)].v=diagonal[vertex_perms[i][0]];
	      (*curr_entry)++;
	      Metis_struct_free(sons[i]);
	      p1[ordering_counter_leaves++] = vertex_perms[i][0];
	    }
	}

      if (inst[curr_inst].type == 0) /* Gremban */
	for(j=0;j<nparts;j++)
	  p1[ordering_counter++] = local_next_unused_vertex+j;
      else /* Toledo or Vaidya */
	for(j=0;j<nparts;j++)
	  p1[ordering_counter++] = local_next_unused_vertex+perm_tmp[j];

    }

  taucs_free(perm_tmp);  

  taucs_free(tmp1);
  taucs_free(weights);
  taucs_free(diags);
  
  /* for(j=0;j<nparts;j++) */
    /* Metis_struct_free(sons[j]); */
  taucs_free(sons);
  
  for (i=0;i<nparts;i++)
    taucs_free(vertex_perms[i]);
  taucs_free(vertex_perms);
  
  if (is_root)
    p1[ordering_counter++] = curr_vertex;

  return(1);
}

#if 0
void print_ccs_mat(taucs_ccs_matrix a)
{
  int i, j;
  taucs_printf("%d %d %d\n",a.n,a.n,a.colptr[a.n]);  
  for(i=0;i<a.n;i++)
    for(j=a.colptr[i];j<a.colptr[i+1];j++)
      taucs_printf("%lg %lg %lg\n",(double)i,(double)a.rowind[j],a.values[j]);
 
}
#endif

int is_perm(int *perm,int n)
{
  int *tmp,i;
  tmp = taucs_calloc(n,sizeof(int));
  for(i=0;i<n;i++)
    {
      assert(perm[i] < n);
      if (tmp[perm[i]])
	{
	  printf("NO WAY!!!\n");exit(345);
	  return(0);
	}
      tmp[perm[i]] = 1;
    }
  taucs_free(tmp);
  return(1);
}

void *taucs_sg_preconditioner_create(taucs_ccs_matrix *A,
				     int **perm,
				     int **invperm,
				     char* ordering,
				     char *gremban_command)
{
#ifdef NOMETIS
  return NULL;
#else
  instruction *inst;
  int preconditioner_n,n;
  int *vertex_perm,*inv_perm;
  int tmp=1,tmp2=0,tmp3=0,i,k;
  Metis_struct *Metis_A;
  taucs_ccs_matrix *symmetric_A;
  int next_unused_vertex,curr_entry;
  int success=1;
  graph *out;
  double *diagonal;
  taucs_ccs_matrix *out1;
  multilevel_args* P;
  int depth;
  taucs_ccs_matrix* PGPT;
  void* snL;
  int *p, *ip;

  double wtime_recursive_create, wtime_supernodal_factor,wtime_factor_llt;

  if (gremban_command[0] == 'r')
    {
      char *p,*p1;
      int k,type;
      
      depth = 0;
      p = gremban_command;
      while((p=(strstr(p,":")))!=NULL)
	{
	  p++;
	  depth++;
	}
      if (depth!=2)
	{
	  printf("Command string should have three parts 'regular:GM/CT/VA:number_of_parts_in_each_level'\n");
	  exit(345);
	}
      p = gremban_command;
      
      p1 = strstr(p,":");
      if (p1)
	*p1=0;

      if (strcmp(p,"regular")!=0)
	{
	  printf("Syntax error in Gremban string. Exiting");
	  exit(345);
	}
      
      p = p1+1;
      p1 = strstr(p,":");
      if (p1)
	*p1=0;
      p1++;
	  
      if (strcmp(p,"GM")==0)
	type = 0;
      else
	if (strcmp(p,"CT")==0)
	  type = 1;
	else
	  if (strcmp(p,"VA")==0)
	    type = 2;
	  else
	    {
	      printf("must choose CT or GM or VA. %s. Exiting\n",p);
	      exit(345);
	    }
	  
      sscanf(p1,"%d",&k);

      if (k < 2)
	{
	  printf("Must divide into at least 2 parts at each level. Exiting\n");
	  exit(345);
	}

      depth = (int)(log(A->n)/log(k))+1;
      inst = (instruction *)taucs_malloc(depth*sizeof(instruction));
      for(i=0;i<depth;i++)
	{
	  inst[i].type = type;
	  inst[i].k = k;
	}
    }
  else
    {
      char *p,*p1,*p2;
      int k;
     
      depth = 1;
      p = gremban_command;
      while((p=(strstr(p,":")))!=NULL)
	{
	  p++;
	  depth++;
	}
      if (depth%2)
	{
	  printf("Command string should have 2 strings for each level\n");
	  exit(345);
	}
      p = gremban_command;

      depth = depth/2;
      if (depth<2)
	{
	  printf("Command string should describe a preconditioner of depth 2 at least\n");
	  exit(345);
	}

      inst = (instruction *)taucs_malloc(depth*sizeof(instruction));
      for(i=0;i<depth;i++)
	{
	  p1 = strstr(p,":");
	  if (p1)
	    *p1=0;
	  p1++;
	  
	  if (strcmp(p,"GM")==0)
	    inst[i].type = 0;
	  else
	    if (strcmp(p,"CT")==0)
	      inst[i].type = 1;
	    else
	      if (strcmp(p,"VA")==0)
		inst[i].type = 2;
	  
	      else
		{
		  printf("must choose CT or GM or VA. Exiting\n");
		  exit(345);
		}
	  
	  p2 = strstr(p1,":");
	  if (p2)
	    *p2=0;
	  
	  sscanf(p1,"%d",&k);
	  inst[i].k = k;
	  p = p2+1;
	}
    }

  n = A->n;
  preconditioner_n = (A->n);

  vertex_perm = (int *)taucs_malloc(n*sizeof(int));
  inv_perm    = (int *)taucs_malloc(n*sizeof(int));
  if ((!vertex_perm)||(!inv_perm))
    {
      taucs_free(vertex_perm);
      taucs_free(inv_perm);
      return(NULL);
    }
  
  for(i=0;i<n;i++)
    {
      vertex_perm[i] = i;
      inv_perm[i] = -1;
    }
  
  for(i=0;i<depth;i++)
    {
      tmp3 += tmp;
      tmp2 += min((A->colptr)[n],(tmp*(tmp-1))/2);
      preconditioner_n += tmp;
      tmp *= inst[i].k;
    }
  
  tmp3 += preconditioner_n;
  tmp2 += (A->colptr)[n] + preconditioner_n;
  out = graph_create(tmp2);
  out->n = 0;
  if (out == NULL)
    {
      taucs_free(vertex_perm);
      taucs_free(inv_perm);
      return(NULL);
    }

  Metis_A = taucs_ccs_matrix_to_Metis_struct(A);
  if (Metis_A == NULL)
    {
      taucs_free(vertex_perm);
      taucs_free(inv_perm);
      free_graph(out);
      return(NULL); 
    }
  
  diagonal = (double *)taucs_malloc(n*sizeof(double));
  if (!diagonal)
    {
      taucs_free(vertex_perm);
      taucs_free(inv_perm);
      free_graph(out);
      Metis_struct_free(Metis_A);
      return(NULL); 
    }
  symmetric_A = taucs_ccs_matrix_to_taucs_ccs_matrix(A,diagonal);
  if (symmetric_A == NULL)
    {
      taucs_free(vertex_perm);
      taucs_free(inv_perm);
      /*taucs_free_graph(out); omer*/
			free_graph(out);
      Metis_struct_free(Metis_A);
      taucs_free(diagonal);
      return(NULL); 
    }

  next_unused_vertex = n+1;

  curr_entry = 0;
  
  p = (int *)taucs_malloc(tmp3*sizeof(int));
  if(!p)
    {
      taucs_free(vertex_perm);
      taucs_free(inv_perm);
      /*taucs_free_graph(out); omer*/
			free_graph(out);
      Metis_struct_free(Metis_A);
      taucs_free(diagonal);
      return(NULL); 
    }

  wtime_recursive_create = taucs_wtime();
  success = create_recursive_preconditioner(out,n,&next_unused_vertex,&curr_entry,Metis_A,vertex_perm,inv_perm,
				  0.0,diagonal,inst,depth-1,0,symmetric_A,ordering,p);
  wtime_recursive_create = taucs_wtime()-wtime_recursive_create;
  taucs_printf("\tRecursive Creation time = % 10.3f seconds\n",wtime_recursive_create);

  taucs_free(diagonal);
  taucs_free(vertex_perm);
  taucs_free(inv_perm);
  /* Metis_struct_free(Metis_A); */
  taucs_ccs_free(symmetric_A);

  if(success == 0)
    return(NULL);

  ip = (int *)taucs_malloc(out->n*sizeof(int));
  *perm = (int *)taucs_malloc(n*sizeof(int));
  *invperm = (int *)taucs_malloc(n*sizeof(int));
  if ((!ip)||(!*perm)||(!*invperm))
    {
      taucs_free(ip);
      taucs_free(*perm);
      taucs_free(*invperm);
      taucs_free(p);
      free_graph(out);
      Metis_struct_free(Metis_A);
      return(NULL); 
    }

  /* is_perm(p,out->n); */
  
  for(i=0;i<out->n;i++)
    ip[p[i]] = i;
  
  for(i=0;i<n;i++)
    (*perm)[i] = p[i];

  for(i=0;i<n;i++)
    (*invperm)[(*perm)[i]] = i;


  out->nent=curr_entry;
  taucs_check_diag_dominant_matrix(out,1);

  k = (out->n)-n;

  out1 = graph_to_ccs_matrix(out);

  if (out1 == NULL)
    return(0);

  /* taucs_ccs_write_ijv(A,"A.ijv"); */
  /* taucs_ccs_write_ijv(out1,"G.ijv"); */
  
  PGPT = taucs_ccs_permute_symmetrically(out1,p,ip);
  
  taucs_ccs_free(out1);
  
  P = (multilevel_args*) taucs_malloc(sizeof(multilevel_args));
  if (!P)
    return(NULL);

  wtime_factor_llt = taucs_wtime();
  taucs_printf("taucs_gremban: factoring, preconditioner has %d rows/cols\n",
	       PGPT->n);

  /*  taucs_ccs_write_ijv( PGPT ,"G.ijv");*/
  /*P->L = taucs_ccs_factor_llt(PGPT,0.0,0);*/
  /*  taucs_ccs_write_ijv( P->L ,"L.ijv");*/

  snL = taucs_ccs_factor_llt_mf(PGPT);
  wtime_factor_llt = taucs_wtime()-wtime_factor_llt;
  taucs_printf("\tFactor LL^t time = % 10.3f seconds\n",wtime_factor_llt);
  if (!snL)
    return(NULL);
  wtime_supernodal_factor = taucs_wtime();
  P->L = taucs_supernodal_factor_to_ccs(snL);
  wtime_supernodal_factor = taucs_wtime()-wtime_supernodal_factor;
  taucs_printf("\tSupernodal-factor-to-ccs factor time = % 10.3f seconds\n",wtime_supernodal_factor);
  taucs_supernodal_factor_free(snL);
  
  taucs_free(p);
  taucs_free(ip);
  
  P->Ztilde = (double*) taucs_malloc((n+k) * sizeof(double));
  P->Rtilde = (double*) taucs_malloc((n+k) * sizeof(double));
  if ((!(P->Ztilde))||(!(P->Rtilde)))
    return(NULL);
  P->n = n;
  P->k = k;
  
  /* printf(">>>%d %d\n",(*perm)[0],(*invperm)[0]); */

  return P;
#endif
}

void taucs_sg_preconditioner_free(void* vP) 
{
  multilevel_args* P = (multilevel_args*) vP;

  taucs_free(P->Rtilde);
  taucs_free(P->Ztilde);
  taucs_ccs_free(P->L);
  taucs_free(P);
}

#endif /* TAUCS_CORE_DOUBLE */

