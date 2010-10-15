/* nb2mesh.c
 *     conversion from xxx.noboite to xxx.mesh(b) or xxx.amdba3
 *     or xxx.(no)boite(b)
 *     to compile use gcc nb2mesh.c libmesh.c -o nb2mesh -lm -O
 *
 * Authored by Pascal J. Frey, Inria-Rocquencourt
 * Copyright (c) Inria, 2000-2003.  All rights reserved.
 * Permission is granted to reproduce, use and distribute 
 * this code for any and all purposes, provided that this
 * notice appears in all copies. 
*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "libmesh.h"

#define min(a,b)       ( ((a) < (b)) ? (a) : (b) )
#define max(a,b)       ( ((b) > (a)) ? (b) : (a) )
#define KA     31
#define KB     57
#define KC     79


static int idir[7] = {0,1,2,3,0,1,2};
#ifndef ubyte
typedef unsigned char  ubyte;
#endif

/* hash table structure for tetras */
typedef struct shtable {
  int     nxt,elt,ind;
} Htable;
typedef Htable * pHtable;

typedef struct stetra {
  int   v[4],ref;
  int   adj[4];
  ubyte voy[4];
} Tetra;
typedef Tetra * pTetra;

typedef struct spoint {
  float  c[3];
  int    new;
} Point;
typedef Point * pPoint;

pHtable hash;
int     nhmax,hnext;


pTetra  tets;
pPoint  points;
double  coefs[10];
int    *tri,*ref,*refv,np,ne,nf,npf,npinit;
short   quiet,surf,dosurf = 0,split4 = 0,ddebug = 0,mtype = 0,addref;
short   tonb = 0;
char    namesurf[256];


/* build hash table */
int hcodeTetra(unsigned int key,int mins,int maxs,int sum,int elt,int ind) {
  pHtable  pht;
  pTetra   pt,pt1;
  pPoint   ppt;
  int      mins1,maxs1,sum1;
  ubyte    i1,i2,i3;

  if ( key >= nhmax ) {
    fprintf(stderr," ## hcodeTetra error %d\n",elt);
    return(0);
  }
  pht = &hash[key];

  /* empty bucket */
  if ( !pht->elt ) {
    pht->elt = elt;
    pht->ind = ind;
    pht->nxt = 0;
    return(1);
  }

  /* search linked elts */
  pt = &tets[elt];
  do {
    pt1 = &tets[pht->elt];

    /* compute key */
    i1 = idir[pht->ind+1];
    i2 = idir[pht->ind+2];
    i3 = idir[pht->ind+3];
    mins1 = min(pt1->v[i1],pt1->v[i2]);
    mins1 = min(mins1,pt1->v[i3]);
    maxs1 = max(pt1->v[i1],pt1->v[i2]);
    maxs1 = max(maxs1,pt1->v[i3]);
    if ( pt1->v[i1] != mins1 && pt1->v[i1] != maxs1 )
      sum1  = KB*pt1->v[i1];
    else if ( pt1->v[i2] != mins1 && pt1->v[i2] != maxs1 )
      sum1  = KB*pt1->v[i2];
    else
      sum1  = KB*pt1->v[i3];
    sum1 += KA*mins1 + KC*maxs1;
 
    /* corresponding face */
    if ( mins1 == mins && maxs1 == maxs && sum1 == sum ) {
      if ( pt1->adj[pht->ind] || pt->adj[ind] ) {
        fprintf(stderr,"  ## adjacency problem. exit.\n");
        fprintf(stderr,"  sum %d  key %d min  %d  max %d   elt %d  ind %d\n",
                sum,key,mins,maxs,elt,ind+1);
        fprintf(stderr," elt  %d: %d %d %d %d\n",
                elt,pt->v[0],pt->v[1],pt->v[2],pt->v[3]);
        fprintf(stderr," vois %d: %d %d %d %d\n",
                pht->elt,pt1->v[0],pt1->v[1],pt1->v[2],pt1->v[3]);
        fprintf(stderr," elt:  vois[%d] = %d\n",
                ind,pt->adj[ind]);
        fprintf(stderr," vois: vois[%d] = %d\n",
                pht->ind,pt1->adj[pht->ind]);
        { 
          FILE *fp;
          int   k;
          
          fprintf(stderr,"Saving  debug.mesh\n");
          fp = fopen("debug.mesh","w");
          fprintf(fp,"MeshVersionFormatted 1\n");
          fprintf(fp,"Dimension\n3\n");
          fprintf(fp,"\nVertices\n%d\n",np);
          for (k=1; k<=np; k++) {
            ppt = &points[k];
            fprintf(fp,"%f %f %f 0\n",
                    ppt->c[0],ppt->c[1],ppt->c[2]);
          }
          fprintf(fp,"\nTetrahedra\n3\n");
          fprintf(fp,"%d %d %d %d 1\n",
                  pt->v[0],pt->v[1],pt->v[2],pt->v[3]);
          fprintf(fp,"%d %d %d %d 2\n",
                  pt1->v[0],pt1->v[1],pt1->v[2],pt1->v[3]);
          pt = &tets[pt1->adj[pht->ind]];
          fprintf(fp,"%d %d %d %d 3\n",
                  pt->v[0],pt->v[1],pt->v[2],pt->v[3]);
          fprintf(fp,"\nEnd\n");
          fclose(fp);
        }
        exit(2);
      }

      /* update neighbor */
      pt->adj[ind]       = pht->elt;
      pt->voy[ind]       = (ubyte)pht->ind;
      pt1->adj[pht->ind] = elt;
      pt1->voy[pht->ind] = (ubyte)ind;
      return(1);
    }
    
    /* link element */
    else if ( !pht->nxt ) {
      pht->nxt = hnext;
      pht      = &hash[hnext];
      if ( !pht ) {
        puts(" ## hash table problem. Exit");
	return(0);
      }
      pht->elt = elt;
      pht->ind = ind;
      hnext    = pht->nxt;
      pht->nxt = 0;

      /* check for size overflow */
      if ( !hnext ) {
        puts("  ## overflow in hash table!");
	return(0);
      }
      return(1);
    }
    pht = &hash[pht->nxt];
  } 
  while (1);

  return(0);
}

int hfaceTetra(unsigned int key,int mins,int maxs,int sum) {
  pHtable  pht;
  pTetra   pt1;
  int      mins1,maxs1,sum1;
  ubyte    i1,i2,i3;

  if ( key >= nhmax )  return(0);
  pht = &hash[key];

  if ( !pht->elt )  return(0);
  do {
    pt1 = &tets[pht->elt];

    /* compute key */
    i1 = idir[pht->ind+1];
    i2 = idir[pht->ind+2];
    i3 = idir[pht->ind+3];
    mins1 = min(pt1->v[i1],pt1->v[i2]);
    mins1 = min(mins1,pt1->v[i3]);
    maxs1 = max(pt1->v[i1],pt1->v[i2]);
    maxs1 = max(maxs1,pt1->v[i3]);
    if ( pt1->v[i1] != mins1 && pt1->v[i1] != maxs1 )
      sum1  = KB*pt1->v[i1];
    else if ( pt1->v[i2] != mins1 && pt1->v[i2] != maxs1 )
      sum1  = KB*pt1->v[i2];
    else
      sum1  = KB*pt1->v[i3];
    sum1 += KA*mins1 + KC*maxs1;
 
    /* corresponding face */
    if ( mins1 == mins && maxs1 == maxs && sum1 == sum )
      return(1);
    
    else if ( !pht->nxt )  return(0);
    pht = &hash[pht->nxt];
  } 
  while (1);

  return(0);
}


/* read noboite file */
int loadNoboite(char *name) {
  FILE    *in;
  pTetra   pt;
  pPoint   ppt;
  int      deb,fin,j,no,icube,npbli,npfixe;
  int      nbele,loele,loelef,nbelef,nbpoi,lopoi,lopoif,nbpoif;
  int      nbsub,losub,nbsubf,losubf,nsde,ideb[1024];
  int      nn,bin,t,k,i;
  char     data[256];
  
  /* attempt to open file */
  if ( strstr(name,".noboiteb") || strstr(name,".boiteb") ) {
    in = fopen(name,"r");
    if ( !in ) {
      fprintf(stderr," File %s not found. Bye.\n",name);
      return(0);
    }
    bin = 1;
    fprintf(stdout," Loading file %s\n",name);
    if ( strstr(name,".boiteb") )  dosurf = 1;
  }
  else if ( strstr(name,".noboite") || strstr(name,".boite") ) {
    in = fopen(name,"r");
    if ( !in ) {
      fprintf(stderr," File %s not found. Bye.\n",name);
      return(0);
    }
    bin = 0;
    fprintf(stdout," Loading file %s\n",name);
    if ( strstr(name,".boite") )  dosurf = 1;
  }
  else {
    sprintf(data,"%s.noboiteb",name);
    if ( !quiet )  fprintf(stdout," Checking %s\n",data);
    in = fopen(data,"r");
    bin = 1;
    if ( !in ) {
      sprintf(data,"%s.noboite",name);
      in = fopen(data,"r");
      bin = 0;
      if ( ! quiet )  fprintf(stdout," Checking %s\n",data);
      if ( !in ) {
        sprintf(data,"%s.boiteb",name);
        in = fopen(data,"r");
        bin = 1;
        if ( ! quiet )  fprintf(stdout," Checking %s\n",data);
        if ( !in ) {
          sprintf(data,"%s.boite",name);
          in = fopen(data,"r");
          bin = 0;
          if ( !quiet )  fprintf(stdout," Checking %s\n",data);
        }
        dosurf = 1;
      }
    }
    if ( !in ) {
      fprintf(stderr," File %s not found. Bye.\n",data);
      return(0);
    }
    fprintf(stdout," Loading file %s\n",data);
  }

  /* read binary file */
  if ( bin == 1 ) {
    fread(&no,sizeof(int),1,in);
    fread(&ne,sizeof(int),1,in);
    fread(&np,sizeof(int),1,in);

    fread(&npfixe,sizeof(int),1,in);
    fread(&icube,sizeof(int),1,in);
    fread(&npbli,sizeof(int),1,in);

    fread(&nbele,sizeof(int),1,in);
    fread(&loele,sizeof(int),1,in);
    fread(&nbelef,sizeof(int),1,in);
    fread(&loelef,sizeof(int),1,in);

    fread(&nbpoi,sizeof(int),1,in);
    fread(&lopoi,sizeof(int),1,in);
    fread(&nbpoif,sizeof(int),1,in);
    fread(&lopoif,sizeof(int),1,in);

    fread(&nbsub,sizeof(int),1,in);
    fread(&losub,sizeof(int),1,in);
    fread(&nbsubf,sizeof(int),1,in);
    fread(&losubf,sizeof(int),1,in);

    fread(&no,sizeof(int),1,in);

    if ( !quiet )
      fprintf(stdout,"  Header: %d %d   %d %d\n",ne,np,npfixe,npbli);
    if ( ne+np == 0 ) {
      fclose(in);
      fprintf(stderr,"  Sorry, no element found. Bye\n");
      return(0);
    }

    /* second record: read tetrahedra */
    tets = (Tetra*)calloc(max(1000,(ne+1)),sizeof(Tetra));
    if ( !tets ) {
      fprintf(stderr,"  ## Not enough memory. Bye.\n");
      exit(2);
    }
    if ( !quiet )  fprintf(stdout,"  Read tetras\n");
    deb = 1;
    fin = loele;
    k   = 1;
    for (j=1; j<=nbele; j++) {
      fread(&no,sizeof(int),1,in);
      for (t=deb; t<=fin; t+=4) {
        pt = &tets[k++];
        nn = fread(&pt->v[0],sizeof(int),4,in);
      }
      fread(&no,sizeof(int),1,in);
      deb += loele;
      fin += loele;
    }
    if ( nbelef != 0 ) {
      fin = deb + loelef - 1;
      fread(&no,sizeof(int),1,in);
      for (t=deb; t<=fin; t+=4) {
        pt = &tets[k++];
        nn = fread(&pt->v[0],sizeof(int),4,in);
      }
      fread(&no,sizeof(int),1,in);
    }

    /* third record: read vertices */
    if ( !quiet )  fprintf(stdout,"  Read vertices\n");
    points = (Point*)malloc(max(1000,np+1)*sizeof(Point));
    assert(points);

    deb = 1;
    fin = lopoi;
    k   = 1;
    for (j=1; j<=nbpoi; j++) {
      fread(&no,sizeof(int),1,in);
      for (t=deb; t<=fin; t+=3) {
        ppt = &points[k++];
        nn = fread(&ppt->c[0],sizeof(float),3,in);
        ppt->new = -1;
      }
      fread(&no,sizeof(int),1,in);
      deb += lopoi;
      fin += lopoi;
    }
    if ( nbpoif != 0 ) {
      fin = deb + lopoif - 1;
      fread(&no,sizeof(int),1,in);
      for (t=deb; t<=fin; t+=3) {
        ppt = &points[k++];
        nn = fread(&ppt->c[0],sizeof(float),3,in);
        ppt->new = -1;
      }
      fread(&no,sizeof(int),1,in);
    }
    
    /* 4th : sub-domains */
    fread(&no,sizeof(int),1,in);
    fread(&nsde,sizeof(int),1,in);
    fread(&no,sizeof(int),1,in);

    /* 5th+ 6th : sub-domains */
    fread(&no,sizeof(int),1,in);
    if ( nsde == 1 )
      fread(&ideb,sizeof(int),3,in);
    else
      fread(&ideb,sizeof(int),3*nsde,in);
    fread(&no,sizeof(int),1,in);
    
    fread(&no,sizeof(int),1,in);
    if ( nsde == 1 )
      for (t=1; t<=3; t++)
        fread(&ideb,sizeof(int),1,in);
    else {
      deb = 1;
      fin = losub;
      for (t=1; t<=fin; t+=3) {
      }
    }
    fread(&no,sizeof(int),1,in);
  }

  /* read ascii file */
  else {
    fscanf(in,"%d %d %d %d %d",&ne,&np,&npfixe,&icube,&npbli);
    fscanf(in,"%d %d %d %d",&nbele,&loele,&nbelef,&loelef);
    fscanf(in,"%d %d %d %d",&nbpoi,&lopoi,&nbpoif,&lopoif);
    fscanf(in,"%d %d %d %d",&nbsub,&losub,&nbsubf,&losubf);

    if ( !quiet )
      fprintf(stdout,"  Header: %d %d   %d %d\n",ne,np,npfixe,npbli);
    if ( ne+np == 0 ) {
      fclose(in);
      fprintf(stderr,"  Sorry, no element found. Bye\n");
      return(0);
    }

    /* second record: read tetrahedra */
    if ( !quiet )  fprintf(stdout,"  Read tetras\n");
    tets = (Tetra *)calloc(max(1000,ne+1),sizeof(Tetra));
    if ( !tets ) {
      fprintf(stderr," ## Not enough memory. Bye.\n");
      exit(2);
    }
    for (j=1; j<=ne; j++) {
      pt = &tets[j];
      fscanf(in,"%d %d %d %d",&pt->v[0],&pt->v[1],&pt->v[2],&pt->v[3]);
    }
    
    /* third record: read vertices */
    if ( !quiet )  fprintf(stdout,"  Read vertices\n");
    points = (Point*)malloc(max(1000,np+1)*sizeof(Point));
    assert(points);
    for (j=1; j<=np; j++) {
      ppt = &points[j];
      fscanf(in,"%f %f %f",&ppt->c[0],&ppt->c[1],&ppt->c[2]);
      ppt->new = -1;
    }
    
    /* 4th : nb subdomains */
    fscanf(in,"%d",&nsde);

    /* 5th : sub-domains */
    if ( nsde == 1 )
      for (k=deb; k<=3; k++)
      fscanf(in,"%d ",&deb);
    else
      for (k=1; k<=3*nsde; k++) {
        fscanf(in,"%d ",&deb);
      }

    for (k=1; k<=ne; k++) {
      pt = &tets[k];
      fscanf(in,"%d",&pt->ref);
    }
    if ( nsde == 1 )
      for (t=1; t<=3; t++)
        fscanf(in,"%d",&ideb);
    else {
      deb = 1;
      fin = losub;
      for (t=1; t<=nbsub; t++) {
        for (i=deb; i<=fin; i++) {
          fscanf(in,"%d",&k);
          pt = &tets[k];
          pt->ref = t;
        }
      }
    }
  }

  fclose(in);
  return(1);
}

int saveNoboite(char *name) {
  FILE    *out;
  pTetra   pt;
  pPoint   ppt;
  int      deb,fin,j,no,icube,npbli,npfixe;
  int      nbele,loele,loelef,nbelef,nbpoi,lopoi,lopoif,nbpoif,
           nerema,isrec;
  int      nbhlo,nbhlof,lohlo,lohlof;
  int      bin,k,l;

  /* attempt to open file */
  if ( strstr(name,".noboiteb") || strstr(name,".boiteb") ) {
    out = fopen(name,"w");
    if ( !out ) {
      fprintf(stderr," Unable to open file %s. Bye.\n",name);
      return(0);
    }
    bin = 1;
    fprintf(stdout," Writing file %s\n",name);
  }
  else if ( strstr(name,".noboite") || strstr(name,".boite") ) {
    out = fopen(name,"w");
    if ( !out ) {
      fprintf(stderr," Unable to open file %s. Bye.\n",name);
      return(0);
    }
    bin = 0;
    fprintf(stdout," Writing file %s\n",name);
  }
  else  return(0);

  /* write binary file */
  nerema = 12 * 16384;
  isrec  = 0;

  /* split record or not */
  if (ne <= nerema) {
    nbele = 1; loele = 4*ne; nbelef = loelef = 0;
    nbpoi = 1; lopoi = 3*np; nbpoif = lopoif = 0;
    nbhlo = 1; lohlo = np;   nbhlof = lohlof = 0;
  }
  else {
    if ( !quiet )  fprintf(stdout,"  splitting records %d",nerema);
    /* nbpoi records with nerema values */
    nbpoi  = 3*np / nerema;   lopoi  = nerema;
    nbpoif = 1;               lopoif = 3*np - nerema*nbpoi;
    nbhlo  = np/nerema;       lohlo  = nerema;
    nbhlof = 1;               lohlof = np - nerema*nbhlo;
    if ( np <= nerema ) {
      nbhlo  = 1;  lohlo = np;
      nbhlof = 0;  lohlof = 0;
    }
    nbele  = 4*ne/nerema;  loele  = nerema; 
    nbelef = 1;            loelef = 4*ne - nerema*nbele;
  }
  npfixe= npf;
  icube = 0;
  npbli = 0;
    
  if ( bin ) {
    no = 68;
    fwrite(&no,sizeof(int),1,out);
    fwrite(&ne,sizeof(int),1,out);
    fwrite(&np,sizeof(int),1,out);

    fwrite(&npfixe,sizeof(int),1,out);
    fwrite(&icube,sizeof(int),1,out);
    fwrite(&npbli,sizeof(int),1,out);

    fwrite(&nbele,sizeof(int),1,out);  fwrite(&loele,sizeof(int),1,out);
    fwrite(&nbelef,sizeof(int),1,out); fwrite(&loelef,sizeof(int),1,out);

    fwrite(&nbpoi,sizeof(int),1,out);  fwrite(&lopoi,sizeof(int),1,out);
    fwrite(&nbpoif,sizeof(int),1,out); fwrite(&lopoif,sizeof(int),1,out);

    fwrite(&nbhlo,sizeof(int),1,out);  fwrite(&lohlo,sizeof(int),1,out);
    fwrite(&nbhlof,sizeof(int),1,out); fwrite(&lohlof,sizeof(int),1,out);

    fwrite(&no,sizeof(int),1,out);
    isrec++;

    /* write simplices */
    if ( !quiet )  fprintf(stdout,"  Writing simplices\n");
    deb = 1;
    fin = loele;
    l = 1;
    for (j=1; j<=nbele; j++) {
      no = (fin-deb+1)*sizeof(int);
      fwrite(&no,sizeof(int),1,out);
      for (k=deb; k<=fin; k+=4) {
        pt = &tets[l++];
        fwrite(&pt->v,sizeof(int),4,out);
      }
      deb += loele;
      fin += loele;
      isrec++;
      fwrite(&no,sizeof(int),1,out);
    }
    if ( deb < 4*ne ) {
      no = (4*ne-deb+1)*sizeof(int);
      fwrite(&no,sizeof(int),1,out);
      for (k=deb; k<=4*ne; k+=4) {
        pt = &tets[l++];
        fwrite(&pt->v,sizeof(int),4,out);
      }
      fwrite(&no,sizeof(int),1,out);
    }
    
    /* write coordinates */
    if ( !quiet )  fprintf(stdout,"  Writing coords\n");
    deb = 1;
    fin = lopoi;
    l   = 1;
    for (j=1; j<=nbpoi; j++) {
      no = (fin-deb+1)*sizeof(float);
      fwrite(&no,sizeof(int),1,out);
      for (k=deb; k<=fin; k+=3) {
        ppt = &points[l++];
        fwrite(&ppt->c,3*sizeof(float),1,out);
      }
      deb += lopoi;
      fin += lopoi;
      isrec++;
      fwrite(&no,sizeof(int),1,out);
    }
    if ( deb < 3*np ) {
      no = (3*np-deb+1)*sizeof(float);
      fwrite(&no,sizeof(int),1,out);
      for (k=deb; k<=3*np; k+=3) {
        ppt = &points[l++];
        fwrite(&ppt->c,3*sizeof(float),1,out);
      }
      isrec++;
      fwrite(&no,sizeof(int),1,out);
    }
    
    /* write coefficients */
    no = 6*sizeof(double);
    fwrite(&no,sizeof(int),1,out);
    fwrite(&coefs,sizeof(double),6,out);
    
    fwrite(&no,sizeof(int),1,out);
    isrec++;
  }
  
  /* write ASCII file */
  else {
    fprintf(out,"%d %d ",ne,np);
    fprintf(out,"%d %d %d ",npfixe,icube,npbli);
    fprintf(out,"%d %d %d %d ",nbele,loele,nbelef,loelef);
    fprintf(out,"%d %d %d %d ",nbpoi,lopoi,nbpoif,lopoif);
    fprintf(out,"%d %d %d %d\n",nbhlo,lohlo,nbhlof,lohlof);
      
    /* write simplices */
    if ( !quiet )  fprintf(stdout,"  Writing simplices\n");
    deb = 1;
    fin = loele;
    l   = 1;
    for (j=1; j<=nbele; j++) {
      for (k=deb; k<=fin; k+=4) {
        pt = &tets[l++];
        fprintf(out,"%d %d %d %d",pt->v[0],pt->v[1],pt->v[2],pt->v[3]);
        if ( k % 12 == 0 )  fprintf(out,"\n");
      }
      if ( k % 12 != 0 )  fprintf(out,"\n");
      deb += loele;
      fin += loele;
    }
    if ( nbelef ) {
      fin = deb + loelef - 1;
      for (k=deb; k<=fin; k+=4) {
        pt = &tets[l++];
        fprintf(out,"%d %d %d %d",pt->v[0],pt->v[1],pt->v[2],pt->v[3]);
        if (k % 12 == 0)  fprintf(out,"\n");
      }
      if (k % 12 != 0) fprintf(out,"\n");
    }
  
    /* write coordinates */
    if ( !quiet )  fprintf(stdout,"  Writing coords\n");
    deb = 1;
    fin = lopoi;
    l   = 1;
    for (j=1; j<=nbpoi; j++) {
      for (k=deb; k<=fin; k+=3) {
        ppt = &points[l++];
        fprintf(out,"%f %f %f",ppt->c[0],ppt->c[1],ppt->c[2]);
        if (k % 9 == 0)  fprintf(out,"\n");
      }
      if (k % 9 != 0)  fprintf(out,"\n");
      deb += lopoi;
      fin += lopoi;
    }
    if ( nbpoif ) {
      fin = deb + lopoif - 1;
      for (k=deb; k<=fin; k+=3) {
        ppt = &points[l++];
        fprintf(out,"%f %f %f",ppt->c[0],ppt->c[1],ppt->c[2]);
        if (k % 9 == 0)  fprintf(out,"\n");
      }
      if (k %9 != 0) fprintf(out,"\n");
    }

    /* write coefficients */
    for (k=1; k<=6; k++) {
      fprintf(out,"%g ",coefs[k]);
    }
    fprintf(out,"\n");
  }
  
  fclose(out);
  return(1);
}


/* load surface description */
int loadSurf(char *name) {
  FILE   *in;
  float   dummyf;
  int     ver,k,msh,degre,nef,nt,nr;
  int     a,b,c,d,e,rf,post,posq,posp;
  char   *ptr,data[256],buf[256];
  
  strcpy(data,name);
  ptr = strstr(data,".noboite");
  if ( ptr ) *ptr = '\0';
  ptr = strstr(data,".boite");
  if ( ptr ) *ptr = '\0';
  strcpy(buf,data);
  ptr = strstr(data,".mesh");
  if ( ptr ) *ptr = '\0';
  strcpy(buf,data);

  msh = 0;
  nf  = nr = npf = 0;
  strcat(data,".faces");
  if ( !quiet )  fprintf(stdout," Checking %s\n",data);
  in = fopen(data,"r");
  if ( in )
    fprintf(stdout," Loading file %s\n",data);
  else {
    strcpy(data,buf);
    strcat(data,".meshb");
    if ( !quiet )  fprintf(stdout," Checking %s\n",data);
    in = ouvrir_mesh(data,"r",&ver);
    msh = 1;
    if ( !in ) {
      strcpy(data,buf);
      strcat(data,".mesh");
      if ( !quiet )  fprintf(stdout," Checking %s\n",data);
      in = ouvrir_mesh(data,"r",&ver);
      if ( !in ) {
        fprintf(stdout,"  No surface file found\n");
        return(0);
      }
    }
    else
    fprintf(stdout," Loading file %s\n",name);
  }

  /* read .faces */
  if ( !msh ) {
    fgets(data,255,in);
    sscanf(data,"%d",&nef);
    if ( !nef ) {
      fprintf(stderr,"  Sorry, no triangles found.\n");
      return(0);
    }
    surf = 1;

    /* read once */
    for (k=1; k<=nef; k++) {
      fscanf(in,"%d",&degre);
      if ( degre < 3 || degre > 4 ) {
	fprintf(stderr,"  Wrong degree %d\n",degre);
	return(0);
      }
      else if ( degre == 3 )
	nf++;
      else if ( degre == 4 )
	nf += 2;
      fgets(data,80,in);
    }

    /* alloc memory */
    if ( !quiet )  fprintf(stdout,"  Triangles   %d\n",nf);
    tri = malloc(max(1000,nf+1)*3*sizeof(int));
    if ( ! tri ) {
      fprintf(stderr,"  Sorry, not enough memory. Bye\n");
      return(0);
    }
    ref = malloc(max(1000,nf+1)*sizeof(int));
    if ( !ref ) {
      fprintf(stderr,"  Sorry, not enough memory. Bye\n");
      return(0);
    }
    
    /* read faces */
    nf = 0;
    nt = 0;
    rewind(in);
    fscanf(in,"%d",&nef);
    for (k=1; k<=nef; k++) {
      fscanf(in,"%d",&degre);
      if ( degre == 3 ) {
	fscanf(in,"%d %d %d %d %d %d %d\n",&a,&b,&c,&rf,&d,&d,&d);
	tri[++nt] = a;
	tri[++nt] = b;
	tri[++nt] = c;
	ref[++nf] = rf;
      }
      else if ( degre == 4 ) {
	fscanf(in,"%d %d %d %d %d %d %d %d %d",&a,&b,&c,&d,&rf,&e,&e,&e,&e);
        tri[++nt] = a;
        tri[++nt] = b;
        tri[++nt] = c;
        ref[++nf] = rf;

	tri[++nt] = a;
        tri[++nt] = c;
        tri[++nt] = d;
        ref[++nf] = rf;
      }
      else
	fgets(data,80,in);
    }
    fclose(in);
    
    /* read points */
    if ( split4) {
      strcpy(data,name);
      ptr = strstr(data,".noboite");
      if ( ptr ) *ptr = '\0';
      ptr = strstr(data,".boite");
      if ( ptr ) *ptr = '\0';
      strcpy(buf,data);
      strcat(buf,".points");
      in = fopen(buf,"r");
      fscanf(in,"%d",&npf);
      fclose(in);
    }
  }

  /* read .mesh */
  else {
    mtype = 1;
    strcpy(namesurf,name);
    posp = chercher_mot_clef(in,Vertices,0);
    fseek(in,posp,SEEK_SET);
    npf    = lire_int(in);
    npinit = npf;
    refv = malloc((npf+1)*sizeof(int));
    if ( !refv ) {
      fprintf(stderr,"  Sorry, not enough memory. Bye\n");
      return(0);
    }
    for (k=1; k<=npf; k++) {
      dummyf  = lire_reel(in);
      dummyf  = lire_reel(in);
      dummyf  = lire_reel(in);
      refv[k] = lire_int(in);
    }
    post = chercher_mot_clef(in,Triangles,0);
    if ( post ) {
      fseek(in,post,SEEK_SET);
      nf = lire_int(in);
    }
    posq = chercher_mot_clef(in,Quadrilaterals,0);
    if ( posq ) {
      fseek(in,posq,SEEK_SET);
      nf += 2*lire_int(in);
    }
    if ( !nf ) {
      fprintf(stderr,"  Sorry, no triangles found.\n");
      return(0);
    }

    /* alloc memory */
    if ( !quiet )  fprintf(stdout,"  Triangles   %d\n",nf);
    tri = malloc(max(1000,nf+1)*3*sizeof(int));
    if ( ! tri ) {
      fprintf(stderr,"  Sorry, not enough memory. Bye\n");
      return(0);
    }
    ref = malloc(max(1000,nf+1)*sizeof(int));
    if ( !ref ) {
      fprintf(stderr,"  Sorry, not enough memory. Bye\n");
      return(0);
    }
    surf = 1;

    /* read mesh */
    nt = nf = 0;
    if ( post ) {
      fseek(in,post,SEEK_SET);
      nef = lire_int(in);
      for (k=1; k<=nef; k++) {
	tri[++nt] = lire_int(in);
	tri[++nt] = lire_int(in);
	tri[++nt] = lire_int(in);
	ref[++nf] = lire_int(in);
      }
    }
    if ( posq ) {
      fseek(in,posq,SEEK_SET);
      nef = lire_int(in);
      for (k=1; k<=nef; k++) {
	a = lire_int(in);
	b = lire_int(in);
	c = lire_int(in);
	d = lire_int(in);
	rf = lire_int(in);

	tri[++nt] = a;
	tri[++nt] = b;
	tri[++nt] = c;
	ref[++nf] = rf;

	tri[++nt] = a;
	tri[++nt] = c;
	tri[++nt] = d;
	ref[++nf] = rf;
      }
    }
    fclose(in);
  }

  return(1);
}


/* hash tetras faces */
int hashFaces() {
  pTetra    pt;
  int       i,i1,i2,i3,k,nt,nbel,mins,maxs;
  int       sum;
  unsigned  int  key;

  if ( hash )   return(1);
  if ( !quiet ) fprintf(stdout,"\n Hashing bdry faces\n");

  /* alloc mem */
  nbel  = 1*ne;
  nhmax = 3*ne;
  hash  = (Htable*)calloc(nhmax+1,sizeof(Htable));
  if ( !hash ) {
    fprintf(stderr," ## Not enough memory to build hash table.\n");
    return(0);
  }

  /* init hash table */
  hnext = nbel;
  for (k=hnext; k<nhmax; k++)
    hash[k].nxt = k+1;

  /* build hash table */
  for (k=1; k<=ne; k++) {
    pt = &tets[k];
    if ( !pt->v[0] )  continue;
    for (i=0; i<4; i++) {
      i1 = idir[i+1];
      i2 = idir[i+2];
      i3 = idir[i+3];

      mins = min(pt->v[i1],pt->v[i2]);
      mins = min(mins,pt->v[i3]);
      maxs = max(pt->v[i1],pt->v[i2]);
      maxs = max(maxs,pt->v[i3]);

      /* compute key */
      if ( pt->v[i1] != mins && pt->v[i1] != maxs )
        sum  = KB*pt->v[i1];
      else if ( pt->v[i2] != mins && pt->v[i2] != maxs )
        sum  = KB*pt->v[i2];
      else
        sum  = KB*pt->v[i3];
      sum += KA*mins + KC*maxs;
      key  = sum % nbel;

      pt->adj[i] = pt->voy[i] = 0;
      if ( !hcodeTetra(key,mins,maxs,sum,k,i) ) {
        fprintf(stderr," ## hashTetra problem %d / %d. Exit.\n",k,ne);
        free(hash);
        return(0);
      }
    }
  }

  nf = nt = 0;
  for (k=1; k<=ne; k++) {
    pt = &tets[k];
    for (i=0; i<4; i++)
      if ( !pt->adj[i] )  nf++;
  }

  if ( !quiet )  printf("  allocate %d faces\n",nf);
  tri = (int*)calloc((nf+1)*3,sizeof(int));
  if ( !tri )  return(0);

  nf = nt = 0;
  for (k=1; k<=ne; k++) {
    pt = &tets[k];
    for (i=0; i<4; i++) {
      if ( !pt->adj[i] ) {
        tri[++nt] = pt->v[idir[i+1]];
        tri[++nt] = pt->v[idir[i+2]];
        tri[++nt] = pt->v[idir[i+3]];
        nf++;
      }
    }
  }

  if ( !quiet )  fprintf(stdout,"   Total bdry faces %d\n",nf);
  return(1);
}

/* hash tetras faces */
int hashTrias() {
  pTetra    pt;
  int       i,i1,i2,i3,k,nbel,mins,maxs;
  int       sum;
  unsigned  int  key;

  if ( hash )   return(1);
  if ( !quiet ) fprintf(stdout,"\n Hashing bdry faces\n");

  /* alloc mem */
  nbel  = 1*ne;
  nhmax = 3*ne;
  hash  = (Htable*)calloc(nhmax+1,sizeof(Htable));
  if ( !hash ) {
    fprintf(stderr," ## Not enough memory to build hash table.\n");
    return(0);
  }

  /* init hash table */
  hnext = nbel;
  for (k=hnext; k<nhmax; k++)
    hash[k].nxt = k+1;

  /* build hash table */
  for (k=1; k<=ne; k++) {
    pt = &tets[k];
    if ( !pt->v[0] )  continue;
    for (i=0; i<4; i++) {
      if ( pt->adj[i] )  continue;
      i1 = idir[i+1];
      i2 = idir[i+2];
      i3 = idir[i+3];

      mins = min(pt->v[i1],pt->v[i2]);
      mins = min(mins,pt->v[i3]);
      maxs = max(pt->v[i1],pt->v[i2]);
      maxs = max(maxs,pt->v[i3]);

      /* compute key */
      if ( pt->v[i1] != mins && pt->v[i1] != maxs )
        sum  = KB*pt->v[i1];
      else if ( pt->v[i2] != mins && pt->v[i2] != maxs )
        sum  = KB*pt->v[i2];
      else
        sum  = KB*pt->v[i3];
      sum += KA*mins + KC*maxs;
      key  = sum % nbel;

      pt->adj[i] = pt->voy[i] = 0;
      if ( !hcodeTetra(key,mins,maxs,sum,k,i) ) {
        fprintf(stderr," ## hashTetra problem %d / %d. Exit.\n",k,ne);
        free(hash);
        return(0);
      }
    }
  }
  return(1);
}


#define EPS   0.0f

/* compute tet quality */
int qualtet(int k) {
  pTetra    pt;
  pPoint    p0,p1,p2,p3;
  double    ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz;
  double    vx,vy,vz,vol;
  
  pt = &tets[k];
  p0 = &points[pt->v[0]];
  p1 = &points[pt->v[1]];
  p2 = &points[pt->v[2]];
  p3 = &points[pt->v[3]];
  
  ax = p0->c[0];
  ay = p0->c[1];
  az = p0->c[2];
  bx = p1->c[0] - ax;  
  by = p1->c[1] - ay;  
  bz = p1->c[2] - az;
  cx = p2->c[0] - ax;  
  cy = p2->c[1] - ay;  
  cz = p2->c[2] - az;
  dx = p3->c[0] - ax;  
  dy = p3->c[1] - ay;  
  dz = p3->c[2] - az;

  vx = cy*dz - dy*cz;
  vy = cz*dx - dz*cx;
  vz = cx*dy - dx*cy;
  vol = bx*vx + by*vy + bz*vz;
  if ( vol <= EPS ) {
    fprintf(stdout,"   tet %d: %d %d %d %d,  vol %E\n",
            k,pt->v[0],pt->v[1],pt->v[2],pt->v[3],vol);
    return(1);
  }
  return(0);
}

void splitTetra() {
  pTetra  pt,pt1;
  pPoint  ppt;
  float   ax,ay,az;
  int     neinit,npc,nnpf,j,k;
  
  /* peculiar treatment */
  npc = 0;
  for (k=1; k<=ne; k++) {
    pt = &tets[k];
    nnpf  = pt->v[0] <= npf;
    nnpf += pt->v[1] <= npf;
    nnpf += pt->v[2] <= npf;
    nnpf += pt->v[3] <= npf;
    if ( nnpf == 4 )  npc++;
  }
  if ( !npc )  return;

  points = (Point*)realloc(points,(np+npc+10)*sizeof(Point));
  assert(points);
  tets = (Tetra*)realloc(tets,(ne+4*npc+10)*sizeof(Tetra));
  if ( !tets )  exit(2);
  neinit = ne;
  npc = 0;
  for (k=1; k<=neinit; k++) {
    pt = &tets[k];
    nnpf  = pt->v[0] <= npf;
    nnpf += pt->v[1] <= npf;
    nnpf += pt->v[2] <= npf;
    nnpf += pt->v[3] <= npf;
    if ( nnpf < 4 )  continue;
    npc++; 
    ax = ay = az = 0.0f;
    for (j=0; j<4; j++) {
      ppt = &points[pt->v[j]];
      ax += ppt->c[0];
      ay += ppt->c[1];
      az += ppt->c[2];
    }

    ppt = &points[++np];
    ppt->c[0] = ax * 0.25;
    ppt->c[1] = ay * 0.25; 
    ppt->c[2] = az * 0.25;
 
    pt1 = &tets[++ne];
    pt1->v[0] = pt->v[1];
    pt1->v[1] = pt->v[3];
    pt1->v[2] = pt->v[2];
    pt1->v[3] = np;
    
    pt1 = &tets[++ne];
    pt1->v[0] = pt->v[0];
    pt1->v[1] = pt->v[2];
    pt1->v[2] = pt->v[3];
    pt1->v[3] = np;
    
    pt1 = &tets[++ne];
    pt1->v[0] = pt->v[0];
    pt1->v[1] = pt->v[3];
    pt1->v[2] = pt->v[1];
    pt1->v[3] = np;

    pt->v[3] = np;
  }
  if ( npc )  fprintf(stdout,"  %d corrected\n",npc);
}  

/* load mesh file */
int loadMesh(char *name) {
  FILE    *in;
  pTetra   pt;
  pPoint   ppt;
  int      k,dim,ver,reff,pos[NbMc],key;

  in = ouvrir_mesh(name,"r",&ver);
  if ( !in ) {
    fprintf(stderr,"\n File %s not found. Bye.\n",name);
    return(0);
  }
  fprintf(stdout," Loading file %s\n",name);
  dosurf = 0;

  /* parse keywords */
  dim = 3;
  for (k=0; k<NbMc; k++)  pos[k] = 0;
  do {
    key = mot_clef_suivant(in);
    if ( key == End ) break;
    pos[key] = ftell(in);
    switch(key) {
    case MeshDimension:
      dim = lire_int(in);
      break;
    case Vertices:
      np = lire_int(in);
      break;
    case Tetrahedra:
      ne = lire_int(in);
      break;
    }
  }
  while (key != End);

  if ( dim != 3 || !np || !ne ) {
    fprintf(stderr,"  Wrong data type.\n");
    return(0);
  }
  
  /* allocation */
  tets = (Tetra*)calloc(max(1000,(ne+1)),sizeof(Tetra));
  assert(tets);
  points = (Point*)malloc(max(1000,np+1)*sizeof(Point));
  assert(points);

  fseek(in,pos[Vertices],SEEK_SET);
  np = lire_int(in);
  refv = malloc((np+1)*sizeof(int));
  if ( !refv ) {
    fprintf(stderr,"Sorry not enough memory. Bye\n");
    return(0);
  }
  for(k=1; k<=np; k++) {
    ppt = &points[k];
    ppt->c[0] = (float)lire_reel(in);
    ppt->c[1] = (float)lire_reel(in);
    ppt->c[2] = (float)lire_reel(in);
    refv[k]   = lire_int(in);
  }

  fseek(in,pos[Tetrahedra],SEEK_SET);
  ne = lire_int(in);
  for(k=1; k<=ne; k++) {
    pt = &tets[k];
    pt->v[0] = lire_int(in);
    pt->v[1] = lire_int(in);
    pt->v[2] = lire_int(in);
    pt->v[3] = lire_int(in);
    reff     = lire_int(in);
  }
  
  coefs[0] = coefs[2] = coefs[4] = 0.0;
  coefs[1] = coefs[3] = coefs[5] = 1.0;
 
  return(1);
}

/* remove dangling faces */
int delDangling() {
  int          f,k,kk,a,b,c,mins,maxs,sum,nfdeb;
  unsigned int key;

  if ( !nf || hash )  return(1);
  fprintf(stdout,"\n Remove dangling faces\n");
  hashTrias();

  nfdeb = nf;
  k = f = 1;
  do {
    a = tri[k+0];
    b = tri[k+1];
    c = tri[k+2];

    mins = min(a,b);
    mins = min(mins,c);
    maxs = max(a,b);
    maxs = max(maxs,c);
    if ( a != mins && a != maxs )
      sum  = KB*a;
    else if ( b != mins && b != maxs )
      sum  = KB*b;
    else
      sum  = KB*c;
    sum += KA*mins + KC*maxs;
    key  = sum % ne;

    if ( !hfaceTetra(key,mins,maxs,sum) ) {
      kk = 3*(nf-1) + 1;
      tri[k+0] = tri[kk+0];
      tri[k+1] = tri[kk+1];
      tri[k+2] = tri[kk+2];
      ref[f]   = ref[nf];
      nf--;
    }
    else {
      k += 3;
      f++;
    }
  }
  while ( k <=3*nf );
  if ( !quiet )  fprintf(stdout,"  %d faces deleted\n",nfdeb-nf);

  free(hash);
  hash = 0;
  return(1);
}
  

/* write mesh file */
int saveMesh(char *name) {
  FILE   *out,*in;
  pTetra  pt;
  pPoint  ppt,p1,p2,p3;
  float   nx,ny,nz;
  int     i,k,l,ver,nnp,nnf,nn,nv,nnv,ntv;
  int    *edg,pos,posnv,nc,ncc,na,naa,cc,aa,bb,rr;
  char   *ptr,namein[256];

  /* output file */
  out = ouvrir_mesh(name,"w",&ver);
  if ( !out ) {
    fprintf(stderr,"\n Unable to open %s. Bye.\n",name);
    return(0);
  }
  fprintf(stdout,"\n Writing %s\n",name);

  /* write mesh header */
  ecrire_commentaire(out,"Mesh generated by nb2mesh (INRIA)");
  formater(out);
  ecrire_mot_clef(out,MeshDimension);
  ecrire_int(out,3);
  formater(out);
  formater(out);

  /* tassage */
  nnp = nnf = 0;
  for (k=1; k<=ne; k++) {
    pt = &tets[k];
    for (i=0; i<4; i++) {
      ppt = &points[pt->v[i]];
      ppt->new = 0;
    }
  }

  for(k=1; k<=np; k++) {
    ppt = &points[k];
    if ( !ppt->new )  ppt->new = ++nnp;
  }

  for (k=1; k<=3*nf; k+=3) {
    p1 = &points[tri[k]];
    p2 = &points[tri[k+1]];
    p3 = &points[tri[k+2]];
    if ( p1->new < 0 || p2->new < 0 || p3->new < 0 )  continue;
    ++nnf;
  }

  /* vertices */
  if ( !quiet && np > 10000 )
    fprintf(stdout,"  Vertices    %8d / %8d\n",nnp,np);
  ecrire_commentaire(out,"Set of mesh vertices");
  ecrire_mot_clef(out,Vertices);
  ecrire_int(out,nnp);
  formater(out);
  for (k=1; k<=np; k++) {
    ppt = &points[k];
    if ( ppt->new < 0 )  continue;
    ecrire_reel(out,ppt->c[0]);
    ecrire_reel(out,ppt->c[1]);
    ecrire_reel(out,ppt->c[2]);
    if ( k <= npinit ) 
      ecrire_int(out,refv[k]);
    else
      ecrire_int(out,0);
    formater(out);
  }
  formater(out);

  /* triangles */
  if ( nnf ) {
    if ( !quiet && np > 10000 )  
      fprintf(stdout,"  Triangles   %8d / %8d\n",nnf,nf);
    ecrire_commentaire(out,"Set of Triangles");
    ecrire_mot_clef(out,Triangles);
    ecrire_int(out,nnf);
    formater(out);
    if ( surf ) {
      l = 0;
      for (k=1; k<=3*nf; k+=3) {
        p1 = &points[tri[k+0]];
        p2 = &points[tri[k+1]];
        p3 = &points[tri[k+2]];
        ++l;
        if ( p1->new < 0 || p2->new < 0 || p3->new < 0 )  continue;
        ecrire_int(out,p1->new);
        ecrire_int(out,p2->new);
        ecrire_int(out,p3->new);
        ecrire_int(out,ref[l]);
        formater(out);
      }
    }
    else {
      for (k=1; k<=3*nf; k+=3) {
        p1 = &points[tri[k+0]];
        p2 = &points[tri[k+1]];
        p3 = &points[tri[k+2]];
        if ( p1->new < 0 || p2->new < 0 || p3->new < 0 )  continue;
        ecrire_int(out,p1->new);
        ecrire_int(out,p2->new);
        ecrire_int(out,p3->new);
        ecrire_int(out,0);
        formater(out);
      }
    }
    formater(out);
  }

  /* specific entities */
  if ( !dosurf && mtype ) {
    strcpy(namein,namesurf);
    ptr = strstr(namein,".noboite");
    if ( ptr )  *ptr = '\0';
    ptr = strstr(namein,".boite");
    if ( ptr )  *ptr = '\0';

    in = ouvrir_mesh(namein,"r",&ver);
    if ( in ) {
      /* corners */
      pos = chercher_mot_clef(in,Corners,0);
      if ( pos ) {
        fseek(in,pos,SEEK_SET);
        nc = lire_int(in);
        if ( !quiet )  fprintf(stdout,"   corners   %6d  ",nc);
        ncc = 0;
        for (k=1; k<=nc; k++) {
          cc = lire_int(in);
          p1 = &points[cc];
          if ( p1->new > 0 ) ncc++;
        }
        if ( !quiet )  fprintf(stdout,"%d\n",ncc);
        ecrire_commentaire(out,"Set of corners");
        ecrire_mot_clef(out,Corners);
        ecrire_int(out,ncc);
        formater(out);
        
        fseek(in,pos,SEEK_SET);
        nc = lire_int(in);
        for (k=1; k<=nc; k++) {
          cc = lire_int(in);
          p1 = &points[cc];
          if ( p1->new > 0 ) {
            ecrire_int(out,p1->new);
            formater(out);
          }
        }
        formater(out);
      }
    
      /* required points */
      pos = chercher_mot_clef(in,RequiredVertices,pos);
      if ( pos ) {
        fseek(in,pos,SEEK_SET);
        nc = lire_int(in);
        if ( !quiet )  fprintf(stdout,"   required %6d ",nc);
        ncc = 0;
        for (k=1; k<=nc; k++) {
          cc = lire_int(in);
          p1 = &points[cc];
          if ( p1->new > 0 ) ncc++;
        }
        if ( !quiet )  fprintf(stdout,"%d ",ncc);
        fseek(in,pos,SEEK_SET);
        nc = lire_int(in);
        ecrire_commentaire(out,"Set of required vertices");
        ecrire_mot_clef(out,RequiredVertices);
        ecrire_int(out,ncc);
        formater(out);
        for (k=1; k<=nc; k++) {
          cc = lire_int(in);
          p1 = &points[cc];
          if ( p1->new > 0 ) {
            ecrire_int(out,p1->new);
            formater(out);
          }
        }
        formater(out);
      }

      /* edges */
      pos = chercher_mot_clef(in,Edges,0);
      if ( pos ) {
        fseek(in,pos,SEEK_SET);
        na = lire_int(in);
        naa = 0;
        for (k=1; k<=na; k++) {
          aa = lire_int(in);
          bb = lire_int(in);
          rr = lire_int(in);
          p1 = &points[aa];
          p2 = &points[bb];
          if ( p1->new > 0 && p2->new > 0 )  naa++;
        }
        if ( naa > 0 ) {
          if ( !quiet )  fprintf(stdout,"   edges     %6d",naa);
          edg = (int*)calloc(na+1,sizeof(int));
          assert(edg);

          fseek(in,pos,SEEK_SET);
          na = lire_int(in);
          ecrire_commentaire(out,"Set of edges");
          ecrire_mot_clef(out,Edges);
          ecrire_int(out,naa);
          formater(out);
          naa = 0;
          for (k=1; k<=na; k++) {
            aa = lire_int(in);
            bb = lire_int(in);
            rr = lire_int(in);
            p1 = &points[aa];
            p2 = &points[bb];
            if ( p1->new > 0 && p2->new > 0 ) {
              ecrire_int(out,p1->new);
              ecrire_int(out,p2->new);
              ecrire_int(out,rr);
              formater(out);
              edg[k] = ++naa;
            }
          }
          formater(out);

          /* ridges */
          pos = chercher_mot_clef(in,Ridges,pos);
          if ( pos ) {
            fseek(in,pos,SEEK_SET);
            na  = lire_int(in);
            naa = 0;
            for (k=1; k<=na; k++) {
              aa = lire_int(in);
              if ( edg[aa] > 0 )  naa++;
            }
            if ( naa > 0 ) {
              if ( !quiet )  fprintf(stdout,"   ridges  %6d / %6d",naa,na);
              fseek(in,pos,SEEK_SET);
              na = lire_int(in);
              ecrire_commentaire(out,"Set of ridges");
              ecrire_mot_clef(out,Ridges);
              ecrire_int(out,naa);
              formater(out);
              for (k=1; k<=na; k++) {
                aa  = lire_int(in);
                if ( edg[aa] > 0 ) {
                  ecrire_int(out,edg[aa]);
                  formater(out);
                }
              }
              formater(out);
            }
          }
          pos = chercher_mot_clef(in,RequiredEdges,pos);
          if ( pos ) {
            fseek(in,pos,SEEK_SET);
            na  = lire_int(in);
            naa = 0;
            for (k=1; k<=na; k++) {
              aa = lire_int(in);
              if ( edg[aa] > 0 )  naa++;
            }

            if ( naa > 0 ) {
              if ( !quiet )  fprintf(stdout,"   req.edges  %6d / %6d",naa,na);
              fseek(in,pos,SEEK_SET);
              na = lire_int(in);
              ecrire_commentaire(out,"Set of required edges");
              ecrire_mot_clef(out,RequiredEdges);
              ecrire_int(out,naa);
              formater(out);
              for (k=1; k<=na; k++) {
                aa  = lire_int(in);
                if ( edg[aa] > 0 ) {
                  ecrire_int(out,edg[aa]);
                  formater(out);
                }
              }
              formater(out);
            }
          }
          fprintf(stdout,"\n");
        }
      }
      
      pos = chercher_mot_clef(in,Normals,0);
      if ( pos ) {
        fseek(in,pos,SEEK_SET);
        nn = lire_int(in);
        if ( !quiet )  fprintf(stdout,"   normals %6d\n",nn);
	ecrire_commentaire(out,"Set of normals");
	ecrire_mot_clef(out,Normals);
	ecrire_int(out,nn);
	formater(out);
        for (k=1; k<=nn; k++) {
          nx = lire_reel(in);
	  ny = lire_reel(in);
	  nz = lire_reel(in);
	  ecrire_reel(out,nx);
	  ecrire_reel(out,ny);
	  ecrire_reel(out,nz);
	  formater(out);
        }
	formater(out);

	posnv = chercher_mot_clef(in,NormalAtVertices,pos);
        if ( posnv ) {
          fseek(in,posnv,SEEK_SET);
          nnv = lire_int(in);
          if ( !quiet )  fprintf(stdout,"   normal at vertices %6d\n",nnv);
	  ecrire_commentaire(out,"Normals at vertices");
	  ecrire_mot_clef(out,NormalAtVertices);
	  ecrire_int(out,nnv);
	  formater(out);
          for (k=1; k<=nnv; k++) {
            nv = lire_int(in);
	    ecrire_int(out,nv);
            nv = lire_int(in);
	    ecrire_int(out,nv);
	    formater(out);
	  }
	  formater(out);
        }
	
	posnv = chercher_mot_clef(in,NormalAtTriangleVertices,pos);
        if ( posnv ) {
	  fseek(in,posnv,SEEK_SET);
          ntv = lire_int(in);
          if ( !quiet )  fprintf(stdout,"   normal at triangle vertices %6d\n",ntv);
	  ecrire_commentaire(out,"Normals at triangle vertices");
	  ecrire_mot_clef(out,NormalAtTriangleVertices);
	  ecrire_int(out,ntv);
	  formater(out);
          for (k=1; k<=ntv; k++) {
            nv = lire_int(in);
	    ecrire_int(out,nv);
            nv = lire_int(in);
	    ecrire_int(out,nv);
            nv = lire_int(in);
	    ecrire_int(out,nv);
	    formater(out);
	  }
	  formater(out);
	}	
      }
    
    }
    fclose(in);
  }

  /* tets */
  if ( !quiet && np > 10000 )
    fprintf(stdout,"  Tetrahedra  %8d\n",ne);
  ecrire_commentaire(out,"Set of tetrahedra");
  ecrire_mot_clef(out,Tetrahedra);
  ecrire_int(out,ne);
  formater(out);

  for (k=1; k<=ne; k++) {
    pt = &tets[k];
    for (i=0; i<4; i++) {
      ppt = &points[pt->v[i]];
      ecrire_int(out,ppt->new);
    }
    /*reft = qualtet(k);
    ecrire_int(out,reft);*/
    ecrire_int(out,pt->ref);
    formater(out);
  }


  ecrire_mot_clef(out,End);
  fclose(out);

  return(1);
}


/* save .amdba3 file */
int saveAmdba(char *name) {
  FILE     *out;
  pTetra    pt;
  pPoint    ppt,p1,p2,p3;
  int       k,l,i,nnp,nnf;

  out = fopen(name,"w");
  if ( !out ) {
    fprintf(stderr,"  Unable to open %s\n",name);
    exit(1);
  }
  fprintf(stdout,"\n Writing %s\n",name);

  /* tassage */
  nnp = nnf = 0;
  for (k=1; k<=ne; k++) {
    pt = &tets[k];
    for (i=0; i<4; i++) {
      ppt = &points[pt->v[i]];
      ppt->new = 0;
    }
  }
  for(k=1; k<=np; k++) {
    ppt = &points[k];
    if ( !ppt->new )  ppt->new = ++nnp;
  }
  for (k=1; k<=3*nf; k+=3) {
    p1 = &points[tri[k]];
    p2 = &points[tri[k+1]];
    p3 = &points[tri[k+2]];
    if ( p1->new < 0 || p2->new < 0 || p3->new < 0 )  continue;
    ++nnf;
  }
  
  fprintf(out,"%d %d %d\n",nnp,ne,nnf);
  for (k=1; k<=np; k++) {
    ppt = &points[k];
    if ( ppt->new < 0 )  continue;
    fprintf(out,"%f %f %f\n",ppt->c[0],ppt->c[1],ppt->c[2]);
    if ( np > 10000 && (k % 1000 == 0) )
      fprintf(stdout,"%10\r",k);
  }

  for (k=1; k<=ne; k++) {
    pt = &tets[k];
    for (i=0; i<4; i++) {
      ppt = &points[pt->v[i]];
      fprintf(out,"%d ",ppt->new);
    }
    fprintf(out,"\n");
    if ( np > 10000 && (k % 1000 == 0) )
      fprintf(stdout,"%10\r",k);
  }

  if ( surf ) {
    l = 0;
    for (k=1; k<=3*nf; k+=3) {
      p1 = &points[tri[k]];
      p2 = &points[tri[k+1]];
      p3 = &points[tri[k+2]];
      ++l;
      if ( p1->new < 0 || p2->new < 0 || p3->new < 0 )  continue;
      fprintf(out,"%d\n",ref[l]);
    }
  } 
  else {
    for (k=1; k<=nnf; k++)
      fprintf(out,"0\n");
  }

  for (k=1; k<=3*nf; k+=3) {
    p1 = &points[tri[k]];
    p2 = &points[tri[k+1]];
    p3 = &points[tri[k+2]];
    if ( p1->new < 0 || p2->new < 0 || p3->new < 0 )  continue;
    fprintf(out,"%d %d %d\n",p1->new,p2->new,p3->new);
    if ( np > 10000 && (k % 1000 == 0) )
      fprintf(stdout,"%10\r",k/3);
  }

  fclose(out);
  return(1);
}


/* save .amdba3 file */
int saveAmFmt(char *name) {
  FILE     *out;
  pTetra    pt;
  pPoint    ppt,p1,p2,p3;
  int       k,i,nnp,nnf;

  out = fopen(name,"w");
  if ( !out ) {
    fprintf(stderr,"  Unable to open %s\n",name);
    exit(1);
  }
  fprintf(stdout,"\n Writing %s\n",name);

  /* tassage */
  nnp = nnf = 0;
  for (k=1; k<=ne; k++) {
    pt = &tets[k];
    for (i=0; i<4; i++) {
      ppt = &points[pt->v[i]];
      ppt->new = 0;
    }
  }
  for(k=1; k<=np; k++) {
    ppt = &points[k];
    if ( !ppt->new )  ppt->new = ++nnp;
  }
  for (k=1; k<=nf; k++) {
    refv[tri[k+0]] = 1;
    refv[tri[k+1]] = 1;
    refv[tri[k+2]] = 1;
    ++nnf;
  }
  
  /* header */
  fprintf(out,"%d %d\n",nnp,ne);
  for (k=1; k<=ne; k++) {
    pt = &tets[k];
    for (i=0; i<4; i++) {
      ppt = &points[pt->v[i]];
      fprintf(out,"%5d ",ppt->new);
    }
    if ( k % 3 == 0 )  fprintf(out,"\n");
  }
  fprintf(out,"\n");

  for (k=1; k<=np; k++) {
    ppt = &points[k];
    if ( ppt->new < 0 )  continue;
    fprintf(out,"%E %E %E\n",ppt->c[0],ppt->c[1],ppt->c[2]);
  }
  fprintf(out,"\n");

  /* refs elements */
  for (k=1; k<=ne; k++) {
    pt = &tets[k];
    fprintf(out,"%5d ",pt->ref);
    if ( k % 10 == 0 )  fprintf(out,"\n");
  }
  fprintf(out,"\n");

  if ( refv ) {
    for (k=1; k<=np; k++) {
      fprintf(out,"%5d ",refv[k]);
      if ( k % 10 == 0 )  fprintf(out,"\n");
    }
  }
  else {
    for (k=1; k<=np; k++) {
      fprintf(out," 1 ",refv[k]);
      if ( k % 10 == 0 )  fprintf(out,"\n");
    }
  }
  fprintf(out,"\n");

  /* refs points */
  
  fclose(out);
  return(1);
}


int addRefTetra() {
  pTetra  pt;
  int    *pile,base,k;
/*
  base = 0;
  dep  = 1;
  pile = (int*)malloc((ne+1)*sizeof(int));
  assert(pile);
 
  do {
    for (k=dep; k<=ne; k++) {
      pt = &tets[k];
      if ( !pt->ref ) break;
    }
    if ( k > ne ) break;

    dep = k + 1;
    base++;
    ipil       = 1;
    pile[ipil] = k;
    while ( ipil ) {
      pt = &tets[ pile[ipil] ];
      pt->ref = base;
      for (i=0; i<4; i++) {
        if ( pt->adj[i] )  continue;
        
      }
    }
    do {
      
    }
    while ();
  }
  while ( dep <= ne );

  free(pile);
  return(1);
*/
}



int main(int argc,char *argv[]) {
  long     i;
  int      ret;
  char     src[256],dest[256],dangfa;
  clock_t  ct0,ct1,ct2,ct3;


  /* defaults */
  ct0 = clock();
  quiet  = 1;
  surf   = 0;
  dosurf = 1;
  dangfa = 0;
  tonb   = 0;
  addref = 0;
  src[0] = dest[0] = '\0';
  ref  = 0;
  refv = 0;
  npinit = 0;

  /* parse args */
  if ( argc < 3 ) {
    fprintf(stdout,"usage: nb2mesh filein fileout [-v] [-s]\n");
    fprintf(stdout,"   filein :  xxx.[no]boite[b]\n");
    fprintf(stdout,"   fileout:  yyy.mesh[b]  or zzz.amdba3\n");
    fprintf(stdout,"   -v     :  verbose mode\n");
    fprintf(stdout,"   -s     :  do not create boundary (surface)\n");
    fprintf(stdout,"   -x     :  split constrained tetras\n");
    fprintf(stdout,"   -df    :  remove dangling faces\n");
    fprintf(stdout,"   -ref   :  add reference to sub-domains\n");
    exit(1);
  }
  else {
    i = 1;
    while ( i < argc ) {
      if ( !strcmp(argv[i],"-v") )
	quiet = 0;
      else if ( !strcmp(argv[i],"-s") )
	dosurf = 0;
      else if ( !strcmp(argv[i],"-x") )
        split4 = 1;
      else if ( !strcmp(argv[i],"-df") )
        dangfa = 1;
      else if ( !strcmp(argv[i],"-ref") )
        addref = 1;
      else if ( !src[0] )
	strcpy(src,argv[i]);
      else if ( !dest[0] ) {
	strcpy(dest,argv[i]);
	if ( !strstr(dest,".mesh") && !strstr(dest,".amdba3") && 
             !strstr(dest,".noboite") && !strstr(dest,".am_fmt") )
	  strcat(dest,".meshb");
        if ( strstr(dest,".noboite") )  tonb = 1;
      }
      i++;
    }
  }

  /* load 3D file */
  ct1 = clock();
  if ( strstr(src,".mesh") )
    ret = loadMesh(src);
  else 
    ret = loadNoboite(src);
  if ( !ret )  exit(1);

  /* load surface */
  if ( !nf || dosurf ) {
    nf  = 0;
    dosurf = 1 - loadSurf(src);
  }
  if ( split4 )  splitTetra();
  ct1 = difftime(clock(),ct1);
  fprintf(stdout,"  Input seconds:     %.2f\n",
          (double)ct1/(double)CLOCKS_PER_SEC);

  /* hash faces */
  if ( dosurf ) {
    ct2 = clock();
    hashFaces();
    ct2 = difftime(clock(),ct2);
    fprintf(stdout,"  Hash seconds:      %.2f\n",
            (double)ct2/(double)CLOCKS_PER_SEC);
  }
  if ( dangfa ) {
    ct2 = clock();
    delDangling();
    ct2 = difftime(clock(),ct2);
    fprintf(stdout,"  Dangling seconds:  %.2f\n",
            (double)ct2/(double)CLOCKS_PER_SEC);
  }
/*
  if ( addref ) {
    hashTrias();
    addRefTetra();    
  }
*/
  /* save file */
  ct3 = clock();
  if ( strstr(dest,".mesh") )
    ret = saveMesh(dest);
  else if ( strstr(dest,".amdba3") )
    ret = saveAmdba(dest);
  else if ( strstr(dest,".am_fmt") )
    ret = saveAmFmt(dest);
  else if ( strstr(dest,".noboite") || strstr(dest,".boite") )
    ret = saveNoboite(dest);
  if ( !ret )  exit(1);
  if ( hash )  free(hash);
  free(points);
  free(tets);
  ct3 = difftime(clock(),ct3);
  fprintf(stdout,"  Output seconds:    %.2f\n",
          (double)ct3/(double)CLOCKS_PER_SEC);

  if ( !quiet ) {
    fprintf(stdout,"\n Statistics\n");
    if ( nf > 0 ) {
      fprintf(stdout,"  Mesh vertices    %8d\n",np);
      fprintf(stdout,"  Mesh triangles   %8d\n",nf);
      fprintf(stdout,"  Mesh tetrahedra  %8d\n",ne);
    }
    else {
      fprintf(stdout,"  Mesh vertices    %8d\n",np);
      fprintf(stdout,"  Mesh tetrahedra  %8d\n",ne);
    }
  }

  ct0 = difftime(clock(),ct0);
  fprintf(stdout,"\n Total running seconds:  %.2f\n",
          (double)ct0/(double)CLOCKS_PER_SEC);

  return(0);
}
