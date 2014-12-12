/*********************************************************************
*	25 - 02 - 97
*	gestion des fichiers image ELF
*
*********************************************************************/

#ifndef CGAL_FILE_XT_H
#define CGAL_FILE_XT_H

#include <stdio.h>
#include <string.h>

#include <CGAL/config.h>

namespace CGAL {
namespace Total {

void permuteLong(char *a)
{
    char tmp;

#ifdef CGAL_LITTLE_ENDIAN 
tmp=a[0]; a[0]=a[3]; a[3]=tmp;
tmp=a[1]; a[1]=a[2]; a[2]=tmp;
#endif
}

void permuteLongTab(long *a,int nb)
{
    int i;

#ifdef CGAL_LITTLE_ENDIAN 
for(i=0;i<nb;i++) permuteLong( (char *) &a[i]);
#endif
}

void permuteShort(char *a)
{
    char tmp;

#ifdef CGAL_LITTLE_ENDIAN 
tmp=a[0]; a[0]=a[1]; a[1]=tmp;
#endif
}

void permuteShortTab(short *a, int nb)
{
    int i;

#ifdef CGAL_LITTLE_ENDIAN 
for(i=0;i<nb;i++) permuteShort( (char *) &a[i]);
#endif
}

void permuteFloat(char *a )
{
	char tmp;
	
#ifdef CGAL_LITTLE_ENDIAN 
tmp=a[0]; a[0]=a[3]; a[3]=tmp;
tmp=a[1]; a[1]=a[2]; a[2]=tmp;
#endif
}


void permuteFloatTab(float *a, int nb)
{
	int i;

#ifdef CGAL_LITTLE_ENDIAN 
for(i=0;i<nb;i++) permuteFloat((char *) &a[i]);
#endif
}


void permuteDouble(char *a )
{
	char tmp;
	
#ifdef CGAL_LITTLE_ENDIAN 
tmp=a[0]; a[0]=a[3]; a[3]=tmp;
tmp=a[1]; a[1]=a[2]; a[2]=tmp;
#endif
}


void permuteDoubleTab(double *a, int nb)
{
	int i;

#ifdef CGAL_LITTLE_ENDIAN 
for(i=0;i<nb;i++) permuteDouble((char *) &a[i]);
#endif
}

/*******************************************************************/
int lire_longueur_trace(const char * nomfich, long * long_trace)
{
    FILE *fin;
    int flag=0;

if ( (fin=fopen(nomfich,"rw"))!=NULL )
 {  
    fseek(fin,4,0);
    if (fread(long_trace,4,1,fin) != 1) {
      flag = -1;
    }
    permuteLong((char *)long_trace);
    fclose(fin);
 }
else flag=-1;
return(flag);
}

/********************************************************************/
int lire_nb_trace(const char * nomfich, long * nb_trace)
{
    FILE *fin;
    int flag=0;

if ( (fin=fopen(nomfich,"rw"))!=NULL )
 {  
    fseek(fin,8,0);
    if (fread(nb_trace,4,1,fin) != 1) {
      flag = -1;
    }
    permuteLong((char *)nb_trace);
    fclose(fin);
 }
else flag=-1;
return(flag);
}

/********************************************************************/
int lire_nb_plan(const char * nomfich, long * nb_plan)
{
    FILE *fin;
    int flag=0;

if ( (fin=fopen(nomfich,"rw"))!=NULL )
 {  
    fseek(fin,12,0);
    if (fread(nb_plan,4,1,fin) != 1) {
      flag = -1;
    }
    permuteLong((char *)nb_plan);
    fclose(fin);
 }
else flag=-1;
return(flag);
}

/********************************************************************/

int lire_nb_octet(const char * nomfich, long * nb_octet)
{
    FILE *fin;
    int flag=0;

if ( (fin=fopen(nomfich,"rw"))!=NULL )
 {  
    fseek(fin,36,0);
    if (fread(nb_octet,4,1,fin) != 1) {
      flag = -1;
    }
    permuteLong((char *)nb_octet);
    fclose(fin);
 }
else flag=-1;
return(flag);
}

/********************************************************************/

int lire_longueur_entete(const char * nomfich, long * long_entete)
{
    FILE *fin;
    int flag=0;

if ( (fin=fopen(nomfich,"rw"))!=NULL )
 {  
    fseek(fin,72,0);
    if (fread(long_entete,4,1,fin) != 1) {
      flag = -1;
    }
    permuteLong((char *)long_entete);
    fclose(fin);
 }
else flag=-1;
return(flag);
}

} // end namespace Total
} // end namespace CGAL


#endif // CGAL_FILE_XT_H
