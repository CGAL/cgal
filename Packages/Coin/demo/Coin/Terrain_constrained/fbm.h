/*****************************************************************

  Prototypes for the fractional Brownian motion algorithm. These
  functions were originally the work of F. Kenton Musgrave.  For
  documentation of the different functions please refer to the book:
  "Texturing and modeling: a procedural approach"
  by David S. Ebert et. al.
  
******************************************************************/

#ifndef _fbm_h
#define _fbm_h

#include  <math.h>

#ifdef __cplusplus
extern "C" {
#endif

#define TRUE    1
#define FALSE   0

typedef struct {
    double x;
    double y;
    double z;
} Vector;
    
float noise3(float vec[]);
double fBm( Vector point, double H, double lacunarity, double octaves, 
	    int init );
#endif

#ifdef __cplusplus
}
#endif




