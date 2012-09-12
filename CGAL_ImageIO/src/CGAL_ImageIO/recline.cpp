// Copyright (c) 2005-2008 ASCLEPIOS Project, INRIA Sophia-Antipolis (France)
// All rights reserved.
//
// This file is part of the ImageIO Library, and as been adapted for
// CGAL (www.cgal.org).
// You can redistribute it and/or  modify it under the terms of the
// GNU Lesser General Public License as published by the Free Software Foundation;
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// These files are provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     :  ASCLEPIOS Project (INRIA Sophia-Antipolis), Laurent Rineau

/*************************************************************************
 * recline.c - tools for recursive filtering of 1D lines
 *
 * $Id$
 *
 * Copyright©INRIA 1999
 *
 * DESCRIPTION: 
 *
 * Recursive filtering of a line (a 1D array)
 * Filter coefficient are static variables.
 *
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * June, 9 1998
 *
 * Copyright Gregoire Malandain, INRIA
 *
 * ADDITIONS, CHANGES
 *
 */



#include <stdio.h>
#include <stdlib.h>
#include <math.h>



#include "recline.h"

static int _VERBOSE_RECLINE_ = 0;

#define EXIT_ON_FAILURE 0
#define EXIT_ON_SUCCESS 1



void printRecursiveCoefficients( RFcoefficientType *RFC )
{
  printf( "denominator:\n" );
  printf( "%f %f %f %f\n", RFC->sd1, RFC->sd2, RFC->sd3, RFC->sd4 );
  printf( "positive numerator:\n" );
  printf( "%f %f %f %f\n", RFC->sp0, RFC->sp1, RFC->sp2, RFC->sp3 );
  printf( "negative numerator:\n" );
  printf( "%f %f %f %f %f\n", RFC->sn0, RFC->sn1, RFC->sn2, RFC->sn3, RFC->sn4 );
  printf( "\n" );
}


RFcoefficientType * InitRecursiveCoefficients( double x, 
					       recursiveFilterType type_filter, 
					       derivativeOrder derivative )
{
  const char *proc="InitRecursiveCoefficients";
  double ex, k1, k2;
  double a0, a1, c0, c1, omega0, omega1, b0, b1;
  double cos0, sin0, cos1, sin1;
  double sumA=0.0, sumC=0.0, aux;

  RFcoefficientType *RFC = NULL;
  RFC = (RFcoefficientType *)malloc( sizeof(RFcoefficientType) );
  if ( RFC == NULL ) {
    if ( _VERBOSE_RECLINE_ != 0 ) 
      fprintf( stderr, "%s: allocation failed\n", proc );
    return( NULL );
  }
  
  RFC->sd1 = RFC->sd2 = RFC->sd3 = RFC->sd4 = 0.0;
  RFC->sp0 = RFC->sp1 = RFC->sp2 = RFC->sp3 = 0.0;
  RFC->sn0 = RFC->sn1 = RFC->sn2 = RFC->sn3 = RFC->sn4 = 0.0;
  
  RFC->type_filter = UNKNOWN_FILTER;
  RFC->derivative  = NODERIVATIVE;
  
  ex = k1 = k2 = 0.0;
  a0 = a1 = c0 = c1 = 0.0;
  b0 = b1 = omega0 = omega1 = 0.0;
  
  /*--- Selon le type de filtrage (filtres de Deriche,
    ou approximation de la gaussienne), x designe
    soit alpha, soit sigma                         ---*/
  
  switch ( type_filter ) {

  case GAUSSIAN_FIDRICH :
    
    if ( x < 0.1 ) {
      if ( _VERBOSE_RECLINE_ != 0 ) {
	fprintf( stderr, "%s: improper value of coefficient (should be >= 0.1).\n", proc );
      }
      free( RFC );
      return( NULL );
    }

    switch ( derivative ) {
    default :
      if ( _VERBOSE_RECLINE_ != 0 ) {
	fprintf( stderr, "%s: improper value of derivative order.\n", proc );
      }
      free( RFC );
      return( NULL );
    case DERIVATIVE_0 :
      a0 =  0.6570033214 / x;
      a1 =  1.978946687  / x;
      c0 = -0.2580640608 / x;
      c1 = -0.2391206463 / x;
      omega0 = 0.6512453378;
      omega1 = 2.05339943;
      b0 = 1.906154352;
      b1 = 1.881305409;
      break;
    case DERIVATIVE_1 :
    case DERIVATIVE_1_CONTOURS :
      a0 = -0.1726729496 / x;
      a1 = -2.003565572  / x;
      c0 =  0.1726730777 / x;
      c1 =  0.4440126835 / x;
      b0 = 1.560644213;
      b1 = 1.594202256;
      omega0 = 0.6995461735;
      omega1 = 2.144671764;
      break;
    case DERIVATIVE_2 :
      a0 = -0.7241334169 / x;
      a1 =  1.688628765  / x;
      c0 =  0.3251949838 / x;
      c1 = -0.7211796018 / x;
      b0 = 1.294951143;
      b1 = 1.427007123;
      omega0 = 0.7789803775;
      omega1 = 2.233566862;
      break;
    case DERIVATIVE_3 :
      a0 =  1.285774106  / x;
      a1 = -0.2896378408 / x;
      c0 = -1.28577129   / x;
      c1 =  0.26249833   / x;
      b0 = 1.01162886;
      b1 = 1.273344739;
      omega0 = 0.9474270928;
      omega1 = 2.337607006;
      break;
    }
    
    omega0 /= x;   sin0 = sin( omega0 );   cos0 = cos( omega0 ); 
    omega1 /= x;   sin1 = sin( omega1 );   cos1 = cos( omega1 ); 
    b0 /= x;
    b1 /= x;

    RFC->sp0  = a0 + c0;
    RFC->sp1  = exp( -b1 ) * (c1 * sin1 - (c0 + 2 * a0) * cos1);
    RFC->sp1 += exp( -b0 ) * (a1 * sin0 - (2 * c0 + a0) * cos0);
    RFC->sp2  = 2.0 * exp( -b0 - b1 ) 
      * ((a0 + c0) * cos1 * cos0 - cos1 * a1 * sin0 - cos0 * c1 * sin1);
    RFC->sp2 += c0 * exp( -2.0 * b0 ) + a0 * exp( -2.0 * b1 );
    RFC->sp3  = exp( -b1 - 2.0 * b0 ) * (c1 * sin1 - c0 * cos1);
    RFC->sp3 += exp( -b0 - 2.0 * b1 ) * (a1 * sin0 - a0 * cos0);
    
    RFC->sd1  = -2.0 * exp( -b1 ) * cos1 - 2.0 * exp( -b0 ) * cos0;
    RFC->sd2  = 4.0 * cos1 * cos0 * exp( -b0 - b1 ) 
      + exp( -2.0 * b1 ) + exp( -2.0 * b0 );
    RFC->sd3 = -2.0 * cos0 * exp( -b0 - 2.0 * b1 ) 
      - 2.0 * cos1 * exp( -b1 - 2.0 * b0 );
    RFC->sd4 = exp( -2.0 * b0 - 2.0 * b1 );
    
    switch ( derivative ) {
    default :
      if ( _VERBOSE_RECLINE_ != 0 ) {
	fprintf( stderr, "%s: improper value of derivative order.\n", proc );
      }
      free( RFC );
      return( NULL );
    case DERIVATIVE_0 :
    case DERIVATIVE_2 :
      RFC->sn1 =   RFC->sp1 - RFC->sd1 * RFC->sp0;
      RFC->sn2 =   RFC->sp2 - RFC->sd2 * RFC->sp0;
      RFC->sn3 =   RFC->sp3 - RFC->sd3 * RFC->sp0;
      RFC->sn4 = - RFC->sd4 * RFC->sp0;
      break;
    case DERIVATIVE_1 :
    case DERIVATIVE_1_CONTOURS :
    case DERIVATIVE_3 :
      RFC->sn1 = - RFC->sp1 + RFC->sd1 * RFC->sp0;
      RFC->sn2 = - RFC->sp2 + RFC->sd2 * RFC->sp0;
      RFC->sn3 = - RFC->sp3 + RFC->sd3 * RFC->sp0;
      RFC->sn4 =   RFC->sd4 * RFC->sp0;
    }
    
    RFC->type_filter = type_filter;
    RFC->derivative  = derivative;
    break;
    
  case GAUSSIAN_DERICHE :
    
    if ( x < 0.1 ) {
      if ( _VERBOSE_RECLINE_ != 0 ) {
	fprintf( stderr, "%s: improper value of coefficient (should be >= 0.1).\n", proc );
      }
      free( RFC );
      return( NULL );
    }

    switch ( derivative ) {
    default :
      if ( _VERBOSE_RECLINE_ != 0 ) {
	fprintf( stderr, "%s: switch to default coefficients (smoothing).\n", proc );
      }
      derivative = DERIVATIVE_0;
    case DERIVATIVE_0 :
      a0     =  1.68;
      omega0 =  0.6318;
      a1     =  3.735;
      b0     =  1.783;
      c0     = -0.6803;
      omega1 =  1.997;
      c1     = -0.2598;
      b1     =  1.723;
      break;
    case DERIVATIVE_1 :
    case DERIVATIVE_1_CONTOURS :
      a0     =  -0.6472;
      omega0 =  0.6719;
      a1     =  -4.531;
      b0     =  1.527;
      c0     =  0.6494;
      omega1 =  2.072;
      c1     =  0.9557;
      b1     =  1.516;
      break;
    case DERIVATIVE_2 :
      a0     = -1.331;
      omega0 =  0.748;
      a1     =  3.661;
      b0     =  1.24;
      c0     =  0.3225;
      omega1 =  2.166;
      c1     = -1.738;
      b1     =  1.314;
    }
	 
    omega0 /= x;   sin0 = sin( omega0 );   cos0 = cos( omega0 ); 
    omega1 /= x;   sin1 = sin( omega1 );   cos1 = cos( omega1 ); 
    b0 /= x;
    b1 /= x;

    /*--- normalisation ---*/
    switch ( derivative ) {
    default :
    case DERIVATIVE_0 :
      sumA  = 2.0 * a1 * exp( b0 ) * cos0 * cos0 - a0 * sin0 * exp( 2.0 * b0 );
      sumA += a0 * sin0 - 2.0 * a1 * exp( b0 );
      sumA /= ( 2.0 * cos0 * exp( b0 ) - exp( 2.0 * b0 ) - 1 ) * sin0;
      sumC  = 2.0 * c1 * exp( b1 ) * cos1 * cos1 - c0 * sin1 * exp( 2.0 * b1 );
      sumC += c0 * sin1 - 2.0 * c1 * exp( b1 );
      sumC /= ( 2.0 * cos1 * exp( b1 ) - exp( 2.0 * b1 ) - 1 ) * sin1;
      break;
    case DERIVATIVE_1 :
      aux   = exp( 4.0 * b0 ) - 4.0 * cos0 * exp( 3.0 * b0 );
      aux  += 2.0 * exp( 2.0 * b0 ) + 4.0 * cos0 * cos0 * exp( 2.0 * b0 );
      aux  += 1.0 - 4.0 * cos0 * exp( b0 );
      sumA  = a0 * cos0 - a1 * sin0 + a1 * sin0 * exp( 2.0 * b0 );
      sumA += a0 * cos0 * exp( 2.0 * b0 ) - 2.0 * a0 * exp( b0 );
      sumA *= exp( b0 ) / aux;
      aux   = exp( 4.0 * b1 ) - 4.0 * cos1 * exp( 3.0 * b1 );
      aux  += 2.0 * exp( 2.0 * b1 ) + 4.0 * cos1 * cos1 * exp( 2.0 * b1 );
      aux  += 1.0 - 4.0 * cos1 * exp( b1 );
      sumC  = c0 * cos1 - c1 * sin1 + c1 * sin1 * exp( 2.0 * b1 );
      sumC += c0 * cos1 * exp( 2.0 * b1 ) - 2.0 * c0 * exp( b1 );
      sumC *= exp( b1 ) / aux;
      /*--- on multiplie les sommes par 2 car on n'a calcule que des demi-sommes 
	et on change le signe car la somme doit etre egale a -1              ---*/
      sumA *= (-2.0);
      sumC *= (-2.0);
      break;
    case DERIVATIVE_1_CONTOURS :
      /*--- la somme de 1 a l'infini est egale a 1 : cela introduit
	un petit biais (reponse un rien superieur a la hauteur du step).
	Avec une somme de 0 a l'infini, c'est pire                       ---*/
      sumA  = a1 * exp( b0 ) - a1 * cos0 * cos0 * exp( b0 );
      sumA += a0 * cos0 * sin0 * exp( b0 ) - a0 * sin0;
      sumA /= sin0 * ( 2.0 * cos0 * exp( b0 ) - exp( 2.0 * b0 ) - 1 );
      sumC  = c1 * exp( b1 ) - c1 * cos1 * cos1 * exp( b1 );
      sumC += c0 * cos1 * sin1 * exp( b1 ) - c0 * sin1;
      sumC /= sin1 * ( 2.0 * cos1 * exp( b1 ) - exp( 2.0 * b1 ) - 1 );
      break;
    case DERIVATIVE_2 :
      aux   = 12.0 * cos0 * exp( 3.0 * b0 ) - 3.0 * exp( 2.0 * b0 );
      aux  += 8.0 * cos0 * cos0 * cos0 * exp( 3.0 * b0 ) 
	- 12.0 * cos0 * cos0 * exp( 4.0 * b0 );
      aux  -= 3.0 * exp( 4.0 * b0 );
      aux  += 6.0 * cos0 * exp( 5.0 * b0 ) -  exp( 6.0 * b0 ) 
	+ 6.0 * cos0 * exp( b0 );
      aux  -= ( 1.0 + 12.0 * cos0 * cos0 * exp( 2.0 * b0 ) );
      sumA  = 4.0 * a0 * sin0 * exp( 3.0 * b0 ) 
	+ a1 * cos0 * cos0 * exp( 4.0 * b0 );
      sumA -= ( 4.0 * a0 * sin0 * exp( b0 ) 
		+ 6.0 * a1 * cos0 * cos0 * exp( 2.0 * b0 ) );
      sumA += 2.0 * a1 * cos0 * cos0 * cos0 * exp( b0 ) 
	- 2.0 * a1 * cos0 * exp( b0 );
      sumA += 2.0 * a1 * cos0 * cos0 * cos0 * exp( 3.0 * b0 ) 
	- 2.0 * a1 * cos0 * exp( 3.0 * b0 );
      sumA += a1 * cos0 * cos0 - a1 * exp( 4.0 * b0 );
      sumA += 2.0 * a0 * sin0 * cos0 * cos0 * exp( b0 ) 
	- 2.0 * a0 * sin0 * cos0 * cos0 * exp( 3.0 * b0 );
      sumA -= ( a0 * sin0 * cos0 * exp( 4.0 * b0 ) + a1 );
      sumA += 6.0 * a1 * exp( 2.0 * b0 ) + a0 * cos0 * sin0;
      sumA *= 2.0 * exp( b0 ) / ( aux * sin0 );
      aux   = 12.0 * cos1 * exp( 3.0 * b1 ) - 3.0 * exp( 2.0 * b1 );
      aux  += 8.0 * cos1 * cos1 * cos1 * exp( 3.0 * b1 ) 
	- 12.0 * cos1 * cos1 * exp( 4.0 * b1 );
      aux  -= 3.0 * exp( 4.0 * b1 );
      aux  += 6.0 * cos1 * exp( 5.0 * b1 ) -  exp( 6.0 * b1 ) 
	+ 6.0 * cos1 * exp( b1 );
      aux  -= ( 1.0 + 12.0 * cos1 * cos1 * exp( 2.0 * b1 ) );
      sumC  = 4.0 * c0 * sin1 * exp( 3.0 * b1 ) 
	+ c1 * cos1 * cos1 * exp( 4.0 * b1 );
      sumC -= ( 4.0 * c0 * sin1 * exp( b1 ) 
		+ 6.0 * c1 * cos1 * cos1 * exp( 2.0 * b1 ) );
      sumC += 2.0 * c1 * cos1 * cos1 * cos1 * exp( b1 ) 
	- 2.0 * c1 * cos1 * exp( b1 );
      sumC += 2.0 * c1 * cos1 * cos1 * cos1 * exp( 3.0 * b1 ) 
	- 2.0 * c1 * cos1 * exp( 3.0 * b1 );
      sumC += c1 * cos1 * cos1 - c1 * exp( 4.0 * b1 );
      sumC += 2.0 * c0 * sin1 * cos1 * cos1 * exp( b1 ) 
	- 2.0 * c0 * sin1 * cos1 * cos1 * exp( 3.0 * b1 );
      sumC -= ( c0 * sin1 * cos1 * exp( 4.0 * b1 ) + c1 );
      sumC += 6.0 * c1 * exp( 2.0 * b1 ) + c0 * cos1 * sin1;
      sumC *= 2.0 * exp( b1 ) / ( aux * sin1 );
      /*--- on divise les sommes par 2 (la somme doit etre egale a 2) ---*/
      sumA /= 2;
      sumC /= 2;
    }
    a0 /= ( sumA + sumC );
    a1 /= ( sumA + sumC );
    c0 /= ( sumA + sumC );
    c1 /= ( sumA + sumC );
    
    /*--- coefficients du calcul recursif ---*/
    RFC->sp0  = a0 + c0;
    RFC->sp1  = exp( -b1 ) * (c1 * sin1 - (c0 + 2 * a0) * cos1);
    RFC->sp1 += exp( -b0 ) * (a1 * sin0 - (2 * c0 + a0) * cos0);
    RFC->sp2  = 2.0 * exp( -b0 - b1 ) 
      * ((a0 + c0) * cos1 * cos0 - cos1 * a1 * sin0 - cos0 * c1 * sin1);
    RFC->sp2 += c0 * exp( -2.0 * b0 ) + a0 * exp( -2.0 * b1 );
    RFC->sp3  = exp( -b1 - 2.0 * b0 ) * (c1 * sin1 - c0 * cos1);
    RFC->sp3 += exp( -b0 - 2.0 * b1 ) * (a1 * sin0 - a0 * cos0);
    
    RFC->sd1  = -2.0 * exp( -b1 ) * cos1 - 2.0 * exp( -b0 ) * cos0;
    RFC->sd2  = 4.0 * cos1 * cos0 * exp( -b0 - b1 ) 
      + exp( -2.0 * b1 ) + exp( -2.0 * b0 );
    RFC->sd3 = -2.0 * cos0 * exp( -b0 - 2.0 * b1 ) 
      - 2.0 * cos1 * exp( -b1 - 2.0 * b0 );
    RFC->sd4 = exp( -2.0 * b0 - 2.0 * b1 );
    
    switch ( derivative ) {
    default :
    case DERIVATIVE_0 :
    case DERIVATIVE_2 :
      RFC->sn1 =   RFC->sp1 - RFC->sd1 * RFC->sp0;
      RFC->sn2 =   RFC->sp2 - RFC->sd2 * RFC->sp0;
      RFC->sn3 =   RFC->sp3 - RFC->sd3 * RFC->sp0;
      RFC->sn4 = - RFC->sd4 * RFC->sp0;
      break;
    case DERIVATIVE_1 :
    case DERIVATIVE_1_CONTOURS :
    case DERIVATIVE_3 :
      RFC->sn1 = - RFC->sp1 + RFC->sd1 * RFC->sp0;
      RFC->sn2 = - RFC->sp2 + RFC->sd2 * RFC->sp0;
      RFC->sn3 = - RFC->sp3 + RFC->sd3 * RFC->sp0;
      RFC->sn4 =   RFC->sd4 * RFC->sp0;
    }
    
    RFC->type_filter = type_filter;
    RFC->derivative  = derivative;
    break;



  default :
    if ( _VERBOSE_RECLINE_ != 0 ) {
      fprintf( stderr, "%s: switch to default recursive filter (Deriche's filters).\n", proc );
    }
    type_filter = ALPHA_DERICHE;
  case ALPHA_DERICHE :

    if ( (x < 0.1) || (x > 1.9) ) {
      if ( _VERBOSE_RECLINE_ != 0 ) {
	fprintf( stderr, "%s: improper value of coefficient (should be >= 0.1 and <= 1.9).\n", proc );
      }
      free( RFC );
      return( NULL );
    }
    ex = exp( (-x) );
    
    switch ( derivative ) {
    default :
      if ( _VERBOSE_RECLINE_ != 0 ) {
	fprintf( stderr, "%s: switch to default coefficients (smoothing).\n", proc );
      }
      derivative = DERIVATIVE_0;
    case DERIVATIVE_0 :
      RFC->sp0 = (1.0 - ex) * (1.0 - ex) / (1.0 + 2.0 * x * ex - ex * ex);
      RFC->sp1 = RFC->sp0 * (x - 1.0) * ex;
      RFC->sn1 = RFC->sp0 * (x + 1.0) * ex;
      RFC->sn2 = (- RFC->sp0) * ex * ex;
      RFC->sd1 = (- 2.0) * ex;
      RFC->sd2 = ex * ex;
      break;
    case DERIVATIVE_1 :
      RFC->sp1 = - (1.0 - ex) * (1.0 - ex) * (1.0 - ex) / (2.0 * (1.0 + ex));
      RFC->sn1 = (- RFC->sp1);
      RFC->sd1 = (- 2.0) * ex;
      RFC->sd2 = ex * ex;	    
      break;
    case DERIVATIVE_1_CONTOURS :
      RFC->sp1 = - (1.0 - ex) * (1.0 - ex);
      RFC->sn1 = (- RFC->sp1);
      RFC->sd1 = (- 2.0) * ex;
      RFC->sd2 = ex * ex;	    
      break;
    case DERIVATIVE_2 :
      k1 = (- 2.0) * (1.0 - ex) * (1.0 - ex) * (1.0 - ex);
      k1 /= (1.0 + ex) * (1.0 + ex) * (1.0 + ex);
      k2 = (1.0 - ex * ex) / (2.0 * ex);
      RFC->sp0 = k1;
      RFC->sp1 = (- k1) * (1.0 + k2) * ex;
      RFC->sn1 = k1 * (1.0 - k2) * ex;
      RFC->sn2 = (- k1) * ex * ex;
      RFC->sd1 = (- 2.0) * ex;
      RFC->sd2 = ex * ex;
      break;
    case DERIVATIVE_3 :
      k1 = (1.0 + x) * ex + (x - 1.0);
      k2 = (1.0 - ex) / k1;
      k1 *= (1.0 - ex) * (1.0 - ex) * (1.0 - ex) * (1.0 - ex);
      k1 /= 2.0 * x * x * ex * ex;
      k1 /= ex + 1.0;
      RFC->sp0 = k1 * x * (k2 + 1.0);
      RFC->sp1 = (- k1) * x * (1.0 + k2 + k2*x) * ex;
      RFC->sn0 = (- RFC->sp0);
      RFC->sn1 = (- RFC->sp1);
      RFC->sd1 = (- 2.0) * ex;
      RFC->sd2 = ex * ex;
    }
    RFC->type_filter = type_filter;
    RFC->derivative  = derivative;
  }

  return( RFC );
}



int RecursiveFilter1D( RFcoefficientType *RFC,
		       double *in, 
		       double *out, 
		       double *work1, 
		       double *work2, 
		       int dim )
{
  const char *proc="RecursiveFilter1D";
  register double rp0, rp1, rp2, rp3;
  register double rd1, rd2, rd3, rd4;
  register double rn0, rn1, rn2, rn3, rn4;
  register int i;
  register double *w0, *w1, *w2, *w3, *w4;
  register double *d0, *d1, *d2, *d3, *d4;

  if ( RFC->type_filter == UNKNOWN_FILTER ) {
    if ( _VERBOSE_RECLINE_ != 0 )
      fprintf( stderr, "%s: unknown type of recursive filter.\n", proc );
    return( EXIT_ON_FAILURE );
  }
  if ( RFC->derivative == NODERIVATIVE ) {
    if ( _VERBOSE_RECLINE_ != 0 )
      fprintf( stderr, "%s: unknown type of derivative.\n", proc );
    return( EXIT_ON_FAILURE );
  }

  rd1 = rd2 = rd3 = rd4 = 0.0;
  rp0 = rp1 = rp2 = rp3 = 0.0;
  rn0 = rn1 = rn2 = rn3 = rn4 = 0.0;
  
  switch( RFC->type_filter ) {
  default :
    if ( _VERBOSE_RECLINE_ != 0 )
      fprintf( stderr, "%s: unknown type of recursive filter.\n", proc );
    return( EXIT_ON_FAILURE );
  case GAUSSIAN_FIDRICH :
  case GAUSSIAN_DERICHE :
    /*--- filtrage generique d'ordre 4 ---*/
    rp0 = RFC->sp0;   rp1 = RFC->sp1;   rp2 = RFC->sp2;   rp3 = RFC->sp3;
    rd1 = RFC->sd1;   rd2 = RFC->sd2;   rd3 = RFC->sd3;   rd4 = RFC->sd4;
    rn1 = RFC->sn1;   rn2 = RFC->sn2;   rn3 = RFC->sn3;   rn4 = RFC->sn4;
    
    /* on positionne les pointeurs 
     */
    w4 = work1;   w3 = w4+1;   w2 = w3+1;   w1 = w2+1;   w0 = w1+1;
    d3 = in+1;    d2 = d3+1;   d1 = d2+1;   d0 = d1+1;
    /*--- calcul de y+ ---*/
    *(w4) = rp0 * *(in);
    *(w3) = rp0 * *(d3) + rp1 * *(in)
          - rd1 * *(w4);   
    *(w2) = rp0 * *(d2) + rp1 * *(d3) + rp2 * *(in)
          - rd1 * *(w3) - rd2 * *(w4);
    *(w1) = rp0 * *(d1) + rp1 * *(d2) + rp2 * *(d3) + rp3 * *(in)
          - rd1 * *(w2) - rd2 * *(w3) - rd3 * *(w4);
    for (i=4; i<dim; i++,w0++,w1++,w2++,w3++,w4++,d0++,d1++,d2++,d3++) 
      *(w0) = rp0 * *(d0) + rp1 * *(d1) + rp2 * *(d2) + rp3 * *(d3)
            - rd1 * *(w1) - rd2 * *(w2) - rd3 * *(w3) - rd4 * *(w4);
    
    /* on positionne les pointeurs 
     */
    w4 = work2+dim-1;   w3 = w4-1;   w2 = w3-1;   w1 = w2-1;   w0 = w1-1;
    d4 = in+dim-1;      d3 = d4-1;   d2 = d3-1;   d1 = d2-1;
    /*--- calcul de y- ---*/
    *(w4) = 0;
    *(w3) = rn1 * *(d4);
    *(w2) = rn1 * *(d3) + rn2 * *(d4) 
          - rd1 * *(w3);
    *(w1) = rn1 * *(d2) + rn2 * *(d3) + rn3 * *(d4) 
          - rd1 * *(w2) - rd2 * *(w3);
    for (i=dim-5; i>=0; i--,w0--,w1--,w2--,w3--,w4--,d1--,d2--,d3--,d4--)
      *(w0) = rn1 * *(d1) + rn2 * *(d2) + rn3 * *(d3) + rn4 * *(d4)
	    - rd1 * *(w1) - rd2 * *(w2) - rd3 * *(w3) - rd4 * *(w4);

    /*--- calcul final ---*/
    w1 = work1;   w2 = work2;   d0 = out;
    for (i=0 ; i<dim ; i++,w1++,w2++,d0++)
      *d0 = *w1 + *w2;
    
    break;

  case ALPHA_DERICHE :
    
    switch( RFC->derivative ) {
    default :
    case DERIVATIVE_0 :
    case DERIVATIVE_2 :

      rp0 = RFC->sp0;   rp1 = RFC->sp1;
      rd1 = RFC->sd1;   rd2 = RFC->sd2;
      rn1 = RFC->sn1;   rn2 = RFC->sn2;
      
      /* on positionne les pointeurs 
       */
      w2 = work1;   w1 = w2+1;   w0 = w1+1;
      d1 = in+1;    d0 = d1+1;
      /*--- calcul de y+ ---*/
      *(w2) = rp0 * *(in);
      *(w1) = rp0 * *(d1) + rp1 * *(in) 
	    - rd1 * *(w2);     
      for (i=2;  i<dim; i++,w0++,w1++,w2++,d0++,d1++)
	*(w0) = rp0 * *(d0) + rp1 * *(d1)
	      - rd1 * *(w1) - rd2 * *(w2);
      
      w2 = work2+dim-1;   w1 = w2-1;   w0 = w1-1;
      d2 = in+dim-1;      d1 = d2-1;
      /*--- calcul de y- ---*/
      *(w2) = 0.0;
      *(w1) = rn1 * *(d2);
      for (i=dim-3; i>=0; i--,w0--,w1--,w2--,d1--,d2--)
	*(w0) = rn1 * *(d1) + rn2 * *(d2)
	      - rd1 * *(w1) - rd2 * *(w2);
      
      /*--- calcul final ---*/
      w1 = work1;   w2 = work2;   d0 = out;
      for (i=0 ; i<dim ; i++,w1++,w2++,d0++)
	*d0 = *w1 + *w2;
      
      break;
      
    case DERIVATIVE_1 :
    case DERIVATIVE_1_CONTOURS :
      rp1 = RFC->sp1;
      rn1 = RFC->sn1;
      rd1 = RFC->sd1;   rd2 = RFC->sd2;
      
      /* on positionne les pointeurs 
       */
      w2 = work1;   w1 = w2+1;   w0 = w1+1;
      d1 = in+1;
      /*--- calcul de y+ ---*/
      *(w2) = 0.0;
      *(w1) = rp1 * *(in);     
      for (i=2;  i<dim; i++,w0++,w1++,w2++,d1++)
	*(w0) = rp1 * *(d1)
	      - rd1 * *(w1) - rd2 * *(w2);
      
      
      w2 = work2+dim-1;   w1 = w2-1;   w0 = w1-1;
      d2 = in+dim-1;      d1 = d2-1;
      /*--- calcul de y- ---*/
      *(w2) = 0.0;
      *(w1) = rn1 * *(d2);
      for (i=dim-3; i>=0; i--,w0--,w1--,w2--,d1--)
	*(w0) = rn1 * *(d1)
	      - rd1 * *(w1) - rd2 * *(w2);
      
      /*--- calcul final ---*/
      w1 = work1;   w2 = work2;   d0 = out;
      for (i=0 ; i<dim ; i++,w1++,w2++,d0++)
	*d0 = *w1 + *w2;
      
      break;

    case DERIVATIVE_3 :
      rp0 = RFC->sp0;   rp1 = RFC->sp1;
      rd1 = RFC->sd1;   rd2 = RFC->sd2;
      rn0 = RFC->sn0;   rn1 = RFC->sn1;
      
      w2 = work1;   w1 = w2+1;   w0 = w1+1;
      d1 = in+1;   d0 = d1+1;
      /*--- calcul de y+ ---*/
      *(w2) = rp0 * *(in);
      *(w1) = rp0 * *(d1) + rp1 * *(in) 
	    - rd1 * *(w2);     
      for (i=2;  i<dim; i++,w0++,w1++,w2++,d0++,d1++)
	*(w0) = rp0 * *(d0) + rp1 * *(d1)
	      - rd1 * *(w1) - rd2 * *(w2);
      
      w2 = work2+dim-1;   w1 = w2-1;   w0 = w1-1;
      d2 = in+dim-1;      d1 = d2-1;   d0 = d1-1;
      /*--- calcul de y- ---*/
      *(w2) = rn0 * *(d2);
      *(w1) = rn0 * *(d1) + rn1 * *(d2) 
	    - rd1 * *(w2);
      for (i=dim-3; i>=0; i--,w0--,w1--,w2--,d0--,d1--)
	*(w0) = rn0 * *(d0) + rn1 * *(d1)
	      - rd1 * *(w1) - rd2 * *(w2);
      
      /*--- calcul final ---*/
      w1 = work1;   w2 = work2;   d0 = out;
      for (i=0 ; i<dim ; i++,w1++,w2++,d0++)
	*d0 = *w1 + *w2;
      
    }
  }
  return( EXIT_ON_SUCCESS );
}

void Recline_verbose ( )
{
  _VERBOSE_RECLINE_ = 1;
}
void Recline_noverbose ( )
{
  _VERBOSE_RECLINE_ = 0;
}
