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
 * recline.h - tools for recursive filtering of 1D lines
 *
 * $Id$
 *
 * Copyright©INRIA 1998
 *
 * DESCRIPTION: 
 *
 * Recursive filtering of a line (a 1D array)
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * http://www.inria.fr/epidaure/personnel/malandain/
 * 
 * CREATION DATE: 
 * June, 9 1998
 *
 * Copyright Gregoire Malandain, INRIA
 *
 * ADDITIONS, CHANGES
 *
 */

/*
 * recursive filtering of a line (a 1D array)
 */

#ifndef _recline_h_
#define _recline_h_



/* The different recursive filter's types.
 *
 * DESCRIPTION:
 *
 * - ALPHA_DERICHE is the first recurvise filter designed
 *   by R. Deriche. See REFERENCES.
 *
 * - with ALPHA_DERICHE's filters, one can either compute
 *   derivatives from order 0 (smoothing) to 3, or extract edges.
 *
 * - GAUSSIAN_DERICHE is a 4th order recursive filter which
 *   approximates the gaussien. See
 *   "Recursively Implementing The Gaussian and Its Derivatives",
 *   R. Deriche, International Conference On Image Processing,
 *   pp 263-267, Singapore, September 1992. Also INRIA research
 *   report.
 *
 * - with GAUSSIAN_DERICHE's filters, one can either compute
 *   derivatives from order 0 (smoothing) to 2, or extract edges.
 *
 * - Extracting edges with ALPHA_DERICHE's filters is faster but
 *   the modulus of the gradient (the estimated height of the step 
 *   edge) depens on the gradient orientation because the filter
 *   is not isotropic. Heights are better estimated with 
 *   GAUSSIAN_DERICHE's filters but they seem not be perfectly
 *   symmetrical.
 *
 * REFERENCES:
 *
 * - "Optimal edge detection using recursive filtering", R. Deriche,
 *   International Journal of Computer Vision, pp 167-187, 1987.
 *
 * - "Recursive filtering and edge tracking: two primary tools
 *    for 3-D edge detection", O. Monga, R. Deriche,
 *   G. Malandain and J.-P. Cocquerez, Image and Vision
 *   Computing 4:9, pp 203-214, August 1991.
 */
typedef enum {
  UNKNOWN_FILTER = 0 /* unknown filter type */,
  ALPHA_DERICHE = 1 /* Deriche's filter (exponential (- alpha |X|)) */,
  GAUSSIAN_DERICHE = 2 /* gaussian approximation (Deriche's coefficients) */,
  GAUSSIAN_FIDRICH = 3 /* gaussian approximation (Fidrich's coefficients) */
} recursiveFilterType;



/* Order of the derivative to be computed.
 *
 * DESCRIPTION:
 *
 * - NODERIVATIVE nothing will be done.
 *
 * - DERIVATIVE_0 means smoothing.
 *
 * - DERIVATIVE_1 first derivative. The normalization
 *   of the filter is made so that the response to the
 *   signal i=x will be 1.
 *
 * - DERIVATIVE_1_CONTOURS first derivative but adapted
 *   to edge detections. The normalization of the filter 
 *   is made so that the response to a step edge is 
 *   the step edge height.
 *
 * - DERIVATIVE_2 second derivative. The normalization
 *   of the filter is made so that the response to the
 *   signal i=x*2/2 will be 1.
 *
 * - DERIVATIVE_3 third derivative. The normalization
 *   of the filter is made so that the response to the
 *   signal i=x*3/6 will be 1.
 */
typedef enum {
  NODERIVATIVE  = -1 /* no derivative (no filtering) */,
  DERIVATIVE_0  = 0 /* smoothing */,
  SMOOTHING     = 0 /* smoothing */,
  DERIVATIVE_1  = 1 /* derivative of order 1 */,
  DERIVATIVE_2  = 2 /* derivative of order 2 */,
  DERIVATIVE_3  = 3 /* derivative of order 3 */,
  DERIVATIVE_1_CONTOURS = 11 /* derivative of order 1, normalization adapted to
				contours. The response to a step-edge is the 
				height of the step. */,
  DERIVATIVE_1_EDGES = 11 /* derivative of order 1, normalization adapted to
				contours. The response to a step-edge is the 
				height of the step. */
} derivativeOrder;



typedef struct {
  /*--- denominateur       ---*/
  double sd1;
  double sd2;
  double sd3;
  double sd4;
  /*--- numerateur positif ---*/
  double sp0;
  double sp1;
  double sp2;
  double sp3;
  /*--- numerateur negatif ---*/
  double sn0;
  double sn1;
  double sn2;
  double sn3;
  double sn4;
  /*--- type de filtre en cours ---*/
  recursiveFilterType type_filter;
  derivativeOrder derivative;
} RFcoefficientType;



/* Initialization of coefficients for recursive filtering.
 *
 * PARAMETERS:
 *
 * - the coefficient is the sigma's value in case of
 *   gaussian filtering, or the alpha's value in case
 *   of Deriche's filters.
 *
 * - the coefficient's value must be larger than 0.1
 *   in case of gaussian filtering, and in the
 *   [0.1,1.9] range in case of Deriche's filters.
 *
 * SEE:
 *
 * - recursiveFilterType
 *
 * - derivativeOrder
 */
extern RFcoefficientType * InitRecursiveCoefficients( double x, /* coefficient's value */
				       recursiveFilterType filterType, /* filter's type */
				       derivativeOrder derivative /* derivative's order */ );



/* 1D recursive filtering along a line.
 *
 * WARNING:
 * Coefficients should already be initialized.
 *
 * SEE:
 *
 * - recursiveFilterType
 *
 * - derivativeOrder
 *
 * RETURN:
 *
 * - 0 in case of error
 *
 * - 1 if successful
 */
extern int RecursiveFilter1D( RFcoefficientType *RFC,
			      double *in, /* input line */ 
			      double *out, /* output line */
			      double *work1, /* first work array */
			      double *work2, /* second work array, 
						could be out if out is different from in */
			      int dim /* lines' length */ );


/* Turn on verbose mode.
 *
 * DESCRIPTION:
 * Some information will be written on stderr when processing.
 */
extern void Recline_verbose ( );

/* Turn off verbose mode.
 *
 * DESCRIPTION:
 * Nothing will be written on stderr when processing.
 */
extern void Recline_noverbose ( );


#endif /* _recline_h_ */
