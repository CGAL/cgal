// ======================================================================
//
// Copyright (c) 1998,2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : src/Fixed_precision_nt.C
// package       : 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Olivier Devillers <Olivivier.Devillers@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#include <CGAL/Fixed_precision_nt.h>

CGAL_BEGIN_NAMESPACE

////// degree 0 constant
static double Fixed_Or2   = 3.0/9007199254740992.0;      
// 3/2^53             semi-static orientation filter
static double Fixed_Ic2   = 3.0/9007199254740992.0;
// 3/2^54             semi-static incircle filter
static double Fixed_Is2   = 6.0 / 18014398509481984.0;        
// 6/2^54             semi-static insphere filter

Oriented_side side_of_oriented_circleC2 (
      const Fixed_precision_nt& x0, const Fixed_precision_nt& y0,
      const Fixed_precision_nt& x1, const Fixed_precision_nt& y1,
      const Fixed_precision_nt& x2, const Fixed_precision_nt& y2,
      const Fixed_precision_nt& x3, const Fixed_precision_nt& y3)
  // relative position of p3 with respect to circle p0p1p2
  // if p0p1p2 is positively oriented,
  // positive side is the interior of the circle
{  
  double X1=x1.to_double()-x0.to_double();
  double Y1=y1.to_double()-y0.to_double();
  double X2=x2.to_double()-x0.to_double();
  double Y2=y2.to_double()-y0.to_double();
  double X3=x3.to_double()-x0.to_double();
  double Y3=y3.to_double()-y0.to_double();
  double R1=X1*X1+Y1*Y1;
  double R2=X2*X2+Y2*Y2;
  double R3=X3*X3+Y3*Y3;
  double M1=X2*Y3-X3*Y2;
  double M2=X1*Y3-X3*Y1;
  double M3=X1*Y2-X2*Y1;
  double det= R1*M1 -R2*M2 +R3*M3;
  if (det >   Fixed_Ic1) return ON_NEGATIVE_SIDE;
  if (det <= -Fixed_Ic1) return ON_POSITIVE_SIDE;
  // static filter failed, error is small
  {
    double error= 
      CGAL_NTS abs(R1*M1) + CGAL_NTS abs(R2*M2) + CGAL_NTS abs(R3*M3) ;
    error *= Fixed_Ic2 ;
    if ( error < Fixed_Bx4 ) error = 0.0;
    if (det >   error) return ON_NEGATIVE_SIDE;
    if (det <  -error) return ON_POSITIVE_SIDE;
    double m1,m2,m3,r1,r2,r3;
    if (error!=0) {
      // dynamic filter failed, error is very small
      Fixed_split(M1,M1,m1,Fixed_S_2_25); 
      Fixed_split(M2,M2,m2,Fixed_S_2_25);
      Fixed_split(M3,M3,m3,Fixed_S_2_25);
      // Minor i is Mi+mi
      Fixed_split(R1,R1,r1,Fixed_S_2_25);
      Fixed_split(R2,R2,r2,Fixed_S_2_25);
      Fixed_split(R3,R3,r3,Fixed_S_2_25);
      // xi^2+yi^2   is Ri+ri
      det  =  R1*M1 -R2*M2 +R3*M3   ;              // most significant
      // exact because necessary less than 2^80 and multiple 2^52
      det += (R1*m1 + r1*M1) - (R2*m2 + r2*M2) + (R3*m3 + r3*M3) ;
      // exact because necessary less than 2^53 and multiple 2^25
      det += r1*m1 -r2*m2 +r3*m3  ;                       // less significant
      if (det>0) return ON_NEGATIVE_SIDE;
      if (det<0) return ON_POSITIVE_SIDE;
    }
    //cocircular case
    if (! Fixed_incircle_perturb) return ON_ORIENTED_BOUNDARY;

    // apply perturbation method
    // first order coefficient is obtained replacing x^2+y^2 by xy
    R1 = X1*Y1;
    R2 = X2*Y2;
    R3 = X3*Y3;
    det= R1*M1 -R2*M2 +R3*M3;
    if (det >   Fixed_Ic1) return ON_NEGATIVE_SIDE;
    if (det <= -Fixed_Ic1) return ON_POSITIVE_SIDE;
    // static filter failed, error is small
    error= 
      (CGAL_NTS abs(R1*M1)+CGAL_NTS abs(R2*M2)+CGAL_NTS abs(R3*M3))
      * Fixed_Ic2 ;
    if ( error < Fixed_Bx4 ) error = 0.0;
    if (det >   error) return ON_NEGATIVE_SIDE;
    if (det <  -error) return ON_POSITIVE_SIDE;
    if (error!=0) {
      // dynamic filter failed, error is very small
      Fixed_split(R1,R1,r1,Fixed_S_2_25);
      Fixed_split(R2,R2,r2,Fixed_S_2_25);
      Fixed_split(R3,R3,r3,Fixed_S_2_25);
      det  =  R1*M1 -R2*M2 +R3*M3   ;
      det += (R1*m1 + r1*M1) - (R2*m2 + r2*M2) + (R3*m3 + r3*M3) ;
      det += r1*m1 -r2*m2 +r3*m3  ;
      if (det>0) return ON_NEGATIVE_SIDE;
      if (det<0) return ON_POSITIVE_SIDE;
    }
    // first order coefficient is null,
    // compute second order coefficient in the pertrubation method
    //  replace x^2+y^2 by y^2
    R1 = Y1*Y1;
    R2 = Y2*Y2;
    R3 = Y3*Y3;
    det= R1*M1 -R2*M2 +R3*M3;
    if (det >   Fixed_Ic1) return ON_NEGATIVE_SIDE;
    if (det <= -Fixed_Ic1) return ON_POSITIVE_SIDE;
    // static filter failed, error is small
    error= 
      (CGAL_NTS abs(R1*M1)+CGAL_NTS abs(R2*M2)+CGAL_NTS abs(R3*M3))
      * Fixed_Ic2 ;
    if ( error < Fixed_Bx4 ) error = 0.0;
    if (det >   error) return ON_NEGATIVE_SIDE;
    if (det <  -error) return ON_POSITIVE_SIDE;
    if (error!=0) {
      // dynamic filter failed, error is very small
      Fixed_split(R1,R1,r1,Fixed_S_2_25);
      Fixed_split(R2,R2,r2,Fixed_S_2_25);
      Fixed_split(R3,R3,r3,Fixed_S_2_25);
      det  =  R1*M1 -R2*M2 +R3*M3   ;
      det += (R1*m1 + r1*M1) - (R2*m2 + r2*M2) + (R3*m3 + r3*M3) ;
      det += r1*m1 -r2*m2 +r3*m3  ;
      if (det>0) return ON_NEGATIVE_SIDE;
      if (det<0) return ON_POSITIVE_SIDE;
    }
    // points are colinear
    return ON_ORIENTED_BOUNDARY;
  }
}

Orientation orientationC3
(   const Fixed_precision_nt& x0, const Fixed_precision_nt& y0, 
    const Fixed_precision_nt& z0,
    const Fixed_precision_nt& x1, const Fixed_precision_nt& y1, 
    const Fixed_precision_nt& z1,
    const Fixed_precision_nt& x2, const Fixed_precision_nt& y2,
    const Fixed_precision_nt& z2,
    const Fixed_precision_nt& x3, const Fixed_precision_nt& y3,
    const Fixed_precision_nt& z3)
{
  double X1=x1.to_double() -x0.to_double();
  double Y1=y1.to_double() -y0.to_double();
  double Z1=z1.to_double() -z0.to_double(); 
  double X2=x2.to_double() -x0.to_double(); 
  double Y2=y2.to_double() -y0.to_double(); 
  double Z2=z2.to_double() -z0.to_double();
  double X3=x3.to_double() -x0.to_double();
  double Y3=y3.to_double() -y0.to_double();
  double Z3=z3.to_double() -z0.to_double();
  double M1=Y2*Z3-Y3*Z2;
  double M2=Y1*Z3-Y3*Z1;
  double M3=Y1*Z2-Y2*Z1;
  double det= X1*M1 - X2*M2 + X3*M3 ;
  if (det >   Fixed_Or1 ) return POSITIVE;
  if (det <  -Fixed_Or1 ) return NEGATIVE;
  /* det is small */
  double error= 
    CGAL_NTS abs(X1*M1) + CGAL_NTS abs(X2*M2) + CGAL_NTS abs(X3*M3);
  error *= Fixed_Or2;
  if ( error < Fixed_Bx3) error = 0.0;
  if (det >   error) return POSITIVE;
  if (det <  -error) return NEGATIVE;
  /* det is very small, exact computation must be done */
  if (error!=0.0) {  // otherwise, exact value 0 already certified
    Fixed_split(M1,M1,Z1,Fixed_S_2_25);
    Fixed_split(M2,M2,Z2,Fixed_S_2_25);
    Fixed_split(M3,M3,Z3,Fixed_S_2_25);
    det  =  X1*M1 - X2*M2 + X3*M3 ;
    det +=  X1*Z1 - X2*Z2 + X3*Z3 ;     /* less significant */
    if (det>0) return POSITIVE;
    if (det<0) return NEGATIVE;
  }
  /* points are coplanar */
  return COPLANAR;
}

Oriented_side side_of_oriented_sphereC3 
(     const Fixed_precision_nt& x0, const Fixed_precision_nt& y0, 
      const Fixed_precision_nt& z0,
      const Fixed_precision_nt& x1, const Fixed_precision_nt& y1,
      const Fixed_precision_nt& z1,
      const Fixed_precision_nt& x2, const Fixed_precision_nt& y2,
      const Fixed_precision_nt& z2,
      const Fixed_precision_nt& x3, const Fixed_precision_nt& y3,
      const Fixed_precision_nt& z3,
      const Fixed_precision_nt& x4, const Fixed_precision_nt& y4,
      const Fixed_precision_nt& z4)
  // relative position of p4 with respect to sphere p0p1p2p3
  // if p0p1p2p3 is positively oriented,
  // positive side is the interior of the sphere
{  double det,error;
  // transform in 4x4 determinant
  // point p4 is inside sphere through oriented tetra p0p1p2p3
  // if the determinant is negative
  double X1=x1.to_double() -x0.to_double() ;
  double Y1=y1.to_double() -y0.to_double() ;
  double Z1=z1.to_double() -z0.to_double() ;
  double X2=x2.to_double() -x0.to_double() ;  //   | X1 X2 X3 X4 |
  double Y2=y2.to_double() -y0.to_double() ;  //   | Y1 Y2 Y3 Y4 |
  double Z2=z2.to_double() -z0.to_double() ;  //   | Z1 Z2 Z3 Z4 |
  double X3=x3.to_double() -x0.to_double() ;  //   | R1 R2 R3 R4 |
  double Y3=y3.to_double() -y0.to_double() ;
  double Z3=z3.to_double() -z0.to_double() ;
  double X4=x4.to_double() -x0.to_double() ;
  double Y4=y4.to_double() -y0.to_double() ;
  double Z4=z4.to_double() -z0.to_double() ;
  double R1= X1*X1+Y1*Y1+Z1*Z1;
  double R2= X2*X2+Y2*Y2+Z2*Z2;
  double R3= X3*X3+Y3*Y3+Z3*Z3;
  double R4= X4*X4+Y4*Y4+Z4*Z4;
  // compute 2x2 minors on the two first lines
  double M12 = X1 * Y2 - X2 * Y1;
  double M23 = X2 * Y3 - X3 * Y2;
  double M34 = X3 * Y4 - X4 * Y3;
  double M13 = X1 * Y3 - X3 * Y1;
  double M14 = X1 * Y4 - X4 * Y1;
  double M24 = X2 * Y4 - X4 * Y2;
  // compute 3x3 minors on the three first lines
  double MM1 =   Z2 * M34 - Z3 * M24 + Z4 * M23 ;
  double MM2 =   Z1 * M34 - Z3 * M14 + Z4 * M13 ;
  double MM3 =   Z1 * M24 - Z2 * M14 + Z4 * M12 ;
  double MM4 =   Z1 * M23 - Z2 * M13 + Z3 * M12 ;
  // compute determinant
  det= (R2 * MM2 - R1 * MM1) + (R4 * MM4 - R3 * MM3) ;
  if (det >=  Fixed_Is1 ) return ON_NEGATIVE_SIDE;
  if (det <= -Fixed_Is1 ) return ON_POSITIVE_SIDE;
  // static filter failed, error is small
  {
    double mm1 =CGAL_NTS abs(Z2*M34) + CGAL_NTS abs(Z3*M24)
      + CGAL_NTS abs(Z4*M23) ;
    double mm2 =CGAL_NTS abs(Z1*M34) + CGAL_NTS abs(Z3*M14) 
      + CGAL_NTS abs(Z4*M13) ;
    double mm3 =CGAL_NTS abs(Z1*M24) + CGAL_NTS abs(Z2*M14) 
      + CGAL_NTS abs(Z4*M12) ;
    double mm4 =CGAL_NTS abs(Z1*M23) + CGAL_NTS abs(Z2*M13) 
      + CGAL_NTS abs(Z3*M12) ;
    error= CGAL_NTS abs(R1)*mm1 + CGAL_NTS abs(R2)*mm2 
           + CGAL_NTS abs(R3)*mm3 + CGAL_NTS abs(R4)*mm4;
    error *= Fixed_Is2;
    if ( error < Fixed_Bx5 ) error = 0.0;
    if (det >   error) return ON_NEGATIVE_SIDE;
    if (det <  -error) return ON_POSITIVE_SIDE;
    if (error!=0) {
      // dynamic filter failed, error is very small
      double a1,a2,a3,a4;
      double b1,b2,b3,b4;
      double c1,c2,c3,c4;
      double d1,d2,d3,d4;
      if (error!=0.0){
	// start exact computation
	Fixed_split(M12,a2,M12,Fixed_S_2_25);
	Fixed_split(M13,a3,M13,Fixed_S_2_25);
	Fixed_split(M14,a4,M14,Fixed_S_2_25);
	Fixed_split(M23,b3,M23,Fixed_S_2_25);
	Fixed_split(M24,b4,M24,Fixed_S_2_25);
	Fixed_split(M34,c4,M34,Fixed_S_2_25);
	// compute 3x3 minors   ci+di is minor i
	c1 =   Z2 *  c4 - Z3 *  b4 + Z4 *  b3 ;      // most significant
	c2 =   Z1 *  c4 - Z3 *  a4 + Z4 *  a3 ;
	c3 =   Z1 *  b4 - Z2 *  a4 + Z4 *  a2 ;
	c4 =   Z1 *  b3 - Z2 *  a3 + Z3 *  a2 ;
	d1 =   Z2 * M34 - Z3 * M24 + Z4 * M23 ;       // less significant
	d2 =   Z1 * M34 - Z3 * M14 + Z4 * M13 ;
	d3 =   Z1 * M24 - Z2 * M14 + Z4 * M12 ;
	d4 =   Z1 * M23 - Z2 * M13 + Z3 * M12 ;
	Fixed_split(d1,a1,d1,Fixed_S_3_25);
	Fixed_split(d2,a2,d2,Fixed_S_3_25);
	Fixed_split(d3,a3,d3,Fixed_S_3_25);
	Fixed_split(d4,a4,d4,Fixed_S_3_25);
	// now ci+ai+di is minor i
	a1 += c1;
	a2 += c2;
	a3 += c3;
	a4 += c4;
	//  now ai+di is minor i
	Fixed_split(a1,c1,a1,Fixed_S_3_51);
	Fixed_split(a2,c2,a2,Fixed_S_3_51);
	Fixed_split(a3,c3,a3,Fixed_S_3_51);
	Fixed_split(a4,c4,a4,Fixed_S_3_51);
	// now ci+ai+di is minor i,  non overlapping 25 bits for each
	Fixed_split(R1,b1,R1,Fixed_S_2_25);
	Fixed_split(R2,b2,R2,Fixed_S_2_25);
	Fixed_split(R3,b3,R3,Fixed_S_2_25);
	Fixed_split(R4,b4,R4,Fixed_S_2_25);
	// now bi+Ri is last line, non overlapping 25 bits for each
	det  = (b2 * c2 - b1 * c1)                          ;      // 1
	det +=                       (b4 * c4 - b3 * c3)    ;      // 1b
	det += (b2 * a2 - b1 * a1) + (b4 * a4 - b3 * a3)    ;      // 2
	det += (R2 * c2 - R1 * c1) + (R4 * c4 - R3 * c3)    ;      // 3
	det += (b2 * d2 - b1 * d1) + (b4 * d4 - b3 * d3)    ;      // 4
	det += (R2 * a2 - R1 * a1) + (R4 * a4 - R3 * a3)    ;      // 5
	det += (R2 * d2 - R1 * d1) + (R4 * d4 - R3 * d3)    ;      // 6
	if (det>0) return ON_NEGATIVE_SIDE;
	if (det<0) return ON_POSITIVE_SIDE;
      }
      // points are cospherical
      
      if (! Fixed_insphere_perturb)return ON_ORIENTED_BOUNDARY;

      // apply perturbation method
      // first order coefficient is obtained replacing x^2+y^2+z^2 by xy
      R1 = X1*Y1;
      R2 = X2*Y2;
      R3 = X3*Y3;
      R4 = X4*Y4;
      // compute determinant
      det= (R2 * MM2 - R1 * MM1) + (R4 * MM4 - R3 * MM3) ;
      if (det >=  Fixed_Is1 ) return ON_NEGATIVE_SIDE;
      if (det <= -Fixed_Is1 ) return ON_POSITIVE_SIDE;
      // static filter failed, error is small
      error= CGAL_NTS abs(R1)*mm1 + CGAL_NTS abs(R2)*mm2 
	     + CGAL_NTS abs(R3)*mm3 + CGAL_NTS abs(R4)*mm4;
      error *= Fixed_Is2;
      if ( error < Fixed_Bx5 ) error = 0.0;
      if (det >   error) return ON_NEGATIVE_SIDE;
      if (det <  -error) return ON_POSITIVE_SIDE;
      if (error!=0.0){
	// dynamic filter failed, error is very small
	// start exact computation
	Fixed_split(R1,b1,R1,Fixed_S_2_25);
	Fixed_split(R2,b2,R2,Fixed_S_2_25);
	Fixed_split(R3,b3,R3,Fixed_S_2_25);
	Fixed_split(R4,b4,R4,Fixed_S_2_25);
	// now bi+Ri is last line, non overlapping 25 bits for each
	det  = (b2 * c2 - b1 * c1)                          ;      // 1
	det +=                       (b4 * c4 - b3 * c3)    ;      // 1b
	det += (b2 * a2 - b1 * a1) + (b4 * a4 - b3 * a3)    ;      // 2
	det += (R2 * c2 - R1 * c1) + (R4 * c4 - R3 * c3)    ;      // 3
	det += (b2 * d2 - b1 * d1) + (b4 * d4 - b3 * d3)    ;      // 4
	det += (R2 * a2 - R1 * a1) + (R4 * a4 - R3 * a3)    ;      // 5
	det += (R2 * d2 - R1 * d1) + (R4 * d4 - R3 * d3)    ;      // 6
	if (det>0) return ON_NEGATIVE_SIDE;
	if (det<0) return ON_POSITIVE_SIDE;
      }
      // first order coefficient is null,
      // compute second order coefficient in the pertrubation method
      //  replace x^2+y^2+z^2 by y^2
      R1 = Y1*Y1;
      R2 = Y2*Y2;
      R3 = Y3*Y3;
      R4 = Y4*Y4;
      // compute determinant
      det= (R2 * MM2 - R1 * MM1) + (R4 * MM4 - R3 * MM3) ;
      if (det >=  Fixed_Is1 ) return ON_NEGATIVE_SIDE;
      if (det <= -Fixed_Is1 ) return ON_POSITIVE_SIDE;
      // static filter failed, error is small
      error= CGAL_NTS abs(R1)*mm1 + CGAL_NTS abs(R2)*mm2
	     + CGAL_NTS abs(R3)*mm3 + CGAL_NTS abs(R4)*mm4;
      error *= Fixed_Is2;
      if ( error < Fixed_Bx5 ) error = 0.0;
      if (det >   error) return ON_NEGATIVE_SIDE;
      if (det <  -error) return ON_POSITIVE_SIDE;
      if (error!=0.0){
	// dynamic filter failed, error is very small
	// start exact computation
	Fixed_split(R1,b1,R1,Fixed_S_2_25);
	Fixed_split(R2,b2,R2,Fixed_S_2_25);
	Fixed_split(R3,b3,R3,Fixed_S_2_25);
	Fixed_split(R4,b4,R4,Fixed_S_2_25);
	// now bi+Ri is last line, non overlapping 25 bits for each
	det  = (b2 * c2 - b1 * c1)                          ;      // 1
	det +=                       (b4 * c4 - b3 * c3)    ;      // 1b
	det += (b2 * a2 - b1 * a1) + (b4 * a4 - b3 * a3)    ;      // 2
	det += (R2 * c2 - R1 * c1) + (R4 * c4 - R3 * c3)    ;      // 3
	det += (b2 * d2 - b1 * d1) + (b4 * d4 - b3 * d3)    ;      // 4
	det += (R2 * a2 - R1 * a1) + (R4 * a4 - R3 * a3)    ;      // 5
	det += (R2 * d2 - R1 * d1) + (R4 * d4 - R3 * d3)    ;      // 6
	if (det>0) return ON_NEGATIVE_SIDE;
	if (det<0) return ON_POSITIVE_SIDE;
      }
      // 2nd order coefficient is null,
      // compute 3rd order coefficient in the pertrubation method
      //  replace x^2+y^2+z^2 by yz
      R1 = Y1*Z1;
      R2 = Y2*Z2;
      R3 = Y3*Z3;
      R4 = Y4*Z4;
      // compute determinant
      det= (R2 * MM2 - R1 * MM1) + (R4 * MM4 - R3 * MM3) ;
      if (det >=  Fixed_Is1 ) return ON_NEGATIVE_SIDE;
      if (det <= -Fixed_Is1 ) return ON_POSITIVE_SIDE;
      // static filter failed, error is small
      error= CGAL_NTS abs(R1)*mm1 + CGAL_NTS abs(R2)*mm2
	     + CGAL_NTS abs(R3)*mm3 + CGAL_NTS abs(R4)*mm4;
      error *= Fixed_Is2;
      if ( error < Fixed_Bx5 ) error = 0.0;
      if (det >   error) return ON_NEGATIVE_SIDE;
      if (det <  -error) return ON_POSITIVE_SIDE;
      if (error!=0.0){
	// dynamic filter failed, error is very small
	// start exact computation
	Fixed_split(R1,b1,R1,Fixed_S_2_25);
	Fixed_split(R2,b2,R2,Fixed_S_2_25);
	Fixed_split(R3,b3,R3,Fixed_S_2_25);
	Fixed_split(R4,b4,R4,Fixed_S_2_25);
	// now bi+Ri is last line, non overlapping 25 bits for each
	det  = (b2 * c2 - b1 * c1)                          ;      // 1
	det +=                       (b4 * c4 - b3 * c3)    ;      // 1b
	det += (b2 * a2 - b1 * a1) + (b4 * a4 - b3 * a3)    ;      // 2
	det += (R2 * c2 - R1 * c1) + (R4 * c4 - R3 * c3)    ;      // 3
	det += (b2 * d2 - b1 * d1) + (b4 * d4 - b3 * d3)    ;      // 4
	det += (R2 * a2 - R1 * a1) + (R4 * a4 - R3 * a3)    ;      // 5
	det += (R2 * d2 - R1 * d1) + (R4 * d4 - R3 * d3)    ;      // 6
	if (det>0) return ON_NEGATIVE_SIDE;
	if (det<0) return ON_POSITIVE_SIDE;
      }
      // 3rd order coefficient is null,
      // compute 4th order coefficient in the pertrubation method
      //  replace x^2+y^2+z^2 by xz
      R1 = X1*Z1;
      R2 = X2*Z2;
      R3 = X3*Z3;
      R4 = X4*Z4;
      // compute determinant
      det= (R2 * MM2 - R1 * MM1) + (R4 * MM4 - R3 * MM3) ;
      if (det >=  Fixed_Is1 ) return ON_NEGATIVE_SIDE;
      if (det <= -Fixed_Is1 ) return ON_POSITIVE_SIDE;
      // static filter failed, error is small
      error= CGAL_NTS abs(R1)*mm1 + CGAL_NTS abs(R2)*mm2
	     + CGAL_NTS abs(R3)*mm3 + CGAL_NTS abs(R4)*mm4;
      error *= Fixed_Is2;
      if ( error < Fixed_Bx5 ) error = 0.0;
      if (det >   error) return ON_NEGATIVE_SIDE;
      if (det <  -error) return ON_POSITIVE_SIDE;
      if (error!=0.0){
	// dynamic filter failed, error is very small
	// start exact computation
	Fixed_split(R1,b1,R1,Fixed_S_2_25);
	Fixed_split(R2,b2,R2,Fixed_S_2_25);
	Fixed_split(R3,b3,R3,Fixed_S_2_25);
	Fixed_split(R4,b4,R4,Fixed_S_2_25);
	// now bi+Ri is last line, non overlapping 25 bits for each
	det  = (b2 * c2 - b1 * c1)                          ;      // 1
	det +=                       (b4 * c4 - b3 * c3)    ;      // 1b
	det += (b2 * a2 - b1 * a1) + (b4 * a4 - b3 * a3)    ;      // 2
	det += (R2 * c2 - R1 * c1) + (R4 * c4 - R3 * c3)    ;      // 3
	det += (b2 * d2 - b1 * d1) + (b4 * d4 - b3 * d3)    ;      // 4
	det += (R2 * a2 - R1 * a1) + (R4 * a4 - R3 * a3)    ;      // 5
	det += (R2 * d2 - R1 * d1) + (R4 * d4 - R3 * d3)    ;      // 6
	if (det>0) return ON_NEGATIVE_SIDE;
	if (det<0) return ON_POSITIVE_SIDE;
      }
      // 4th order coefficient is null,
      // compute 5th order coefficient in the pertrubation method
      //  replace x^2+y^2+z^2 by z^2
      R1 = Z1*Z1;
      R2 = Z2*Z2;
      R3 = Z3*Z3;
      R4 = Z4*Z4;
      // compute determinant
      det= (R2 * MM2 - R1 * MM1) + (R4 * MM4 - R3 * MM3) ;
      if (det >=  Fixed_Is1 ) return ON_NEGATIVE_SIDE;
      if (det <= -Fixed_Is1 ) return ON_POSITIVE_SIDE;
      // static filter failed, error is small
      error= CGAL_NTS abs(R1)*mm1 + CGAL_NTS abs(R2)*mm2
	     + CGAL_NTS abs(R3)*mm3 + CGAL_NTS abs(R4)*mm4;
      error *= Fixed_Is2;
      if ( error < Fixed_Bx5 ) error = 0.0;
      if (det >   error) return ON_NEGATIVE_SIDE;
      if (det <  -error) return ON_POSITIVE_SIDE;
      if (error!=0.0){
	// dynamic filter failed, error is very small
	// start exact computation
	Fixed_split(R1,b1,R1,Fixed_S_2_25);
	Fixed_split(R2,b2,R2,Fixed_S_2_25);
	Fixed_split(R3,b3,R3,Fixed_S_2_25);
	Fixed_split(R4,b4,R4,Fixed_S_2_25);
	// now bi+Ri is last line, non overlapping 25 bits for each
	det  = (b2 * c2 - b1 * c1)                          ;      // 1
	det +=                       (b4 * c4 - b3 * c3)    ;      // 1b
	det += (b2 * a2 - b1 * a1) + (b4 * a4 - b3 * a3)    ;      // 2
	det += (R2 * c2 - R1 * c1) + (R4 * c4 - R3 * c3)    ;      // 3
	det += (b2 * d2 - b1 * d1) + (b4 * d4 - b3 * d3)    ;      // 4
	det += (R2 * a2 - R1 * a1) + (R4 * a4 - R3 * a3)    ;      // 5
	det += (R2 * d2 - R1 * d1) + (R4 * d4 - R3 * d3)    ;      // 6
	if (det>0) return ON_NEGATIVE_SIDE;
	if (det<0) return ON_POSITIVE_SIDE;
      }
      // 5th order coefficient is null,
      // compute 6th order coefficient in the pertrubation method
      //  replace x^2+y^2+z^2 by x^2
      R1 = X1*X1;
      R2 = X2*X2;
      R3 = X3*X3;
      R4 = X4*X4;
      // compute determinant
      det= (R2 * MM2 - R1 * MM1) + (R4 * MM4 - R3 * MM3) ;
      if (det >=  Fixed_Is1 ) return ON_NEGATIVE_SIDE;
      if (det <= -Fixed_Is1 ) return ON_POSITIVE_SIDE;
      // static filter failed, error is small
      error= CGAL_NTS abs(R1)*mm1 + CGAL_NTS abs(R2)*mm2 
	     + CGAL_NTS abs(R3)*mm3 + CGAL_NTS abs(R4)*mm4;
      error *= Fixed_Is2;
      if ( error < Fixed_Bx5 ) error = 0.0;
      if (det >   error) return ON_NEGATIVE_SIDE;
      if (det <  -error) return ON_POSITIVE_SIDE;
      if (error!=0.0){
	// dynamic filter failed, error is very small
	// start exact computation
	Fixed_split(R1,b1,R1,Fixed_S_2_25);
	Fixed_split(R2,b2,R2,Fixed_S_2_25);
	Fixed_split(R3,b3,R3,Fixed_S_2_25);
	Fixed_split(R4,b4,R4,Fixed_S_2_25);
	// now bi+Ri is last line, non overlapping 25 bits for each
	det  = (b2 * c2 - b1 * c1)                          ;      // 1
	det +=                       (b4 * c4 - b3 * c3)    ;      // 1b
	det += (b2 * a2 - b1 * a1) + (b4 * a4 - b3 * a3)    ;      // 2
	det += (R2 * c2 - R1 * c1) + (R4 * c4 - R3 * c3)    ;      // 3
	det += (b2 * d2 - b1 * d1) + (b4 * d4 - b3 * d3)    ;      // 4
	det += (R2 * a2 - R1 * a1) + (R4 * a4 - R3 * a3)    ;      // 5
	det += (R2 * d2 - R1 * d1) + (R4 * d4 - R3 * d3)    ;      // 6
	if (det>0) return ON_NEGATIVE_SIDE;
	if (det<0) return ON_POSITIVE_SIDE;
      }
      // 6th order coefficient is null,
      // points are coplanar
      return ON_ORIENTED_BOUNDARY;
    }
  }
  return ON_ORIENTED_BOUNDARY;
 // code never goes here, it is just to avoid compilation warning
}

CGAL_END_NAMESPACE
