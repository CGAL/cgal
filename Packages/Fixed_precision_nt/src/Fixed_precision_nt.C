// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : src/Fixed_precision_nt.C
// revision      : $Revision$
// revision_date : $Date$
// package       : Fixed number type
// author(s)     : Olivier Devillers <Olivivier.Devillers@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================


// Fixed implement 24 bits integers using float
// The user has to initiate a multiplicator factor to be applied to
// original data.
// if CGAL is defined, Fixed are encapsulated in class Fixed_precision_nt


//  for compilation outside CGAL do something :
//  s/Fixed_/Fixed_/g
//  s!CGAL/Fixed_precision_nt.h!Fixed.h!


#include <CGAL/Fixed_precision_nt.h>
#include <CGAL/Interval_arithmetic/_FPU.h>
/*

// Fixed.h for use outside CGAL below
// there is no class Fixed, float are used directly and must be
// rounded (on user responsability)


// intialisation, bounds
bool Fixed_init(float);
  // return true if init succeed false otherwise
  // Precondition : b positive
  // Precondition : non yet initialized
  // parameter b is the upper bound on absolute value on all entries
float Fixed_unit_value();
float Fixed_upper_bound();

// round to the nearest legal fixed
void Fixed_round(float&);

// geometric predicates
int Fixed_orientation(float x0, float y0, 
		      float x1, float y1, 
		      float x2, float y2);
int Fixed_orientation(  float x0, float y0, float z0,
			float x1, float y1, float z1,
			float x2, float y2, float z2,
			float x3, float y3, float z3);
int Fixed_insphere (
      float x0, float y0,
      float x1, float y1,
      float x2, float y2,
      float x3, float y3);
  // relative position of p3 with respect to circle p0p1p2
  // if p0p1p2 is positively oriented,
  // positive side is the interior of the circle
  // perturbation mode apply to cocircular and non colinear 4uples
void  Fixed_perturb_incircle();
void  Fixed_unperturb_incircle();
bool  Fixed_is_perturbed_incircle();
int Fixed_insphere (
      float x0, float y0, float z0,
      float x1, float y1, float z1,
      float x2, float y2, float z2,
      float x3, float y3, float z3,
      float x4, float y4, float z4);
  // relative position of p4 with respect to sphere p0p1p2p3
  // if p0p1p2p3 is positively oriented,
  // positive side is the interior of the sphere
  // perturbation mode apply to cocircular and non coplanar 5uples
void  Fixed_perturb_insphere();
void  Fixed_unperturb_insphere();
bool  Fixed_is_perturbed_insphere();
*/


#ifdef CGAL_FIXED_PRECISION_NT_H  // we are using CGAL

#define CGAL_Fixed_public       inline
#include <CGAL/assertions.h>
#include <CGAL/double.h>
#include <CGAL/number_utils.h>
#define CGAL_Fixed_abs abs
#define CGAL_Fixed_rounding_precondition_msg( v ) \
    CGAL_warning_msg((( v >=-Fixed_B24) && ( v <=Fixed_B24)),\
		     "Warning : Fixed_precision_nt overflow")


CGAL_BEGIN_NAMESPACE



#else //CGAL_FIXED_PRECISION_NT_H

#include <iostream.h>
#define CGAL_Fixed_public
extern "C" double fabs(double);
#define CGAL_Fixed_abs fabs
#define CGAL_Fixed_rounding_precondition_msg( v ) \
   if( ( v <-Fixed_B24) || ( v >Fixed_B24)) \
    cerr<<"Warning : Fixed overflow |"<<f<<"| > "<<Fixed_B24<<endl
typedef int bool;
const int false=0;
const int true =0;

#endif// else CGAL_FIXED_PRECISION_NT_H





// ======================================================================
//--------- static variables declaration and initialization
// ======================================================================

// static marks to decide if perturbation are used
static bool Fixed_incircle_perturb= false;
static bool Fixed_insphere_perturb= false;

////// degree 0 constant
static double Fixed_Or2   = 3.0/9007199254740992.0;      
// 3/2^53             semi-static orientation filter
static double Fixed_Ic2   = 3.0/9007199254740992.0;
// 3/2^54             semi-static incircle filter
static double Fixed_Is2   = 6.0 / 18014398509481984.0;        
// 6/2^54             semi-static insphere filter
////// degree 1 constant
static double Fixed_B0=0.0;              
//   precision on coordinates
// assumed to be 1.0 to fix value below
static double Fixed_B24   = 16777216.0;          
// 2^24               bound on coordinates
static double Fixed_SNAP  = 6755399441055744.0;
// 3 * 2^51           for rounding of coordinates
////// degree 2 constant
static double Fixed_S_2_25= 453347182355485940514816.0;
// 3 * 2^77           split to bit B0^2 * 2^25
////// degree 3 constant
static double Fixed_S_3_25= Fixed_S_2_25;           
// 3 * 2^77           split to bit B0^3 * 2^25
static double Fixed_S_3_51= 30423614405477505635920876929024.0;           
// 3 * 2^103          split to bit B0^3 * 2^51
static double Fixed_Or1   = 37748736.0;
// 9*2^22             static error for orientation filter
static double Fixed_Bx3   = 0.5;     
// 2^(-1)             half degree 3 unit for semi-static filter
////// degree 4 constant
static double Fixed_Ic1   = 2533274790395904.0;    
// 9*2^48;           static error for incircle filter
static double Fixed_Bx4   = 0.5;     
// 2^(-1)             half degree 4 unit for semi-static filter
////// degree 5 constant
static double Fixed_Is1   = 1416709944860893564108800.0;    
// 75*2^74;           static error for insphere filter
static double Fixed_Bx5   = 0.5;     
// 2^(-1)             half degree 5 unit for semi-static filter

// ======================================================================
//--------- static access and parametrization functions
// ======================================================================


CGAL_Fixed_public float Fixed_unit_value(){return Fixed_B0;}
CGAL_Fixed_public float Fixed_upper_bound(){return Fixed_B24;}
CGAL_Fixed_public void  Fixed_perturb_incircle() 
{Fixed_incircle_perturb=true;}
CGAL_Fixed_public void  Fixed_unperturb_incircle() 
{Fixed_incircle_perturb=false;}
CGAL_Fixed_public bool  Fixed_is_perturbed_incircle() 
{return Fixed_incircle_perturb;}
CGAL_Fixed_public void  Fixed_perturb_insphere() 
{Fixed_insphere_perturb=true;}
CGAL_Fixed_public void  Fixed_unperturb_insphere() 
{Fixed_insphere_perturb=false;}
CGAL_Fixed_public bool  Fixed_is_perturbed_insphere() 
{return Fixed_insphere_perturb;}

// ======================================================================
//--------- initialize the upper bound on fixed
// ======================================================================

CGAL_Fixed_public bool Fixed_init(float b)
  // return true if init succeed false otherwise
  // Precondition : b positive
  // Precondition : non yet initialized
  // parameter b is the upper bound on absolute value on all entries
{
  bool result = true;
  if (b<=0) b=-b;
  if (Fixed_B0!=0) result =  false;
  float D = ((double) b) * 4503599627370497.0;  //2^52 + 1
  D = (((double) b)+D )-D ;  // get the power of two closer to b
  if (D<b) D *=2;
  Fixed_B0 = D/Fixed_B24; // B24 = 2^24
  Fixed_B24 = D;
  Fixed_SNAP  *= Fixed_B0;
  Fixed_S_2_25*= Fixed_B0*Fixed_B0;
  Fixed_S_3_25*= Fixed_B0*Fixed_B0*Fixed_B0;
  Fixed_S_3_51*= Fixed_B0*Fixed_B0*Fixed_B0;
  Fixed_Or1   *= Fixed_B0*Fixed_B0*Fixed_B0;
  Fixed_Bx3   *= Fixed_B0*Fixed_B0*Fixed_B0;
  Fixed_Ic1   *= Fixed_B0*Fixed_B0*Fixed_B0*Fixed_B0;
  Fixed_Bx4   *= Fixed_B0*Fixed_B0*Fixed_B0*Fixed_B0;
  Fixed_Is1   *= 
    Fixed_B0*Fixed_B0*Fixed_B0*Fixed_B0*Fixed_B0;
  Fixed_Bx5   *= 
    Fixed_B0*Fixed_B0*Fixed_B0*Fixed_B0*Fixed_B0;
  return result;
};

// ======================================================================
//--------- ( round to the nearest legal fixed)
// ======================================================================

CGAL_Fixed_public void Fixed_round(float& f)
{
  CGAL_Fixed_rounding_precondition_msg(f);
  f = ( f+ Fixed_SNAP ) - Fixed_SNAP ;
}

// ======================================================================
//--------- splitting numbers for exact computation
// ======================================================================




inline void Fixed_split
(double a, double &a1, double &a0, const double &format)
  // split a in two numbers. a=a1+a0, a1 get the most significant digit
  // format specify the range of the input, and how to split
  // if format is S_i_j it means that a is of degree i
  // (in terms of original coordinates) the splitting bit is B0^i*2^j
{ a1 = a+format; a1-= format; a0 = a-a1; }


// ======================================================================
//--------- geometric predicates
// ======================================================================


CGAL_Fixed_public
int Fixed_cmp_dist(float x0, float y0, 
		   float x1, float y1, 
		   float x2, float y2)
{
  double X1 = (double)x1-(double)x0;
  double Y1 = (double)y1-(double)y0;
  double X2 = (double)x2-(double)x0;
  double Y2 = (double)y2-(double)y0;
  X1 = (X2*X2+Y2*Y2)-(X1*X1+Y1*Y1);         // exact
  return (X1>0) ? 1 : (X1==0) ? 0 : -1;
}

CGAL_Fixed_public
int Fixed_cmp_dist(float x0, float y0, float z0, 
		   float x1, float y1, float z1,
		   float x2, float y2, float z2)
{
  double X1 = (double)x1-(double)x0;
  double Y1 = (double)y1-(double)y0;
  double Z1 = (double)z1-(double)z0;
  double X2 = (double)x2-(double)x0;
  double Y2 = (double)y2-(double)y0;
  double Z2 = (double)z2-(double)z0;
  X1 = (X2*X2+Y2*Y2+Z2*Z2)-(X1*X1+Y1*Y1+Z1*Z1);         // exact
  return (X1>0) ? 1 : (X1==0) ? 0 : -1;
}

CGAL_Fixed_public
int Fixed_orientation(float x0, float y0, 
		      float x1, float y1, 
		      float x2, float y2)
{
    /*  points are assumed to be distincts */
    double det = ( (double)(x1) - (double)(x0))
                *( (double)(y2) - (double)(y0))
                -( (double)(x2) - (double)(x0))
                *( (double)(y1) - (double)(y0));
    /* stay inside double precision, thus it is exact */
    if (det>0) return 1;
    if (det)   return -1;
    /* the points are collinear */
    return 0;
}


CGAL_Fixed_public
int Fixed_orientation(  float x0, float y0, float z0,
			float x1, float y1, float z1,
			float x2, float y2, float z2,
			float x3, float y3, float z3)
{
  double X1=(double)(x1) -(double)(x0);
  double Y1=(double)(y1) -(double)(y0);
  double Z1=(double)(z1) -(double)(z0); 
  double X2=(double)(x2) -(double)(x0); 
  double Y2=(double)(y2) -(double)(y0); 
  double Z2=(double)(z2) -(double)(z0);
  double X3=(double)(x3) -(double)(x0);
  double Y3=(double)(y3) -(double)(y0);
  double Z3=(double)(z3) -(double)(z0);
  double M1=Y2*Z3-Y3*Z2;
  double M2=Y1*Z3-Y3*Z1;
  double M3=Y1*Z2-Y2*Z1;
  double det= X1*M1 - X2*M2 + X3*M3 ;
  if (det >   Fixed_Or1 ) return 1;
  if (det <  -Fixed_Or1 ) return -1;
  /* det is small */
  double error= 
    CGAL_Fixed_abs(X1*M1) + CGAL_Fixed_abs(X2*M2) + CGAL_Fixed_abs(X3*M3);
  error *= Fixed_Or2;
  if ( error < Fixed_Bx3) error = 0.0;
  if (det >   error) return 1;
  if (det <  -error) return -1;
  /* det is very small, exact computation must be done */
  if (error!=0.0) {  // otherwise, exact value 0 already certified
    Fixed_split(M1,M1,Z1,Fixed_S_2_25);
    Fixed_split(M2,M2,Z2,Fixed_S_2_25);
    Fixed_split(M3,M3,Z3,Fixed_S_2_25);
    det  =  X1*M1 - X2*M2 + X3*M3 ;
    det +=  X1*Z1 - X2*Z2 + X3*Z3 ;     /* less significant */
    if (det>0) return 1;
    if (det<0) return -1;
  }
  /* points are coplanar */
  return 0;
}


CGAL_Fixed_public
int Fixed_insphere (
      float x0, float y0,
      float x1, float y1,
      float x2, float y2,
      float x3, float y3)

  // relative position of p3 with respect to circle p0p1p2
  // if p0p1p2 is positively oriented,
  // positive side is the interior of the circle

{  
  double X1=(double)(x1)-(double)(x0);
  double Y1=(double)(y1)-(double)(y0);
  double X2=(double)(x2)-(double)(x0);
  double Y2=(double)(y2)-(double)(y0);
  double X3=(double)(x3)-(double)(x0);
  double Y3=(double)(y3)-(double)(y0);
  double R1=X1*X1+Y1*Y1;
  double R2=X2*X2+Y2*Y2;
  double R3=X3*X3+Y3*Y3;
  double M1=X2*Y3-X3*Y2;
  double M2=X1*Y3-X3*Y1;
  double M3=X1*Y2-X2*Y1;
  double det= R1*M1 -R2*M2 +R3*M3;
  if (det >   Fixed_Ic1) return -1;
  if (det <= -Fixed_Ic1) return 1;
  // static filter failed, error is small
  {
    double error= 
      CGAL_Fixed_abs(R1*M1) + CGAL_Fixed_abs(R2*M2) + CGAL_Fixed_abs(R3*M3) ;
    error *= Fixed_Ic2 ;
    if ( error < Fixed_Bx4 ) error = 0.0;
    if (det >   error) return -1;
    if (det <  -error) return 1;
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
      if (det>0) return -1;
      if (det<0) return 1;
    }
    //cocircular case
    if (! Fixed_incircle_perturb) return 0;

    // apply perturbation method
    // first order coefficient is obtained replacing x^2+y^2 by xy
    R1 = X1*Y1;
    R2 = X2*Y2;
    R3 = X3*Y3;
    det= R1*M1 -R2*M2 +R3*M3;
    if (det >   Fixed_Ic1) return -1;
    if (det <= -Fixed_Ic1) return 1;
    // static filter failed, error is small
    error= 
      (CGAL_Fixed_abs(R1*M1)+CGAL_Fixed_abs(R2*M2)+CGAL_Fixed_abs(R3*M3))
      * Fixed_Ic2 ;
    if ( error < Fixed_Bx4 ) error = 0.0;
    if (det >   error) return -1;
    if (det <  -error) return 1;
    if (error!=0) {
      // dynamic filter failed, error is very small
      Fixed_split(R1,R1,r1,Fixed_S_2_25);
      Fixed_split(R2,R2,r2,Fixed_S_2_25);
      Fixed_split(R3,R3,r3,Fixed_S_2_25);
      det  =  R1*M1 -R2*M2 +R3*M3   ;
      det += (R1*m1 + r1*M1) - (R2*m2 + r2*M2) + (R3*m3 + r3*M3) ;
      det += r1*m1 -r2*m2 +r3*m3  ;
      if (det>0) return -1;
      if (det<0) return 1;
    }
    // first order coefficient is null,
    // compute second order coefficient in the pertrubation method
    //  replace x^2+y^2 by y^2
    R1 = Y1*Y1;
    R2 = Y2*Y2;
    R3 = Y3*Y3;
    det= R1*M1 -R2*M2 +R3*M3;
    if (det >   Fixed_Ic1) return -1;
    if (det <= -Fixed_Ic1) return 1;
    // static filter failed, error is small
    error= 
      (CGAL_Fixed_abs(R1*M1)+CGAL_Fixed_abs(R2*M2)+CGAL_Fixed_abs(R3*M3))
      * Fixed_Ic2 ;
    if ( error < Fixed_Bx4 ) error = 0.0;
    if (det >   error) return -1;
    if (det <  -error) return 1;
    if (error!=0) {
      // dynamic filter failed, error is very small
      Fixed_split(R1,R1,r1,Fixed_S_2_25);
      Fixed_split(R2,R2,r2,Fixed_S_2_25);
      Fixed_split(R3,R3,r3,Fixed_S_2_25);
      det  =  R1*M1 -R2*M2 +R3*M3   ;
      det += (R1*m1 + r1*M1) - (R2*m2 + r2*M2) + (R3*m3 + r3*M3) ;
      det += r1*m1 -r2*m2 +r3*m3  ;
      if (det>0) return -1;
      if (det<0) return 1;
    }
    // points are colinear
    return 0;
  }
}



CGAL_Fixed_public
int Fixed_insphere (
      float x0, float y0, float z0,
      float x1, float y1, float z1,
      float x2, float y2, float z2,
      float x3, float y3, float z3,
      float x4, float y4, float z4)

  // relative position of p4 with respect to sphere p0p1p2p3
  // if p0p1p2p3 is positively oriented,
  // positive side is the interior of the sphere

{
  double det,error;
  // transform in 4x4 determinant
  // point p4 is inside sphere through oriented tetra p0p1p2p3
  // if the determinant is negative
  double X1=(double)(x1) -(double)(x0) ;
  double Y1=(double)(y1) -(double)(y0) ;
  double Z1=(double)(z1) -(double)(z0) ;
  double X2=(double)(x2) -(double)(x0) ;  //   | X1 X2 X3 X4 |
  double Y2=(double)(y2) -(double)(y0) ;  //   | Y1 Y2 Y3 Y4 |
  double Z2=(double)(z2) -(double)(z0) ;  //   | Z1 Z2 Z3 Z4 |
  double X3=(double)(x3) -(double)(x0) ;  //   | R1 R2 R3 R4 |
  double Y3=(double)(y3) -(double)(y0) ;
  double Z3=(double)(z3) -(double)(z0) ;
  double X4=(double)(x4) -(double)(x0) ;
  double Y4=(double)(y4) -(double)(y0) ;
  double Z4=(double)(z4) -(double)(z0) ;
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
  if (det >=  Fixed_Is1 ) return -1;
  if (det <= -Fixed_Is1 ) return 1;
  // static filter failed, error is small
  {
    double mm1 =CGAL_Fixed_abs(Z2*M34) + CGAL_Fixed_abs(Z3*M24) 
      + CGAL_Fixed_abs(Z4*M23) ;
    double mm2 =CGAL_Fixed_abs(Z1*M34) + CGAL_Fixed_abs(Z3*M14) 
      + CGAL_Fixed_abs(Z4*M13) ;
    double mm3 =CGAL_Fixed_abs(Z1*M24) + CGAL_Fixed_abs(Z2*M14) 
      + CGAL_Fixed_abs(Z4*M12) ;
    double mm4 =CGAL_Fixed_abs(Z1*M23) + CGAL_Fixed_abs(Z2*M13) 
      + CGAL_Fixed_abs(Z3*M12) ;
    error= CGAL_Fixed_abs(R1)*mm1 + CGAL_Fixed_abs(R2)*mm2 
           + CGAL_Fixed_abs(R3)*mm3 + CGAL_Fixed_abs(R4)*mm4;
    error *= Fixed_Is2;
    if ( error < Fixed_Bx5 ) error = 0.0;
    if (det >   error) return -1;
    if (det <  -error) return 1;
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
	if (det>0) return -1;
	if (det<0) return 1;
      }
      // points are cospherical
      
      if (! Fixed_insphere_perturb)return 0;

      // apply perturbation method
      // first order coefficient is obtained replacing x^2+y^2+z^2 by xy
      R1 = X1*Y1;
      R2 = X2*Y2;
      R3 = X3*Y3;
      R4 = X4*Y4;
      // compute determinant
      det= (R2 * MM2 - R1 * MM1) + (R4 * MM4 - R3 * MM3) ;
      if (det >=  Fixed_Is1 ) return -1;
      if (det <= -Fixed_Is1 ) return 1;
      // static filter failed, error is small
      error= CGAL_Fixed_abs(R1)*mm1 + CGAL_Fixed_abs(R2)*mm2 
	     + CGAL_Fixed_abs(R3)*mm3 + CGAL_Fixed_abs(R4)*mm4;
      error *= Fixed_Is2;
      if ( error < Fixed_Bx5 ) error = 0.0;
      if (det >   error) return -1;
      if (det <  -error) return 1;
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
	if (det>0) return -1;
	if (det<0) return 1;
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
      if (det >=  Fixed_Is1 ) return -1;
      if (det <= -Fixed_Is1 ) return 1;
      // static filter failed, error is small
      error= CGAL_Fixed_abs(R1)*mm1 + CGAL_Fixed_abs(R2)*mm2
	     + CGAL_Fixed_abs(R3)*mm3 + CGAL_Fixed_abs(R4)*mm4;
      error *= Fixed_Is2;
      if ( error < Fixed_Bx5 ) error = 0.0;
      if (det >   error) return -1;
      if (det <  -error) return 1;
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
	if (det>0) return -1;
	if (det<0) return 1;
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
      if (det >=  Fixed_Is1 ) return -1;
      if (det <= -Fixed_Is1 ) return 1;
      // static filter failed, error is small
      error= CGAL_Fixed_abs(R1)*mm1 + CGAL_Fixed_abs(R2)*mm2
	     + CGAL_Fixed_abs(R3)*mm3 + CGAL_Fixed_abs(R4)*mm4;
      error *= Fixed_Is2;
      if ( error < Fixed_Bx5 ) error = 0.0;
      if (det >   error) return -1;
      if (det <  -error) return 1;
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
	if (det>0) return -1;
	if (det<0) return 1;
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
      if (det >=  Fixed_Is1 ) return -1;
      if (det <= -Fixed_Is1 ) return 1;
      // static filter failed, error is small
      error= CGAL_Fixed_abs(R1)*mm1 + CGAL_Fixed_abs(R2)*mm2
	     + CGAL_Fixed_abs(R3)*mm3 + CGAL_Fixed_abs(R4)*mm4;
      error *= Fixed_Is2;
      if ( error < Fixed_Bx5 ) error = 0.0;
      if (det >   error) return -1;
      if (det <  -error) return 1;
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
	if (det>0) return -1;
	if (det<0) return 1;
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
      if (det >=  Fixed_Is1 ) return -1;
      if (det <= -Fixed_Is1 ) return 1;
      // static filter failed, error is small
      error= CGAL_Fixed_abs(R1)*mm1 + CGAL_Fixed_abs(R2)*mm2
	     + CGAL_Fixed_abs(R3)*mm3 + CGAL_Fixed_abs(R4)*mm4;
      error *= Fixed_Is2;
      if ( error < Fixed_Bx5 ) error = 0.0;
      if (det >   error) return -1;
      if (det <  -error) return 1;
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
	if (det>0) return -1;
	if (det<0) return 1;
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
      if (det >=  Fixed_Is1 ) return -1;
      if (det <= -Fixed_Is1 ) return 1;
      // static filter failed, error is small
      error= CGAL_Fixed_abs(R1)*mm1 + CGAL_Fixed_abs(R2)*mm2 
	     + CGAL_Fixed_abs(R3)*mm3 + CGAL_Fixed_abs(R4)*mm4;
      error *= Fixed_Is2;
      if ( error < Fixed_Bx5 ) error = 0.0;
      if (det >   error) return -1;
      if (det <  -error) return 1;
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
	if (det>0) return -1;
	if (det<0) return 1;
      }
      // 6th order coefficient is null,
      // points are coplanar
      return 0;
    }
  }
  return 0; // code never goes here, it is just to avoid compilation warning
}

#ifdef CGAL_FIXED_PRECISION_NT_H
// ======================================================================
//--------- static access and parametrization functions
// ======================================================================
float Fixed_precision_nt::unit_value()
{return Fixed_unit_value();}
float Fixed_precision_nt::upper_bound()
{return Fixed_upper_bound();}
void  Fixed_precision_nt::perturb_incircle() 
{Fixed_perturb_incircle();}
void  Fixed_precision_nt::unperturb_incircle() 
{Fixed_unperturb_incircle();}
bool  Fixed_precision_nt::is_perturbed_incircle() 
{return Fixed_is_perturbed_incircle();}
void  Fixed_precision_nt::perturb_insphere() 
{Fixed_perturb_insphere();}
void  Fixed_precision_nt::unperturb_insphere() 
{Fixed_unperturb_insphere();}
bool  Fixed_precision_nt::is_perturbed_insphere() 
{return Fixed_is_perturbed_insphere();}

extern void force_ieee_double_precision();

bool Fixed_precision_nt::init(float b) 
{ 
#ifdef __i386               // processor Intel 386
  force_ieee_double_precision();
#endif
  return Fixed_init(b);
}

// ======================================================================
//--------- geometric predicates
// ======================================================================

//template <>
Comparison_result
cmp_dist_to_pointC2(
  const Fixed_precision_nt x0, const Fixed_precision_nt y0,
  const Fixed_precision_nt x1, const Fixed_precision_nt y1,
  const Fixed_precision_nt x2, const Fixed_precision_nt y2)
{
  return (Comparison_result) sign(Fixed_cmp_dist(
			     x0.to_float(),  y0.to_float(), 
			     x1.to_float(),  y1.to_float(), 
			     x2.to_float(),  y2.to_float()));
}

//template <>
Comparison_result
cmp_dist_to_pointC3(
  const Fixed_precision_nt x0, const Fixed_precision_nt y0,
  const Fixed_precision_nt z0,
  const Fixed_precision_nt x1, const Fixed_precision_nt y1,
  const Fixed_precision_nt z1,
  const Fixed_precision_nt x2, const Fixed_precision_nt y2,
  const Fixed_precision_nt z2)
{
  return (Comparison_result) sign(Fixed_cmp_dist(
			     x0.to_float(),  y0.to_float(), z0.to_float(),
			     x1.to_float(),  y1.to_float(), z1.to_float(),
			     x2.to_float(),  y2.to_float(), z2.to_float()));
}

//template <>
Orientation orientationC2
(const Fixed_precision_nt x0, const Fixed_precision_nt y0, 
 const Fixed_precision_nt x1, const Fixed_precision_nt y1, 
 const Fixed_precision_nt x2, const Fixed_precision_nt y2)
{
  return (Orientation) sign (Fixed_orientation( 
                             x0.to_float(),  y0.to_float(), 
			     x1.to_float(),  y1.to_float(), 
			     x2.to_float(),  y2.to_float()));
}

//template <>
Orientation orientationC3
(   const Fixed_precision_nt x0, const Fixed_precision_nt y0, 
    const Fixed_precision_nt z0,
    const Fixed_precision_nt x1, const Fixed_precision_nt y1, 
    const Fixed_precision_nt z1,
    const Fixed_precision_nt x2, const Fixed_precision_nt y2,
    const Fixed_precision_nt z2,
    const Fixed_precision_nt x3, const Fixed_precision_nt y3,
    const Fixed_precision_nt z3)
{
  return (Orientation) sign (Fixed_orientation(  
                              x0.to_float(),  y0.to_float(),  z0.to_float(),
			      x1.to_float(),  y1.to_float(),  z1.to_float(),
			      x2.to_float(),  y2.to_float(),  z2.to_float(),
			      x3.to_float(),  y3.to_float(),  z3.to_float()));
}

//template <>
Oriented_side side_of_oriented_circleC2 (
      const Fixed_precision_nt x0, const Fixed_precision_nt y0,
      const Fixed_precision_nt x1, const Fixed_precision_nt y1,
      const Fixed_precision_nt x2, const Fixed_precision_nt y2,
      const Fixed_precision_nt x3, const Fixed_precision_nt y3)
  // relative position of p3 with respect to circle p0p1p2
  // if p0p1p2 is positively oriented,
  // positive side is the interior of the circle
{
  return (Oriented_side) sign (Fixed_insphere(  
                                     x0.to_float(),  y0.to_float(),
				     x1.to_float(),  y1.to_float(),
				     x2.to_float(),  y2.to_float(),
				     x3.to_float(),  y3.to_float()));
}

//template <>
Oriented_side side_of_oriented_sphereC3 
(     const Fixed_precision_nt x0, const Fixed_precision_nt y0, 
      const Fixed_precision_nt z0,
      const Fixed_precision_nt x1, const Fixed_precision_nt y1,
      const Fixed_precision_nt z1,
      const Fixed_precision_nt x2, const Fixed_precision_nt y2,
      const Fixed_precision_nt z2,
      const Fixed_precision_nt x3, const Fixed_precision_nt y3,
      const Fixed_precision_nt z3,
      const Fixed_precision_nt x4, const Fixed_precision_nt y4,
      const Fixed_precision_nt z4)
  // relative position of p4 with respect to sphere p0p1p2p3
  // if p0p1p2p3 is positively oriented,
  // positive side is the interior of the sphere
{
  return (Oriented_side) sign (Fixed_insphere(   
                            x0.to_float(),  y0.to_float(),  z0.to_float(),
			    x1.to_float(),  y1.to_float(),  z1.to_float(),
			    x2.to_float(),  y2.to_float(),  z2.to_float(),
			    x3.to_float(),  y3.to_float(),  z3.to_float(),
			    x4.to_float(),  y4.to_float(),  z4.to_float()));
}


// ======================================================================
//--------- constructors             CGAL only
// ======================================================================
Fixed_precision_nt::Fixed_precision_nt():_value(0.0){}

Fixed_precision_nt::Fixed_precision_nt
(const Fixed_precision_nt& f):_value(f._value){}

Fixed_precision_nt::Fixed_precision_nt(double f):_value(f)
{ Fixed_round(_value); }

Fixed_precision_nt::Fixed_precision_nt(int f):_value(f)
{ Fixed_round(_value); }



// ======================================================================
//--------- NT requirement           CGAL only
// ======================================================================


Fixed_precision_nt Fixed_precision_nt::operator=
(const Fixed_precision_nt& f)
{_value = f._value;return *this;}

// only numbers between -B24 and B24 are authorized
// except overflow only valid fixed are constructed
bool is_valid(Fixed_precision_nt f)
{
  if ((CGAL::to_double(f)>Fixed_B24)||(CGAL::to_double(f)<-Fixed_B24))
    return false;    // bigger than largest authorized value
  double v = CGAL::to_double(f) / Fixed_B0; // should be integer
  v -= abs(v); // should be 0
  if (v!=0) return false;
  return true;
}
bool is_finite(Fixed_precision_nt f)
{return is_valid(f);}

bool  operator==(Fixed_precision_nt a, Fixed_precision_nt b)
{   return (CGAL::to_double(a) == CGAL::to_double(b) );}
bool  operator!=(Fixed_precision_nt a, Fixed_precision_nt b)
{   return (CGAL::to_double(a) != CGAL::to_double(b) );}
bool  operator<(Fixed_precision_nt a, Fixed_precision_nt b)
{   return (CGAL::to_double(a) < CGAL::to_double(b) );}
bool  operator>(Fixed_precision_nt a, Fixed_precision_nt b)
{   return (CGAL::to_double(a) > CGAL::to_double(b) );}
bool  operator<=(Fixed_precision_nt a, Fixed_precision_nt b)
{   return (CGAL::to_double(a) <= CGAL::to_double(b) );}
bool  operator>=(Fixed_precision_nt a, Fixed_precision_nt b)
{   return (CGAL::to_double(a) >= CGAL::to_double(b) );}

Fixed_precision_nt  operator+
(Fixed_precision_nt a, Fixed_precision_nt b)
{   return Fixed_precision_nt(CGAL::to_double(a) + CGAL::to_double(b) );}
Fixed_precision_nt  operator-
(Fixed_precision_nt a, Fixed_precision_nt b)
{   return Fixed_precision_nt(CGAL::to_double(a) - CGAL::to_double(b) );}
Fixed_precision_nt  operator*
(Fixed_precision_nt a, Fixed_precision_nt b)
{   return Fixed_precision_nt(CGAL::to_double(a) * CGAL::to_double(b) );}
Fixed_precision_nt  operator-( Fixed_precision_nt b)
{   return Fixed_precision_nt( - CGAL::to_double(b) );}
Fixed_precision_nt Fixed_precision_nt::operator+=
(const Fixed_precision_nt& f)
{*this = *this+f;return *this;}
Fixed_precision_nt Fixed_precision_nt::operator-=
(const Fixed_precision_nt& f)
{*this = *this-f;return *this;}
Fixed_precision_nt Fixed_precision_nt::operator*=
(const Fixed_precision_nt& f)
{*this = *this*f;return *this;}
Fixed_precision_nt Fixed_precision_nt::operator/=
(const Fixed_precision_nt& f)
{(CGAL::to_double(f) == 0.0) ? *this = *this/f
                     : *this = Fixed_precision_nt::upper_bound();
 return *this;}
Fixed_precision_nt  operator/
(Fixed_precision_nt a, Fixed_precision_nt b)
{   return Fixed_precision_nt( 
      (CGAL::to_double(b)) ? CGAL::to_double(a) / CGAL::to_double(b) 
      : Fixed_precision_nt::upper_bound() );}

// ======================================================================
//--------- non official NT requirement IO       CGAL only
// ======================================================================

std::ostream &operator<<(std::ostream &os, Fixed_precision_nt a)
{ return os << CGAL::to_double(a); }
std::istream &operator>>(std::istream &is, Fixed_precision_nt a)
{ float f;  is >>f; a=Fixed_precision_nt(f); return is; }

#undef CGAL_Fixed_public
#undef CGAL_Fixed_abs
#undef Fixed_to_double
#undef CGAL_Fixed_rounding_precondition_msg


CGAL_END_NAMESPACE

#endif  //CGAL_FIXED_PRECISION_NT_H
