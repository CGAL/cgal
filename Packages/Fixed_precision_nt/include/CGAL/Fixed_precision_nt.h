// ======================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
// file          : include/CGAL/Fixed_precision_nt.h
// package       : 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Olivier Devillers <Olivivier.Devillers@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================


#ifndef CGAL_FIXED_PRECISION_NT_H
#define CGAL_FIXED_PRECISION_NT_H

#include <CGAL/basic.h>
#include <iostream>
#include <CGAL/tags.h>
#include <CGAL/enum.h>
#include <CGAL/assertions.h>
#include <CGAL/double.h>
#include <CGAL/number_utils.h>
#include <CGAL/Interval_base.h>

CGAL_BEGIN_NAMESPACE

// Fixed_precision_nt implement 24 bits integers using float
// The user has to initiate a multiplicator factor to be applied to
// original data.
// ======================================================================
//--------- static variables declaration and initialization
// ======================================================================

// static marks to decide if perturbation are used
static bool Fixed_incircle_perturb= false;
static bool Fixed_insphere_perturb= false;

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




class Fixed_precision_nt
{
private:
  float _value;          // value of the number
  
public:
  // constructors
  Fixed_precision_nt() {_value=0;}
  Fixed_precision_nt(const Fixed_precision_nt&f) {_value=f._value;}
  Fixed_precision_nt(double f) {_value=f; round(); }
  Fixed_precision_nt(int f) {_value=f; round(); }
  //access functions
  double to_double() const {return (double)_value;}
  float  to_float() const {return _value;}
  // NT requirements
  Fixed_precision_nt operator= (const Fixed_precision_nt&f)
                        {_value= f._value;  return *this;}
  Fixed_precision_nt operator+=(const Fixed_precision_nt&f)
                        {_value+=f._value;  return *this;}
  Fixed_precision_nt operator-=(const Fixed_precision_nt&f)
                        {_value-=f._value;  return *this;}
  Fixed_precision_nt operator*=(const Fixed_precision_nt&f)
                        {_value*=f._value;  return *this;}
  Fixed_precision_nt operator/=(const Fixed_precision_nt&f)
                        {_value/=f._value; return *this;}
  // access and parametrization of static members
  inline static bool init(float f);
  void round() {_value = ( _value+ Fixed_SNAP ) - Fixed_SNAP ;}
  static float unit_value() {return Fixed_B0;}
  static float upper_bound(){return Fixed_B24;}
  static void perturb_incircle(){Fixed_incircle_perturb=true;}
  static void unperturb_incircle(){Fixed_incircle_perturb=false;}
  static bool is_perturbed_incircle(){return Fixed_incircle_perturb;}
  static void perturb_insphere(){Fixed_insphere_perturb=true;}
  static void unperturb_insphere(){Fixed_insphere_perturb=false;}
  static bool is_perturbed_insphere(){return Fixed_insphere_perturb;}
};

inline double to_double(Fixed_precision_nt a){ return a.to_double(); }
inline bool is_finite(Fixed_precision_nt) { return true; }
inline bool is_valid(Fixed_precision_nt) { return true; }

inline
Interval_base
to_interval (Fixed_precision_nt a)
{
  return a.to_double();
}

inline Fixed_precision_nt  operator- (Fixed_precision_nt a)
{   return Fixed_precision_nt( - (a.to_double()) );}
inline Fixed_precision_nt  operator+
(Fixed_precision_nt a, Fixed_precision_nt b)
{   return Fixed_precision_nt(a.to_double() + b.to_double() );}
inline Fixed_precision_nt  operator-
(Fixed_precision_nt a, Fixed_precision_nt b)
{   return Fixed_precision_nt(a.to_double() - b.to_double() );}
inline Fixed_precision_nt  operator*
(Fixed_precision_nt a, Fixed_precision_nt b)
{   return Fixed_precision_nt(a.to_double() * b.to_double() );}
inline Fixed_precision_nt  operator/
(Fixed_precision_nt a, Fixed_precision_nt b)
{   return Fixed_precision_nt(a.to_double() / b.to_double() );}

inline bool  operator<
(Fixed_precision_nt a, Fixed_precision_nt b)
{   return (a.to_double() <  b.to_double() );}
inline bool  operator<=
(Fixed_precision_nt a, Fixed_precision_nt b)
{   return (a.to_double() <= b.to_double() );}
inline bool  operator>
(Fixed_precision_nt a, Fixed_precision_nt b)
{   return (a.to_double() >  b.to_double() );}
inline bool  operator>=
(Fixed_precision_nt a, Fixed_precision_nt b)
{   return (a.to_double() >= b.to_double() );}
inline bool  operator==
(Fixed_precision_nt a, Fixed_precision_nt b)
{   return (a.to_double() == b.to_double() );}
inline bool  operator!=
(Fixed_precision_nt a, Fixed_precision_nt b)
{   return (a.to_double() != b.to_double() );}


inline std::ostream &operator<<(std::ostream &os, const Fixed_precision_nt& a)
{ return os << a.to_double(); }
inline std::istream &operator>>(std::istream &is, Fixed_precision_nt& a)
{ float f;  is >>f; a=Fixed_precision_nt(f); return is; }

// ======================================================================
//--------- non official mysterious NT requirement
// ======================================================================

inline Number_tag number_type_tag(Fixed_precision_nt)
{ return Number_tag(); }

inline io_Operator io_tag(Fixed_precision_nt)
{ return io_Operator(); }

// ======================================================================
//--------- initialize the upper bound on fixed
// ======================================================================
inline
bool Fixed_precision_nt::init(float b)
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

inline Comparison_result
cmp_dist_to_pointC2(
  const Fixed_precision_nt& x0, const Fixed_precision_nt& y0,
  const Fixed_precision_nt& x1, const Fixed_precision_nt& y1,
  const Fixed_precision_nt& x2, const Fixed_precision_nt& y2)
{
  return CGAL_NTS compare(
	CGAL_NTS square(x0.to_double()-x1.to_double()) 
                + CGAL_NTS square(y0.to_double()-y1.to_double()),
        CGAL_NTS square(x0.to_double()-x2.to_double()) 
                + CGAL_NTS square(y0.to_double()-y2.to_double()));
}


inline Orientation orientationC2
(const Fixed_precision_nt& x0, const Fixed_precision_nt& y0, 
 const Fixed_precision_nt& x1, const Fixed_precision_nt& y1, 
 const Fixed_precision_nt& x2, const Fixed_precision_nt& y2)
{
    /*  points are assumed to be distincts */
    double det = ( x1.to_double() - x0.to_double())
                *( y2.to_double() - y0.to_double())
                -( x2.to_double() - x0.to_double())
                *( y1.to_double() - y0.to_double());
    /* stay inside double precision, thus it is exact */
    if (det>0) return LEFTTURN;
    if (det)   return RIGHTTURN;
    /* the points are collinear */
    return COLLINEAR;
}


Oriented_side side_of_oriented_circleC2 (
      const Fixed_precision_nt& x0, const Fixed_precision_nt& y0,
      const Fixed_precision_nt& x1, const Fixed_precision_nt& y1,
      const Fixed_precision_nt& x2, const Fixed_precision_nt& y2,
      const Fixed_precision_nt& x3, const Fixed_precision_nt& y3);
  // relative position of p3 with respect to circle p0p1p2
  // if p0p1p2 is positively oriented,
  // positive side is the interior of the circle

inline Comparison_result
cmp_dist_to_pointC3(
  const Fixed_precision_nt& x0, const Fixed_precision_nt& y0,
        const Fixed_precision_nt& z0,
  const Fixed_precision_nt& x1, const Fixed_precision_nt& y1,
        const Fixed_precision_nt& z1,
  const Fixed_precision_nt& x2, const Fixed_precision_nt& y2,
        const Fixed_precision_nt& z2)
{
  return CGAL_NTS compare(
	CGAL_NTS square(x0.to_double()-x1.to_double()) 
                + CGAL_NTS square(y0.to_double()-y1.to_double())
                + CGAL_NTS square(z0.to_double()-z1.to_double()),
        CGAL_NTS square(x0.to_double()-x2.to_double()) 
                + CGAL_NTS square(y0.to_double()-y2.to_double())
                + CGAL_NTS square(z0.to_double()-z2.to_double()));
}

Orientation orientationC3
(   const Fixed_precision_nt& x0, const Fixed_precision_nt& y0, 
    const Fixed_precision_nt& z0,
    const Fixed_precision_nt& x1, const Fixed_precision_nt& y1, 
    const Fixed_precision_nt& z1,
    const Fixed_precision_nt& x2, const Fixed_precision_nt& y2,
    const Fixed_precision_nt& z2,
    const Fixed_precision_nt& x3, const Fixed_precision_nt& y3,
    const Fixed_precision_nt& z3);

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
      const Fixed_precision_nt& z4);
  // relative position of p4 with respect to sphere p0p1p2p3
  // if p0p1p2p3 is positively oriented,
  // positive side is the interior of the sphere

CGAL_END_NAMESPACE

#endif  //CGAL_FIXED_PRECISION_NT_H
