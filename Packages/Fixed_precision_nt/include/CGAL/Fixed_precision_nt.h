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
// file          : include/CGAL/Fixed_precision_nt.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Fixed number type
// author(s)     : Olivier Devillers <Olivivier.Devillers@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_FIXED_PRECISION_NT_H
#define CGAL_FIXED_PRECISION_NT_H

#include <CGAL/config.h>
#include <iostream>
#ifndef CGAL_TAGS_H
#include <CGAL/tags.h>
#endif // CGAL_TAGS_H
#include <CGAL/enum.h>

CGAL_BEGIN_NAMESPACE

// Fixed_precision_nt implement 24 bits integers using float
// The user has to initiate a multiplicator factor to be applied to
// original data.

class Fixed_precision_nt
{
private:
  float _value;          // value of the number
  
public:
  // constructors
  Fixed_precision_nt();
  Fixed_precision_nt(const Fixed_precision_nt&);
  Fixed_precision_nt(double);
  Fixed_precision_nt(int);
  //access functions
  double to_double() const;
  float  to_float() const;
  // NT requirements
  Fixed_precision_nt operator= (const Fixed_precision_nt&);
  Fixed_precision_nt operator+=(const Fixed_precision_nt&);
  Fixed_precision_nt operator-=(const Fixed_precision_nt&);
  Fixed_precision_nt operator*=(const Fixed_precision_nt&);
  Fixed_precision_nt operator/=(const Fixed_precision_nt&);
  // access and parametrization of static members
  static bool init(float);
  static float unit_value();
  static float upper_bound();
  static void perturb_incircle();
  static void unperturb_incircle();
  static bool is_perturbed_incircle();
  static void perturb_insphere();
  static void unperturb_insphere();
  static bool is_perturbed_insphere();
};

// ======================================================================
//--------- access fonctions
// ======================================================================

inline double Fixed_precision_nt::to_double() const{return (double) _value;}

inline float  Fixed_precision_nt::to_float()  const{return _value;}

inline double to_double(const Fixed_precision_nt f) {return f.to_double();}

// ======================================================================
//--------- geometric predicates
// ======================================================================


//template <>
Comparison_result
cmp_dist_to_pointC2(
  const Fixed_precision_nt x0, const Fixed_precision_nt y0,
  const Fixed_precision_nt x1, const Fixed_precision_nt y1,
  const Fixed_precision_nt x2, const Fixed_precision_nt y2);

//template <>
Comparison_result
cmp_dist_to_pointC3(
  const Fixed_precision_nt x0, const Fixed_precision_nt y0,
  const Fixed_precision_nt z0,
  const Fixed_precision_nt x1, const Fixed_precision_nt y1,
  const Fixed_precision_nt z1,
  const Fixed_precision_nt x2, const Fixed_precision_nt y2,
  const Fixed_precision_nt z2);

//template <>
Orientation orientationC2(
  const Fixed_precision_nt x0, const Fixed_precision_nt y0,
  const Fixed_precision_nt x1, const Fixed_precision_nt y1,
  const Fixed_precision_nt x2, const Fixed_precision_nt y2);

//template <>
Oriented_side side_of_oriented_circleC2 (
  const Fixed_precision_nt x0, const Fixed_precision_nt y0,
  const Fixed_precision_nt x1, const Fixed_precision_nt y1,
  const Fixed_precision_nt x2, const Fixed_precision_nt y2,
  const Fixed_precision_nt x3, const Fixed_precision_nt y3);

//template <>
Orientation orientationC3(  
  const Fixed_precision_nt x0, const Fixed_precision_nt y0,
  const Fixed_precision_nt z0,
  const Fixed_precision_nt x1, const Fixed_precision_nt y1,
  const Fixed_precision_nt z1,
  const Fixed_precision_nt x2, const Fixed_precision_nt y2,
  const Fixed_precision_nt z2,
  const Fixed_precision_nt x3, const Fixed_precision_nt y3,
  const Fixed_precision_nt z3);

//template <>
Oriented_side side_of_oriented_sphereC3 (
  const Fixed_precision_nt x0, const Fixed_precision_nt y0,
  const Fixed_precision_nt z0,
  const Fixed_precision_nt x1, const Fixed_precision_nt y1,
  const Fixed_precision_nt z1,
  const Fixed_precision_nt x2, const Fixed_precision_nt y2,
  const Fixed_precision_nt z2,
  const Fixed_precision_nt x3, const Fixed_precision_nt y3,
  const Fixed_precision_nt z3,
  const Fixed_precision_nt x4, const Fixed_precision_nt y4,
  const Fixed_precision_nt z4);



// ======================================================================
//--------- NT requirement
// ======================================================================

// only numbers between -B24 and B24 are authorized
// except overflow only valid fixed are constructed
bool is_valid(Fixed_precision_nt);
bool is_finite(Fixed_precision_nt);

bool  operator==(Fixed_precision_nt a, Fixed_precision_nt b);
bool  operator!=(Fixed_precision_nt a, Fixed_precision_nt b);
bool  operator<(Fixed_precision_nt a, Fixed_precision_nt b);
bool  operator>(Fixed_precision_nt a, Fixed_precision_nt b);
bool  operator<=(Fixed_precision_nt a, Fixed_precision_nt b);
bool  operator>=(Fixed_precision_nt a, Fixed_precision_nt b);

Fixed_precision_nt  
operator+(Fixed_precision_nt a, Fixed_precision_nt b);
Fixed_precision_nt  
operator-(Fixed_precision_nt a, Fixed_precision_nt b);
Fixed_precision_nt  
operator*(Fixed_precision_nt a, Fixed_precision_nt b);
Fixed_precision_nt  
operator-( Fixed_precision_nt b);
Fixed_precision_nt  
operator/(Fixed_precision_nt a, Fixed_precision_nt b);

// ======================================================================
//--------- non official NT requirement IO
// ======================================================================

std::ostream &operator<<(std::ostream &os, Fixed_precision_nt a);
std::istream &operator>>(std::istream &is, Fixed_precision_nt a);

// ======================================================================
//--------- non official mysterious NT requirement
// ======================================================================

inline Number_tag number_type_tag(Fixed_precision_nt)
{ return Number_tag(); }

inline io_Operator io_tag(Fixed_precision_nt)
{ return io_Operator(); }

CGAL_END_NAMESPACE

#ifdef CGAL_INTERVAL_ARITHMETIC_H
#include <CGAL/Interval_arithmetic/IA_Fixed.h>
#endif // CGAL_INTERVAL_ARITHMETIC_H


#endif // CGAL_FIXED_PRECISION_NT_H.

