// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// source        : test_kernel_programs.fw
// file          : predicates_on_points_rat_leda_2.h
// revision      : 3.8
// revision_date : 08 Oct 2000 
// author(s)     : Stefan Schirra
//
// maintainer    : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de> 
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_PREDICATES_ON_POINTS_RAT_LEDA_2_H
#define CGAL_PREDICATES_ON_POINTS_RAT_LEDA_2_H

#include <CGAL/rat_leda.h>
#include <CGAL/predicates_on_points_2.h>

CGAL_BEGIN_NAMESPACE


template <>
inline
Orientation
orientation<use_rat_leda_kernel>(const Point_2<use_rat_leda_kernel>& p,
                                 const Point_2<use_rat_leda_kernel>& q,
                                 const Point_2<use_rat_leda_kernel>& r)
{
/*
  typedef use_rat_leda_kernel::Point_2_base  RPoint_2;
  return static_cast<Orientation>( ::orientation((const RPoint_2&)p,
                                                 (const RPoint_2&)q,
                                                 (const RPoint_2&)r) );
*/
  return static_cast<Orientation>(
      ::orientation((const use_rat_leda_kernel::Point_2_base&)p,
                    (const use_rat_leda_kernel::Point_2_base&)q,
                    (const use_rat_leda_kernel::Point_2_base&)r) );
}

// #ifdef CGAL_CFG_NO_KOENIG_LOOKUP
template <>
inline
bool
collinear<use_rat_leda_kernel>(const Point_2<use_rat_leda_kernel>& p,
                               const Point_2<use_rat_leda_kernel>& q,
                               const Point_2<use_rat_leda_kernel>& r)
{
/*
  typedef use_rat_leda_kernel::Point_2_base  RPoint_2;
  return ::collinear((const RPoint_2&)p,
                     (const RPoint_2&)q,
                     (const RPoint_2&)r );
*/
  return
      ::collinear((const use_rat_leda_kernel::Point_2_base&)p,
                  (const use_rat_leda_kernel::Point_2_base&)q,
                  (const use_rat_leda_kernel::Point_2_base&)r);
}
// #endif // CGAL_CFG_NO_KOENIG_LOOKUP

CGAL_END_NAMESPACE

#endif // CGAL_PREDICATES_ON_POINTS_RAT_LEDA_2_H
