// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : include/CGAL/kernel_to_kernel.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 
#ifndef CGAL_KERNEL_TO_KERNEL_H
#define CGAL_KERNEL_TO_KERNEL_H

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <LEDA/rat_point.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class NumberType>
struct Cartesian_double_to_Homogeneous
{
  typedef Point_2< Homogeneous< NumberType> >    Point2;
  typedef Segment_2< Homogeneous< NumberType> >  Segment;

  Point2
  operator()(  const Point_2<Cartesian<double> >& p) const
  { return Point2( NumberType(p.x()), NumberType(p.y()) ); }

  Segment
  operator()(  const Segment_2<Cartesian<double> >& s) const
  {
    return Segment( Point2( NumberType(s.source().x()),
                            NumberType(s.source().y()) ),
                    Point2( NumberType(s.target().x()),
                            NumberType(s.target().y()) ) );
  }
};

#ifdef CGAL_USE_LEDA
struct Cartesian_double_to_H_double_int
{
  typedef Point_2< Homogeneous< double> >    Point2;
  typedef Segment_2< Homogeneous< double> >  Segment;

  Segment
  operator()(  const Segment_2< Cartesian< double> >& s) const
  {
    leda_rat_point rs =  leda_point(s.source().x(), s.source().y());
    leda_rat_point rt =  leda_point(s.target().x(), s.target().y());

    return Segment(
      Point2(::to_double(rs.X()),::to_double(rs.Y()),::to_double(rs.W())),
      Point2(::to_double(rt.X()),::to_double(rt.Y()),::to_double(rt.W())) );
  }
};

struct Cartesian_float_to_H_double_int
{
  typedef Point_2< Homogeneous< double> >    Point2;
  typedef Segment_2< Homogeneous< double> >  Segment;

  Segment
  operator()(  const Segment_2< Cartesian< float> >& s) const
  {
    leda_rat_point rs =  leda_point(s.source().x(), s.source().y());
    leda_rat_point rt =  leda_point(s.target().x(), s.target().y());

    return Segment(
      Point2(::to_double(rs.X()),::to_double(rs.Y()),::to_double(rs.W())),
      Point2(::to_double(rt.X()),::to_double(rt.Y()),::to_double(rt.W())) );
  }
};
#endif // CGAL_USE_LEDA

CGAL_END_NAMESPACE

#endif // CGAL_KERNEL_TO_KERNEL_H
