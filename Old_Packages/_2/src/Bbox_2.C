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
// release_date  : 2000, July 30
// 
// source        : Bbox_2.fw
// file          : Bbox_2.C
// package       : _2 (3.6)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 3.6
// revision_date : 30 Jul 2000 
// author(s)     : Andreas Fabri
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_BBOX_2_H
#include <CGAL/Bbox_2.h>
#endif // CGAL_BBOX_2_H

CGAL_BEGIN_NAMESPACE

Bbox_2::Bbox_2()
{ new ( static_cast< void*>(ptr)) Fourtuple<double>(); }

Bbox_2::Bbox_2(double x_min, double y_min,
               double x_max, double y_max)
{
  new ( static_cast< void*>(ptr)) Fourtuple<double>(x_min, y_min,
                                                    x_max, y_max);
}

bool Bbox_2::operator==(const Bbox_2 &b) const
{
  return    xmin() == b.xmin() && xmax() == b.xmax()
         && ymin() == b.ymin() && ymax() == b.ymax();
}

bool Bbox_2::operator!=(const Bbox_2 &b) const
{
  return ! (b == *this);
}

CGAL_END_NAMESPACE


