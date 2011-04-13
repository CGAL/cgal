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
// file          : Bbox_2.C
// package       : _2
// revision      : $Revision$
// revision_date : $Date$
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


