// Copyright (c) 2016 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Efi Fogel <efifogel@gmail.com>

#ifndef CGAL_POLYCURVE_2_IOSTREAM_H
#define CGAL_POLYCURVE_2_IOSTREAM_H

#include <CGAL/license/Arrangement_on_surface_2.h>


#include <iostream>

#include <CGAL/basic.h>
#include <CGAL/Arr_geometry_traits/Polycurve_2.h>
#include <CGAL/Arr_segment_traits_2.h>

namespace CGAL {
namespace internal {

/*! Output operator for a polyline. */
template <typename OutputStream, typename SubcurveType_2, typename PointType_2>
void write_polyline(OutputStream& os,
                    const Polycurve_2<SubcurveType_2, PointType_2>& xcv)
{
  typedef SubcurveType_2                                Subcurve_type_2;
  typedef PointType_2                                   Point_type_2;
  typedef Polycurve_2<Subcurve_type_2, Point_type_2>    Curve_2;

  os << xcv.number_of_subcurves();       // export the number of points.

  // Write the points.
  typename Curve_2::Point_const_iterator it = xcv.points_begin();
  while (it != xcv.points_end()) os << " " << *it++;
}

/*! Output operator for a polyline of type CGAL::Arr_segment_2<Kernel_>. */
template <typename OutputStream, typename Kernel_, typename PointType_2>
void write_polycurve(OutputStream& os,
                     const Polycurve_2<CGAL::Arr_segment_2<Kernel_>, PointType_2>& xcv)
{ write_polyline(os, xcv); }

/*! Output operator for a polyline of type CGAL::Segment_2<Kernel_>. */
template <typename OutputStream, typename Kernel_, typename PointType_2>
void write_polycurve(OutputStream& os,
                     const Polycurve_2<CGAL::Segment_2<Kernel_>, PointType_2>& xcv)
{ write_polyline(os, xcv); }

/*! Output operator for a polycurve. */
template <typename OutputStream, typename SubcurveType_2, typename PointType_2>
void write_polycurve(OutputStream& os,
                     const Polycurve_2<SubcurveType_2, PointType_2>& xcv)
{
  typedef SubcurveType_2                                Subcurve_type_2;
  typedef PointType_2                                   Point_type_2;
  typedef Polycurve_2<Subcurve_type_2, Point_type_2>    Curve_2;

  os << xcv.number_of_subcurves();       // export the number of subcurves.

  // Write the subcurves.
  typename Curve_2::Subcurve_const_iterator it = xcv.subcurves_begin();
  while (it != xcv.subcurves_end()) os << " " << *it++;
}

/*! Output operator for a polycurve. */
template <typename OutputStream, typename SubcurveType_2, typename PointType_2>
OutputStream& operator<<(OutputStream& os,
                         const Polycurve_2<SubcurveType_2, PointType_2>& xcv)
{
  write_polycurve(os, xcv);
  return os;
}

/*! Input operator for a polyline. */
template <typename InputStream, typename SubcurveType_2, typename PointType_2>
void read_polyline(InputStream& is,
                   Polycurve_2<SubcurveType_2, PointType_2>& xcv)
{
  typedef SubcurveType_2                                Subcurve_type_2;
  typedef PointType_2                                   Point_type_2;
  typedef Polycurve_2<Subcurve_type_2, Point_type_2>    Curve_2;

  std::size_t num;              // read the number of points.
  is >> num;
  if (0 == num) return;

  Point_type_2 ps;
  is >> ps;
  if (1 == num) {
    xcv = Curve_2(Subcurve_type_2(ps, ps));
    return;
  }

  // Read the points.
  Point_type_2 pt;
  Point_type_2* ps_p = &ps;
  Point_type_2* pt_p = &pt;
  std::list<Subcurve_type_2> subcurves;
  for (std::size_t i = 1; i < num; ++i) {
    is >> *pt_p;
    subcurves.push_back(Subcurve_type_2(*ps_p, *pt_p));
    std::swap(ps_p, pt_p);
  }

  xcv = Curve_2(subcurves.begin(), subcurves.end());    // create the polycurve
}

/*! Input operator for a polyline of type CGAL::Arr_segment_2<Kernel_>. */
template <typename InputStream, typename Kernel_, typename PointType_2>
void read_polycurve(InputStream& is,
                    Polycurve_2<CGAL::Arr_segment_2<Kernel_>, PointType_2>& xcv)
{ read_polyline(is, xcv); }

/*! Input operator for a polyline of type CGAL::Segment_2<Kernel_>. */
template <typename InputStream, typename Kernel_, typename PointType_2>
void read_polycurve(InputStream& is,
                    Polycurve_2<CGAL::Segment_2<Kernel_>, PointType_2>& xcv)
{ read_polyline(is, xcv); }

/*! Input operator for a polycurve. */
template <typename InputStream, typename SubcurveType_2, typename PointType_2>
void read_polycurve(InputStream& is,
                    Polycurve_2<SubcurveType_2, PointType_2>& xcv)
{
  typedef SubcurveType_2                                Subcurve_type_2;
  typedef PointType_2                                   Point_type_2;
  typedef Polycurve_2<Subcurve_type_2, Point_type_2>    Curve_2;

  std::size_t num;                      // read the number of subcurves.
  is >> num;

  // Read the subcurves.
  std::list<Subcurve_type_2> subcurves;
  for (std::size_t i = 0; i < num; ++i) {
    Subcurve_type_2 subcurve;
    is >> subcurve;
    subcurves.push_back(subcurve);
  }

  xcv = Curve_2(subcurves.begin(), subcurves.end());    // create the polycurve
}

/*! Input operator for a polycurve. */
template <typename InputStream, typename SubcurveType_2, typename PointType_2>
InputStream& operator>>(InputStream& is,
                        Polycurve_2<SubcurveType_2, PointType_2>& xcv)
{
  read_polycurve(is, xcv);
  return is;
}

}
}

#endif
