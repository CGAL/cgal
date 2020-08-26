// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein        <wein@post.tau.ac.il>
//                 Efi Fogel       <efif@post.tau.ac.il>

#ifndef CGAL_BSO_INTERNAL_FUNCTIONS_H
#define CGAL_BSO_INTERNAL_FUNCTIONS_H

#include <CGAL/license/Boolean_set_operations_2.h>


#include <CGAL/Boolean_set_operations_2/Gps_default_traits.h>
#include <iterator>

namespace CGAL {

template <typename Traits>
const typename Traits::Polygon_2&
convert_polygon (const typename Traits::Polygon_2& polygon,
                 Traits&)
{
  return polygon;
}

template <typename Traits>
const typename Traits::Polygon_with_holes_2&
convert_polygon (const typename Traits::Polygon_with_holes_2& polygon,
                 Traits&)
{
  return polygon;
}

template <typename Kernel>
General_polygon_2<typename Gps_polyline_traits<Kernel>::Base>
convert_polygon (const Polygon_2<Kernel>& polygon,
                 typename Gps_polyline_traits<Kernel>::Traits&)
{
  using Polyline_traits = typename Gps_polyline_traits<Kernel>::Base;
  return General_polygon_2<Polyline_traits>
    (Polyline_traits::make_curve_2
     (CGAL::make_range (polygon.vertices_begin(),
                        polygon.vertices_end()),
      true));
}

template <typename Kernel>
General_polygon_with_holes_2<General_polygon_2
                             <typename Gps_polyline_traits<Kernel>::Base> >
convert_polygon (const Polygon_with_holes_2<Kernel>& pwh,
                 typename Gps_polyline_traits<Kernel>::Traits&)
{
  using Polyline_traits = typename Gps_polyline_traits<Kernel>::Base;

  General_polygon_with_holes_2<General_polygon_2<Polyline_traits> > out
    (General_polygon_2<Polyline_traits>
     (Polyline_traits::make_curve_2
      (CGAL::make_range (pwh.outer_boundary().vertices_begin(),
                         pwh.outer_boundary().vertices_end()),
       true)));

  for (const Polygon_2<Kernel>& h : pwh.holes())
    out.add_hole
      (General_polygon_2<Polyline_traits>
       (Polyline_traits::make_curve_2
        (CGAL::make_range (h.vertices_begin(),
                           h.vertices_end()),
         true)));

  return out;
}

template <typename OutputIterator, typename Traits>
OutputIterator convert_polygon_back (OutputIterator output,
                                     const typename Traits::Polygon_2&,
                                     Traits&)
{
  return output;
}

template <typename OutputIterator, typename Traits>
OutputIterator convert_polygon_back (OutputIterator output,
                                     const typename Traits::Polygon_with_holes_2&,
                                     Traits&)
{
  return output;
}

template <typename OutputIterator, typename Kernel>
struct Polygon_converter
{
  using Polyline_traits = typename Gps_polyline_traits<Kernel>::Traits;

  using Input_type = typename Polyline_traits::Polygon_with_holes_2;
  using Input_polygon = typename Input_type::General_polygon_2;

  using Output_type = Polygon_with_holes_2<Kernel>;
  using Output_polygon = typename Output_type::Polygon_2;

  OutputIterator& output;
  Polygon_converter (OutputIterator& output) : output (output) { }

  void operator() (const Input_type& pwh) const
  {
    Output_polygon outer_boundary;
    convert_polygon (pwh.outer_boundary(), outer_boundary);

    std::vector<Output_polygon> holes;

    for (const Input_polygon& h : pwh.holes())
    {
      holes.emplace_back ();
      convert_polygon (h, holes.back());
    }

    *(output ++) = Output_type (outer_boundary, holes.begin(), holes.end());
  }

private:

  void convert_polygon (const Input_polygon& input,
                        Output_polygon& out) const
  {
    for (typename Input_polygon::Curve_const_iterator
           it = input.curves_begin(); it != input.curves_end(); ++ it)
      for (typename Input_polygon::X_monotone_curve_2::Point_const_iterator
             it2 = it->points_begin(); it2 != it->points_end(); ++ it2)
        out.push_back (*it2);
  }
};

template <typename OutputIterator, typename Kernel>
struct Polygon_converter_output_iterator
  : boost::function_output_iterator<Polygon_converter<OutputIterator, Kernel> >
{
  using Base = boost::function_output_iterator<Polygon_converter<OutputIterator, Kernel> >;

  OutputIterator& output;

  Polygon_converter_output_iterator (OutputIterator& output)
    : Base (output), output (output)
  {

  }

  operator OutputIterator() const { return output; }
};

template <typename OutputIterator, typename Kernel>
Polygon_converter_output_iterator<OutputIterator, Kernel>
convert_polygon_back (OutputIterator& output,
                      const Polygon_2<Kernel>&,
                      typename Gps_polyline_traits<Kernel>::Traits&)
{
  return Polygon_converter_output_iterator<OutputIterator, Kernel>(output);
}

template <typename OutputIterator, typename Kernel>
Polygon_converter_output_iterator<OutputIterator, Kernel>
convert_polygon_back (OutputIterator& output,
                      const Polygon_with_holes_2<Kernel>&,
                      typename Gps_polyline_traits<Kernel>::Traits&)
{
  return Polygon_converter_output_iterator<OutputIterator, Kernel>(output);
}

/// \name _do_intersect() functions.
//@{

template <class Pgn1, class Pgn2, class Traits>
inline bool _do_intersect(const Pgn1& pgn1, const Pgn2& pgn2, Traits& tr)
{
  General_polygon_set_2<Traits> gps(tr);
  gps.insert(convert_polygon(pgn1, tr));
  return (gps.do_intersect(convert_polygon(pgn2, tr)));
}

template <class Pgn1, class Pgn2>
inline bool _do_intersect(const Pgn1& pgn1, const Pgn2& pgn2)
{
  typename Gps_default_traits<Pgn1>::Traits    tr;
  return _do_intersect(pgn1, pgn2, tr);
}

//@}
/// \name _oriented_side() functions.
//@{

template <class Obj, class Pgn, class Traits>
inline
Oriented_side _oriented_side(const Obj& obj, const Pgn& pgn, Traits& tr)
{
  General_polygon_set_2<Traits> gps(tr);
  gps.insert(convert_polygon(pgn, tr));
  return (gps.oriented_side(obj));
}

template <class Obj, class Pgn>
inline Oriented_side _oriented_side(const Obj& obj, const Pgn& pgn)
{
  typename Gps_default_traits<Pgn>::Traits    tr;
  return _oriented_side(obj, pgn, tr);
}

//@}
/// \name _intersection() functions.
//@{

template <class Pgn1, class Pgn2, class OutputIterator, class Traits>
inline OutputIterator _intersection(const Pgn1& pgn1, const Pgn2& pgn2,
                                    OutputIterator out, Traits& tr)
{
  General_polygon_set_2<Traits> gps(tr);
  gps.insert(convert_polygon(pgn1, tr));
  gps.intersection(convert_polygon(pgn2, tr));
  return (gps.polygons_with_holes(convert_polygon_back(out, pgn1, tr)));
}

template <class Pgn1, class Pgn2, class OutputIterator>
inline OutputIterator _intersection(const Pgn1& pgn1, const Pgn2& pgn2,
                                    OutputIterator out)
{
  typename Gps_default_traits<Pgn1>::Traits    tr;
  return (_intersection(pgn1, pgn2, out, tr));
}

//@}
/// \name _join() functions.
//@{

template <class Traits>
inline bool _is_empty (const typename Traits:: Polygon_2& pgn, Traits& tr)
{
  typedef typename Traits::Curve_const_iterator Curve_const_iterator;
  const std::pair<Curve_const_iterator, Curve_const_iterator>& itr_pair =
    tr.construct_curves_2_object()(pgn);
  return (itr_pair.first == itr_pair.second);
}

template <class Traits>
inline bool _is_empty (const typename Traits::Polygon_with_holes_2&, Traits&)
{
  return false;
}

template <typename Kernel>
inline bool _is_empty (const Polygon_2<Kernel>& polygon,
                       typename Gps_polyline_traits<Kernel>::Traits&)
{
  return (polygon.size() == 0);
}

template <typename Kernel>
inline bool _is_empty (const Polygon_with_holes_2<Kernel>& pwh,
                       typename Gps_polyline_traits<Kernel>::Traits&)
{
  return (pwh.outer_boundary().size() == 0);
}

template <class Pgn1, class Pgn2, class Pwh, class Traits>
inline bool _join(const Pgn1& pgn1, const Pgn2& pgn2,
                  Pwh& res, Traits& tr)
{
  if (_is_empty(pgn1, tr) || _is_empty(pgn2, tr))
    return false;

  General_polygon_set_2<Traits> gps(tr);
  gps.insert(convert_polygon(pgn1, tr));
  gps.join(convert_polygon(pgn2, tr));
  if (gps.number_of_polygons_with_holes() == 1)
  {
    Oneset_iterator<Pwh> oi (res);
    gps.polygons_with_holes(convert_polygon_back(oi, pgn1, tr));
    return true;
  }

  // the polygon doesnt intersect, the original pgn1, pgn2 contain the union
  return false;
}

template <class Pgn1, class Pgn2, class Pwh>
inline bool _join(const Pgn1& pgn1, const Pgn2& pgn2, Pwh& res)
{
  typename Gps_default_traits<Pgn1>::Traits  tr;
  return _join(pgn1, pgn2, res, tr);
}

//@}
/// \name _difference() functions.
//@{

template <class Pgn1, class Pgn2, class OutputIterator, class Traits>
inline OutputIterator _difference(const Pgn1& pgn1, const Pgn2& pgn2,
                                  OutputIterator oi, Traits& tr)
{
  General_polygon_set_2<Traits> gps(tr);
  gps.insert(convert_polygon(pgn1, tr));
  gps.difference(convert_polygon(pgn2, tr));
  return gps.polygons_with_holes(convert_polygon_back(oi, pgn1, tr));
}

template <class Pgn1, class Pgn2, class OutputIterator>
inline OutputIterator _difference(const Pgn1& pgn1, const Pgn2& pgn2,
                                  OutputIterator oi)
{
  typename Gps_default_traits<Pgn1>::Traits  tr;
  return _difference(pgn1, pgn2, oi, tr);
}

//@}
/// \name _symmetric_difference() functions.
//@{

template <class Pgn1, class Pgn2, class OutputIterator, class Traits>
inline OutputIterator _symmetric_difference(const Pgn1& pgn1, const Pgn2& pgn2,
                                            OutputIterator oi, Traits& tr)
{
  General_polygon_set_2<Traits> gps(tr);
  gps.insert(convert_polygon(pgn1, tr));
  gps.symmetric_difference(convert_polygon(pgn2, tr));
  return gps.polygons_with_holes(convert_polygon_back(oi, pgn1, tr));
}

template <class Pgn1, class Pgn2, class OutputIterator>
inline OutputIterator _symmetric_difference(const Pgn1& pgn1, const Pgn2& pgn2,
                                            OutputIterator oi)
{
  typename Gps_default_traits<Pgn1>::Traits    tr;
  return _symmetric_difference(pgn1, pgn2, oi, tr);
}

//@}
/// \name _complement() functions.
//@{

template <class Pgn, class Pwh, class Traits>
void _complement(const Pgn& pgn, Pwh& res,
                 Traits& tr)
{
  General_polygon_set_2<Traits> gps(tr);
  gps.insert(convert_polygon(pgn, tr));
  gps.complement();
  Oneset_iterator<Pwh> oi(res);
  gps.polygons_with_holes(convert_polygon_back(oi, pgn, tr));
}

template <class Pgn, class Pwh>
void _complement(const Pgn& pgn, Pwh& res)
{
  typename Gps_default_traits<Pgn>::Traits    tr;
  _complement(pgn, res, tr);
}

//@}

} //namespace CGAL

#endif
