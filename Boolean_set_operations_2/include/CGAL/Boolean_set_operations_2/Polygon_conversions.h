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
// Author(s): Simon Giraudot  <simon.giraudot@geometryfactory.com>

#ifndef CGAL_BSO_POLYGON_CONVERSIONS_H
#define CGAL_BSO_POLYGON_CONVERSIONS_H

#include <boost/range/join.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <CGAL/license/Boolean_set_operations_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/General_polygon_2.h>
#include <CGAL/General_polygon_with_holes_2.h>
#include <CGAL/Boolean_set_operations_2/Gps_default_traits.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Single.h>
#include <CGAL/Iterator_range.h>

namespace CGAL {

template <typename InputIterator>
struct Is_Kernel_Polygon_2_iterator
{
  static constexpr bool value = false;
};

template <typename Kernel>
struct Is_Kernel_Polygon_2_iterator<Polygon_2<Kernel> >
{
  static constexpr bool value = true;
};

template <typename Kernel>
struct Is_Kernel_Polygon_2_iterator<Polygon_with_holes_2<Kernel> >
{
  static constexpr bool value = true;
};

template <typename InputIterator>
using Enable_if_Polygon_2_iterator =
  typename std::enable_if<Is_Kernel_Polygon_2_iterator<InputIterator>::value>::type;

template <typename InputIterator>
using Disable_if_Polygon_2_iterator =
  typename std::enable_if<!Is_Kernel_Polygon_2_iterator<InputIterator>::value>::type;

// Convert Polygon_2 to General_polygon_2<Polyline_traits>
template <typename ArrTraits, typename Kernel, typename Container>
General_polygon_2<ArrTraits>
convert_polygon(const Polygon_2<Kernel, Container>& polygon)
{
  return General_polygon_2<ArrTraits>
    (ArrTraits::make_curve_2
     (CGAL::make_range(polygon.vertices_begin(), polygon.vertices_end()), true));
}

// Convert Polygon_2 to General_polygon_2<Polyline_traits>
template <typename ArrTraits, typename Kernel, typename Container>
General_polygon_2<ArrTraits>
convert_polygon(const Polygon_2<Kernel, Container>& polygon,
                const ArrTraits& traits)
{
  auto ctr = traits.construct_curve_2_object();
  if (polygon.is_empty()) return General_polygon_2<ArrTraits>();
  return ctr(boost::range::join(CGAL::make_range(polygon.vertices_begin(),
                                                 polygon.vertices_end()),
                                CGAL::make_single(*polygon.vertices_begin())));
}

// Convert Polygon_with_holes_2 to General_polygon_with_holes_2<Polyline_traits>
template <typename ArrTraits, typename Kernel, typename Container>
General_polygon_with_holes_2<General_polygon_2<ArrTraits> >
convert_polygon(const Polygon_with_holes_2<Kernel, Container>& pwh) {
  General_polygon_with_holes_2<General_polygon_2<ArrTraits> >
    out(convert_polygon<ArrTraits>(pwh.outer_boundary()));
  for (const Polygon_2<Kernel, Container>& h : pwh.holes())
    out.add_hole(convert_polygon<ArrTraits>(h));
  return out;
}

// Convert Polygon_with_holes_2 to General_polygon_with_holes_2<Polyline_traits>
template <typename ArrTraits, typename Kernel, typename Container>
General_polygon_with_holes_2<General_polygon_2<ArrTraits> >
convert_polygon(const Polygon_with_holes_2<Kernel, Container>& pwh,
                const ArrTraits& traits) {
  typedef General_polygon_2<ArrTraits>          General_pgn_2;
  typedef Polygon_2<Kernel, Container>          Pgn_2;
  auto converter = [&](const Pgn_2& pgn)->General_pgn_2 {
    return convert_polygon(pgn, traits);
  };
  return General_polygon_with_holes_2<General_polygon_2<ArrTraits>>
    (convert_polygon(pwh.outer_boundary(), traits),
     boost::make_transform_iterator(pwh.holes().begin(), converter),
     boost::make_transform_iterator(pwh.holes().end(), converter));
}

// The following should go away... Besides
// 1. the mapping should be based on Container as well, and
// 2. Theoretically it could differ for Polygon_with_holes.
template <typename Kernel>
struct Gps_polyline_traits {
  typedef Polygon_2<Kernel>                             Pgn;
  typedef typename Gps_default_traits<Pgn>::Arr_traits  Segment_traits;
  typedef Arr_polyline_traits_2<Segment_traits>         Polyline_traits;
  typedef Gps_traits_2<Polyline_traits>                 Traits;
};

// Convert CGAL::Polygon_2 to General_polygon_2<Polyline_traits>
template <typename InputIterator, typename Kernel>
boost::transform_iterator
  <std::function<General_polygon_2
                 <typename Gps_polyline_traits<Kernel>::Polyline_traits>
                 (typename std::iterator_traits<InputIterator>::reference)>,
   InputIterator>
convert_polygon_iterator(InputIterator it,
                         typename Gps_polyline_traits<Kernel>::Traits& tr)
{
  using Polyline_traits = typename Gps_polyline_traits<Kernel>::Polyline_traits;
  using Input_type = typename std::iterator_traits<InputIterator>::reference;
  using Return_type = General_polygon_2<Polyline_traits>;
  using Function_type = std::function<Return_type(Input_type)>;

  Function_type func =
    std::bind(convert_polygon<Kernel>, std::placeholders::_1, std::ref(tr));

  return boost::transform_iterator<Function_type, InputIterator>(it, func);
}

// Polygon converter unary function
// Converts General_polygon_with_holes_2<Polyline_traits> to Polygon_with_holes_2
template <typename OutputIterator, typename Kernel>
struct Polygon_converter {
  using Polyline_traits = typename Gps_polyline_traits<Kernel>::Traits;
  using Input_type = typename Polyline_traits::Polygon_with_holes_2;
  using Input_polygon = typename Input_type::General_polygon_2;
  using Output_type = Polygon_with_holes_2<Kernel>;
  using Output_polygon = typename Output_type::Polygon_2;

  OutputIterator& output;
  Polygon_converter(OutputIterator& output) : output(output) {}

  void operator()(const Input_type& pwh) const {
    Output_polygon outer_boundary;
    convert_polygon(pwh.outer_boundary(), outer_boundary);

    std::vector<Output_polygon> holes;

    for (const Input_polygon& h : pwh.holes()) {
      holes.emplace_back();
      convert_polygon(h, holes.back());
    }

    *(output++) = Output_type(outer_boundary, holes.begin(), holes.end());
  }

private:
  void convert_polygon(const Input_polygon& input, Output_polygon& out) const {
    for (auto cit = input.curves_begin(); cit != input.curves_end(); ++cit) {
      auto end = cit->points_end();
      // Skip last point, which is a duplication of the first point
      --end;
      for (auto pit = cit->points_begin(); pit != end; ++pit)
        out.push_back(*pit);
    }
  }
};

// Function output iterator wrapping OutputIterator with conversion to
// OutputIterator
template <typename OutputIterator, typename Kernel>
struct Polygon_converter_output_iterator :
  boost::function_output_iterator<Polygon_converter<OutputIterator, Kernel> >
{
  using Converter = Polygon_converter<OutputIterator, Kernel>;
  using Base = boost::function_output_iterator<Converter>;

  OutputIterator& output;
  Polygon_converter_output_iterator(OutputIterator& output) :
    Base(output),
    output(output)
  {}

  operator OutputIterator() const { return output; }
};

// Convert General_polygon_2<Polyline_traits> to Polygon_2
template <typename Kernel, typename Container, typename ArrTraits>
Polygon_2<Kernel, Container>
convert_polygon_back(const General_polygon_2<ArrTraits>& gpgn) {
  Polygon_2<Kernel, Container> pgn;
  for (auto cit = gpgn.curves_begin(); cit != gpgn.curves_end(); ++cit) {
    auto end = cit->points_end();
    --end;      // skip last point, which is a duplication of the first point
    for (auto pit = cit->points_begin(); pit != end; ++pit) pgn.push_back(*pit);
  }
  return pgn;
}

// Convert General_polygon_with_holes_2<Polyline_traits> to Polygon_with_holes_2
template <typename Kernel, typename Container, typename ArrTraits>
Polygon_with_holes_2<Kernel, Container>
convert_polygon_back(const General_polygon_with_holes_2
                       <General_polygon_2<ArrTraits> >& gpwh)
{
  typedef Polygon_2<Kernel, Container>                  Pgn_2;
  typedef General_polygon_2<ArrTraits>                  General_pgn_2;
  auto converter = [](const General_pgn_2& gpgn)->Pgn_2 {
    return convert_polygon_back<Kernel, Container>(gpgn);
  };
  return Polygon_with_holes_2<Kernel, Container>
    (convert_polygon_back<Kernel, Container>(gpwh.outer_boundary()),
     boost::make_transform_iterator(gpwh.holes().begin(), converter),
     boost::make_transform_iterator(gpwh.holes().end(), converter));
}

// Converts General_polygon_with_holes_2<Polyline_traits> to
// Polygon_with_holes_2
template <typename Kernel, typename OutputIterator>
Polygon_converter_output_iterator<OutputIterator, Kernel>
convert_polygon_back(OutputIterator& output)
{ return Polygon_converter_output_iterator<OutputIterator, Kernel>(output); }

// Utility for checking if polygon remains the same after being
// converted and back
template <typename Kernel>
Polygon_2<Kernel> test_conversion(const Polygon_2<Kernel>& polygon)
{
  typedef Polygon_2<Kernel>                             Pgn;
  typedef typename Gps_default_traits<Pgn>::Arr_traits  Segment_traits;
  typedef Arr_polyline_traits_2<Segment_traits>         Polyline_traits;
  typedef Gps_traits_2<Polyline_traits>                 Traits;

  using Traits_polygon = typename Traits::Polygon_2;
  using Traits_polygon_with_holes = typename Traits::Polygon_with_holes_2;
  Polyline_traits traits;

  Traits_polygon polygon2 = convert_polygon(polygon, traits);
  Traits_polygon_with_holes polygon3(polygon2);

  Polygon_with_holes_2<Kernel> out;
  Oneset_iterator<Polygon_with_holes_2<Kernel> > iterator(out);
  auto converter = convert_polygon_back<Kernel>(iterator);
  *converter++ = polygon3;

  return out.outer_boundary();
}

}

#endif // CGAL_BSO_POLYGON_CONVERSIONS_H
