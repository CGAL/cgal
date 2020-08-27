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
// Author(s)     : Simon Giraudot <simon.giraudot@geometryfactory.com>

#ifndef CGAL_BSO_POLYGON_CONVERSIONS_H
#define CGAL_BSO_POLYGON_CONVERSIONS_H

#include <CGAL/license/Boolean_set_operations_2.h>

namespace CGAL
{

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
using Enable_if_Polygon_2_iterator
= typename std::enable_if<Is_Kernel_Polygon_2_iterator<InputIterator>::value>::type;

template <typename InputIterator>
using Disable_if_Polygon_2_iterator
= typename std::enable_if<!Is_Kernel_Polygon_2_iterator<InputIterator>::value>::type;

// Default Polygon conversion (identity)
template <typename Traits>
const typename Traits::Polygon_2&
convert_polygon (const typename Traits::Polygon_2& polygon,
                 Traits&)
{
  return polygon;
}

// Default Polygon_with_holes conversion (identity)
template <typename Traits>
const typename Traits::Polygon_with_holes_2&
convert_polygon (const typename Traits::Polygon_with_holes_2& polygon,
                 Traits&)
{
  return polygon;
}

// Convert CGAL::Polygon_2 to General_polygon_2<Polyline_traits>
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

// Convert CGAL::Polygon_with_holes_2 to General_polygon_with_holes_2<Polyline_traits>
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

// Convert CGAL::Polygon_2 to General_polygon_2<Polyline_traits>
template <typename InputIterator, typename Kernel>
boost::transform_iterator
<std::function
 <General_polygon_2<typename Gps_polyline_traits<Kernel>::Base>
  (typename std::iterator_traits<InputIterator>::reference)>,
 InputIterator>
convert_polygon_iterator (InputIterator it,
                          typename Gps_polyline_traits<Kernel>::Traits& tr)
{
  using Polyline_traits = typename Gps_polyline_traits<Kernel>::Base;
  using Input_type = typename std::iterator_traits<InputIterator>::reference;
  using Return_type = General_polygon_2<Polyline_traits>;
  using Function_type = std::function<Return_type(Input_type)>;

  Function_type func = std::bind (convert_polygon<Kernel>, std::placeholders::_1, std::ref(tr));

  return boost::transform_iterator<Function_type, InputIterator>
    (it, func);
}


// Default Identity back conversion for Polygon_2
template <typename OutputIterator, typename Traits>
OutputIterator convert_polygon_back (OutputIterator output,
                                     const typename Traits::Polygon_2&, // used for indirection
                                     Traits&)
{
  return output;
}

// Default Identity back conversion for Polygon_with_holes_2
template <typename OutputIterator, typename Traits>
OutputIterator convert_polygon_back (OutputIterator output,
                                     const typename Traits::Polygon_with_holes_2&,  // used for indirection
                                     Traits&)
{
  return output;
}

// Polygon converter unary function
// Converts General_polygon_with_holes_2<Polyline_traits> to CGAL::Polygon_with_holes_2
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
    {
      // Skip last point which is the repetition of the first point
      typename Input_polygon::X_monotone_curve_2::Point_const_iterator
        end = it->points_end();
      -- end;
      for (typename Input_polygon::X_monotone_curve_2::Point_const_iterator
             it2 = it->points_begin(); it2 != end; ++ it2)
        out.push_back (*it2);
    }
  }
};


// Function output iterator wrapping OutputIterator with conversion to OutputIterator
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

// Converts General_polygon_with_holes_2<Polyline_traits> to CGAL::Polygon_with_holes_2
// (indirection with Polygon_2)
template <typename OutputIterator, typename Kernel>
Polygon_converter_output_iterator<OutputIterator, Kernel>
convert_polygon_back (OutputIterator& output,
                      const Polygon_2<Kernel>&,  // used for indirection
                      typename Gps_polyline_traits<Kernel>::Traits&)
{
  return Polygon_converter_output_iterator<OutputIterator, Kernel>(output);
}

// Converts General_polygon_with_holes_2<Polyline_traits> to CGAL::Polygon_with_holes_2
// (indirection with Polygon_with_holes_2)
template <typename OutputIterator, typename Kernel>
Polygon_converter_output_iterator<OutputIterator, Kernel>
convert_polygon_back (OutputIterator& output,
                      const Polygon_with_holes_2<Kernel>&, // used for indirection
                      typename Gps_polyline_traits<Kernel>::Traits&)
{
  return Polygon_converter_output_iterator<OutputIterator, Kernel>(output);
}

// Utility for checking if polygon remains the same after being
// converted and back
template <typename Kernel>
Polygon_2<Kernel>
test_conversion (const Polygon_2<Kernel>& polygon)
{
  using Polyline_traits = typename Gps_polyline_traits<Kernel>::Traits;
  using Traits_polygon = typename Polyline_traits::Polygon_with_holes_2;
  Polyline_traits traits;

  auto polygon2 = convert_polygon(polygon, traits);
  Traits_polygon polygon3 (polygon2);

  Polygon_with_holes_2<Kernel> out;

  Oneset_iterator<Polygon_with_holes_2<Kernel> > iterator(out);
  auto converter = convert_polygon_back (iterator, polygon, traits);
  *(converter ++) = polygon3;

  return out.outer_boundary();
}

}

#endif // CGAL_BSO_POLYGON_CONVERSIONS_H
