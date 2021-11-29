// Copyright (c) 2020 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sébastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_SIMPLIFY_POLYLINE_H
#define CGAL_POLYGON_MESH_PROCESSING_SIMPLIFY_POLYLINE_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <type_traits>


namespace CGAL {
namespace Polygon_mesh_processing {
namespace experimental {

enum Polyline_simplification_algorithms { DOUGLAS_PEUCKER, ITERATIVE };

template <typename PointRangeIn, typename PointRangeOut,
          typename NamedParametersIn, typename NamedParametersOut>
void simplify_polyline(const PointRangeIn& input,
                             PointRangeOut& output,
                       const double max_squared_frechet_distance,
                       const NamedParametersIn& np_in,
                       const NamedParametersOut& np_out)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  static_assert(std::is_same<
    typename std::iterator_traits<typename PointRangeIn::const_iterator>::value_type,
    typename std::iterator_traits<typename PointRangeOut::const_iterator>::value_type >::value, "");

  typedef typename GetPointMap<PointRangeIn, NamedParametersIn>::type Point_map_in;
  typedef typename GetPointMap<PointRangeOut, NamedParametersOut>::type Point_map_out;
  typedef typename Point_set_processing_3::GetK<PointRangeIn, NamedParametersIn>::Kernel Kernel;

  Point_map_in in_pm = choose_parameter<Point_map_in>(get_parameter(np_in, internal_np::point_map));
  Point_map_out out_pm = choose_parameter<Point_map_out>(get_parameter(np_out, internal_np::point_map));

  const Polyline_simplification_algorithms algorithm =
    choose_parameter(get_parameter(np_in, internal_np::algorithm), DOUGLAS_PEUCKER);

  switch(algorithm)
  {
    case ITERATIVE:
    {
      const bool is_closed = input.front()==input.back();

      std::size_t nb_points = is_closed ? input.size()-1 : input.size();

      // skip points in the input range that do not contains any information
      if (nb_points<=2)
      {
        output.reserve(input.size());
        for (const auto& p : input)
        {
          output.push_back(p);
          put(out_pm, output.back(), get(in_pm, p));
        }
        return;
      }

      auto is_valid_approx = [&input, &in_pm, max_squared_frechet_distance](
        std::size_t b, std::size_t e,
        const typename Kernel::Line_3& line)
      {
        typename Kernel::Compare_squared_distance_3 compare_squared_distance;
        for (std::size_t i=b+1; i<e; ++i)
        {

          if (compare_squared_distance(get(in_pm, input[i]), line, max_squared_frechet_distance) == LARGER)
            return false;
        }
        return true;
      };

      std::size_t bi=0;
      while(bi!=nb_points)
      {
        std::size_t ei=bi+2;

        output.push_back(input[bi]);
        put(out_pm, output.back(), get(in_pm, input[bi]));

        while(ei<nb_points)
        {
          typename Kernel::Line_3 sl(get(in_pm,input[bi]), get(in_pm,input[ei]));

          if (is_valid_approx(bi,ei,sl))
            ++ei; // we skip ei-1
          else
          {
            bi=ei-1; // ei-1 shall not be skipt
            break;
          }
        }
        if(ei>=nb_points) break;
      }
      output.push_back(input[nb_points-1]);
      put(out_pm, output.back(), get(in_pm, input[nb_points-1]));
      if (is_closed)
      {
        output.push_back(input.back());
        put(out_pm, output.back(), get(in_pm, input.back()));
      }
      return;
    }
    case DOUGLAS_PEUCKER:
    {
      const bool is_closed = input.front()==input.back();

      std::size_t nb_points = is_closed ? input.size()-1 : input.size();

      if (nb_points<=2)
      {
        output.reserve(input.size());
        for (const auto& p : input)
        {
          output.push_back(p);
          put(out_pm, output.back(), get(in_pm, p));
        }
        return;
      }

      std::vector< std::pair<std::size_t, std::size_t> > ranges;
      ranges.push_back(std::make_pair(0, nb_points-1));

      std::vector<bool> kept(input.size(), false);
      if (is_closed) kept[nb_points]=true;
      while( !ranges.empty() )
      {
        std::size_t rb, re;
        std::tie(rb, re) = ranges.back();
        ranges.pop_back();
        kept[rb]=true;
        kept[re]=true;
        if (rb+1==re) continue;

        typename Kernel::Line_3 line(get(in_pm, input[rb]), get(in_pm, input[re]));
        double max_d = max_squared_frechet_distance;
        std::size_t max_i = 0;
        for (std::size_t i=rb; i<re; ++i)
        {
          double d = squared_distance(line, get(in_pm, input[i]));
          if (d > max_d)
          {
            max_d = d;
            max_i = i;
          }
        }
        if (max_i != 0)
        {
          ranges.push_back( std::make_pair(max_i, re) );
          ranges.push_back( std::make_pair(rb, max_i) );
        }
      }
      std::size_t nb_kept=0;
      for (std::size_t i=0; i<input.size(); ++i)
        if (kept[i]) ++nb_kept;

      output.reserve(nb_kept);
      for (std::size_t i=0; i<input.size(); ++i)
        if (kept[i])
        {
          output.push_back(input[i]);
          put(out_pm, output.back(), get(in_pm, input[i]));
        }
      //TODO if is_closed-==true, shall we add en extra step to see if we can remove output.front() and output[output.size()-2] (inital endpoints)
    }
  }
}


template <typename PointRangeIn, typename PointRangeOut>
void simplify_polyline(const PointRangeIn& input,
                             PointRangeOut& output,
                       const double max_squared_frechet_distance)
{
  simplify_polyline(input, output, max_squared_frechet_distance,
                    parameters::all_default(), parameters::all_default());
}

template <typename PointRangeIn, typename PointRangeOut, typename NamedParametersIn>
void simplify_polyline(const PointRangeIn& input,
                             PointRangeOut& output,
                       const double max_squared_frechet_distance,
                       const NamedParametersIn& np_in)
{
  simplify_polyline(input, output, max_squared_frechet_distance, np_in, parameters::all_default());
}

} } } // end of CGAL::Polygon_mesh_processing::experimental namespace


#endif // CGAL_POLYGON_MESH_PROCESSING_SIMPLIFY_POLYLINE_H
