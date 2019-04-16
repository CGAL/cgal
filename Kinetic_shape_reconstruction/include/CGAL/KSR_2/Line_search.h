// Copyright (c) 2019 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSR_2_LINE_SEARCH_H
#define CGAL_KSR_2_LINE_SEARCH_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <memory>

#include <CGAL/KSR/utils.h>
#include <CGAL/KSR_2/Vertex.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_primitive.h>

//#define USE_AABB_TREE

namespace CGAL
{

namespace KSR_2
{

template <typename GeomTraits>
class Line_search
{
public:
  typedef GeomTraits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Line_2 Line_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Segment_3 Segment_3;

  typedef std::pair<Segment_2, KSR::size_t> Segment_with_id;

  typedef std::vector<Segment_with_id> Segments_with_ids;
  typedef typename Segments_with_ids::iterator iterator;

  static Segment_3 to_3d (const Segment_2& segment) 
  {
    return Segment_3 (Point_3 (segment.source().x(), segment.source().y(), 0),
                      Point_3 (segment.target().x(), segment.target().y(), 0));
  }

  struct Segment_3_of_segment_iterator_property_map
  {
    typedef iterator key_type;
    typedef Segment_3 value_type;
    typedef value_type reference;
    typedef boost::readable_property_map_tag category;

    inline friend reference get (const Segment_3_of_segment_iterator_property_map&, key_type it)
    {
      return to_3d (it->first);
    }
    
  };
  struct Source_of_segment_iterator_property_map
  {
    typedef iterator key_type;
    typedef Point_3 value_type;
    typedef value_type reference;
    typedef boost::readable_property_map_tag category;

    inline friend reference get (const Source_of_segment_iterator_property_map&, key_type it)
    {
      return reference (it->first.source().x(), it->first.source().y(), 0);
    }
    
  };
  
  typedef CGAL::AABB_primitive<iterator,
                               Segment_3_of_segment_iterator_property_map,
                               Source_of_segment_iterator_property_map,
                               Tag_false,
                               Tag_false> Primitive;
  typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Tree;
  
  typedef typename Tree::template Intersection_and_primitive_id<Segment_3>::Type Intersection_type;
  typedef boost::optional<Intersection_type> Segment_intersection;
  
  
private:

  Segments_with_ids m_input;
  std::unique_ptr<Tree> m_tree;

public:

  Line_search() { }

  void add_line (KSR::size_t line_idx, const Segment_2& segment, FT discretization_step)
  {
    FT size = CGAL::approximate_sqrt (segment.squared_length());
    std::size_t nb_items = std::size_t(size / discretization_step) + 1;

    m_input.reserve (m_input.size() + nb_items);

    for (std::size_t i = 0; i < nb_items; ++ i)
    {
      Segment_2 seg (CGAL::barycenter (segment.source(), i / FT(nb_items),
                                       segment.target()),
                     CGAL::barycenter (segment.source(), (i+1) / FT(nb_items),
                                       segment.target()));
      m_input.push_back (std::make_pair (seg, line_idx));
    }

  }

  void build()
  {
    m_tree = std::make_unique<Tree>(m_input.begin(), m_input.end());
  }

  void compute_intersected_lines (const Segment_2& segment, std::vector<std::pair<KSR::size_t, Point_2> >& intersected_lines)
  {
    Segment_3 query = to_3d (segment);
    m_tree->all_intersections
      (query,
       boost::make_function_output_iterator
       ([&](const Segment_intersection& intersection) -> void
        {
          if (intersection)
          {
            const Point_3* p = boost::get<Point_3>(&(intersection->first));
            if(p)
              intersected_lines.push_back (std::make_pair (intersection->second->second, Point_2 (p->x(), p->y())));
          }
        }));
  }
};


}} // namespace CGAL::KSR_2


#endif // CGAL_KSR_2_LINE_SEARCH_H
