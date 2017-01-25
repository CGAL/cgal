// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein        <wein@post.tau.ac.il>

#ifndef CGAL_DXF_BSOP_READER_H
#define CGAL_DXF_BSOP_READER_H

#include <CGAL/license/Boolean_set_operations_2.h>


#include <CGAL/IO/Dxf_reader.h>
#include <iostream>

#include <list>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/squared_distance_2.h>

namespace CGAL {

template <class Kernel_>
class Dxf_bsop_reader
{
private:

  typedef Kernel_                                   Kernel;
  typedef typename Kernel::FT                       FT;
  typedef typename Kernel::Point_2                  Point_2;
  typedef typename Kernel::Circle_2                 Circle_2;

  typedef std::list<std::pair<Point_2, double> >    Dxf_polygon_2;
  typedef std::list<Dxf_polygon_2>                  Dxf_polygons_list;

  typedef std::pair<Point_2, FT>                    Dxf_circle_2;
  typedef std::list<Dxf_circle_2>                   Dxf_circles_list;
  typedef CGAL::Dxf_reader<Kernel>                  Dxf_reader;

  typedef CGAL::Gps_circle_segment_traits_2<Kernel> Traits_2;
  typedef typename Traits_2::Point_2                Arc_point_2;
  typedef typename Traits_2::Curve_2                Curve_2;
  typedef typename Traits_2::X_monotone_curve_2     X_monotone_curve_2;

public:

  typedef typename Traits_2::Polygon_2              Circ_polygon_2;
  typedef typename Traits_2::Polygon_with_holes_2   Circ_polygon_with_holes_2;

private:

  typedef std::vector<Circ_polygon_2>               Circ_pgn_vec;
  typedef std::vector<Circ_polygon_with_holes_2>    Circ_pgn_with_holes_vec;

  typedef CGAL::General_polygon_set_2<Traits_2>     Circ_pgn_set_2;

public:

  /*!
   * Read a set of circular polygons from a DXF file.
   * \param in_file A file stream in DXF format.
   * \param pgns An output iterator of Circ_polygon_2 objects.
   * \param pgns_with_holes An output iterator of Circ_polygon_with_holes_2
   *                        objects.
   * \param simplify Should we simplify the polygons (in case the file
   *                 contains polygons that are not simple).
   * \return A pair of past-the-end iterators for both ranges.
   */ 
  template <class PolygonOutputIterator,
            class PolygonWithHolesOutputIterator>
  std::pair<PolygonOutputIterator, PolygonWithHolesOutputIterator>
  operator() (std::ifstream& in_file,
              PolygonOutputIterator pgns,
              PolygonWithHolesOutputIterator pgns_with_holes,
              bool simplify = true)
  {
    // Read polygons and circles from the DXF file.
    Circ_pgn_set_2    gps;

    Dxf_polygons_list polygons;
    Dxf_circles_list  circles;
    Dxf_reader        reader;

    reader (in_file, polygons, circles);

    // Convert the circles, such that each circle becomes a simple polygon
    // with two edges (the upper and the lower half-circles).
    Traits_2                             traits;
    typename Traits_2::Make_x_monotone_2 make_x_monotone = 
                                             traits.make_x_monotone_2_object();
    typename Dxf_circles_list::iterator  circ_it;
    CGAL::Object                         obj_vec[3];
    CGAL::Object                        *obj_begin = (obj_vec + 0);
    CGAL::Object                        *obj_end;
    X_monotone_curve_2                   cv1, cv2;

    for (circ_it = circles.begin(); circ_it != circles.end(); ++circ_it)
    {
      // Break the circle into two x-monotone circular arcs.
      const Dxf_circle_2&  dxf_circ = *circ_it;
      Curve_2              circ (dxf_circ.first, dxf_circ.second);

      obj_end = make_x_monotone (circ, obj_begin);
      CGAL_assertion(obj_end - obj_begin == 2);

      CGAL::assign(cv1, obj_vec[0]);
      CGAL::assign(cv2, obj_vec[1]);

      // Generate the corresponding polygon.
      Circ_polygon_2       pgn;
      pgn.push_back (cv1);
      pgn.push_back (cv2);
      *pgns = pgn;
      ++pgns;
    }

    circles.clear();

    // Convert the DXF polygons into circular polygons.
    Kernel                                ker;
    typename Kernel::Equal_2              equal = ker.equal_2_object();
    typename Dxf_polygons_list::iterator  pgn_it;
    typename Dxf_polygon_2::iterator      curr, next;
    Point_2                               ps, pt;
    Circle_2                              supp_circ;
    std::size_t                           n_subarcs;
    std::size_t                           i;

    for (pgn_it = polygons.begin(); pgn_it != polygons.end(); ++pgn_it)
    {
      Circ_polygon_2        pgn;

      for (curr = pgn_it->begin(); curr != pgn_it->end(); ++curr)
      {
        // Get the current vertex ps and the next vertex pt (in a circular
        // sense).
        next = curr;
        ++next;

        ps = curr->first;

        if (next != pgn_it->end())
          pt = next->first;
        else
          pt = pgn_it->begin()->first;

        // Check whether a "bulge" is defined between the two vertices.
        if (curr->second) 
        {
          // A non-zero bulge: ps and pt are connected by a circular arc.
          // We compute the center and the squared radius of its supporting
          // circle.
          CGAL_assertion (! equal (ps, pt));

          const FT bulge = curr->second;
          const FT common = (1 - CGAL::square(bulge)) / (4*bulge);
          const FT x_coord = ((ps.x() + pt.x())/2) + common*(ps.y() - pt.y());
          const FT y_coord = ((ps.y() + pt.y())/2) + common*(pt.x() - ps.x());
          const FT sqr_bulge = CGAL::square(bulge);
          const FT sqr_rad = CGAL::squared_distance(ps, pt) * 
                             (1/sqr_bulge + 2 + sqr_bulge) / 16; 

          // Construct the arc: A positive bulge means the arc is
          // counterclockwise oriented and a negative bulge means that the arc
          // is clockwise oriented.
          if(bulge > 0)
            supp_circ = Circle_2 (Point_2 (x_coord, y_coord), sqr_rad,
                                  CGAL::COUNTERCLOCKWISE);
          else
            supp_circ = Circle_2 (Point_2 (x_coord, y_coord), sqr_rad,
                                  CGAL::CLOCKWISE);

          Curve_2 circ_arc (supp_circ, 
                            Arc_point_2 (ps.x(), ps.y()),
                            Arc_point_2 (pt.x(), pt.y()));

          // Break the arc into x-monotone subarcs (there can be at most
          // three subarcs) and add them to the polygon.
          obj_end = make_x_monotone (circ_arc, obj_begin);
          n_subarcs = (obj_end - obj_begin);
          CGAL_assertion (n_subarcs <= 3);

          for (i = 0; i < n_subarcs; i++)
          {
            if (CGAL::assign (cv1, obj_vec[i]))
              pgn.push_back (cv1);
          }
        }
        else
        {
          // A zero bulge: ps and pt are connected by a straight line segment.
          if (! equal (ps, pt))
            pgn.push_back (X_monotone_curve_2 (ps, pt));
        }
      }

      // Make sure that the polygon is counterclockwise oriented.
      if (pgn.orientation() == CGAL::CLOCKWISE)
        pgn.reverse_orientation();

      // Perform polygon simplification if necessary.
      if (!simplify)
      {
        *pgns = pgn;
        ++pgns;
      }
      else
      {
        Circ_polygon_with_holes_2  pgn_with_holes;
        gps.simplify (pgn, pgn_with_holes);

        if (pgn_with_holes.number_of_holes() == 0)
        {
          // A simple polygon after all ...
          *pgns = pgn_with_holes.outer_boundary();
          ++pgns;
        }
        else
        {
          *pgns_with_holes = pgn_with_holes;
          ++pgns_with_holes;
        }
      }
    }

    polygons.clear();

    // Return a pair of past-the-end iterators.
    return (std::make_pair (pgns, pgns_with_holes));
  }
};

} //namespace CGAL

#endif
