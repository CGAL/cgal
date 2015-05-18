// Copyright (c) 2013 INRIA Sophia-Antipolis (France),
//               2014-2015 GeometryFactory (France).
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
// Author(s) : Jane Tournois, Raul Gallegos, Pierre Alliez
//

#ifndef CGAL_LLOYD_MOVE_2_H
#define CGAL_LLOYD_MOVE_2_H

#include <CGAL/Mesh_2/Uniform_sizing_field_2.h>
#include <CGAL/Constrained_voronoi_diagram_2.h>

namespace CGAL
{
namespace Mesh_2
{
  template<typename CDT,
           typename SizingField
           = Uniform_sizing_field_2<typename CDT::Geom_traits> >
  class Lloyd_move_2
  {
    typedef typename CDT::Vertex_handle          Vertex_handle;
    typedef typename CDT::Geom_traits::Point_2   Point_2;
    typedef typename CDT::Geom_traits::Segment_2  Segment;
    typedef typename CDT::Geom_traits::Triangle_2 Triangle;
    typedef typename CDT::Geom_traits::FT         FT;
    typedef typename CDT::Geom_traits::Vector_2   Vector_2;

    typedef typename CGAL::Constrained_voronoi_diagram_2<CDT> CVD;
    typedef typename CVD::Cvd_cell               Cvd_cell;

  public:
    typedef SizingField Sizing_field;

  public:
    Vector_2 operator()(Vertex_handle v,
          const CDT& cdt,
          const Sizing_field& sizing_field = Sizing_field()) const
    {
      if(cdt.are_there_incident_constraints(v))
        return CGAL::NULL_VECTOR;

      Point_2 p = v->point();
      Vector_2 move = CGAL::NULL_VECTOR;
      FT sum_masses(0);

      Cvd_cell cell = CGAL::dual(cdt, v);
      if(cell.is_infinite() || cell.is_empty())
        return CGAL::NULL_VECTOR; //center of mass is at infinity!

      CGAL_assertion(cell.number_of_vertices() > 2);

      typename Cvd_cell::segment_iterator sit = cell.segments_begin();
      for( ; sit != cell.segments_end(); ++sit)
      {
        Segment s = *sit;
        Triangle tri(p, s.source(), s.target());
        Point_2 tri_centroid = CGAL::centroid(tri);

        // Compute mass
        FT density = density_2d(tri_centroid, sizing_field);
        FT abs_area = CGAL::abs(tri.area());
        FT mass = abs_area * density;

        move = move + mass * Vector_2(p, tri_centroid);
        sum_masses += mass;
      }

      CGAL_assertion(sum_masses != 0.0);
      return move / sum_masses;
    }

  private:
    FT density_2d(const Point_2& p,
                  const Sizing_field& sizing_field) const
    {
      FT s = sizing_field(p);
      CGAL_assertion( s > 0. );

      // 1 / s^(d+2)
      return ( 1/(s*s*s*s) );
    }

#ifdef CGAL_MESH_2_OPTIMIZER_VERBOSE
  public:
    static std::string name() { return std::string("Lloyd"); }
#endif

  };

} //end namespace Mesh_2
} //end namespace CGAL

#endif //CGAL_LLOYD_MOVE_2_H
