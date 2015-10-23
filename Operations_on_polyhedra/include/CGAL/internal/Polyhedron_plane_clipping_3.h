// Copyright (c) 2015 GeometryFactory (France).
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
// Author(s)     : Sebastien Loriot

#ifndef CGAL_INTERNAL_POLYHEDRON_PLANE_CLIPPING_3_H
#define CGAL_INTERNAL_POLYHEDRON_PLANE_CLIPPING_3_H

#include <CGAL/corefinement_operations.h>
#include <CGAL/iterator.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/array.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_2_filtered_projection_traits_3.h>
#include <CGAL/convex_hull_3.h>

namespace CGAL{
namespace corefinement{

namespace internal{
 template <class HDS,class T>
 class Builder_from_T_2 : public CGAL::Modifier_base<HDS> {
   typedef std::map<typename T::Vertex_handle,unsigned> Vertex_map;

   const T& t;
   template <class Builder>
   static unsigned get_vertex_index( Vertex_map& vertex_map,
                                     typename T::Vertex_handle vh,
                                     Builder& builder,
                                     unsigned& vindex)
   {
     std::pair<typename Vertex_map::iterator,bool>
       res=vertex_map.insert(std::make_pair(vh,vindex));
     if (res.second){
       builder.add_vertex(vh->point());
         ++vindex;
     }
     return res.first->second;
   }

   public:
   Builder_from_T_2(const T& t_):t(t_)
   {
     CGAL_assertion(t.dimension()==2);
   }
   void operator()( HDS& hds) {
     // Postcondition: `hds' is a valid polyhedral surface.
     CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
     Vertex_map vertex_map;
     //start the surface
     B.begin_surface( t.number_of_vertices(), t.number_of_faces());
     unsigned vindex=0;
     for (typename T::Finite_faces_iterator it=t.finite_faces_begin();
                                            it!=t.finite_faces_end();++it)
     {
       unsigned i0=get_vertex_index(vertex_map,it->vertex(0),B,vindex);
       unsigned i1=get_vertex_index(vertex_map,it->vertex(1),B,vindex);
       unsigned i2=get_vertex_index(vertex_map,it->vertex(2),B,vindex);
       B.begin_facet();
       B.add_vertex_to_facet( i0 );
       B.add_vertex_to_facet( i1 );
       B.add_vertex_to_facet( i2 );
       B.end_facet();
     }
     B.end_surface();
   }
 };
} // end of namespace internal

template <class Polyhedron, class Plane_3>
Polyhedron clip_to_bbox(const Bbox_3& bbox, const Plane_3& plane)
{
  typedef typename Polyhedron::Traits::Kernel Kernel;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Segment_3 Segment_3;
  cpp11::array<typename Kernel::Point_3,8> corners= {{
    Point_3(bbox.xmin(),bbox.ymin(),bbox.zmin()),
    Point_3(bbox.xmin(),bbox.ymax(),bbox.zmin()),
    Point_3(bbox.xmax(),bbox.ymax(),bbox.zmin()),
    Point_3(bbox.xmax(),bbox.ymin(),bbox.zmin()),
    Point_3(bbox.xmin(),bbox.ymin(),bbox.zmax()),
    Point_3(bbox.xmin(),bbox.ymax(),bbox.zmax()),
    Point_3(bbox.xmax(),bbox.ymax(),bbox.zmax()),
    Point_3(bbox.xmax(),bbox.ymin(),bbox.zmax())
  }};

  cpp11::array<CGAL::Oriented_side,8> orientations = {{
    plane.oriented_side(corners[0]),
    plane.oriented_side(corners[1]),
    plane.oriented_side(corners[2]),
    plane.oriented_side(corners[3]),
    plane.oriented_side(corners[4]),
    plane.oriented_side(corners[5]),
    plane.oriented_side(corners[6]),
    plane.oriented_side(corners[7])
  }};

  std::vector<Point_3> intersection_points;
  // first look for intersections at corners
  for (int i=0; i<8; ++i)
    if (orientations[i]==ON_ORIENTED_BOUNDARY)
      intersection_points.push_back(corners[i]);
  // second look for intersections on edges
  cpp11::array<int,24> edge_indices = {{ // 2 *12 edges
    0,1, 1,2, 2,3, 3,0, // bottom face edges
    4,5, 5,6, 6,7, 7,4, // top face edges
    0,4, 1,5, 2,6, 3,7
  }};

  for (int i=0; i<12; ++i)
  {
    int i1=edge_indices[2*i], i2=edge_indices[2*i+1];
    if (orientations[i1]==ON_ORIENTED_BOUNDARY) continue;
    if (orientations[i2]==ON_ORIENTED_BOUNDARY) continue;
    if (orientations[i1]!=orientations[i2])
      intersection_points.push_back(
        get<Point_3>(
          *CGAL::intersection(plane, Segment_3(corners[i1], corners[i2]) )
        )
      );
  }

  Polyhedron P;

  //if less that 3 points there will be nothing to clipping.
  if (intersection_points.size()<3) return P;

  //triangulate the set of intersection points (I know it's overkill)
  typedef CGAL::Triangulation_2_filtered_projection_traits_3<Kernel>   P_traits;
  typedef CGAL::Delaunay_triangulation_2<P_traits> DT;
  DT dt(P_traits(plane.orthogonal_vector()));
  dt.insert(intersection_points.begin(),
            intersection_points.end());

  //now create the polyhedron from the triangulation
  internal::Builder_from_T_2< typename Polyhedron::HalfedgeDS,DT > builder(dt);
  P.delegate(builder);

  return P;
}

template <class Polyhedron>
Polyhedron* clip_polyhedron(Polyhedron& P, Polyhedron& clipping_polyhedron)
{
  std::pair <Polyhedron*,int> result;
  typedef CGAL::Polyhedron_corefinement<Polyhedron> Corefinement;
  Corefinement coref;
  CGAL::Emptyset_iterator emptyset_iterator;
  coref(P, clipping_polyhedron, emptyset_iterator,
        &result, Corefinement::Intersection_tag);

  return result.first;
}

template <class Polyhedron, class Plane_3>
Polyhedron* clip_polyhedron(const Polyhedron& P, const Plane_3& p)
{
  CGAL::Bbox_3 bbox = CGAL::bbox_3(P.points_begin(), P.points_end());
  //extend the bbox a bit to avoid border cases
  double xd=(bbox.xmax()-bbox.xmin())/100;
  double yd=(bbox.ymax()-bbox.ymin())/100;
  double zd=(bbox.zmax()-bbox.zmin())/100;
  bbox=CGAL::Bbox_3(bbox.xmin()-xd, bbox.ymin()-yd, bbox.zmin()-zd,
                    bbox.xmax()+xd, bbox.ymax()+yd, bbox.zmax()+zd);
  Polyhedron clipping_polyhedron=clip_to_bbox<Polyhedron>(bbox, p);

  if (clipping_polyhedron.empty()) //no intersection, result is all or nothing
  {
    if (p.oriented_side(*P.points_begin())==ON_POSITIVE_SIDE)
      return new Polyhedron();
    else
      return new Polyhedron(P);
  }
  Polyhedron copy(P);
  return clip_polyhedron(copy, clipping_polyhedron);
}

} } // CGAL::corefinement


#endif // CGAL_INTERNAL_POLYHEDRON_PLANE_CLIPPING_3_H
