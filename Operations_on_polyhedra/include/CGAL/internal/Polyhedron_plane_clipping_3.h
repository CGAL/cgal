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
#include <CGAL/Triangulation_2_projection_traits_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
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
  typedef CGAL::Triangulation_2_projection_traits_3<Kernel>   P_traits;
  typedef CGAL::Delaunay_triangulation_2<P_traits> DT;
  DT dt(P_traits(plane.orthogonal_vector()));
  dt.insert(intersection_points.begin(),
            intersection_points.end());

  // tangency with the bbox -> no intersection
  if (dt.dimension()!=2) return P;

  //now create the polyhedron from the triangulation
  internal::Builder_from_T_2< typename Polyhedron::HalfedgeDS,DT > builder(dt);
  P.delegate(builder);

  return P;
}

template <class Polyhedron, class Plane_3>
Polyhedron clip_bbox(const Bbox_3& bbox, const Plane_3& plane)
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

  std::vector<Point_3> points;
  // first look for intersections at corners
  for (int i=0; i<8; ++i)
    if (orientations[i]==ON_ORIENTED_BOUNDARY)
      points.push_back(corners[i]);
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
      points.push_back(
        get<Point_3>(
          *CGAL::intersection(plane, Segment_3(corners[i1], corners[i2]) )
        )
      );
  }

  Polyhedron P;

  //if less that 3 points there will be nothing to clipping.
  if (points.size()<3) return P;

  for (int i=0; i<8; ++i)
    if (orientations[i]==ON_NEGATIVE_SIDE)
      points.push_back(corners[i]);

  // take the convex hull of the points on the negative side+intersection points
  // overkill...
  CGAL::convex_hull_3(points.begin(), points.end(), P);

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
  if(P.empty()) return new Polyhedron();
  CGAL::Bbox_3 bbox( CGAL::bbox_3(P.points_begin(), P.points_end()) );
  //extend the bbox a bit to avoid border cases
  double xd=(bbox.xmax()-bbox.xmin())/100;
  double yd=(bbox.ymax()-bbox.ymin())/100;
  double zd=(bbox.zmax()-bbox.zmin())/100;
  bbox=CGAL::Bbox_3(bbox.xmin()-xd, bbox.ymin()-yd, bbox.zmin()-zd,
                    bbox.xmax()+xd, bbox.ymax()+yd, bbox.zmax()+zd);
  Polyhedron clipping_polyhedron=clip_bbox<Polyhedron>(bbox, p);

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

template <class Polyhedron, class Plane_3>
std::pair<Polyhedron*,Polyhedron*> split_polyhedron(const Polyhedron& P, const Plane_3& p)
{
 std::pair<Polyhedron*, Polyhedron*> res;
  if(P.empty()) {res.first=res.second= new Polyhedron(); return res;}
  CGAL::Bbox_3 bbox( CGAL::bbox_3(P.points_begin(), P.points_end()) );
  //extend the bbox a bit to avoid border cases
  double xd=(bbox.xmax()-bbox.xmin())/100;
  double yd=(bbox.ymax()-bbox.ymin())/100;
  double zd=(bbox.zmax()-bbox.zmin())/100;
  bbox=CGAL::Bbox_3(bbox.xmin()-xd, bbox.ymin()-yd, bbox.zmin()-zd,
                    bbox.xmax()+xd, bbox.ymax()+yd, bbox.zmax()+zd);
  //First Polyhedron
  Polyhedron clipping_polyhedron=clip_bbox<Polyhedron>(bbox, p);


  if (clipping_polyhedron.empty()) //no intersection, result is all or nothing
  {
    if (p.oriented_side(*P.points_begin())==ON_POSITIVE_SIDE)
    {
      res.first = new Polyhedron(); res.second = new Polyhedron(P);
    }
    else
    {
      res.first = new Polyhedron(P); res.second = new Polyhedron();
    }
    return res;
  }

  Polyhedron copy(P);
  res.first = clip_polyhedron(copy, clipping_polyhedron);
  //Second Polyhedron
  Plane_3 p_op = p.opposite();
  clipping_polyhedron=clip_bbox<Polyhedron>(bbox, p_op);

  copy = P;
  res.second = clip_polyhedron(copy, clipping_polyhedron);
  return res;

}

namespace internal{

template<class Polyhedron>
struct Edge_is_marked4coref{
  std::set<typename Polyhedron::Halfedge_handle>& marked_halfedges;
  typedef bool value_type;
  typedef value_type reference;
  typedef std::pair<typename Polyhedron::Halfedge_handle,Polyhedron*> key_type;
  typedef boost::read_write_property_map_tag category;

  Edge_is_marked4coref(std::set<typename Polyhedron::Halfedge_handle>& mh)
  : marked_halfedges(mh)
  {}

  friend reference get(Edge_is_marked4coref& map,const key_type& key) {
    return map.marked_halfedges.count(key.first)!=0;
  }
  friend void put(Edge_is_marked4coref& map,key_type key,value_type v) {
    if (v) map.marked_halfedges.insert(key.first);
    else  map.marked_halfedges.erase(key.first);
  }
};

template<class Polyhedron>
struct Edge_is_marked{
  const std::set<typename Polyhedron::Halfedge_handle>* marked_halfedges;
  typedef bool value_type;
  typedef value_type reference;
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor key_type;
  typedef boost::readable_property_map_tag category;

  Edge_is_marked(){}
  Edge_is_marked(const std::set<typename Polyhedron::Halfedge_handle>& mh)
  : marked_halfedges(&mh)
  {}

  friend reference get(const Edge_is_marked& map,const key_type& key) {
    return map.marked_halfedges->count(key.halfedge())!=0;
  }
};

} //end of internal namespace

template <class Polyhedron, class Plane_3>
void inplace_clip_open_polyhedron(Polyhedron& P, const Plane_3& p)
{
  CGAL::Bbox_3 bbox( CGAL::bbox_3(P.points_begin(), P.points_end()) );
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
      P.clear();
    return;
  }

  // set for marking edges of P intersected by the clipping plane
  std::set<typename Polyhedron::Halfedge_handle> marked_halfedges;
  internal::Edge_is_marked4coref<Polyhedron> cr_edge_is_marked(marked_halfedges);
  typedef typename Polyhedron::Traits::Kernel K;
  typedef CGAL::Node_visitor_refine_polyhedra<
            Polyhedron,K,internal::Edge_is_marked4coref<Polyhedron> >
                  Split_visitor;
    Split_visitor visitor(NULL, true, cr_edge_is_marked);
    CGAL::Intersection_of_Polyhedra_3<Polyhedron,K,Split_visitor>
      polyline_intersections(visitor);
    CGAL::Emptyset_iterator emptyset_iterator;
    // corefinement P and clipping_polyhedron
    polyline_intersections(P,clipping_polyhedron,emptyset_iterator);

    // extract connected components bounded by marked edges
    internal::Edge_is_marked<Polyhedron> edge_is_marked(marked_halfedges);
    namespace PMP=Polygon_mesh_processing;
    std::map<typename Polyhedron::Face_handle, std::size_t> face_ccs;
    std::size_t nb_cc=PMP::connected_components(P,
      boost::make_assoc_property_map(face_ccs),
      PMP::parameters::edge_is_constrained_map(edge_is_marked)
        .face_index_map(get(boost::face_external_index,P))
    );

    // remove cc on the positive side of the plane
    std::vector<bool> cc_handled(nb_cc, false);
    std::vector<std::size_t> ccs_to_remove;
    BOOST_FOREACH(typename Polyhedron::Face_handle f, faces(P))
    {
      std::size_t cc_id=face_ccs[f];
      if (cc_handled[cc_id]) continue;

      //look for a vertex not on the intersection
      typename Polyhedron::Halfedge_handle h=f->halfedge();
      for(int i=0;i<3;++i){
        bool no_marked_edge=true;
        BOOST_FOREACH(typename Polyhedron::Halfedge_handle h, halfedges_around_target(h, P))
          if ( marked_halfedges.count(h) )
            no_marked_edge=false;
        if (no_marked_edge){
          if ( p.oriented_side(h->vertex()->point())==ON_POSITIVE_SIDE )
            ccs_to_remove.push_back(cc_id);
          cc_handled[cc_id]=true;
          if (--nb_cc==0) break;
          break;
        }
        h=h->next();
      }
    }

  //now remove the faces on the positive side
  PMP::remove_connected_components(P,
    ccs_to_remove,
    boost::make_assoc_property_map(face_ccs),
    PMP::parameters::vertex_index_map(get(boost::vertex_external_index,P))
  );
}

} } // CGAL::corefinement


#endif // CGAL_INTERNAL_POLYHEDRON_PLANE_CLIPPING_3_H
