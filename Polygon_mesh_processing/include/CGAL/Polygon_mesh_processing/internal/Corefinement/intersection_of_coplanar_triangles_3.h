// Copyright (c) 2016 GeometryFactory (France).
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

#ifndef CGAL_PMP_INTERNAL_COREFINEMENT_INTERSECTION_OF_COPLANAR_TRIANGLES_3_H
#define CGAL_PMP_INTERNAL_COREFINEMENT_INTERSECTION_OF_COPLANAR_TRIANGLES_3_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>


#include <CGAL/Polygon_mesh_processing/internal/Corefinement/Intersection_type.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/property_map.h>

#include <bitset>

namespace CGAL{
namespace Corefinement{

template <class TriangleMesh, class VertexPointMap>
struct Intersect_coplanar_faces_3{
 // typedefs
  typedef typename boost::property_traits<VertexPointMap>::value_type Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel Input_kernel;
  typedef CGAL::Exact_predicates_exact_constructions_kernel Exact_kernel;

  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;

  typedef Coplanar_intersection<TriangleMesh, Exact_kernel> Inter_pt_info;
// data members
  const TriangleMesh &tm1, &tm2;
  const VertexPointMap &vpm1, &vpm2;
// constructor
  Intersect_coplanar_faces_3(const TriangleMesh& tm1_,
                             const TriangleMesh& tm2_,
                             const VertexPointMap& vpm1_,
                             const VertexPointMap& vpm2_)
  : tm1(tm1_), tm2(tm2_), vpm1(vpm1_), vpm2(vpm2_)
  {}

  typename Exact_kernel::Point_3
  to_exact(const Point& p)
  {
    return typename Exact_kernel::Point_3(p.x(), p.y(), p.z());
  }

// function construction Inter_pt_info objects
  //build the intersection point as the target vertex of h1 in the face of h2
  Inter_pt_info operator()(halfedge_descriptor h1,halfedge_descriptor h2)
  {
    Inter_pt_info res;
    res.type_1=ON_VERTEX;
    res.type_2=ON_FACE;
    res.info_1=h1,
    res.info_2=h2;
    res.point=to_exact(get(vpm1,target(h1,tm1)));
    return res;
  }


  //constructor for intersection of edges. prev and curr are two points on an edge of the first facet (preserving the
  //orientation of the facet). This edge is intersected by h2 from the second facet.
  //
  //The rational is the following: we first check whether curr and prev are on the same edge. I so we create
  //an intersection point between two edges. Otherwise, the point is a vertex of the second facet included into
  //the first facet.
  //
  //(V,F) : point initialy constructed
  //(V,E) : (V,F) updated by get_orientation_and_update_info_2 (i.e lies on one edge)
  //(V,V) : (V,E) updated by get_orientation_and_update_info_2 (i.e lies on two edges)
  //(E,E) : created in the following function when prev and curr lie on the same edge
  //(E,V) : (E,E) updated by get_orientation_and_update_info_2 (always done as lies on two edges)
  //(E,F) : impossible
  //(F,V) : detected when curr and prev and not on the same edge
  //(F,E) : impossible
  //(F,F) : impossible
  //
  Inter_pt_info
  operator()(Inter_pt_info ipt_prev, Inter_pt_info ipt_curr,
             halfedge_descriptor h1, halfedge_descriptor h2)
  {
    Inter_pt_info res;
    res.type_2=ON_EDGE;
    res.info_2=h2;

    if (ipt_prev.type_1==ON_VERTEX && next(ipt_prev.info_1, tm1) == ipt_curr.info_1){
      CGAL_assertion(ipt_curr.type_1!=ON_FACE);
      res.type_1=ON_EDGE;
      res.info_1=ipt_curr.info_1;
    }
    else{
      if(ipt_curr.type_1==ON_VERTEX && ipt_prev.info_1 == ipt_curr.info_1){
        CGAL_assertion(ipt_prev.type_1!=ON_FACE);
        res.type_1=ON_EDGE;
        res.info_1=ipt_curr.info_1;
      }
      else{
        if (ipt_curr.type_1==ON_EDGE && ipt_prev.type_1==ON_EDGE &&  ipt_curr.info_1==ipt_prev.info_1){
          res.type_1=ON_EDGE;
          res.info_1=ipt_curr.info_1;
        }
        else{
          //ipt_curr and ipt_prev are not on the same edge of the first facet.
          //The intersection point to be computed is a VERTEX of the second facet
          res.type_1=ON_FACE;
          res.info_1=h1;
          res.type_2=ON_VERTEX;

          //this is used to select the correct endpoint of the edge of the second facet
          typename Exact_kernel::Collinear_3 is_collinear = Exact_kernel().collinear_3_object();
          if ( !is_collinear(ipt_prev.point,ipt_curr.point,to_exact(get(vpm2,target(res.info_2,tm2)) ) ) ){
            res.info_2=prev(res.info_2,tm2);
            CGAL_assertion( is_collinear(ipt_prev.point,ipt_curr.point,to_exact(get(vpm2,target(res.info_2,tm2))) ) );
          }
          res.point = to_exact( get(vpm2, target(res.info_2,tm2)) );
          return res;
        }
      }
    }
    //
    //handle degenerate case when two edges overlap
    //at least one of the two vertex has already been found as a vertex of a facet. Here we set it for the second point
    if(ipt_prev.type_2!=ON_FACE && ipt_curr.type_2!=ON_FACE && (ipt_prev.type_1==ON_VERTEX ||
       ipt_prev.type_2==ON_VERTEX) && (ipt_curr.type_1==ON_VERTEX || ipt_curr.type_2==ON_VERTEX))
    {
      typename Exact_kernel::Collinear_3 is_collinear = Exact_kernel().collinear_3_object();
      if ( is_collinear(ipt_prev.point,
                        ipt_curr.point,
                        to_exact(get(vpm2, source(res.info_2,tm2))) ) )
      {
        res.info_2=prev(res.info_2,tm2);
        res.type_2=ON_VERTEX;
        res.point = to_exact( get(vpm2, target(res.info_2,tm2)) );
        return res;
      }
      if ( is_collinear(ipt_prev.point,
                        ipt_curr.point,
                        to_exact(get(vpm2, target(res.info_2,tm2)) ) ) )
      {
        res.type_2=ON_VERTEX;
        res.point = to_exact( get(vpm2, target(res.info_2,tm2)) );
        return res;
      }
    }

    //handle regular intersection of two edges
    typename Exact_kernel::Construct_line_line_intersection_point_3 intersect;
    res.point = intersect(to_exact(get(vpm2, target(res.info_2,tm2))),
                          to_exact(get(vpm2, source(res.info_2,tm2))),
                          ipt_prev.point, ipt_curr.point );
    return res;
  }

//
  CGAL::Orientation
  get_orientation_and_update_info_2(halfedge_descriptor h2,Inter_pt_info& ipt)
  {
    typename Exact_kernel::Coplanar_orientation_3 orient=
      Exact_kernel().coplanar_orientation_3_object();

    Orientation res = orient(to_exact(get(vpm2,source(h2,tm2))),
                             to_exact(get(vpm2,target(h2,tm2))),
                             to_exact(get(vpm2,target(next(h2,tm2),tm2))),
                             ipt.point);

    if ( (ipt.type_1==ON_VERTEX || ipt.type_1==ON_EDGE) && res==COLLINEAR){
      if (ipt.type_2==ON_FACE){ //detect a case (VERTEX,EDGE)
        ipt.type_2=ON_EDGE;
        ipt.info_2=h2;
      }
      else{
        //detect a case (VERTEX,VERTEX) or (EDGE,VERTEX)
        CGAL_assertion(ipt.type_2==ON_EDGE);
        ipt.type_2=ON_VERTEX;
        if (next(ipt.info_2,tm2)!=h2){
          CGAL_assertion(next(h2,tm2)==ipt.info_2);
          ipt.info_2=h2;
        }
      }
    }
    return res;
  }

  void cutoff_face(
    halfedge_descriptor h2,
    std::list<Inter_pt_info>& inter_pts,
    halfedge_descriptor h1)
  {
    if ( inter_pts.empty() ) return;
    typedef typename std::list<Inter_pt_info>::iterator Iterator;

    std::map<Inter_pt_info*,Orientation> orientations;
    BOOST_FOREACH(Inter_pt_info& ipt, inter_pts)
      orientations[ &ipt ]=get_orientation_and_update_info_2(h2,ipt);

    CGAL_assertion_code(int pt_added=0;)

    Inter_pt_info* prev = &(*boost::prior(inter_pts.end()));
    bool inter_pts_size_g_2 = inter_pts.size() > 2;
    Iterator stop = inter_pts_size_g_2 ? inter_pts.end() : boost::prior(inter_pts.end());
    for (Iterator it=inter_pts.begin();it!=stop;++it)
    {
      Inter_pt_info* curr=&(*it);
      if (!inter_pts_size_g_2) std::swap(prev,curr);
      Orientation or_prev=orientations[prev],or_curr=orientations[curr];
      if ( (or_prev==POSITIVE && or_curr==NEGATIVE) || (or_prev==NEGATIVE && or_curr==POSITIVE) )
      {
        Iterator it_curr = inter_pts_size_g_2 ? it:boost::next(it);
        prev=&(* inter_pts.insert( it_curr,operator()(*prev,*curr,h1,h2) ) );
        orientations[prev]=COLLINEAR;
        CGAL_assertion_code(++pt_added;)
      }
      prev=&(*it);
    }

    CGAL_assertion(pt_added<3);
    Iterator it=inter_pts.begin();
    std::size_t nb_interpt=inter_pts.size();
    //this boolean allows to reverse order of intersection points in case there were 3 remaining intersection points
    //and the point in the middle was removed. In that case the order must be reversed to preserve the orientations
    //of the last edge:
    // A---X---B  --> AB  to be consistent with the other cases this should be BA!
    // X---B---A  --> BA
    // B---A---X  --> BA
    //

    bool should_revert_list=false;

    while(it!=inter_pts.end())
    {
      if (orientations[&(*it)]==NEGATIVE){
        inter_pts.erase(it++);
        if (--nb_interpt == 2 && it!=inter_pts.end() && boost::next(it)==inter_pts.end()) should_revert_list=true;
      }
      else
        ++it;
    }
    if (should_revert_list && nb_interpt==2) inter_pts.reverse();
  }
};

template <class TriangleMesh, class VertexPointMap, class Exact_kernel>
void intersection_coplanar_faces(
  typename boost::graph_traits<TriangleMesh>::face_descriptor f1,
  typename boost::graph_traits<TriangleMesh>::face_descriptor f2,
  const TriangleMesh& tm1,
  const TriangleMesh& tm2,
  const VertexPointMap& vpm1,
  const VertexPointMap& vpm2,
  std::list< Coplanar_intersection<TriangleMesh, Exact_kernel> >& inter_pts)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;

  halfedge_descriptor h1=halfedge(f1,tm1), h2=halfedge(f2,tm2);

  Intersect_coplanar_faces_3<TriangleMesh, VertexPointMap>
    intersect_cpln(tm1, tm2, vpm1, vpm2);

  // We will add in `inter_pts` the initial triangle of h1
  inter_pts.push_back( intersect_cpln(h1,h2) );
  inter_pts.push_back( intersect_cpln(next(h1,tm1),h2) );
  inter_pts.push_back( intersect_cpln(next(next(h1,tm1),tm1),h2) );

  // We now cut the initial triangle with the halfspaces defined by each
  // oriented edge of the second triangle face
  intersect_cpln.cutoff_face(h2,inter_pts,h1);
  intersect_cpln.cutoff_face(next(h2,tm2),inter_pts,h1);
  intersect_cpln.cutoff_face(next(next(h2,tm2),tm2),inter_pts,h1);
}

} } //namespace CGAL::Corefinement


#endif //CGAL_PMP_INTERNAL_COREFINEMENT_INTERSECTION_OF_COPLANAR_TRIANGLES_3_H
