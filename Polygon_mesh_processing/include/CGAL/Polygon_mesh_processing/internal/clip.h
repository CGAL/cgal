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

#ifndef CGAL_POLYGON_MESH_PROCESSING_CLIP_H
#define CGAL_POLYGON_MESH_PROCESSING_CLIP_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>


#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <CGAL/convex_hull_3.h>

namespace CGAL{
namespace Polygon_mesh_processing {

namespace internal
{
template <class TriangleMesh,
          class Ecm,
          class NamedParameters1,
          class NamedParameters2>
bool
clip_open_impl(      TriangleMesh& tm,
                     TriangleMesh& clipper,
               Ecm ecm,
               const NamedParameters1& np_tm,
               const NamedParameters2& np_c)
{
  // first corefine the meshes
  corefine(tm, clipper, np_tm, np_c);

  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters1>::type Vpm;
  typedef typename GetFaceIndexMap<TriangleMesh,
                                   NamedParameters1>::type Fid_map;
  typedef typename GetVertexIndexMap<TriangleMesh,
                                     NamedParameters1>::type Vid_map;

  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters2>::type GeomTraits;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::face_descriptor face_descriptor;

  Fid_map fid_map = boost::choose_param(get_param(np_tm, internal_np::face_index),
                                        get_property_map(boost::face_index, tm));
  Vid_map vid_map = boost::choose_param(get_param(np_tm, internal_np::vertex_index),
                                        get_property_map(boost::vertex_index, tm));
  Vpm vpm1 = boost::choose_param(get_param(np_tm, internal_np::vertex_point),
                                 get_property_map(vertex_point, tm));
  Vpm vpm_c = boost::choose_param(get_param(np_c, internal_np::vertex_point),
                                  get_property_map(vertex_point, clipper));

  // init indices if needed
  helpers::init_face_indices(tm, fid_map);
  helpers::init_vertex_indices(tm, vid_map);

  // set the connected component id of each face
  std::vector<std::size_t> face_cc(num_faces(tm), std::size_t(-1));
  std::size_t nb_cc =
    connected_components(tm,
                         bind_property_maps(fid_map, make_property_map(face_cc)),
                         parameters::face_index_map(fid_map).
                         edge_is_constrained_map(ecm));


  boost::dynamic_bitset<> cc_not_handled(nb_cc);
  cc_not_handled.set();
  std::vector <std::size_t> ccs_to_remove;
  /// \todo clipper has been modified, this is not robust if inexact constructions are used
  CGAL::Side_of_triangle_mesh<TriangleMesh, GeomTraits, Vpm>
    side_of(clipper, vpm_c);
  BOOST_FOREACH(face_descriptor f, faces(tm))
  {
    std::size_t cc_id = face_cc[ get(fid_map, f) ];
    if ( !cc_not_handled.test(cc_id) ) continue;

    halfedge_descriptor h=halfedge(f, tm);
    for(int i=0;i<3;++i)
    {
      bool no_marked_edge=true;
      BOOST_FOREACH(halfedge_descriptor h2, halfedges_around_target(h, tm))
        if ( get(ecm, edge(h2, tm)) ){
          no_marked_edge=false;
          break;
        }
      if (no_marked_edge){
        if ( side_of( get(vpm1, target(h, tm) ) ) == ON_UNBOUNDED_SIDE )
          ccs_to_remove.push_back(cc_id);
        cc_not_handled.reset(cc_id);
        break;
      }
      h=next(h, tm);
    }
    if (!cc_not_handled.any()) break;
  }

  if (cc_not_handled.any())
  {
    ///\todo handle using barycenters? won't work for coplanar faces
  }
  //now remove the cc
  remove_connected_components(tm,
    ccs_to_remove,
    bind_property_maps(fid_map, make_property_map(face_cc)),
    np_tm);

  return true;
}

/// \TODO move this to property_map.h?
template <class Set>
struct Constrained_edge_map
{
  typedef boost::read_write_property_map_tag      category;
  typedef bool                                  value_type;
  typedef bool                                  reference;
  typedef typename Set::value_type              key_type;

  Constrained_edge_map()
    : edge_set(NULL)
  {}

  Constrained_edge_map(Set& set)
    : edge_set(&set)
  {}

  friend bool get(const Constrained_edge_map<Set>& map, key_type k)
  {
    return map.edge_set->count(k)!=0;
  }

  friend void put(Constrained_edge_map<Set>& map, key_type k, bool b)
  {
    if (b)
      map.edge_set->insert(k);
    else
      map.edge_set->erase(k);
  }
private:
  Set* edge_set;
};

template <class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2>
bool
clip_open_impl(      TriangleMesh& tm,
                     TriangleMesh& clipper,
               boost::param_not_found,
               const NamedParameters1& np_tm,
               const NamedParameters2& np_c)
{
  typedef typename boost::graph_traits<TriangleMesh>
    ::edge_descriptor edge_descriptor;
  boost::unordered_set<edge_descriptor> constrained_edges;
  Constrained_edge_map<boost::unordered_set<edge_descriptor> >
    cst_map(constrained_edges);

  return clip_open_impl(tm, clipper,
                        cst_map,
                        np_tm.edge_is_constrained_map(cst_map),
                        np_c);
}

} // end of internal namespace

#ifndef DOXYGEN_RUNNING

///\todo clipper const!
/// requires face_index_map, vertex_index_map for np_tm
/// requires face_index_map for np_c
/// if edge_is_constrained_map is not provided in np_tm a default one is
/// provided using boost::unordered_set<edge_descriptor>
template <class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2>
bool
clip(      TriangleMesh& tm,
     /*const*/ TriangleMesh& clipper,
           bool close,
     const NamedParameters1& np_tm,
     const NamedParameters2& np_c)
{
  if (close && is_closed(tm))
    return corefine_and_compute_intersection(tm, clipper, tm, np_tm, np_c);

  return internal::clip_open_impl(tm, clipper,
      get_param(np_tm, internal_np::edge_is_constrained), np_tm, np_c);
}

/// \todo document me
template <class Plane_3,
          class TriangleMesh,
          class NamedParameters>
Oriented_side
clip_to_bbox(const Plane_3& plane,
             const Bbox_3& bbox,
                   TriangleMesh& tm_out,
             const NamedParameters& np )
{
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type Geom_traits;
  typedef typename Geom_traits::Point_3 Point_3;
  typedef typename Geom_traits::Segment_3 Segment_3;
  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters>::type Vpm;

  Vpm vpm_out = boost::choose_param(get_param(np, internal_np::vertex_point),
                                    get_property_map(boost::vertex_point, tm_out));


  cpp11::array<Point_3,8> corners= {{
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

  // look for intersections on edges
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
        boost::get<Point_3>(
          *intersection(plane, Segment_3(corners[i1], corners[i2]) )
        )
      );
  }


  Oriented_side last_os = ON_ORIENTED_BOUNDARY;
  for (int i=0; i<8; ++i)
    if (orientations[i]!=ON_ORIENTED_BOUNDARY)
    {
      if (last_os==ON_ORIENTED_BOUNDARY)
        last_os=orientations[i];
      else
      {
        if(last_os!=orientations[i])
        {
          last_os=ON_ORIENTED_BOUNDARY;
          break;
        }
      }
    }

  // the intersection is the full bbox
  if (last_os!=ON_ORIENTED_BOUNDARY)
    return last_os;

  //add points on negative side and on the plane
  for (int i=0; i<8; ++i)
    if (orientations[i]!=ON_POSITIVE_SIDE)
      points.push_back(corners[i]);

  // take the convex hull of the points on the negative side+intersection points
  // overkill...
  Polyhedron_3<Geom_traits> P;
  CGAL::convex_hull_3(points.begin(), points.end(), P);
  copy_face_graph(P, tm_out,
                  Emptyset_iterator(), Emptyset_iterator(), Emptyset_iterator(),
                  get(vertex_point, P), vpm_out);
  return ON_ORIENTED_BOUNDARY;
}


// convenience overload
template <class TriangleMesh,
          class NamedParameters1>
bool
clip(      TriangleMesh& tm,
     /*const*/ TriangleMesh& clipper,
           bool close,
     const NamedParameters1& np_tm)
{
  return clip(tm, clipper, close, np_tm, parameters::all_default());
}

// convenience overload
template <class TriangleMesh>
bool
clip(      TriangleMesh& tm,
     /*const*/ TriangleMesh& clipper,
           bool close)
{
  return clip(tm, clipper, close, parameters::all_default());
}

// works only with the default point map, for more complex use cases, use
// clip_to_bbox first and the other overload of clip with two meshes
/// \todo document me
template <class TriangleMesh,
          class Plane_3>
void clip(      TriangleMesh& tm,
          const Plane_3& plane,
          bool close)
{
  if( boost::begin(faces(tm))==boost::end(faces(tm)) ) return;
  CGAL::Bbox_3 bbox = ::CGAL::Polygon_mesh_processing::bbox(tm);
  //extend the bbox a bit to avoid border cases
  double xd=(bbox.xmax()-bbox.xmin())/100;
  double yd=(bbox.ymax()-bbox.ymin())/100;
  double zd=(bbox.zmax()-bbox.zmin())/100;
  bbox=CGAL::Bbox_3(bbox.xmin()-xd, bbox.ymin()-yd, bbox.zmin()-zd,
                    bbox.xmax()+xd, bbox.ymax()+yd, bbox.zmax()+zd);
  TriangleMesh clipper;
  Oriented_side os = clip_to_bbox(plane, bbox, clipper, parameters::all_default());

  switch(os)
  {
    case ON_NEGATIVE_SIDE:
      return; // nothing to clip, the full mesh is on the negative side
    case ON_POSITIVE_SIDE:
      clear(tm); // clear the mesh that is fully on the positive side
      return;
    default:
      clip(tm, clipper, close);
  }
}

#endif // !DOXYGEN_RUNNING

} } //end of namespace CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_CLIP_H
