// Copyright (c) 2018-2023 GeometryFactory (France).
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

#ifndef CGAL_POLYGON_MESH_PROCESSING_REMESH_PLANAR_PATCHES_H
#define CGAL_POLYGON_MESH_PROCESSING_REMESH_PLANAR_PATCHES_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/border.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/properties.h>

#include <boost/dynamic_bitset.hpp>
#include <boost/iterator/function_output_iterator.hpp>
#include <boost/container/small_vector.hpp>

#include <algorithm>
#include <unordered_map>

namespace CGAL{

namespace Polygon_mesh_processing {

namespace Planar_segmentation{

template <class PM>
struct Face_map
{
  typedef typename boost::property_traits<PM>::value_type key_type;
  typedef typename boost::property_traits<PM>::key_type value_type;
  typedef value_type reference;
  typedef boost::writable_property_map_tag category;

  Face_map(PM pm, const std::vector<std::size_t>& face_ids)
    : pm(pm)
    , face_ids(face_ids)
  {}

  friend
  void put(Face_map m, key_type k, value_type v)
  {
    put(m.pm, v, static_cast<key_type>(m.face_ids[k]));
  }

  PM pm;
  const std::vector<std::size_t>& face_ids;
};

template <class PM>
struct Vertex_map
{
  typedef typename boost::property_traits<PM>::value_type key_type;
  typedef typename boost::property_traits<PM>::key_type value_type;
  typedef value_type reference;
  typedef boost::writable_property_map_tag category;

  Vertex_map(PM pm) : pm(pm) {}

  friend
  void put(Vertex_map m, key_type k, value_type v)
  {
    put(m.pm, v, k);
  }

  PM pm;
};

template <class TriangleMesh>
struct Default_visitor
{
  void operator()(TriangleMesh&) const {}
};

template <class PolygonMeshOut, class VertexCornerMapOut>
struct Face_index_tracker_base
{
  typedef boost::graph_traits<PolygonMeshOut> GT;

  Face_index_tracker_base(VertexCornerMapOut vertex_corner_map)
    : m_v2v2_map(vertex_corner_map)
  {}

  Vertex_map<VertexCornerMapOut> v2v_map()
  {
    return m_v2v2_map;
  }

  Vertex_map<VertexCornerMapOut> m_v2v2_map;
};

template <class PolygonMeshOut>
struct Face_index_tracker_base<PolygonMeshOut, internal_np::Param_not_found>
{
  typedef boost::graph_traits<PolygonMeshOut> GT;
  using Vmap = Constant_property_map<std::size_t,
                                     typename boost::graph_traits<PolygonMeshOut>::vertex_descriptor>;

  Face_index_tracker_base(internal_np::Param_not_found) {}
  Vmap v2v_map() { return Vmap(); }
};

template <class PolygonMeshOut, class VertexCornerMapOut, class FacePatchMapOut>
struct Face_index_tracker
  : public Face_index_tracker_base<PolygonMeshOut, VertexCornerMapOut>
{
  typedef boost::graph_traits<PolygonMeshOut> GT;

  Face_index_tracker(VertexCornerMapOut vertex_corner_map, FacePatchMapOut face_patch_map)
    : Face_index_tracker_base<PolygonMeshOut, VertexCornerMapOut>(vertex_corner_map)
    , m_f2f_map(face_patch_map, face_ids)
  {}

  std::vector<std::size_t> face_ids;

  void register_faces_of_cc(std::size_t nb_faces, std::size_t i)
  {
    face_ids.insert(face_ids.end(), nb_faces, i);
  }

  Face_map<FacePatchMapOut>
  f2f_map()
  {
    return m_f2f_map;
  }

  Face_map<FacePatchMapOut> m_f2f_map;
};


template <class PolygonMeshOut, class VertexCornerMapOut>
struct Face_index_tracker<PolygonMeshOut, VertexCornerMapOut, internal_np::Param_not_found>
  : public Face_index_tracker_base<PolygonMeshOut, VertexCornerMapOut>
{
  using Fmap = Constant_property_map<std::size_t,
                                     typename boost::graph_traits<PolygonMeshOut>::face_descriptor>;

  Face_index_tracker(VertexCornerMapOut vertex_corner_map,internal_np::Param_not_found)
    : Face_index_tracker_base<PolygonMeshOut, VertexCornerMapOut>(vertex_corner_map)
  {}

  void register_faces_of_cc(std::size_t /*nb_triangles*/, std::size_t /*in_patch_id*/) {}

  Fmap f2f_map() { return Fmap(); }
};

template <class Vector_3, class PatchNormalMap>
void init_face_normals(std::vector<Vector_3>& face_normals,
                       std::size_t nb_patches,
                       PatchNormalMap patch_normal_map)
{
  face_normals.resize(nb_patches);
  for (std::size_t i=0; i<nb_patches; ++i)
    face_normals[i] = get(patch_normal_map, i);
}

template <class Vector_3>
void init_face_normals(std::vector<Vector_3>& face_normals,
                       std::size_t nb_patches,
                       ::CGAL::internal_np::Param_not_found)
{
  face_normals.assign(nb_patches, NULL_VECTOR);
}

inline std::size_t init_id()
{
  return std::size_t(-2);
}

inline std::size_t default_id()
{
  return std::size_t(-1);
}

inline bool is_init_id(std::size_t i)
{
  return i == init_id();
}

inline bool is_corner_id(std::size_t i)
{
  return i < init_id();
}

template <typename Vector_3>
bool is_vector_positive(const Vector_3& normal)
{
  if (normal.x()==0)
  {
    if (normal.y()==0)
      return normal.z() > 0;
    else
      return normal.y() > 0;
  }
  else
    return normal.x() > 0;
}

struct FaceInfo2
{
  FaceInfo2():m_in_domain(-1){}
  int m_in_domain;

  void set_in_domain()   { m_in_domain=1; }
  void set_out_domain()  { m_in_domain=0; }
  bool visited() const   { return m_in_domain!=-1; }
  bool in_domain() const { return m_in_domain==1; }
};

template <typename Kernel,
          typename TriangleMesh,
          typename VertexPointMap,
          typename edge_descriptor>
bool is_edge_between_coplanar_faces(edge_descriptor e,
                                    const TriangleMesh& tm,
                                    double coplanar_cos_threshold,
                                    const VertexPointMap& vpm)
{
  typedef typename boost::property_traits<VertexPointMap>::reference Point_ref_3;
  if (is_border(e, tm)) return false;
  typename boost::graph_traits<TriangleMesh>::halfedge_descriptor
    h =  halfedge(e, tm);
  Point_ref_3 p = get(vpm, source(h, tm) );
  Point_ref_3 q = get(vpm, target(h, tm) );
  Point_ref_3 r = get(vpm, target(next(h, tm), tm) );
  Point_ref_3 s = get(vpm, target(next(opposite(h, tm), tm), tm) );

  if (coplanar_cos_threshold==-1)
    return coplanar(p, q, r, s);
  else
  {
    typename Kernel::Compare_dihedral_angle_3 pred;
    return pred(p, q, r, s, typename Kernel::FT(coplanar_cos_threshold)) == CGAL::LARGER;
  }
}

template <typename Kernel,
          typename TriangleMesh,
          typename VertexPointMap,
          typename halfedge_descriptor,
          typename EdgeIsConstrainedMap>
bool is_target_vertex_a_corner(halfedge_descriptor h,
                               EdgeIsConstrainedMap edge_is_constrained,
                               const TriangleMesh& tm,
                               double coplanar_cos_threshold,
                               const VertexPointMap& vpm)
{
  typedef typename Kernel::Point_3 Point_3;
  typedef typename boost::graph_traits<TriangleMesh> graph_traits;

  halfedge_descriptor h2 = graph_traits::null_halfedge();
  for(halfedge_descriptor h_loop : halfedges_around_target(h, tm))
  {
    if (h_loop==h) continue;
    if (get(edge_is_constrained, edge(h_loop, tm)))
    {
      if (h2 != graph_traits::null_halfedge()) return true;
      h2=h_loop;
    }
  }

  // handle case when the graph of constraints does not contains only cycle
  // (for example when there is a tangency between surfaces and is shared)
  if (h2 == graph_traits::null_halfedge()) return true;

  const Point_3& p = get(vpm, source(h, tm));
  const Point_3& q = get(vpm, target(h, tm));
  const Point_3& r = get(vpm, source(h2, tm));

  if (coplanar_cos_threshold==-1)
    return !collinear(p, q, r);
  else
  {
    typename Kernel::Compare_angle_3 pred;
    return pred(p, q, r, typename Kernel::FT(coplanar_cos_threshold))==CGAL::SMALLER;
  }
}

template <typename Kernel,
          typename TriangleMesh,
          typename EdgeIsConstrainedMap,
          typename VertexPointMap>
void
mark_constrained_edges(
  const TriangleMesh& tm,
  EdgeIsConstrainedMap edge_is_constrained,
  double coplanar_cos_threshold,
  const VertexPointMap& vpm)
{
  for(typename boost::graph_traits<TriangleMesh>::edge_descriptor e : edges(tm))
  {
    if (!get(edge_is_constrained,e))
      if (!is_edge_between_coplanar_faces<Kernel>(e, tm, coplanar_cos_threshold, vpm))
        put(edge_is_constrained, e, true);
  }
}

template <typename Kernel,
          typename TriangleMesh,
          typename VertexPointMap,
          typename EdgeIsConstrainedMap,
          typename VertexCornerIdMap>
std::size_t
mark_corner_vertices(
  const TriangleMesh& tm,
  EdgeIsConstrainedMap& edge_is_constrained,
  VertexCornerIdMap& vertex_corner_id,
  double coplanar_cos_threshold,
  const VertexPointMap& vpm)
{
  typedef boost::graph_traits<TriangleMesh> graph_traits;
  std::size_t corner_id = 0;
  for(typename graph_traits::edge_descriptor e : edges(tm))
  {
    if (get(edge_is_constrained, e))
    {
      put(vertex_corner_id, source(e, tm), init_id());
      put(vertex_corner_id, target(e, tm), init_id());
    }
  }
  for(typename graph_traits::edge_descriptor e : edges(tm))
  {
    if (!get(edge_is_constrained, e)) continue;
    typename graph_traits::halfedge_descriptor h = halfedge(e, tm);

    if (is_init_id(get(vertex_corner_id, target(h, tm))))
    {
      if (is_target_vertex_a_corner<Kernel>(h, edge_is_constrained, tm, coplanar_cos_threshold, vpm))
        put(vertex_corner_id, target(h, tm), corner_id++);
      else
        put(vertex_corner_id, target(h, tm), default_id());
    }
    if (is_init_id(get(vertex_corner_id, source(h, tm))))
    {
      if (is_target_vertex_a_corner<Kernel>(opposite(h, tm), edge_is_constrained, tm, coplanar_cos_threshold, vpm))
        put(vertex_corner_id, source(h, tm), corner_id++);
      else
        put(vertex_corner_id, source(h, tm), default_id());
    }
  }

  return corner_id;
}

template <typename CDT>
void mark_face_triangles(CDT& cdt)
{
  //look for a triangle inside the domain of the face
  typename CDT::Face_handle fh = cdt.infinite_face();
  fh->info().set_out_domain();
  std::vector<typename CDT::Edge> queue;
  for (int i=0; i<3; ++i)
    queue.push_back(typename CDT::Edge(fh, i) );
  while(true)
  {
    typename CDT::Edge e = queue.back();
    queue.pop_back();
    e=cdt.mirror_edge(e);
    if (e.first->info().visited()) continue;
    if (cdt.is_constrained(e))
    {
      queue.clear();
      queue.push_back(e);
      break;
    }
    else
    {
      for(int i=1; i<3; ++i)
      {
        typename CDT::Edge candidate(e.first, (e.second+i)%3);
        if (!candidate.first->neighbor(candidate.second)->info().visited())
          queue.push_back( candidate );
      }
      e.first->info().set_out_domain();
    }
  }
  // now extract triangles inside the face
  while(!queue.empty())
  {
    typename CDT::Edge e = queue.back();
    queue.pop_back();
    if (e.first->info().visited()) continue;
    e.first->info().set_in_domain();

    for(int i=1; i<3; ++i)
    {
      typename CDT::Edge candidate(e.first, (e.second+i)%3);
      if (!cdt.is_constrained(candidate) &&
          !candidate.first->neighbor(candidate.second)->info().visited())
      {
        queue.push_back( cdt.mirror_edge(candidate) );
      }
    }
  }
}

template < typename GT,
           typename Vb = Triangulation_vertex_base_2<GT> >
class Triangulation_vertex_with_corner_id_2
  : public Vb
{
  std::size_t _id = -1;
public:
  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other          Vb2;
    typedef Triangulation_vertex_with_corner_id_2<GT, Vb2>   Other;
  };

  template<class ... T>
  Triangulation_vertex_with_corner_id_2(const T& ... t)
    : Vb(t...)
  {}

  const std::size_t& corner_id() const { return _id; }
        std::size_t& corner_id()       { return _id; }
};


template <typename Kernel>
bool add_triangle_faces(const std::vector< std::pair<std::size_t, std::size_t> >& csts,
                        typename Kernel::Vector_3 normal,
                        const std::vector<typename Kernel::Point_3>& corners,
                        std::vector<boost::container::small_vector<std::size_t, 3> >& out_faces)
{
  typedef Projection_traits_3<Kernel>                            P_traits;
  typedef Triangulation_vertex_with_corner_id_2<P_traits>              Vb;
  typedef Triangulation_face_base_with_info_2<FaceInfo2,P_traits>     Fbb;
  typedef Constrained_triangulation_face_base_2<P_traits,Fbb>          Fb;
  typedef Triangulation_data_structure_2<Vb,Fb>                       TDS;
  typedef No_constraint_intersection_requiring_constructions_tag     Itag;
  typedef Constrained_Delaunay_triangulation_2<P_traits, TDS, Itag>   CDT;
  typedef typename CDT::Vertex_handle                       Vertex_handle;
  typedef typename CDT::Face_handle                           Face_handle;

  typedef typename Kernel::Point_3 Point_3;

  std::size_t expected_nb_pts = csts.size()/2;
  std::vector<std::size_t> corner_ids;
  corner_ids.reserve(expected_nb_pts);

  typedef std::pair<std::size_t, std::size_t> Id_pair;
  for(const Id_pair& p : csts)
  {
    CGAL_assertion(p.first<corners.size());
    corner_ids.push_back(p.first);
  }

  bool reverse_face_orientation = is_vector_positive(normal);
  if (reverse_face_orientation)
    normal=-normal;

  // create cdt and insert points
  P_traits p_traits(normal);
  CDT cdt(p_traits);

  // now do the point insert and info set
  typedef typename Pointer_property_map<Point_3>::const_type Pmap;
  typedef Spatial_sort_traits_adapter_2<P_traits,Pmap> Search_traits;

  spatial_sort(corner_ids.begin(), corner_ids.end(),
               Search_traits(make_property_map(corners),p_traits));

  Vertex_handle v_hint;
  Face_handle hint;
  for (std::size_t corner_id : corner_ids)
  {
    v_hint = cdt.insert(corners[corner_id], hint);
    if (v_hint->corner_id()!=std::size_t(-1) && v_hint->corner_id()!=corner_id)
      return false; // handle case of points being identical upon projection
    v_hint->corner_id()=corner_id;
    hint=v_hint->face();
  }

  // note that nbv might be different from points.size() in case of hole
  // tangent to the principal CCB
  CGAL_assertion_code(std::size_t nbv=cdt.number_of_vertices();)

  // insert constrained edges
  std::unordered_map<std::size_t, typename CDT::Vertex_handle> vertex_map;
  for(typename CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin(),
                                             end = cdt.finite_vertices_end(); vit!=end; ++vit)
  {
    vertex_map[vit->corner_id()]=vit;
  }

  std::vector< std::pair<typename CDT::Vertex_handle, typename CDT::Vertex_handle> > local_csts;
  local_csts.reserve(csts.size());
  try{
    for(const Id_pair& p : csts)
    {
      CGAL_assertion(vertex_map.count(p.first)!=0 && vertex_map.count(p.second)!=0);
      cdt.insert_constraint(vertex_map[p.first], vertex_map[p.second]);
    }
  }catch(typename CDT::Intersection_of_constraints_exception&)
  {
    // intersection of constraints probably due to the projection
    return false;
  }
  CGAL_assertion(cdt.number_of_vertices() == nbv);

  if (cdt.dimension()!=2) return false;

  mark_face_triangles(cdt);

  for (typename CDT::Finite_faces_iterator fit=cdt.finite_faces_begin(),
                                           end=cdt.finite_faces_end(); fit!=end; ++fit)
  {
    if (!fit->info().in_domain()) continue;
    if (cdt.is_infinite(fit)) return false;

    if (reverse_face_orientation)
      out_faces.push_back( {fit->vertex(1)->corner_id(),
                            fit->vertex(0)->corner_id(),
                            fit->vertex(2)->corner_id()} );
    else
      out_faces.push_back( {fit->vertex(0)->corner_id(),
                            fit->vertex(1)->corner_id(),
                            fit->vertex(2)->corner_id()} );
  }

  return true;
}

template <typename Kernel,
          typename TriangleMesh,
          typename VertexCornerIdMap,
          typename EdgeIsConstrainedMap,
          typename FaceCCIdMap,
          typename VertexPointMap>
std::pair<std::size_t, std::size_t>
tag_corners_and_constrained_edges(const TriangleMesh& tm,
                                  double coplanar_cos_threshold,
                                  VertexCornerIdMap& vertex_corner_id,
                                  EdgeIsConstrainedMap& edge_is_constrained,
                                  FaceCCIdMap& face_cc_ids,
                                  const VertexPointMap& vpm)
{
  typedef typename boost::graph_traits<TriangleMesh> graph_traits;
  // mark constrained edges
  mark_constrained_edges<Kernel>(tm, edge_is_constrained, coplanar_cos_threshold, vpm);

  // mark connected components (cc) delimited by constrained edges
  std::size_t nb_cc = Polygon_mesh_processing::connected_components(
    tm, face_cc_ids, parameters::edge_is_constrained_map(edge_is_constrained));

  if (coplanar_cos_threshold!=-1)
  {
    for(typename graph_traits::edge_descriptor e : edges(tm))
    {
      if (get(edge_is_constrained, e) && !is_border(e, tm))
      {
        typename graph_traits::halfedge_descriptor h = halfedge(e, tm);
        if ( get(face_cc_ids, face(h, tm))==get(face_cc_ids, face(opposite(h, tm), tm)) )
          put(edge_is_constrained, e, false);
      }
    }
  }

  std::size_t nb_corners =
    mark_corner_vertices<Kernel>(tm, edge_is_constrained, vertex_corner_id, coplanar_cos_threshold, vpm);

  return std::make_pair(nb_corners, nb_cc);
}

template <typename Kernel,
          typename TriangleMesh,
          typename VertexCornerIdMap,
          typename EdgeIsConstrainedMap,
          typename FaceCCIdMap,
          typename VertexPointMap,
          typename IndexTracking>
bool decimate_impl(const TriangleMesh& tm,
                   std::pair<std::size_t, std::size_t>& nb_corners_and_nb_cc,
                   VertexCornerIdMap& vertex_corner_id,
                   EdgeIsConstrainedMap& edge_is_constrained,
                   FaceCCIdMap& face_cc_ids,
                   const VertexPointMap& vpm,
                   bool do_not_triangulate_faces,
                   std::vector< typename Kernel::Point_3 >& corners,
                   std::vector< boost::container::small_vector<std::size_t,3> >& out_faces,
                   IndexTracking& f_id_tracker,
                   std::vector< typename Kernel::Vector_3 >& face_normals)
{
  typedef typename boost::graph_traits<TriangleMesh> graph_traits;
  typedef typename graph_traits::halfedge_descriptor halfedge_descriptor;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits::face_descriptor face_descriptor;
  typedef std::pair<std::size_t, std::size_t> Id_pair;
  typedef typename boost::property_traits<FaceCCIdMap>::value_type PID;

  // compute the new mesh
  std::vector< std::vector< boost::container::small_vector<std::size_t,3> > > faces_per_cc(nb_corners_and_nb_cc.second);
  boost::dynamic_bitset<> cc_to_handle(nb_corners_and_nb_cc.second);
  cc_to_handle.set();

  bool all_patches_successfully_remeshed = true;
  do
  {
    std::vector< std::vector<Id_pair> > face_boundaries(nb_corners_and_nb_cc.second);
    std::vector<bool> face_boundaries_valid(nb_corners_and_nb_cc.second, true);

    std::vector<vertex_descriptor> corner_id_to_vd(nb_corners_and_nb_cc.first, graph_traits::null_vertex());
    std::vector<bool> duplicated_corners(nb_corners_and_nb_cc.first, false);
    auto check_corner = [&corner_id_to_vd, &duplicated_corners](std::size_t corner_id, vertex_descriptor vd)
    {
      if (corner_id_to_vd[corner_id]!=graph_traits::null_vertex() && corner_id_to_vd[corner_id]!=vd)
        duplicated_corners[corner_id]=true;
    };

    // collect maximal constrained edges per cc
    for(halfedge_descriptor h : halfedges(tm))
    {
      if (!get(edge_is_constrained, edge(h, tm)) || is_border(h, tm)) continue;

      std::size_t i1 = get(vertex_corner_id, source(h, tm));
      if ( is_corner_id(i1) )
      {
        check_corner(i1, source(h, tm));
        halfedge_descriptor h_init = h;
        std::size_t cc_id = get(face_cc_ids, face(h_init, tm));
        if (!cc_to_handle.test(cc_id)) continue;
        do{
          std::size_t i2 = get(vertex_corner_id, target(h_init, tm));
          if ( is_corner_id(i2) )
          {
            check_corner(i2, target(h_init, tm));
            face_boundaries[ cc_id ].push_back( Id_pair(i1,i2) );
            if (face_normals[ cc_id ] == NULL_VECTOR)
            {
              face_normals[ cc_id ] = normal(get(vpm, source(h, tm)),
                                             get(vpm, target(h, tm)),
                                             get(vpm, target(next(h, tm), tm)));
            }
            break;
          }

          do{
            h_init=opposite(next(h_init, tm), tm);
          } while( !get(edge_is_constrained, edge(h_init, tm)) );
          h_init=opposite(h_init, tm);
        }
        while(true);
      }
    }

    for (std::size_t cc_id = cc_to_handle.find_first();
                         cc_id < cc_to_handle.npos;
                         cc_id = cc_to_handle.find_next(cc_id))
    {
      std::vector< boost::container::small_vector<std::size_t,3> >& cc_faces = faces_per_cc[cc_id];
      cc_faces.clear();

      std::vector< Id_pair >& csts = face_boundaries[cc_id];

      if (!face_boundaries_valid[cc_id]) continue;

      // do not remesh a patch containing duplicated vertices
      for (auto c : csts)
        if (duplicated_corners[c.first] || duplicated_corners[c.second])
        {
          csts.clear(); // this will trigger the copy of the current patch rather than a remeshing
          break;
        }
#ifdef CGAL_DEBUG_DECIMATION
      std::cout << "csts.size() " <<  csts.size() << "\n";
#endif

      if (csts.size()==3)
      {
        cc_faces.push_back( { csts[0].first,
                              csts[0].second,
                              csts[0].first==csts[1].first ||
                              csts[0].second==csts[1].first ?
                              csts[1].second:csts[1].first} );
        cc_to_handle.set(cc_id, 0);
      }
      else
      {
        if (csts.size() > 3 && do_not_triangulate_faces)
        {
          // TODO this is not optimal at all since we already have the set of contraints,
          //      we could work on the graph on constraint and recover only the orientation
          //      of the edge. To be done if someone find it too slow.
          std::vector<halfedge_descriptor> hborders;
          CGAL::Face_filtered_graph<TriangleMesh> ffg(tm, static_cast<PID>(cc_id), face_cc_ids);
          extract_boundary_cycles(ffg, std::back_inserter(hborders));

          if (hborders.size()==1)
          {
            cc_faces.resize(1);
            for (halfedge_descriptor h : halfedges_around_face(hborders[0], ffg))
            {
              std::size_t cid = get(vertex_corner_id, target(h, tm));
              if (is_corner_id(cid))
                cc_faces.back().push_back(cid);
            }
            std::reverse(cc_faces.back().begin(), cc_faces.back().end());
            cc_to_handle.set(cc_id, 0);
            continue;
          }
        }

        if (csts.size() > 3 && add_triangle_faces<Kernel>(csts, face_normals[cc_id], corners, cc_faces))
          cc_to_handle.set(cc_id, 0);
        else
        {
          //TODO: shall we try to plug pseudo-cdt?
#ifdef CGAL_DEBUG_DECIMATION
          static int fail_case_id=0;
          std::cout << "  DEBUG: Failed to remesh a patch, case #" << fail_case_id  << std::endl;
          std::ofstream debug("failed_remesh_"+std::to_string(fail_case_id)+".polylines.txt");
          debug << std::setprecision(17);
          for (auto c : csts)
            debug << "2 " << corners[c.first] << " " << corners[c.second] << "\n";
          debug.close();
          std::cout << "  normal used is " << face_normals[cc_id] << "\n";
          debug.open("normal"+std::to_string(fail_case_id)+".polylines.txt");
          debug << "2 " << corners[csts[0].first] << " " << corners[csts[0].first]+face_normals[cc_id] << "\n";
          debug.close();
          ++fail_case_id;
#endif
          all_patches_successfully_remeshed = false;
          // make all vertices of the patch a corner
          CGAL::Face_filtered_graph<TriangleMesh> ffg(tm, static_cast<PID>(cc_id), face_cc_ids);
          std::vector<vertex_descriptor> new_corners;
          for (vertex_descriptor v : vertices(ffg))
          {
            std::size_t i = get(vertex_corner_id, v);
            if ( !is_corner_id(i) )
            {
              i = nb_corners_and_nb_cc.first++;
              put(vertex_corner_id, v, i);
              corners.push_back(get(vpm, v));
              new_corners.push_back(v);
            }
          }
          // add all the faces of the current patch
          for (face_descriptor f : faces(ffg))
          {
            halfedge_descriptor h = halfedge(f, tm);
            cc_faces.push_back({ get(vertex_corner_id, source(h,tm)),
                                 get(vertex_corner_id, target(h,tm)),
                                  get(vertex_corner_id, target(next(h,tm), tm)) });
          }
          // reset flag for neighbor connected components only if interface has changed
          for (vertex_descriptor v : new_corners)
          {
            for (halfedge_descriptor h : halfedges_around_target(halfedge(v, tm), tm))
            {
              if (!is_border(h, tm))
              {
                std::size_t other_cc_id = get(face_cc_ids, face(h, tm));
                cc_to_handle.set(other_cc_id, 1);
                face_boundaries_valid[ other_cc_id ]=false;
              }
            }
          }
          cc_to_handle.set(cc_id, 0);
        }
      }
    }
  }
  while(cc_to_handle.any());

  std::size_t cc_id=0;
  for (const std::vector<boost::container::small_vector<std::size_t,3>>& cc_trs : faces_per_cc)
  {
    out_faces.insert(out_faces.end(), cc_trs.begin(), cc_trs.end());
    f_id_tracker.register_faces_of_cc(cc_trs.size(), cc_id++);
  }

  return all_patches_successfully_remeshed;
}

template <typename Kernel,
          typename TriangleMeshIn,
          typename PolygonMeshOut,
          typename VertexCornerIdMap,
          typename EdgeIsConstrainedMap,
          typename FaceCCIdMap,
          typename VertexPointMapIn,
          typename VertexPointMapOut,
          typename VertexCornerMapOut,
          typename FacePatchMapOut,
          typename Visitor>
bool decimate_impl(const TriangleMeshIn& tm_in,
                   PolygonMeshOut& pm_out,
                   std::pair<std::size_t, std::size_t> nb_corners_and_nb_cc,
                   VertexCornerIdMap& vertex_corner_id,
                   EdgeIsConstrainedMap& edge_is_constrained,
                   FaceCCIdMap& face_cc_ids,
                   const VertexPointMapIn& vpm_in,
                   const VertexPointMapOut& vpm_out,
                   bool do_not_triangulate_faces,
                   VertexCornerMapOut vcorner_map_out,
                   FacePatchMapOut fpatch_map_out,
                   Visitor& visitor,
                   std::vector<typename Kernel::Vector_3>& face_normals)
{
  typedef typename boost::graph_traits<TriangleMeshIn> graph_traits;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename Kernel::Point_3 Point_3;

  Face_index_tracker<PolygonMeshOut, VertexCornerMapOut, FacePatchMapOut>
    f_id_tracker(vcorner_map_out, fpatch_map_out);

  //collect corners
  std::vector< Point_3 > corners(nb_corners_and_nb_cc.first);
  for(vertex_descriptor v : vertices(tm_in))
  {
    std::size_t i = get(vertex_corner_id, v);
    if ( is_corner_id(i) )
    {
      corners[i]=get(vpm_in, v);
    }
  }

  std::vector< boost::container::small_vector<std::size_t,3> > faces;
  bool remeshing_failed = decimate_impl<Kernel>(tm_in,
                                                nb_corners_and_nb_cc,
                                                vertex_corner_id,
                                                edge_is_constrained,
                                                face_cc_ids,
                                                vpm_in,
                                                do_not_triangulate_faces,
                                                corners,
                                                faces,
                                                f_id_tracker,
                                                face_normals);

  if (!is_polygon_soup_a_polygon_mesh(faces))
  {
    return false;
  }

  visitor(pm_out);
  polygon_soup_to_polygon_mesh(corners, faces, pm_out,
                               parameters::point_to_vertex_map(f_id_tracker.v2v_map()).
                                           polygon_to_face_map(f_id_tracker.f2f_map()),
                               parameters::vertex_point_map(vpm_out));
  return remeshing_failed;
}

template <typename vertex_descriptor,
          typename Point_3,
          typename OutputIterator>
void extract_meshes_containing_a_point(
  const Point_3& pt,
  const std::map<Point_3, std::multimap<std::size_t, vertex_descriptor> >& point_to_vertex_maps,
  OutputIterator out)
{
  typedef std::pair<const std::size_t, vertex_descriptor> Pair_type;
  for(const Pair_type& p : point_to_vertex_maps.find(pt)->second)
    *out++=p.first;
}

template <typename TriangleMesh,
          typename Point_3,
          typename vertex_descriptor,
          typename EdgeIsConstrainedMap,
          typename VertexIsSharedMap,
          typename VertexPointMap>
void mark_boundary_of_shared_patches_as_constrained_edges(
  std::vector<TriangleMesh*>& mesh_ptrs,
  std::map<Point_3, std::multimap<std::size_t, vertex_descriptor> >& point_to_vertex_maps,
  std::vector<EdgeIsConstrainedMap>& edge_is_constrained_maps,
  std::vector<VertexIsSharedMap>& vertex_shared_maps,
  const std::vector<VertexPointMap>& vpms)
{
  typedef boost::graph_traits<TriangleMesh> graph_traits;
  typedef typename graph_traits::edge_descriptor edge_descriptor;
  typedef typename graph_traits::halfedge_descriptor halfedge_descriptor;

  std::size_t mesh_id = 0;
  for(TriangleMesh* tm_ptr : mesh_ptrs)
  {
    TriangleMesh& tm=*tm_ptr;
    EdgeIsConstrainedMap& edge_is_constrained = edge_is_constrained_maps[mesh_id];
    VertexIsSharedMap& is_vertex_shared = vertex_shared_maps[mesh_id];

    for(edge_descriptor e : edges(tm))
    {
      if (is_border(e, tm)) continue; //border edges will be automatically marked as constrained

      halfedge_descriptor h = halfedge(e, tm);
      vertex_descriptor src = source(h, tm), tgt = target(h, tm);
      if (get(is_vertex_shared, src) && get(is_vertex_shared, tgt))
      {
        //extract the set of meshes having both vertices
        std::set<std::size_t> src_set, tgt_set, inter_set;
        extract_meshes_containing_a_point(get(vpms[mesh_id], src),
                                          point_to_vertex_maps,
                                          std::inserter(src_set, src_set.begin()));
        extract_meshes_containing_a_point(get(vpms[mesh_id], tgt),
                                          point_to_vertex_maps,
                                          std::inserter(tgt_set, tgt_set.begin()));

        std::set_intersection(src_set.begin(), src_set.end(),
                              tgt_set.begin(), tgt_set.end(),
                              std::inserter(inter_set, inter_set.begin()));

        std::multimap<std::size_t, vertex_descriptor>& mesh_to_vertex_src =
          point_to_vertex_maps[get(vpms[mesh_id], src)];
        std::multimap<std::size_t, vertex_descriptor>& mesh_to_vertex_tgt =
          point_to_vertex_maps[get(vpms[mesh_id], tgt)];

        std::set<Point_3> incident_face_points;
        incident_face_points.insert(get(vpms[mesh_id], target(next(h, tm), tm)));
        h=opposite(h, tm);
        incident_face_points.insert(get(vpms[mesh_id], target(next(h, tm), tm)));

        // we mark as constrained edge, any edge that is shared between more than 2 meshes
        // such that at least one of the two incident faces to the edge are not present in
        // all the meshes containing the edge
        for(std::size_t other_mesh_id : inter_set)
        {
          TriangleMesh* other_tm_ptr = mesh_ptrs[other_mesh_id];
          if (other_tm_ptr==&tm) continue;

          std::vector<vertex_descriptor> srcs, tgts;
          auto it = mesh_to_vertex_src.find(other_mesh_id);
          while (it!=mesh_to_vertex_src.end() && it->first==other_mesh_id)
            srcs.push_back(it++->second);
          it = mesh_to_vertex_tgt.find(other_mesh_id);
          while (it!=mesh_to_vertex_tgt.end() && it->first==other_mesh_id)
            tgts.push_back(it++->second);

          for (vertex_descriptor other_src : srcs)
          {
            for (vertex_descriptor other_tgt : tgts)
            {
              std::pair<halfedge_descriptor, bool> hres = halfedge(other_src, other_tgt, *other_tm_ptr);
              if (hres.second)
              {
                if (is_border_edge(hres.first, *other_tm_ptr))
                {
                  put(edge_is_constrained, e, true);
                  break;
                }
                if (incident_face_points.count(
                      get(vpms[other_mesh_id], target(next(hres.first, *other_tm_ptr), *other_tm_ptr)))==0)
                {
                  put(edge_is_constrained, e, true);
                  break;
                }
                hres.first=opposite(hres.first, *other_tm_ptr);
                if (incident_face_points.count(
                      get(vpms[other_mesh_id], target(next(hres.first, *other_tm_ptr), *other_tm_ptr)))==0)
                {
                  put(edge_is_constrained, e, true);
                  break;
                }
              }
            }
          }
        }
      }
    }
    ++mesh_id;
  }
}

template<typename Point_3,
         typename vertex_descriptor,
         typename VertexIsSharedMap>
void propagate_corner_status(
  std::vector<VertexIsSharedMap>& vertex_corner_id_maps,
  std::map<Point_3, std::multimap<std::size_t, vertex_descriptor> >& point_to_vertex_maps,
  std::vector< std::pair<std::size_t, std::size_t> >& nb_corners_and_nb_cc_all)
{
  typedef std::pair<const Point_3, std::multimap<std::size_t, vertex_descriptor> > Pair_type;
  for(Pair_type& p : point_to_vertex_maps)
  {
    // if one vertex is a corner, all should be
    typedef std::pair<const std::size_t, vertex_descriptor> Map_pair_type;
    bool is_corner=false;
    for (Map_pair_type& mp : p.second)
    {
      std::size_t mesh_id = mp.first;
      if ( is_corner_id( get(vertex_corner_id_maps[mesh_id], mp.second) ))
      {
        is_corner=true;
        break;
      }
    }
    if (is_corner)
    {
      for(Map_pair_type& mp : p.second)
      {
        std::size_t mesh_id = mp.first;
        if ( !is_corner_id(get(vertex_corner_id_maps[mesh_id], mp.second)) )
        {
          put(vertex_corner_id_maps[mesh_id], mp.second,
            nb_corners_and_nb_cc_all[mesh_id].first++);
        }
      }
    }
  }
}

template <typename Kernel,
          typename TriangleMeshRange,
          typename MeshMap,
          typename VertexPointMap,
          typename TagFunction>
bool decimate_meshes_with_common_interfaces_impl(TriangleMeshRange& meshes,
                                                 MeshMap mesh_map,
                                                 const TagFunction& tag_function,
                                                 const std::vector<VertexPointMap>& vpms,
                                                 bool do_not_triangulate_faces)
{
  typedef typename boost::property_traits<MeshMap>::value_type Triangle_mesh;
  typedef typename std::iterator_traits<typename TriangleMeshRange::iterator>::value_type Mesh_descriptor;
  typedef typename boost::property_traits<VertexPointMap>::value_type Point_3;
  typedef typename boost::graph_traits<Triangle_mesh> graph_traits;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits::edge_descriptor edge_descriptor;
  typedef typename graph_traits::face_descriptor face_descriptor;
  typedef typename graph_traits::halfedge_descriptor halfedge_descriptor;

  // declare and init all property maps
  typedef typename boost::property_map<Triangle_mesh, CGAL::dynamic_face_property_t<std::size_t> >::type Face_cc_map;
  typedef typename boost::property_map<Triangle_mesh, CGAL::dynamic_edge_property_t<bool> >::type Edge_is_constrained_map;
  typedef typename boost::property_map<Triangle_mesh, CGAL::dynamic_vertex_property_t<std::size_t> >::type Vertex_corner_id_map;
  typedef typename boost::property_map<Triangle_mesh, CGAL::dynamic_vertex_property_t<bool> >::type Vertex_is_shared_map;

  std::vector<Vertex_is_shared_map > vertex_shared_maps;
  std::vector<Edge_is_constrained_map> edge_is_constrained_maps;
  std::vector<Vertex_corner_id_map> vertex_corner_id_maps;
  std::vector<Face_cc_map> face_cc_ids_maps;
  const std::size_t nb_meshes = meshes.size();
  std::vector<bool> mesh_has_non_manifold_vertices(nb_meshes, false);
  vertex_corner_id_maps.reserve(nb_meshes);
  vertex_shared_maps.reserve(nb_meshes);
  edge_is_constrained_maps.reserve(nb_meshes);
  face_cc_ids_maps.reserve(nb_meshes);

  std::vector<Triangle_mesh*> mesh_ptrs;
  mesh_ptrs.reserve(nb_meshes);
  for(Mesh_descriptor& md : meshes)
    mesh_ptrs.push_back( &(mesh_map[md]) );

  auto has_non_manifold_vertices = [](const Triangle_mesh& tm)
  {
    std::vector<halfedge_descriptor> nmvs;
    non_manifold_vertices(tm, std::back_inserter(nmvs));
    return !nmvs.empty();
  };

  std::size_t mesh_id=0;
  for(Triangle_mesh* tm_ptr : mesh_ptrs)
  {
    Triangle_mesh& tm = *tm_ptr;
    mesh_has_non_manifold_vertices[mesh_id] = has_non_manifold_vertices(tm);
    vertex_shared_maps.push_back( get(CGAL::dynamic_vertex_property_t<bool>(), tm) );
    for(vertex_descriptor v : vertices(tm))
      put(vertex_shared_maps.back(), v, false);

    edge_is_constrained_maps.push_back( get(CGAL::dynamic_edge_property_t<bool>(), tm) );
    for(edge_descriptor e : edges(tm))
      put(edge_is_constrained_maps.back(), e, false);

    vertex_corner_id_maps.push_back( get(CGAL::dynamic_vertex_property_t<std::size_t>(), tm) );
    for(vertex_descriptor v : vertices(tm))
      put(vertex_corner_id_maps.back(), v, Planar_segmentation::init_id());

    face_cc_ids_maps.push_back( get(CGAL::dynamic_face_property_t<std::size_t>(), tm) );
    for(face_descriptor f : faces(tm))
      put(face_cc_ids_maps.back(), f, -1);
    ++mesh_id;
  }

  std::map<Point_3, std::multimap<std::size_t, vertex_descriptor> > point_to_vertex_maps;

  //start by detecting and marking all shared vertices
  mesh_id = 0;
  for(Triangle_mesh* tm_ptr : mesh_ptrs)
  {
    Triangle_mesh& tm = *tm_ptr;

    for(vertex_descriptor v : vertices(tm))
    {
      std::multimap<std::size_t, vertex_descriptor>& mesh_id_to_vertex =
        point_to_vertex_maps[get(vpms[mesh_id], v)];
      if (!mesh_id_to_vertex.empty())
        put(vertex_shared_maps[mesh_id], v, true);
      if (mesh_id_to_vertex.size()==1)
      {
        std::pair<std::size_t, vertex_descriptor> other=*mesh_id_to_vertex.begin();
        put(vertex_shared_maps[other.first], other.second, true);
      }
      mesh_id_to_vertex.insert( std::make_pair(mesh_id, v) );
    }

    ++mesh_id;
  }

  //then detect edge on the boundary of shared patches and mark them as constrained
  mark_boundary_of_shared_patches_as_constrained_edges(mesh_ptrs, point_to_vertex_maps, edge_is_constrained_maps, vertex_shared_maps, vpms);

  // first tag corners and constrained edges
  std::vector< std::pair<std::size_t, std::size_t> > nb_corners_and_nb_cc_all(nb_meshes);
  mesh_id=0;
  for(Triangle_mesh* tm_ptr : mesh_ptrs)
  {
    Triangle_mesh& tm = *tm_ptr;

    //reset face cc ids as it was set by coplanarity_segmentation_with_pca
    for(face_descriptor f : faces(tm))
      put(face_cc_ids_maps[mesh_id], f, -1);

    if (!mesh_has_non_manifold_vertices[mesh_id])
      nb_corners_and_nb_cc_all[mesh_id] = tag_function(tm,
                                                       vertex_corner_id_maps[mesh_id],
                                                       edge_is_constrained_maps[mesh_id],
                                                       face_cc_ids_maps[mesh_id],
                                                       vpms[mesh_id]);
    else
    {
      nb_corners_and_nb_cc_all[mesh_id]={0,1};
      for (vertex_descriptor vd : vertices(tm))
      {
        if (get(vertex_shared_maps[mesh_id], vd))
        {
          put(vertex_corner_id_maps[mesh_id], vd, true);
          ++nb_corners_and_nb_cc_all[mesh_id].first;
        }
      }
    }
    ++mesh_id;
  }

  // extra step to propagate is_corner to all meshes to make sure shared vertices are kept
  propagate_corner_status(vertex_corner_id_maps, point_to_vertex_maps, nb_corners_and_nb_cc_all);

  // TODO: make identical patches normal identical (up to the sign). Needed only in the approximate case

// now call the decimation
  // storage of all new triangles and all corners
  std::vector< std::vector< Point_3 > > all_corners(nb_meshes);
  std::vector< std::vector< boost::container::small_vector<std::size_t,3> > > all_faces(nb_meshes);
  bool res = true;
  std::vector<bool> to_be_processed(nb_meshes, true);
  bool loop_again;
//  bool no_remeshing_issue = true;
  do{
    loop_again = false;
    for(std::size_t mesh_id=0; mesh_id<nb_meshes; ++mesh_id)
    {
      if (!to_be_processed[mesh_id]) continue;
#ifdef CGAL_DEBUG_DECIMATION
      std::cout << "Handling mesh #" << mesh_id << "\n";
      if (mesh_has_non_manifold_vertices[mesh_id])
        std::cout << "  mesh has non-manifold vertices, mesh is kept as is.\n";
#endif
      if (mesh_has_non_manifold_vertices[mesh_id])
      {
        to_be_processed[mesh_id]=false;
        continue;
      }
      all_faces[mesh_id].clear();
      Triangle_mesh& tm = *mesh_ptrs[mesh_id];

      //collect corners
      std::vector< Point_3 >& corners = all_corners[mesh_id];
      if (corners.empty())
      {
        corners.resize(nb_corners_and_nb_cc_all[mesh_id].first);
        for(vertex_descriptor v : vertices(tm))
        {
          std::size_t i = get(vertex_corner_id_maps[mesh_id], v);
          if ( is_corner_id(i) )
            corners[i]=get(vpms[mesh_id], v);
        }
      }
      std::size_t ncid=corners.size();

      typedef internal_np::Param_not_found PNF;
      PNF pnf;
      Face_index_tracker<Triangle_mesh, PNF, PNF> tracker(pnf, pnf);
      std::vector< typename Kernel::Vector_3 > face_normals(nb_corners_and_nb_cc_all[mesh_id].second, NULL_VECTOR);
      bool all_patches_successfully_remeshed =
        decimate_impl<Kernel>(tm,
                              nb_corners_and_nb_cc_all[mesh_id],
                              vertex_corner_id_maps[mesh_id],
                              edge_is_constrained_maps[mesh_id],
                              face_cc_ids_maps[mesh_id],
                              vpms[mesh_id],
                              do_not_triangulate_faces,
                              corners,
                              all_faces[mesh_id],
                              tracker,
                              face_normals) &&
        is_polygon_soup_a_polygon_mesh(all_faces[mesh_id]);
#ifdef CGAL_DEBUG_DECIMATION
      std::cout << "all_patches_successfully_remeshed? " << all_patches_successfully_remeshed << "\n";
#endif
      if (!all_patches_successfully_remeshed)
      {
//        no_remeshing_issue=false;
        // iterate over points newly marked as corners
        std::set<std::size_t> mesh_ids;
        for (std::size_t cid=ncid; cid<corners.size(); ++cid)
        {
          typedef std::pair<const std::size_t, vertex_descriptor> Map_pair_type;
          auto find_res = point_to_vertex_maps.find(corners[cid]);
          assert(find_res != point_to_vertex_maps.end());
          for(Map_pair_type& mp : find_res->second)
          {
            std::size_t other_mesh_id = mp.first;
            if ( other_mesh_id!=mesh_id  && !is_corner_id(get(vertex_corner_id_maps[other_mesh_id], mp.second)))
            {
              mesh_ids.insert(other_mesh_id);
              put(vertex_corner_id_maps[other_mesh_id], mp.second,
                nb_corners_and_nb_cc_all[other_mesh_id].first++);
              if (!all_corners[other_mesh_id].empty())
                all_corners[other_mesh_id].push_back(corners[cid]);
            }
          }
        }
        for (std::size_t mid : mesh_ids)
          if (!to_be_processed[mid])
          {
#ifdef CGAL_DEBUG_DECIMATION
            if (!loop_again)
              std::cout << "setting for another loop\n";
#endif
            loop_again=true;
            to_be_processed[mid] = true;
          }
      }
      to_be_processed[mesh_id] = false;
    }
  }
  while(loop_again);

  // now create the new meshes:
  for(std::size_t mesh_id=0; mesh_id<nb_meshes; ++mesh_id)
  {
    Triangle_mesh& tm = *mesh_ptrs[mesh_id];
    if (mesh_has_non_manifold_vertices[mesh_id])
    {
//      no_remeshing_issue = false;
      continue;
    }
    CGAL_assertion(is_polygon_soup_a_polygon_mesh(all_faces[mesh_id]));

    //clear(tm);
    tm.clear_without_removing_property_maps();
    polygon_soup_to_polygon_mesh(all_corners[mesh_id], all_faces[mesh_id],
                                 tm, parameters::default_values(), parameters::vertex_point_map(vpms[mesh_id]));
  }

  return res;
}

template <typename Kernel,
          typename TriangleMeshRange,
          typename MeshMap,
          typename VertexPointMap>
bool decimate_meshes_with_common_interfaces_impl(TriangleMeshRange& meshes,
                                                 MeshMap mesh_map,
                                                 double coplanar_cos_threshold,
                                                 const std::vector<VertexPointMap>& vpms,
                                                 bool do_not_triangulate_faces)
{
  typedef typename boost::property_traits<MeshMap>::value_type Triangle_mesh;
  auto tag_function = [coplanar_cos_threshold](Triangle_mesh& tm,
                                               typename boost::property_map<Triangle_mesh, CGAL::dynamic_vertex_property_t<std::size_t> >::type vertex_corner_id,
                                               typename boost::property_map<Triangle_mesh, CGAL::dynamic_edge_property_t<bool> >::type edge_is_constrained,
                                               typename boost::property_map<Triangle_mesh, CGAL::dynamic_face_property_t<std::size_t> >::type face_cc_ids,
                                               VertexPointMap vpm)
  {
    return tag_corners_and_constrained_edges<Kernel>(tm,
                                                     coplanar_cos_threshold,
                                                     vertex_corner_id,
                                                     edge_is_constrained,
                                                     face_cc_ids,
                                                     vpm);
  };

  return decimate_meshes_with_common_interfaces_impl<Kernel>(meshes,
                                                             mesh_map,
                                                             tag_function,
                                                             vpms,
                                                             do_not_triangulate_faces);

}

} //end of namespace Planar_segmentation


/*!
 *  \ingroup PMP_meshing_grp
 *  generates a new triangle mesh `pm_out` with the minimal number of triangles while preserving the shape of `tm_in`.
 *  In practice, this means that connected components of edge-connected faces belonging to the same plane are
 *  first extracted (each such connected component is called a *patch*). Then, the connected components of
 *  vertex-connected patch border edges belonging to the same line are extracted. Endpoints of such components and
 *  vertices incident to more than two patches (or two patches + one mesh boundary) are called *corners*.
 *  `pm_out` contains the 2D constrained Delaunay triangulation of each patch using only corner vertices
 *  on the boundary of the patch.
 *
 *  \warning if `tm_in` contains a non-manifold vertex, `pm_out` will be empty. Those vertices must be
 *           duplicated with `duplicate_non_manifold_vertices()` to get an output.
 *
 *  \tparam TriangleMeshIn a model of `HalfedgeListGraph` and `FaceListGraph`
 *  \tparam PolygonMeshOut a model of `MutableFaceGraph`
 *  \tparam NamedParametersIn a sequence of \ref bgl_namedparameters "Named Parameters"
 *  \tparam NamedParametersOut a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 *  \param tm_in input triangle mesh
 *  \param pm_out output polygon mesh
 *  \param np_in an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
 *
 *  \cgalNamedParamsBegin
 *    \cgalParamNBegin{vertex_point_map}
 *      \cgalParamDescription{a property map associating points to the vertices of `tm_in`}
 *      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMeshIn>::%vertex_descriptor`
 *                     as key type and `GeomTraits::Point_3` as value type, `GeomTraits` being the type of the parameter `geom_traits`}
 *      \cgalParamDefault{`boost::get(CGAL::vertex_point, tm_in)`}
 *      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                      must be available in `TriangleMeshIn`.}
 *    \cgalParamNEnd
 *    \cgalParamNBegin{geom_traits}
 *      \cgalParamDescription{an instance of a geometric traits class}
 *      \cgalParamType{a class model of `Kernel`}
 *      \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *      \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *    \cgalParamNEnd
 *    \cgalParamNBegin{edge_is_constrained_map}
 *      \cgalParamDescription{a property map filled by this function and that will contain `true` if an edge is on the border of a patch and `false` otherwise.}
 *      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMeshIn>::%edge_descriptor`
 *                     as key type and `bool` as value type}
 *      \cgalParamDefault{None}
 *    \cgalParamNEnd
 *    \cgalParamNBegin{face_patch_map}
 *      \cgalParamDescription{a property map filled by this function and that will contain for each face the id
 *                            of its patch in the range `[0, number of patches - 1]`}
 *      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMeshIn>::%face_descriptor`
 *                     as key type and `std::size_t` as value type}
 *      \cgalParamDefault{None}
 *    \cgalParamNEnd
 *    \cgalParamNBegin{vertex_corner_map}
 *      \cgalParamDescription{a property map filled by this function and that will contain for each vertex that is a corner
 *                            an id in the range `[0, number of corners - 1]`, and `std::size_t(-1)` otherwise.}
 *      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMeshIn>::%vertex_descriptor`
 *                     as key type and `std::size_t` as value type}
 *      \cgalParamDefault{None}
 *   \cgalParamNEnd
 *    \cgalParamNBegin{cosine_of_maximum_angle}
 *      \cgalParamDescription{The maximum angle, given as a cosine,
 *                            (i) between the normals of the supporting planes of adjacent faces such that they are considered coplanar, and
 *                            (ii) for the smallest angle between the supporting line of a segment and an adjacent segment such that they are considered collinear.}
 *      \cgalParamType{`FT` type from the `geom_traits` parameter}
 *      \cgalParamDefault{1, which means exact coplanarity and collinearity}
 *      \cgalParamExtra{The value must be in the interval `[0,1]`}
 *   \cgalParamNEnd
 *  \cgalNamedParamsEnd
 *
 *  \param np_out an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
 *
 *  \cgalNamedParamsBegin
 *    \cgalParamNBegin{do_not_triangulate_faces}
 *      \cgalParamDescription{if `true`, faces of `pm_out` will not be triangulated, but the one with more than one connected component of the boundary.}
 *      \cgalParamType{`bool`}
 *      \cgalParamDefault{false}
 *    \cgalParamNEnd
 *    \cgalParamNBegin{vertex_point_map}
 *      \cgalParamDescription{a property map associating points to the vertices of `pm_out`}
 *      \cgalParamType{a class model of `WritablePropertyMap` with `boost::graph_traits<PolygonMeshOut>::%vertex_descriptor`
 *                     as key type and `GeomTraits::Point_3` as value type, `GeomTraits` being the type of the parameter `geom_traits`}
 *      \cgalParamDefault{`boost::get(CGAL::vertex_point, pm_out)`}
 *      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                      must be available in `PolygonMeshOut`.}
 *    \cgalParamNEnd
 *    \cgalParamNBegin{face_patch_map}
 *      \cgalParamDescription{a property map filled by this function and that will contain for each face the id
 *                            of its patch in the range `[0, number of patches - 1]`,
 *                            the patch id of two identical patches in the input and output meshes being equal.}
 *      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMeshOut>::%face_descriptor`
 *                     as key type and `std::size_t` as value type}
 *      \cgalParamDefault{None}
 *    \cgalParamNEnd
 *    \cgalParamNBegin{vertex_corner_map}
 *      \cgalParamDescription{a property map filled by this function and that will contain for each vertex its corner
 *                            an id in the range `[0, number of corners - 1]`}
 *      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMeshOut>::%vertex_descriptor`
 *                     as key type and `std::size_t` as value type}
 *      \cgalParamDefault{None}
 *    \cgalParamNEnd
 *    \cgalParamNBegin{visitor}
 *      \cgalParamDescription{a callable with `visitor(pm_out)` being called once `tm_in` is no longer needed
 *                            and before `pm_out` starts being built. It should be used in the case when `tm_in` and `pm_out` are the same mesh,
 *                            so that `pm_out` can be cleared before being filled.}
 *      \cgalParamType{`visitor(pm_out)` must be a valid expression.}
 *      \cgalParamDefault{None}
 *    \cgalParamNEnd
 *  \cgalNamedParamsEnd
 */
template <typename TriangleMeshIn,
          typename PolygonMeshOut,
          typename NamedParametersIn = parameters::Default_named_parameters,
          typename NamedParametersOut = parameters::Default_named_parameters>
void remesh_planar_patches(const TriangleMeshIn& tm_in,
                                 PolygonMeshOut& pm_out,
                           const NamedParametersIn& np_in = parameters::default_values(),
                           const NamedParametersOut& np_out = parameters::default_values())
{
  typedef typename GetGeomTraits<TriangleMeshIn, NamedParametersIn>::type  Traits;
  typedef typename GetVertexPointMap <TriangleMeshIn, NamedParametersIn>::const_type VPM_in;
  typedef typename GetVertexPointMap <TriangleMeshIn, NamedParametersOut>::type VPM_out;
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename boost::graph_traits<TriangleMeshIn> graph_traits;
  typedef typename graph_traits::edge_descriptor edge_descriptor;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits::face_descriptor face_descriptor;

  double coplanar_cos_threshold = - choose_parameter(get_parameter(np_in, internal_np::cosine_of_maximum_angle), 1);
  CGAL_precondition(coplanar_cos_threshold<=0 && coplanar_cos_threshold>=-1);

  // initialize property maps (fill user provided or user internal ones)
  typedef typename boost::property_map<TriangleMeshIn,
                                       dynamic_edge_property_t<bool> >::const_type Default_ECM;
  typedef typename boost::property_map<TriangleMeshIn,
                                       dynamic_vertex_property_t<std::size_t> >::const_type Default_VCM;
  typedef typename boost::property_map<TriangleMeshIn,
                                       dynamic_face_property_t<std::size_t> >::const_type Default_FCM;

  typename internal_np::Lookup_named_param_def< internal_np::edge_is_constrained_t,
                                                NamedParametersIn,
                                                Default_ECM>::type
    edge_is_constrained = choose_parameter<Default_ECM>(get_parameter(np_in, internal_np::edge_is_constrained),
                                                        dynamic_edge_property_t<bool>(), tm_in);

  typename internal_np::Lookup_named_param_def< internal_np::vertex_corner_map_t,
                                                NamedParametersIn,
                                                Default_VCM>::type
    vertex_corner_id = choose_parameter<Default_VCM>(get_parameter(np_in, internal_np::vertex_corner_map),
                                                     dynamic_vertex_property_t<std::size_t>(), tm_in);

  typename internal_np::Lookup_named_param_def< internal_np::face_patch_t,
                                                NamedParametersIn,
                                                Default_FCM>::type
    face_cc_ids = choose_parameter<Default_FCM>(get_parameter(np_in, internal_np::face_patch),
                                                dynamic_face_property_t<std::size_t>(), tm_in);

  for(edge_descriptor e : edges(tm_in)) put(edge_is_constrained, e, false);
  for(vertex_descriptor v : vertices(tm_in)) put(vertex_corner_id, v, Planar_segmentation::default_id());
  for(face_descriptor f : faces(tm_in)) put(face_cc_ids, f, -1);

  VPM_in vpm_in = choose_parameter(get_parameter(np_in, internal_np::vertex_point),
                                   get_const_property_map(vertex_point, tm_in));

  VPM_out vpm_out = choose_parameter(get_parameter(np_out, internal_np::vertex_point),
                                     get_property_map(vertex_point, pm_out));

  typename internal_np::Lookup_named_param_def< internal_np::visitor_t,
                                                NamedParametersOut,
                                                Planar_segmentation::Default_visitor<PolygonMeshOut>>::type
    visitor = choose_parameter<Planar_segmentation::Default_visitor<PolygonMeshOut>>(get_parameter(np_out, internal_np::visitor));

  bool do_not_triangulate_faces = choose_parameter(get_parameter(np_out, internal_np::do_not_triangulate_faces), false);

  std::pair<std::size_t, std::size_t> nb_corners_and_nb_cc =
    Planar_segmentation::tag_corners_and_constrained_edges<Traits>(tm_in, coplanar_cos_threshold, vertex_corner_id, edge_is_constrained, face_cc_ids, vpm_in);
  std::vector< typename Traits::Vector_3 > face_normals(nb_corners_and_nb_cc.second, NULL_VECTOR);
  Planar_segmentation::decimate_impl<Traits>(tm_in, pm_out,
                                             nb_corners_and_nb_cc,
                                             vertex_corner_id,
                                             edge_is_constrained,
                                             face_cc_ids,
                                             vpm_in, vpm_out,
                                             do_not_triangulate_faces,
                                             get_parameter(np_out, internal_np::vertex_corner_map),
                                             get_parameter(np_out, internal_np::face_patch),
                                             visitor,
                                             face_normals);
}

/*!
 *  \ingroup PMP_meshing_grp
 *  generates a new triangle mesh `pm_out` with the minimal number of triangles from a partition of `tm_in`.
 *  The terminology used here and the global idea is very similar to what is done by `remesh_planar_patches()`
 *  except that here the partition into patches and corner identification is provided by the user.
 *  It allows to have a remeshing of almost coplanar regions, detected for example using the region growing algorithm
 *  with the functions `region_growing_of_planes_on_faces()` and `detect_corners_of_regions()`.
 *  If a patch cannot be triangulated, it is left untouched in the output and all its vertices become corners
 *  so that the output is still a valid conformal triangle mesh.
 *  \returns `true` if all patches could be triangulated and `false` otherwise.
 *
 *  \tparam TriangleMeshIn a model of `HalfedgeListGraph` and `FaceListGraph`
 *  \tparam PolygonMeshOut a model of `MutableFaceGraph`
 *
 *  \tparam FacePatchMap a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMeshIn>::%face_descriptor`
 *                       as key type and `std::size_t` as value type
 *  \tparam EdgeIsConstrainedMap a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMeshIn>::%edge_descriptor`
 *                               as key type and `bool` as value type
 *  \tparam VertexCornerMap a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMeshIn>::%vertex_descriptor`
 *                          as key type and `std::size_t` as value type
 *  \tparam NamedParametersIn a sequence of \ref bgl_namedparameters "Named Parameters"
 *  \tparam NamedParametersOut a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 *  \param tm_in input triangle mesh
 *  \param pm_out output polygon mesh
 *  \param nb_patches the number of patches in the partition
 *  \param nb_corners the number of corners
 *  \param face_patch_map a property map that contains for each face the id of its patch in the range `[0, nb_patches]`
 *  \param vertex_corner_map a property map that contains for each vertex that is a corner an id in the range `[0, nb_corners - 1]`,
                             and `std::size_t(-1)` otherwise.
 *  \param ecm a property map that contains `true` if an edge is on the border of a patch and `false` otherwise.
 *  \param np_in an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
 *
 *  \cgalNamedParamsBegin
 *    \cgalParamNBegin{patch_normal_map}
 *      \cgalParamDescription{a property map providing for each patch the normal of the supporting plane of the patch (used to triangulate it)}
 *      \cgalParamType{a class model of `ReadPropertyMap` with the value type of `FacePatchMap` as key and
 *                     `GeomTraits::Vector_3` as value type, `GeomTraits` being the type of the parameter `geom_traits`}
 *      \cgalParamDefault{If not provided, patch normals will be estimated using corners of the patches}
 *    \cgalParamNEnd
 *    \cgalParamNBegin{vertex_point_map}
 *      \cgalParamDescription{a property map associating points to the vertices of `tm_in`}
 *      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMeshIn>::%vertex_descriptor`
 *                     as key type and `GeomTraits::Point_3` as value type, `GeomTraits` being the type of the parameter `geom_traits`}
 *      \cgalParamDefault{`boost::get(CGAL::vertex_point, tm_in)`}
 *      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                      must be available in `TriangleMeshIn`.}
 *    \cgalParamNEnd
 *    \cgalParamNBegin{geom_traits}
 *      \cgalParamDescription{an instance of a geometric traits class}
 *      \cgalParamType{a class model of `Kernel`}
 *      \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *      \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *    \cgalParamNEnd
 *  \cgalNamedParamsEnd
 *
 *  \param np_out an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
 *
 *  \cgalNamedParamsBegin
 *    \cgalParamNBegin{do_not_triangulate_faces}
 *      \cgalParamDescription{if `true`, faces of `pm_out` will not be triangulated, but the one with more than one connected component of the boundary.}
 *      \cgalParamType{`bool`}
 *      \cgalParamDefault{false}
 *    \cgalParamNEnd
 *    \cgalParamNBegin{vertex_point_map}
 *      \cgalParamDescription{a property map associating points to the vertices of `pm_out`}
 *      \cgalParamType{a class model of `WritablePropertyMap` with `boost::graph_traits<PolygonMeshOut>::%vertex_descriptor`
 *                     as key type and `GeomTraits::Point_3` as value type, `GeomTraits` being the type of the parameter `geom_traits`}
 *      \cgalParamDefault{`boost::get(CGAL::vertex_point, pm_out)`}
 *      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                      must be available in `PolygonMeshOut`.}
 *    \cgalParamNEnd
 *    \cgalParamNBegin{face_patch_map}
 *      \cgalParamDescription{a property map filled by this function and that will contain for each face the id
 *                            of its patch in the range `[0, number of patches - 1]`}
 *      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMeshOut>::%face_descriptor`
 *                     as key type and `std::size_t` as value type}
 *      \cgalParamDefault{None}
 *    \cgalParamNEnd
 *    \cgalParamNBegin{vertex_corner_map}
 *      \cgalParamDescription{a property map filled by this function and that will contain for each vertex its corner
 *                            an id in the range `[0, number of corners - 1]`}
 *      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMeshOut>::%vertex_descriptor`
 *                     as key type and `std::size_t` as value type}
 *      \cgalParamDefault{None}
 *    \cgalParamNEnd
 *    \cgalParamNBegin{visitor}
 *      \cgalParamDescription{a callable with `visitor(pm_out)` being called once `tm_in` is no longer needed
 *                            and before `pm_out` starts being built. It should be used in the case when `tm_in` and `pm_out` are the same mesh,
 *                            so that `pm_out` can be cleared before being filled.}
 *      \cgalParamType{`visitor(pm_out)` must be a valid expression.}
 *      \cgalParamDefault{None}
 *    \cgalParamNEnd
 *  \cgalNamedParamsEnd
 */
template <typename TriangleMeshIn,
          typename PolygonMeshOut,
          typename FacePatchMap,
          typename EdgeIsConstrainedMap,
          typename VertexCornerMap,
          typename NamedParametersIn = parameters::Default_named_parameters,
          typename NamedParametersOut = parameters::Default_named_parameters>
bool remesh_almost_planar_patches(const TriangleMeshIn& tm_in,
                                        PolygonMeshOut& pm_out,
                                   std::size_t nb_patches,
                                   std::size_t nb_corners,
                                   FacePatchMap face_patch_map,
                                   VertexCornerMap vertex_corner_map,
                                   EdgeIsConstrainedMap ecm,
                                   const NamedParametersIn& np_in = parameters::default_values(),
                                   const NamedParametersOut& np_out = parameters::default_values())
{
  typedef typename GetGeomTraits<TriangleMeshIn, NamedParametersIn>::type  Traits;
  typedef typename GetVertexPointMap <TriangleMeshIn, NamedParametersIn>::const_type VPM_in;
  typedef typename GetVertexPointMap <TriangleMeshIn, NamedParametersOut>::type VPM_out;
  using parameters::choose_parameter;
  using parameters::get_parameter;

  VPM_in vpm_in = choose_parameter(get_parameter(np_in, internal_np::vertex_point),
                                   get_const_property_map(vertex_point, tm_in));

  VPM_out vpm_out = choose_parameter(get_parameter(np_out, internal_np::vertex_point),
                                     get_property_map(vertex_point, pm_out));

  typename internal_np::Lookup_named_param_def< internal_np::visitor_t,
                                                NamedParametersOut,
                                                Planar_segmentation::Default_visitor<PolygonMeshOut>>::type
    visitor = choose_parameter<Planar_segmentation::Default_visitor<PolygonMeshOut>>(get_parameter(np_out, internal_np::visitor));

  bool do_not_triangulate_faces = choose_parameter(get_parameter(np_out, internal_np::do_not_triangulate_faces), false);

  std::vector< typename Traits::Vector_3 > face_normals;
  Planar_segmentation::init_face_normals(face_normals, nb_patches, get_parameter(np_in, internal_np::patch_normal_map));
  return Planar_segmentation::decimate_impl<Traits>(tm_in, pm_out,
                                                    std::make_pair(nb_corners, nb_patches),
                                                    vertex_corner_map, ecm, face_patch_map, vpm_in, vpm_out,
                                                    do_not_triangulate_faces,
                                                    get_parameter(np_out, internal_np::vertex_corner_map),
                                                    get_parameter(np_out, internal_np::face_patch),
                                                    visitor,
                                                    face_normals);
}

#ifndef DOXYGEN_RUNNING
// MeshMap must be a mutable lvalue pmap with Triangle_mesh as value_type
template <typename TriangleMeshRange, typename MeshMap>
bool decimate_meshes_with_common_interfaces(TriangleMeshRange& meshes, double coplanar_cos_threshold, MeshMap mesh_map, bool do_not_triangulate_faces=false)
{
  CGAL_assertion(coplanar_cos_threshold<0);
  typedef typename boost::property_traits<MeshMap>::value_type Triangle_mesh;
  typedef typename std::iterator_traits<typename TriangleMeshRange::iterator>::value_type Mesh_descriptor;
  typedef typename boost::property_map<Triangle_mesh, boost::vertex_point_t>::type VPM;
  typedef typename boost::property_traits<VPM>::value_type Point_3;
  typedef typename Kernel_traits<Point_3>::type Kernel;

  // todo turn into a range of named parameter
  std::vector<typename boost::property_map<Triangle_mesh, boost::vertex_point_t>::type > vpms;
  vpms.reserve(meshes.size());

  for(Mesh_descriptor& md : meshes)
    vpms.push_back( get(boost::vertex_point, mesh_map[md]) );
  return Planar_segmentation::decimate_meshes_with_common_interfaces_impl<Kernel>(meshes, mesh_map, coplanar_cos_threshold, vpms, do_not_triangulate_faces);
}

template <class TriangleMesh>
bool decimate_meshes_with_common_interfaces(std::vector<TriangleMesh>& meshes, double coplanar_cos_threshold=-1)
{
  return decimate_meshes_with_common_interfaces(meshes, coplanar_cos_threshold, CGAL::Identity_property_map<TriangleMesh>());
}
#endif

} } // end of CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_REMESH_PLANAR_PATCHES_H
