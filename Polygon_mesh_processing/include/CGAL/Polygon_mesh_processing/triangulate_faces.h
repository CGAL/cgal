// Copyright (c) 2010-2011  GeometryFactory Sarl (France).
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
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_POLYGON_MESH_PROCESSING_TRIANGULATE_FACES_H
#define CGAL_POLYGON_MESH_PROCESSING_TRIANGULATE_FACES_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/Euler_operations.h>

#ifndef CGAL_TRIANGULATE_FACES_DO_NOT_USE_CDT2
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_2_projection_traits_3.h>
#else
#include <CGAL/use.h>
#endif

#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Eigen_diagonalize_traits.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <boost/range/size.hpp>
#include <boost/foreach.hpp>

#include <queue>
#include <vector>
#include <utility>
#include <CGAL/array.h>

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {

template <class PM
          , typename VertexPointMap
          , typename Kernel>
class Triangulate_modifier
{
  typedef Kernel Traits;

  typedef typename boost::graph_traits<PM>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PM>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<PM>::edge_descriptor edge_descriptor;
  typedef typename Kernel::Point_3 Point;

  struct Face_info {
    typename boost::graph_traits<PM>::halfedge_descriptor e[3];
    bool is_external;
  };

  typedef typename boost::property_traits<VertexPointMap>::reference Point_ref;
  VertexPointMap _vpmap;
  Traits _traits;

public:
  Triangulate_modifier(VertexPointMap vpmap, const Traits& traits = Traits())
    : _vpmap(vpmap), _traits(traits)
  {
  }

  template <class Face_handle>
  bool is_external(Face_handle fh) const {
    return fh->info().is_external;
  }

  bool triangulate_face(face_descriptor f, PM& pmesh, bool use_cdt)
  {
    Emptyset_iterator out;
    return triangulate_face(f, pmesh, use_cdt, out);
  }

  // Default for out is Emptyset_iterator(),
  // and out is only used for triangulate_with_cdt_2().
  template <typename OutputIterator>
  bool triangulate_face(face_descriptor f, PM& pmesh, bool use_cdt, OutputIterator out)
  {
    typedef typename Traits::FT FT;

    typename Traits::Vector_3 normal =
      Polygon_mesh_processing::compute_face_normal(
        f, pmesh, CGAL::Polygon_mesh_processing::parameters::geom_traits(_traits)
                                                            .vertex_point_map(_vpmap));

    if(normal == typename Traits::Vector_3(0,0,0))
      return false;

    std::size_t original_size = CGAL::halfedges_around_face(halfedge(f, pmesh), pmesh).size();
    if(original_size == 4)
    {
      halfedge_descriptor v0, v1, v2, v3;
      v0 = halfedge(f, pmesh);
      Point_ref p0 = get(_vpmap, target(v0, pmesh));
      v1 = next(v0, pmesh);
      Point_ref p1 = get(_vpmap, target(v1, pmesh));
      v2 = next(v1, pmesh);
      Point_ref p2 = get(_vpmap, target(v2, pmesh));
      v3 = next(v2, pmesh);
      Point_ref p3 = get(_vpmap, target(v3, pmesh));

      /* Chooses the diagonal that will split the quad in two triangles that maximize
       * the scalar product of of the un-normalized normals of the two triangles.
       * The lengths of the un-normalized normals (computed using cross-products of two vectors)
       *  are proportional to the area of the triangles.
       * Maximize the scalar product of the two normals will avoid skinny triangles,
       * and will also taken into account the cosine of the angle between the two normals.
       * In particular, if the two triangles are oriented in different directions,
       * the scalar product will be negative.
       */
      FT p1p3 = CGAL::cross_product(p2-p1,p3-p2) * CGAL::cross_product(p0-p3,p1-p0);
      FT p0p2 = CGAL::cross_product(p1-p0,p1-p2) * CGAL::cross_product(p3-p2,p3-p0);
      if(p0p2>p1p3)
      {
        // CGAL::Euler::split_face(v0, v2, pmesh);
        halfedge_descriptor hnew = halfedge(add_edge(pmesh), pmesh);
        face_descriptor fnew = add_face(pmesh);
        CGAL::internal::insert_tip(hnew, v2, pmesh);
        CGAL::internal::insert_tip(opposite(hnew, pmesh), v0, pmesh);
        set_face(hnew, face(v0, pmesh), pmesh);
        CGAL::internal::set_face_in_face_loop(opposite(hnew, pmesh), fnew, pmesh);
        set_halfedge(face(hnew, pmesh), hnew, pmesh);
        set_halfedge(face(opposite(hnew, pmesh), pmesh), opposite(hnew, pmesh), pmesh);
        *out++ = fnew;
      }
      else
      {
        // CGAL::Euler::split_face(v1, v3, pmesh);
        halfedge_descriptor hnew = halfedge(add_edge(pmesh), pmesh);
        face_descriptor fnew = add_face(pmesh);
        CGAL::internal::insert_tip(hnew, v3, pmesh);
        CGAL::internal::insert_tip(opposite(hnew, pmesh), v1, pmesh);
        set_face(hnew, face(v1, pmesh), pmesh);
        CGAL::internal::set_face_in_face_loop(opposite(hnew, pmesh), fnew, pmesh);
        set_halfedge(face(hnew, pmesh), hnew, pmesh);
        set_halfedge(face(opposite(hnew, pmesh), pmesh), opposite(hnew, pmesh), pmesh);
        *out++ = fnew;
      }
    }
    else
    {
#ifndef CGAL_TRIANGULATE_FACES_DO_NOT_USE_CDT2
      if (use_cdt)
      {
        typedef CGAL::Triangulation_2_projection_traits_3<Traits>   P_traits;
        typedef CGAL::Triangulation_vertex_base_with_info_2<halfedge_descriptor,
                                                            P_traits>        Vb;
        typedef CGAL::Triangulation_face_base_with_info_2<Face_info,
                                                          P_traits>          Fb1;
        typedef CGAL::Constrained_triangulation_face_base_2<P_traits, Fb1>   Fb;
        typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                  TDS;
        typedef CGAL::Exact_intersections_tag                                Itag;
        typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits,
                                                           TDS,
                                                           Itag>             CDT;
        P_traits cdt_traits(normal);
        CDT cdt(cdt_traits);
        return triangulate_face_with_CDT(f, pmesh, cdt, out);
      }
#else
      CGAL_USE(use_cdt);
#endif
      // Don't use out if no cdt.
      return triangulate_face_with_hole_filling(f, pmesh);
    }
    return true;
  }

  template <class CDT, typename OutputIterator>
  bool triangulate_face_with_CDT(face_descriptor f, PM& pmesh, CDT& cdt, OutputIterator out)
  {
    std::size_t original_size = CGAL::halfedges_around_face(halfedge(f, pmesh), pmesh).size();

    // Halfedge_around_facet_circulator
    typedef typename CDT::Vertex_handle Tr_Vertex_handle;
    halfedge_descriptor start = halfedge(f, pmesh);
    halfedge_descriptor h = start;
    Tr_Vertex_handle previous, first;
    do
    {
      Tr_Vertex_handle vh = cdt.insert(get(_vpmap, target(h, pmesh)));
      if (first == Tr_Vertex_handle()) {
        first = vh;
      }
      vh->info() = h;
      if(previous != Tr_Vertex_handle() && previous != vh) {
        cdt.insert_constraint(previous, vh);
      }
      previous = vh;
      h = next(h, pmesh);

    } while( h != start );
    cdt.insert_constraint(previous, first);

    // sets mark is_external
    for(typename CDT::All_faces_iterator fit = cdt.all_faces_begin(),
          end = cdt.all_faces_end();
        fit != end; ++fit)
    {
      fit->info().is_external = false;
    }
    std::queue<typename CDT::Face_handle> face_queue;
    face_queue.push(cdt.infinite_vertex()->face());
    while(! face_queue.empty() )
    {
      typename CDT::Face_handle fh = face_queue.front();
      face_queue.pop();

      if(fh->info().is_external)
        continue;

      fh->info().is_external = true;
      for(int i = 0; i <3; ++i)
      {
        if(!cdt.is_constrained(typename CDT::Edge(fh, i)))
        {
          face_queue.push(fh->neighbor(i));
        }
      }
    }

    if(cdt.dimension() != 2 ||
       cdt.number_of_vertices() != original_size)
      return false;


    // then modify the polyhedron
    // make_hole. (see comment in function body)
    this->make_hole(halfedge(f, pmesh), pmesh);

    for(typename CDT::Finite_edges_iterator eit = cdt.finite_edges_begin(),
          end = cdt.finite_edges_end();
        eit != end; ++eit)
    {
      typename CDT::Face_handle fh = eit->first;
      const int index = eit->second;
      typename CDT::Face_handle opposite_fh = fh->neighbor(eit->second);
      const int opposite_index = opposite_fh->index(fh);

      const Tr_Vertex_handle va = fh->vertex(cdt. cw(index));
      const Tr_Vertex_handle vb = fh->vertex(cdt.ccw(index));

      if( ! (is_external(fh) && is_external(opposite_fh))//not both fh are external
          && ! cdt.is_constrained(*eit) )                  //and edge is not constrained
      {
        // strictly internal edge
        halfedge_descriptor hnew = halfedge(add_edge(pmesh), pmesh),
          hnewopp = opposite(hnew, pmesh);

        fh->info().e[index] = hnew;
        opposite_fh->info().e[opposite_index] = hnewopp;

        set_target(hnew,    target(va->info(), pmesh), pmesh);
        set_target(hnewopp, target(vb->info(), pmesh), pmesh);
      }
      if( cdt.is_constrained(*eit) ) //edge is constrained
      {
        if(!is_external(fh)) {
          fh->info().e[index] = va->info();
        }
        if(!is_external(opposite_fh)) {
          opposite_fh->info().e[opposite_index] = vb->info();
        }
      }
    }
    for(typename CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(),
          end = cdt.finite_faces_end();
        fit != end; ++fit)
    {
      if(!is_external(fit))
      {
        halfedge_descriptor h0 = fit->info().e[0];
        halfedge_descriptor h1 = fit->info().e[1];
        halfedge_descriptor h2 = fit->info().e[2];
        CGAL_assertion(h0 != halfedge_descriptor());
        CGAL_assertion(h1 != halfedge_descriptor());
        CGAL_assertion(h2 != halfedge_descriptor());

        set_next(h0, h1, pmesh);
        set_next(h1, h2, pmesh);
        set_next(h2, h0, pmesh);

        // Euler::fill_hole(h0, pmesh);
        face_descriptor new_face = add_face(pmesh);
        BOOST_FOREACH(halfedge_descriptor hd, CGAL::halfedges_around_face(h0, pmesh)) {
          set_face(hd, new_face, pmesh);
        }
        set_halfedge(new_face, h0, pmesh);
        *out++ = new_face;
      }
    }
    return true;
  }

  bool triangulate_face_with_hole_filling(face_descriptor f, PM& pmesh)
  {
    namespace PMP = CGAL::Polygon_mesh_processing;

    // gather halfedges around the face
    std::vector<Point> hole_points;
    std::vector<vertex_descriptor> border_vertices;
    CGAL_assertion(CGAL::halfedges_around_face(halfedge(f, pmesh), pmesh).size() > 0);
    BOOST_FOREACH(halfedge_descriptor h, CGAL::halfedges_around_face(halfedge(f, pmesh), pmesh))
    {
      vertex_descriptor v = source(h, pmesh);
      hole_points.push_back( get(_vpmap, v) );
      border_vertices.push_back(v);
    }

    // use hole filling
    typedef CGAL::Triple<int, int, int> Face_indices;
    std::vector<Face_indices> patch;
    PMP::triangulate_hole_polyline(hole_points, std::back_inserter(patch),
                                   PMP::parameters::geom_traits(_traits));

    if(patch.empty())
      return false;

    // triangulate the hole
    std::map< std::pair<int, int> , halfedge_descriptor > halfedge_map;
    int i=0;
    BOOST_FOREACH(halfedge_descriptor h, CGAL::halfedges_around_face(halfedge(f, pmesh), pmesh))
    {
      int j = std::size_t(i+1) == hole_points.size() ? 0 : i+1;
      halfedge_map[ std::make_pair(i, j) ] = h;
      ++i;
    }

    bool first = true;
    std::vector<halfedge_descriptor> hedges;
    hedges.reserve(4);
    BOOST_FOREACH(const Face_indices& triangle, patch)
    {
      if (first)
        first=false;
      else
        f=add_face(pmesh);

      cpp11::array<int, 4> indices =
        make_array( triangle.first,
                    triangle.second,
                    triangle.third,
                    triangle.first );
      for (int i=0; i<3; ++i)
      {
        typename std::map< std::pair<int, int> , halfedge_descriptor >::iterator insert_res =
          halfedge_map.insert(
            std::make_pair( std::make_pair(indices[i], indices[i+1]),
                            boost::graph_traits<PM>::null_halfedge() ) ).first;
        if (insert_res->second == boost::graph_traits<PM>::null_halfedge())
        {
          halfedge_descriptor nh = halfedge(add_edge(pmesh), pmesh);
          insert_res->second=nh;
          halfedge_map[std::make_pair(indices[i+1], indices[i])]=opposite(nh, pmesh);
        }
        hedges.push_back(insert_res->second);
      }
      hedges.push_back(hedges.front());
      for(int i=0; i<3;++i)
      {
        set_next(hedges[i], hedges[i+1], pmesh);
        set_face(hedges[i], f, pmesh);
        set_target(hedges[i], border_vertices[indices[i+1]], pmesh);
      }
      set_halfedge(f, hedges[0], pmesh);
      hedges.clear();
    }
    return true;
  }

  template<typename FaceRange>
  bool operator()(FaceRange face_range, PM& pmesh, bool use_cdt)
  {
   bool result = true;
    // One need to store facet handles into a vector, because the list of
    // facets of the polyhedron will be modified during the loop, and
    // that invalidates the range [facets_begin(), facets_end()[.
    std::vector<face_descriptor> facets;
    facets.reserve(std::distance(boost::begin(face_range), boost::end(face_range)));

    //only consider non-triangular faces
    BOOST_FOREACH(face_descriptor fit, face_range)
      if ( next( next( halfedge(fit, pmesh), pmesh), pmesh)
        !=       prev( halfedge(fit, pmesh), pmesh) )
        facets.push_back(fit);

    // Iterates on the vector of face descriptors
    BOOST_FOREACH(face_descriptor f, facets)
    {
     if(!this->triangulate_face(f, pmesh, use_cdt))
       result = false;
    }
    return result;
  }

  void make_hole(halfedge_descriptor h, PM& pmesh)
  {
    //we are not using Euler::make_hole because it has a precondition
    //that the hole is not made on the boundary of the mesh
    //here we allow making a hole on the boundary, and the pair(s) of
    //halfedges that become border-border are fixed by the connectivity
    //setting made in operator()
    CGAL_assertion(!is_border(h, pmesh));
    face_descriptor fd = face(h, pmesh);

    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(h, pmesh))
    {
      CGAL::internal::set_border(hd, pmesh);
    }
    remove_face(fd, pmesh);
  }


}; // end class Triangulate_modifier

}//end namespace internal

/**
* \ingroup PMP_meshing_grp
* triangulates a single face of a polygon mesh. This function depends on the package \ref PkgTriangulation2
* @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
* @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
*
* @param f face to be triangulated
* @param pmesh the polygon mesh to which the face to be triangulated belongs
* @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
*
*
* \cgalNamedParamsBegin
*    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
*   If this parameter is omitted, an internal property map for
*   `CGAL::vertex_point_t` must be available in `PolygonMesh`\cgalParamEnd
*    \cgalParamBegin{geom_traits} a geometric traits class instance \cgalParamEnd
* \cgalNamedParamsEnd
*
* @return `true` if the face has been triangulated.
*/
template<typename PolygonMesh, typename NamedParameters>
bool triangulate_face(typename boost::graph_traits<PolygonMesh>::face_descriptor f,
                      PolygonMesh& pmesh,
                      const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  //VertexPointMap
  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VPMap;
  VPMap vpmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_property_map(vertex_point, pmesh));

  //Kernel
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type Kernel;
  Kernel traits = choose_parameter(get_parameter(np, internal_np::geom_traits), Kernel());

  //Option
  bool use_cdt = choose_parameter(get_parameter(np, internal_np::use_delaunay_triangulation), true);

  typedef typename internal_np::Lookup_named_param_def<
    internal_np::output_iterator_t,
    NamedParameters,
    Emptyset_iterator>::type Output_iterator;

  Output_iterator out = choose_parameter(
    get_parameter(np, internal_np::output_iterator), Emptyset_iterator());

  internal::Triangulate_modifier<PolygonMesh, VPMap, Kernel> modifier(vpmap, traits);
  return modifier.triangulate_face(f, pmesh, use_cdt, out);
}

template<typename PolygonMesh>
bool triangulate_face(typename boost::graph_traits<PolygonMesh>::face_descriptor f,
                      PolygonMesh& pmesh)
{
  return triangulate_face(f, pmesh, CGAL::Polygon_mesh_processing::parameters::all_default());
}

/**
* \ingroup PMP_meshing_grp
* triangulates given faces of a polygon mesh. This function depends on the package \ref PkgTriangulation2
*
* @tparam FaceRange range of `boost::graph_traits<PolygonMesh>::%face_descriptor`,
          model of `Range`.
          Its iterator type is `InputIterator`.
* @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
* @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
*
* @param face_range the range of faces to be triangulated
* @param pmesh the polygon mesh to be triangulated
* @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
*   If this parameter is omitted, an internal property map for
*   `CGAL::vertex_point_t` must be available in `PolygonMesh`\cgalParamEnd
*    \cgalParamBegin{geom_traits} a geometric traits class instance \cgalParamEnd
* \cgalNamedParamsEnd
*
* @return `true` if all the faces have been triangulated.
* @see triangulate_face()
*/
template <typename FaceRange, typename PolygonMesh, typename NamedParameters>
bool triangulate_faces(FaceRange face_range,
                       PolygonMesh& pmesh,
                       const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  //VertexPointMap
  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VPMap;
  VPMap vpmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_property_map(vertex_point, pmesh));

  //Kernel
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type Kernel;
  Kernel traits = choose_parameter(get_parameter(np, internal_np::geom_traits), Kernel());

  //Option
  bool use_cdt = choose_parameter(get_parameter(np, internal_np::use_delaunay_triangulation), true);

  internal::Triangulate_modifier<PolygonMesh, VPMap, Kernel> modifier(vpmap, traits);
  return modifier(face_range, pmesh, use_cdt);
}

template <typename FaceRange, typename PolygonMesh>
bool triangulate_faces(FaceRange face_range, PolygonMesh& pmesh)
{
  return triangulate_faces(face_range, pmesh, CGAL::Polygon_mesh_processing::parameters::all_default());
}

/**
* \ingroup PMP_meshing_grp
* triangulates all faces of a polygon mesh. This function depends on the package \ref PkgTriangulation2
* @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
* @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
*
* @param pmesh the polygon mesh to be triangulated
* @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
*   If this parameter is omitted, an internal property map for
*   `CGAL::vertex_point_t` must be available in `PolygonMesh`\cgalParamEnd
*    \cgalParamBegin{geom_traits} a geometric traits class instance \cgalParamEnd
* \cgalNamedParamsEnd
*
* @return `true` if all the faces have been triangulated.
* @see triangulate_face()
*/
template <typename PolygonMesh, typename NamedParameters>
bool triangulate_faces(PolygonMesh& pmesh,
                       const NamedParameters& np)
{
  return triangulate_faces(faces(pmesh), pmesh, np);
}

template <typename PolygonMesh>
bool triangulate_faces(PolygonMesh& pmesh)
{
  return triangulate_faces(faces(pmesh), pmesh, CGAL::Polygon_mesh_processing::parameters::all_default());
}

/// \cond SKIP_IN_MANUAL
template<typename GeomTraits>
bool is_planar_2(
  const std::vector<typename GeomTraits::Point_3>& polyline_3d,
  const typename GeomTraits::Plane_3& fitted_plane,
  GeomTraits traits) {

  // Typedefs.
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Compute_squared_length_3 Squared_length_3;
  typedef typename GeomTraits::Compute_squared_distance_3 Squared_distance_3;
  typedef typename GeomTraits::Compute_scalar_product_3 Scalar_product_3;
  typedef typename GeomTraits::Collinear_3 Collinear_3;

  const Squared_length_3 squared_length_3 = traits.compute_squared_length_3_object();
  const Squared_distance_3 squared_distance_3 = traits.compute_squared_distance_3_object();
  const Scalar_product_3 scalar_product_3 = traits.compute_scalar_product_3_object();
  const Collinear_3 collinear_3 = traits.collinear_3_object();

  CGAL_precondition(
    polyline_3d.size() >= 3);

  // Distance criteria.

  // Tolerance.
  const FT dist_tol = FT(1) / FT(100000);

  // Compute max distance.
  FT max_dist = -FT(1);
  for (std::size_t i = 0; i < polyline_3d.size(); ++i) {
    const Point_3& p = polyline_3d[i];
    const FT dist = static_cast<FT>(
      CGAL::sqrt(CGAL::to_double(squared_distance_3(p, fitted_plane))));
    max_dist = (CGAL::max)(dist, max_dist);
  }
  CGAL_assertion(max_dist != -FT(1));

  // Angle criteria.

  // Tolerance.
  const FT angle_tol = FT(5); // degrees
  const FT normal_tol = static_cast<FT>(
    std::cos(CGAL::to_double(
      (angle_tol * static_cast<FT>(CGAL_PI)) / FT(180))));

  // Compute fitted plane normal.
  Vector_3 normal = fitted_plane.orthogonal_vector();
  FT normal_length = static_cast<FT>(
    CGAL::sqrt(CGAL::to_double(squared_length_3(normal))));
  CGAL_assertion(normal_length > FT(0));
  const Vector_3 ref_normal = normal / normal_length;

  // Compute average normal of the hole.
  FT x = FT(0), y = FT(0), z = FT(0);
  std::size_t num_normals = 0;
  const Point_3& ref_point = polyline_3d[0];
  for (std::size_t i = 1; i < polyline_3d.size() - 1; ++i) {
    const std::size_t ip = i + 1;

    const Point_3& p1 = ref_point;
    const Point_3& p2 = polyline_3d[i];
    const Point_3& p3 = polyline_3d[ip];

    // Skip in case we have collinear points.
    if (collinear_3(p1, p2, p3))
      continue;

    normal = CGAL::normal(p1, p2, p3);
    normal_length = static_cast<FT>(
      CGAL::sqrt(CGAL::to_double(squared_length_3(normal))));
    CGAL_assertion(normal_length > FT(0));
    normal /= normal_length;

    x += normal.x(); y += normal.y(); z += normal.z();
    ++num_normals;
  }
  CGAL_assertion(num_normals >= 1);
  x /= static_cast<FT>(num_normals);
  y /= static_cast<FT>(num_normals);
  z /= static_cast<FT>(num_normals);

  normal = Vector_3(x, y, z);
  normal_length = static_cast<FT>(
    CGAL::sqrt(CGAL::to_double(squared_length_3(normal))));
  const Vector_3 avg_normal = normal / normal_length;

  const FT cos_value =
    CGAL::abs(scalar_product_3(avg_normal, ref_normal));

  // Check planarity.
  const bool is_planar = (
    (  max_dist <= dist_tol ) &&
    ( cos_value >= normal_tol ));
  return is_planar;
}
/// \endcond

/*!
  \ingroup hole_filling_grp
  \brief triangulates a planar hole in a polygon mesh.

  If the hole is planar, this function uses the 2D constrained Delaunay triangulation
  in order to close the hole. The constraints are the border edges of the hole.

  The hole must not contain any non-manifold vertex, nor self-intersections.
  The patch generated does not introduce non-manifold edges nor degenerate triangles.
  If a hole cannot be triangulated, `pmesh` is not modified and nothing is recorded in `out`.

  \tparam PolygonMesh a model of `MutableFaceGraph`
  \tparam OutputIterator a model of `OutputIterator` holding
  `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces.
  \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"

  \param pmesh polygon mesh which has the hole
  \param border_halfedge a border halfedge incident to the hole
  \param out iterator over patch faces
  \param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
        If this parameter is omitted, an internal property map for
        `CGAL::vertex_point_t` should be available in `PolygonMesh`
        \cgalParamEnd
    \cgalParamBegin{geom_traits} a geometric traits class instance \cgalParamEnd
  \cgalNamedParamsEnd

  \return `out`
*/
template<
typename PolygonMesh,
typename OutputIterator,
typename NamedParameters>
OutputIterator triangulate_hole_with_cdt_2(
  PolygonMesh& pmesh,
  typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
  OutputIterator out,
  const NamedParameters& np) {

  typedef Halfedge_around_face_circulator<PolygonMesh> Hedge_around_face_circulator;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type VPM;
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type Kernel;

  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Plane_3 Plane_3;

  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;

  const VPM vpm = parameters::choose_parameter(
    parameters::get_parameter(np, internal_np::vertex_point),
    get_const_property_map(boost::vertex_point, pmesh));
  const Kernel traits = parameters::choose_parameter(
    parameters::get_parameter(np, internal_np::geom_traits), Kernel());

  Hedge_around_face_circulator circ(border_halfedge, pmesh);
  Hedge_around_face_circulator done(circ);

  std::vector<Point_3> polyline_3d;
  do {
    polyline_3d.push_back(get(vpm, target(*circ, pmesh)));
  } while (++circ != done);

  // Plane fitting.
  typedef Exact_predicates_inexact_constructions_kernel Local_kernel;
  typedef typename Local_kernel::FT Local_FT;
  typedef typename Local_kernel::Point_3 Local_point_3;
  typedef typename Local_kernel::Plane_3 Local_plane_3;
  typedef Cartesian_converter<Kernel, Local_kernel> To_local_converter;

  const To_local_converter to_local_converter;
  std::vector<Local_point_3> points;
  points.reserve(polyline_3d.size());
  for (std::size_t i = 0; i < polyline_3d.size(); ++i)
    points.push_back(to_local_converter(polyline_3d[i]));

  Local_plane_3 fitted_plane;
  Local_point_3 fitted_centroid;
  linear_least_squares_fitting_3(
    points.begin(), points.end(), fitted_plane, fitted_centroid,
    CGAL::Dimension_tag<0>(), Local_kernel(),
    CGAL::Eigen_diagonalize_traits<Local_FT, 3>());

  const Plane_3 plane = Plane_3(
    static_cast<FT>(fitted_plane.a()),
    static_cast<FT>(fitted_plane.b()),
    static_cast<FT>(fitted_plane.c()),
    static_cast<FT>(fitted_plane.d()));

  // Checking simplicity and planarity.
  std::vector<Point_2> polyline_2d;
  polyline_2d.reserve(polyline_3d.size());
  for (std::size_t i = 0; i < polyline_3d.size(); ++i)
    polyline_2d.push_back(plane.to_2d(polyline_3d[i]));

  const bool is_simple =
    is_simple_2(polyline_2d.begin(), polyline_2d.end(), traits);
  const bool is_planar =
    is_planar_2(polyline_3d, plane, traits);

  const bool is_planar_hole = is_simple && is_planar;
  if (!is_planar_hole)
    return triangulate_hole(
      pmesh, border_halfedge, out, np);

  face_descriptor new_face = add_face(pmesh);
  set_halfedge(new_face, border_halfedge, pmesh);
  do {
    set_face(*circ, new_face, pmesh);
  } while (++circ != done);

  // Triangulating.
  const bool use_cdt = true;
  internal::Triangulate_modifier<PolygonMesh, VPM, Kernel> modifier(vpm, traits);
  const bool success_with_cdt_2 =
    modifier.triangulate_face(new_face, pmesh, use_cdt, out);
  CGAL_assertion(success_with_cdt_2);
  return out;
}

} // end namespace Polygon_mesh_processing

} // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_POLYGON_MESH_PROCESSING_TRIANGULATE_FACES_H
