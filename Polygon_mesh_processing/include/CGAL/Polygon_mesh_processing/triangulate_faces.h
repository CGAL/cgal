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

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_2_projection_traits_3.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <boost/range/size.hpp>
#include <boost/foreach.hpp>

#include <queue>
#include <vector>
#include <utility>

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

  typedef CGAL::Triangulation_2_projection_traits_3<Traits>   P_traits;

  typedef CGAL::Triangulation_vertex_base_with_info_2<halfedge_descriptor,
                                                      P_traits>        Vb;

  struct Face_info {
    typename boost::graph_traits<PM>::halfedge_descriptor e[3];
    bool is_external;
  };

  typedef CGAL::Triangulation_face_base_with_info_2<Face_info,
                                                    P_traits>          Fb1;
  typedef CGAL::Constrained_triangulation_face_base_2<P_traits, Fb1>   Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                  TDS;
  typedef CGAL::Exact_intersections_tag                                Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits,
                                                     TDS,
                                                     Itag>             CDT;

  typedef typename boost::property_traits<VertexPointMap>::reference Point_ref;
  VertexPointMap _vpmap;

public:
  Triangulate_modifier(VertexPointMap vpmap)
    : _vpmap(vpmap)
  {
  }

  bool is_external(typename CDT::Face_handle fh) const {
    return fh->info().is_external;
  }

  bool triangulate_face(face_descriptor f, PM& pmesh)
  {
    typedef typename Traits::FT FT;
    typename Traits::Vector_3 normal =
      Polygon_mesh_processing::compute_face_normal(f, pmesh);
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
        CGAL::Euler::split_face(v0, v2, pmesh);
      }
      else
      {
        CGAL::Euler::split_face(v1, v3, pmesh);
      }
    }
    else
    {
      P_traits cdt_traits(normal);
      CDT cdt(cdt_traits);

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

          Euler::fill_hole(h0, pmesh);
        }
      }
    }
    return true;
  }

  template<typename FaceRange>
  bool operator()(FaceRange face_range, PM& pmesh)
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
     if(!this->triangulate_face(f, pmesh))
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
* triangulates a single face of a polygon mesh. This function depends on the package \ref PkgTriangulation2Summary
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
*   `CGAL::vertex_point_t` should be available in `PolygonMesh`\cgalParamEnd
*    \cgalParamBegin{geom_traits} a geometric traits class instance \cgalParamEnd
* \cgalNamedParamsEnd
*
* @return `true` if the face has been triangulated.
*/
template<typename PolygonMesh, typename NamedParameters>
bool  triangulate_face(typename boost::graph_traits<PolygonMesh>::face_descriptor f,
                      PolygonMesh& pmesh,
                      const NamedParameters& np)
{
  using boost::choose_param;
  using boost::get_param;

  //VertexPointMap
  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VPMap;
  VPMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                             get_property_map(vertex_point, pmesh));
  //Kernel
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type Kernel;

  internal::Triangulate_modifier<PolygonMesh, VPMap, Kernel> modifier(vpmap);
  return modifier.triangulate_face(f, pmesh);
}

template<typename PolygonMesh>
bool triangulate_face(typename boost::graph_traits<PolygonMesh>::face_descriptor f,
                      PolygonMesh& pmesh)
{
  return triangulate_face(f, pmesh, CGAL::Polygon_mesh_processing::parameters::all_default());
}

/**
* \ingroup PMP_meshing_grp
* triangulates given faces of a polygon mesh. This function depends on the package \ref PkgTriangulation2Summary
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
*   `CGAL::vertex_point_t` should be available in `PolygonMesh`\cgalParamEnd
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
  using boost::choose_param;
  using boost::get_param;

  //VertexPointMap
  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VPMap;
  VPMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                             get_property_map(vertex_point, pmesh));
  //Kernel
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type Kernel;

  internal::Triangulate_modifier<PolygonMesh, VPMap, Kernel> modifier(vpmap);
  return modifier(face_range, pmesh);
}

template <typename FaceRange, typename PolygonMesh>
bool triangulate_faces(FaceRange face_range, PolygonMesh& pmesh)
{
  return triangulate_faces(face_range, pmesh, CGAL::Polygon_mesh_processing::parameters::all_default());
}

/**
* \ingroup PMP_meshing_grp
* triangulates all faces of a polygon mesh. This function depends on the package \ref PkgTriangulation2Summary
* @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
* @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
*
* @param pmesh the polygon mesh to be triangulated
* @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
*   If this parameter is omitted, an internal property map for
*   `CGAL::vertex_point_t` should be available in `PolygonMesh`\cgalParamEnd
*    \cgalParamBegin{geom_traits} a geometric traits class instance \cgalParamEnd
* \cgalNamedParamsEnd
*
* @return `true` if all the faces have been triangulated.
* @see triangulate_face()
*/
template <typename PolygonMesh, typename NamedParameters>
bool  triangulate_faces(PolygonMesh& pmesh,
                       const NamedParameters& np)
{
  return triangulate_faces(faces(pmesh), pmesh, np);
}

template <typename PolygonMesh>
bool triangulate_faces(PolygonMesh& pmesh)
{
  return triangulate_faces(faces(pmesh), pmesh, CGAL::Polygon_mesh_processing::parameters::all_default());
}

} // end namespace Polygon_mesh_processing

} // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_POLYGON_MESH_PROCESSING_TRIANGULATE_FACES_H
