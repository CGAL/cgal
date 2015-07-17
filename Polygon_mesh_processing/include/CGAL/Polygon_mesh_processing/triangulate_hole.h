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
// Author(s)     : Ilker O. Yaz

#ifndef CGAL_POLYGON_MESH_PROCESSING_TRIANGULATE_HOLE_H
#define CGAL_POLYGON_MESH_PROCESSING_TRIANGULATE_HOLE_H

#include <CGAL/Polygon_mesh_processing/internal/Hole_filling/Triangulate_hole_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/internal/Hole_filling/Triangulate_hole_polyline.h>
#include <CGAL/Polygon_mesh_processing/refine.h>
#include <CGAL/Polygon_mesh_processing/fair.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <CGAL/boost/graph/helpers.h>

#include <boost/tuple/tuple.hpp>

namespace CGAL {

namespace Polygon_mesh_processing {

  /*!
  \ingroup PkgPolygonMeshProcessing
  triangulates a hole in a polygon mesh.
  The hole must not contain any non-manifold vertex.
  The patch generated does not introduce non-manifold edges nor degenerate triangles.
  If a hole cannot be triangulated, `pmesh` is not modified and nothing is recorded in `out`.

  @tparam PolygonMesh a model of `MutableFaceGraph`
          that has an internal property map for `CGAL::vertex_point_t`
  @tparam OutputIterator a model of `OutputIterator`
    holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces.
  @tparam NamedParameters a sequence of \ref namedparameters

  @param pmesh polygon mesh containing the hole
  @param border_halfedge a border halfedge incident to the hole
  @param out iterator over patch faces
  @param np optional sequence of \ref namedparameters among the ones listed below

  \cgalNamedParamsBegin
     \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh` \cgalParamEnd
     \cgalParamBegin{use_delaunay_triangulation} if `true`, use the Delaunay triangulation facet search space \cgalParamEnd
     \cgalParamBegin{geom_traits} a geometric traits class instance \cgalParamEnd
  \cgalNamedParamsEnd

  @return `out`

  \todo handle islands
  @todo Replace border_halfedge by a range of border halfedges.
        The first one would describe the hole,
        the other ones would describe the islands.
  @todo Then, insert the holes vertices in the set of possibilities
        for connecting vertices together
  @todo handle the case where an island is reduced to a point
  */  
  template<typename PolygonMesh,
           typename OutputIterator,
           typename NamedParameters>
  OutputIterator
  triangulate_hole(PolygonMesh& pmesh,
              typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
              OutputIterator out,
              const NamedParameters& np)
  {
    using boost::choose_param;
    using boost::get_param;

    bool use_dt3 =
#ifdef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
      false;
#else
      choose_param(get_param(np, use_delaunay_triangulation), true);
#endif

    CGAL_precondition(face(border_halfedge, pmesh) == boost::graph_traits<PolygonMesh>::null_face());

    return internal::triangulate_hole_polygon_mesh(pmesh,
      border_halfedge,
      out,
      choose_param(get_param(np, vertex_point), get(CGAL::vertex_point, pmesh)),
      use_dt3,
      choose_param(get_param(np, geom_traits), typename GetGeomTraits<PolygonMesh,NamedParameters>::type()))
      .first;
  }

  template<typename PolygonMesh, typename OutputIterator>
  OutputIterator
    triangulate_hole(PolygonMesh& pmesh,
      typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
      OutputIterator out)
  {
    return triangulate_hole(pmesh, border_halfedge, out,
      CGAL::Polygon_mesh_processing::parameters::all_default());
  }

  template<typename PM, typename VertexRange>
  void test_in_edges(const PM& pmesh, const VertexRange& patch)
  {
    BOOST_FOREACH(typename boost::graph_traits<PM>::vertex_descriptor v0, patch)
    {
      typename boost::graph_traits<PM>::in_edge_iterator e, e_end;
      for (boost::tie(e, e_end) = in_edges(v0, pmesh); e != e_end; e++)
      {
        typename boost::graph_traits<PM>::halfedge_descriptor he = halfedge(*e, pmesh);
        if (is_border(he, pmesh)) { continue; }

        CGAL_assertion(v0 == target(he, pmesh) || v0 == source(he, pmesh));
      }
    }
  }

  /*!
  \ingroup PkgPolygonMeshProcessing
  @brief triangulates and refines a hole in a polygon mesh.

  @tparam PolygonMesh must be model of `MutableFaceGraph`
          that has an internal property map for `CGAL::vertex_point_t`
  @tparam FacetOutputIterator model of `OutputIterator`
     holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces.
  @tparam VertexOutputIterator model of `OutputIterator`
     holding `boost::graph_traits<PolygonMesh>::%vertex_descriptor` for patch vertices.
  @tparam NamedParameters a sequence of \ref namedparameters

  @param pmesh polygon mesh which has the hole
  @param border_halfedge a border halfedge incident to the hole
  @param face_out output iterator over patch faces
  @param vertex_out output iterator over patch vertices without including the boundary
  @param np optional sequence of \ref namedparameters among the ones listed below

  \cgalNamedParamsBegin
     \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh` \cgalParamEnd
     \cgalParamBegin{density_control_factor} factor to control density of the ouput mesh, where larger values cause denser refinements, as in `refine()` \cgalParamEnd
     \cgalParamBegin{use_delaunay_triangulation} if `true`, use the Delaunay triangulation facet search space \cgalParamEnd
     \cgalParamBegin{geom_traits} a geometric traits class instance \cgalParamEnd
  \cgalNamedParamsEnd

  @return pair of `face_out` and `vertex_out`

  \sa CGAL::Polygon_mesh_processing::triangulate_hole()
  \sa CGAL::Polygon_mesh_processing::refine()

  \todo handle islands
  */
  template<typename PolygonMesh,
           typename FaceOutputIterator,
           typename VertexOutputIterator,
           typename NamedParameters>
  std::pair<FaceOutputIterator, VertexOutputIterator>
    triangulate_and_refine_hole(PolygonMesh& pmesh,
      typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
      FaceOutputIterator face_out,
      VertexOutputIterator vertex_out,
      const NamedParameters& np)
  {
    std::vector<typename boost::graph_traits<PolygonMesh>::face_descriptor> patch;
    triangulate_hole(pmesh, border_halfedge, std::back_inserter(patch), np);
    face_out = std::copy(patch.begin(), patch.end(), face_out);

    test_in_edges(pmesh, vertices(pmesh));

    return refine(pmesh, patch, face_out, vertex_out, np);
  }

  template<typename PolygonMesh,
           typename FaceOutputIterator,
           typename VertexOutputIterator>
  std::pair<FaceOutputIterator, VertexOutputIterator>
    triangulate_and_refine_hole(PolygonMesh& pmesh,
       typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
       FaceOutputIterator face_out,
       VertexOutputIterator vertex_out)
  {
    return triangulate_and_refine_hole(pmesh, border_halfedge,
      face_out, vertex_out,
      CGAL::Polygon_mesh_processing::parameters::all_default());
  }

  /*!
  \ingroup PkgPolygonMeshProcessing
  @brief triangulates, refines and fairs a hole in a polygon mesh.

  @tparam PolygonMesh a model of `MutableFaceGraph`
          that has an internal property map for `CGAL::vertex_point_t`
  @tparam FaceOutputIterator model of `OutputIterator`
      holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces
  @tparam VertexOutputIterator model of `OutputIterator`
      holding `boost::graph_traits<PolygonMesh>::%vertex_descriptor` for patch vertices
  @tparam NamedParameters a sequence of \ref namedparameters

  @param pmesh polygon mesh which has the hole
  @param border_halfedge a border halfedge incident to the hole
  @param face_out output iterator over patch faces
  @param vertex_out output iterator over patch vertices without including the boundary
  @param np optional sequence of \ref namedparameters among the ones listed below

  \cgalNamedParamsBegin
     \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh` \cgalParamEnd
     \cgalParamBegin{use_delaunay_triangulation} if `true`, use the Delaunay triangulation facet search space \cgalParamEnd
     \cgalParamBegin{density_control_factor} factor to control density of the ouput mesh, where larger values cause denser refinements, as in `refine()` \cgalParamEnd
     \cgalParamBegin{fairing_continuity} tangential continuity of the output surface patch \cgalParamEnd
     \cgalParamBegin{sparse_linear_solver} an instance of the sparse linear solver used for fairing \cgalParamEnd
     \cgalParamBegin{geom_traits} a geometric traits class instance \cgalParamEnd
  \cgalNamedParamsEnd

  @return tuple of
  - `bool`: `true` if fairing is successful
  - `face_out`
  - `vertex_out`

  \sa CGAL::Polygon_mesh_processing::triangulate_hole()
  \sa CGAL::Polygon_mesh_processing::refine()
  \sa CGAL::Polygon_mesh_processing::fair()

  \todo handle islands
  */
  template<typename PolygonMesh,
           typename FaceOutputIterator,
           typename VertexOutputIterator,
           typename NamedParameters>
  CGAL::cpp11::tuple<bool, FaceOutputIterator, VertexOutputIterator>
  triangulate_refine_and_fair_hole(PolygonMesh& pmesh,
    typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
    FaceOutputIterator face_out,
    VertexOutputIterator vertex_out,
    const NamedParameters& np)
  {
    std::vector<typename boost::graph_traits<PolygonMesh>::vertex_descriptor> patch;

    CGAL_assertion(CGAL::is_triangle_mesh(pmesh));

    face_out = triangulate_and_refine_hole
      (pmesh, border_halfedge, face_out, std::back_inserter(patch), np).first;

    CGAL_assertion(CGAL::is_triangle_mesh(pmesh));

    test_in_edges(pmesh, patch);

    bool fair_success = fair(pmesh, patch, np);

    vertex_out = std::copy(patch.begin(), patch.end(), vertex_out);
    return CGAL::cpp11::make_tuple(fair_success, face_out, vertex_out);
  }

  template<typename PolygonMesh,
           typename FaceOutputIterator,
           typename VertexOutputIterator>
  CGAL::cpp11::tuple<bool, FaceOutputIterator, VertexOutputIterator>
  triangulate_refine_and_fair_hole(PolygonMesh& pmesh,
        typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
        FaceOutputIterator face_out,
        VertexOutputIterator vertex_out)
  {
    return triangulate_refine_and_fair_hole(pmesh, border_halfedge,
      face_out, vertex_out,
      CGAL::Polygon_mesh_processing::parameters::all_default());
  }

  /*!
  \ingroup PkgPolygonMeshProcessing
  creates triangles to fill the hole defined by points in the range `points`.
  Triangles are recorded into `out` using the indices of the input points in the range `points`.
  Note that no degenerate triangles will be produced.
  If no triangulation can be found, then nothing is recorded in `out`.

  The point range `third_points` indicates for each pair of consecutive points in the range `points`,
  the third point of the facet this segment is incident to.

  Note that the ranges `points` and `third_points` may or may not contain duplicated first point at the end of sequence.

  @tparam OutputIteratorValueType value type of `OutputIterator`
    having a constructor `OutputIteratorValueType(int p0, int p1, int p2)` available.
    It defaults to `value_type_traits<OutputIterator>::%type`, and can be omitted when the default is fine.
  @tparam PointRange range of points, model of `Range`.
    Its iterator type is `InputIterator`.
  @tparam OutputIterator model of `OutputIterator`
     holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces
  @tparam NamedParameters a sequence of \ref namedparameters

  @param points the range of input points
  @param third_points the range of third points
  @param out iterator over output patch triangles
  @param np optional sequence of \ref namedparameters among the ones listed below

  \cgalNamedParamsBegin
     \cgalParamBegin{use_delaunay_triangulation} if `true`, use the Delaunay triangulation facet search space \cgalParamEnd
     \cgalParamBegin{geom_traits} a geometric traits class instance \cgalParamEnd
  \cgalNamedParamsEnd

  \todo handle islands
  */
  template </*typename OutputIteratorValueType,*/
            typename PointRange,
            typename OutputIterator,
            typename NamedParameters>
  OutputIterator
  triangulate_hole_polyline(const PointRange& points,
                            const PointRange& third_points,
                            OutputIterator out,
                            const NamedParameters& np)
  {
    using boost::choose_param;
    using boost::get_param;

    bool use_dt3 =
#ifdef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
      false;
#else
      choose_param(get_param(np, use_delaunay_triangulation), true);
#endif

    typedef CGAL::internal::Weight_min_max_dihedral_and_area      Weight;
    typedef CGAL::internal::Weight_calculator<Weight,
                  CGAL::internal::Is_not_degenerate_triangle>  WC;
    typedef std::vector<std::pair<int, int> > Holes;
    typedef std::back_insert_iterator<Holes>  Holes_out;

    Holes holes;//just to check there is no holes left

    typedef typename value_type_traits<OutputIterator>::type OutputIteratorValueType;
    CGAL::internal::Tracer_polyline_incomplete<OutputIteratorValueType, OutputIterator, Holes_out>
      tracer(out, Holes_out(holes));
    
    typedef typename PointRange::iterator InIterator;
    typedef typename std::iterator_traits<InIterator>::value_type Point;

    triangulate_hole_polyline(points, third_points, tracer, WC(),
      use_dt3,
      choose_param(get_param(np, CGAL::geom_traits),
        typename CGAL::Kernel_traits<Point>::Kernel()));

    CGAL_assertion(holes.empty());
    return tracer.out;
  }

  template <typename PointRange,
            typename OutputIterator>
  OutputIterator
  triangulate_hole_polyline(const PointRange& points,
                            const PointRange& third_points,
                            OutputIterator out)
  {
    return triangulate_hole_polyline(points, third_points, out,
      CGAL::Polygon_mesh_processing::parameters::all_default());
  }

  /*!
  \ingroup PkgPolygonMeshProcessing
  same as above but the range of third points is omitted. They are not
  taken into account in the cost computation that leads the hole filling. 
*/
  template <typename PointRange,
            typename OutputIterator,
            typename NamedParameters>
  OutputIterator
  triangulate_hole_polyline(const PointRange& points,
                            OutputIterator out,
                            const NamedParameters& np)
  {
    return triangulate_hole_polyline(points, PointRange(), out, np);
  }

  template <typename PointRange,
            typename OutputIterator>
  OutputIterator
  triangulate_hole_polyline(const PointRange& points,
                            OutputIterator out)
  {
    return triangulate_hole_polyline(points, out,
      CGAL::Polygon_mesh_processing::parameters::all_default());
  }

} //end namespace Polygon_mesh_processing

} //end namespace CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_TRIANGULATE_HOLE_H
