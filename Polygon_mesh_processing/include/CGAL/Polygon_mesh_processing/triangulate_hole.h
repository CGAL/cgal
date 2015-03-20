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
#include <CGAL/Default.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>

#include <boost/tuple/tuple.hpp>

namespace CGAL {

namespace Polygon_mesh_processing {

  /*!
  \ingroup PkgPolygonMeshProcessing
  triangulates a hole in a polygon mesh.
  The hole must not contain any non-manifold vertex.
  The patch generated does not introduce non-manifold edges nor degenerate triangles.
  If a hole cannot be triangulated, `pmesh` is not modified and nothing is put in `out`.

  @tparam PolygonMesh a model of `MutableFaceGraph`
  @tparam OutputIterator a model of `OutputIterator`
    holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces.

  @param pmesh polygon mesh containing the hole
  @param border_halfedge a border halfedge incident to the hole
  @param out iterator over patch faces
  @param use_delaunay_triangulation if `true`, use the Delaunay triangulation facet search space

  @return `out`

  \todo handle islands
  @todo Replace border_halfedge by a range of border halfedges.
        The first one would describe the hole,
        the other ones would describe the islands.
  @todo Then, insert the holes vertices in the set of possibilities
        for connecting vertices together
  @todo handle the case where an island is reduced to a point
  \todo SUBMISSION: VertexPointMap
  */
  template<typename PolygonMesh,
           typename OutputIterator,
           class P, class T, class R>
  OutputIterator
    triangulate_hole(PolygonMesh& pmesh,
      typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
      OutputIterator out,
      const pmp_bgl_named_params<P, T, R>& p)
  {
    using boost::choose_param;
    using boost::get_param;

    bool use_dt3 =
#ifdef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
      false;
#else
      choose_param(get_param(p, use_delaunay_triangulation), true);
#endif

    CGAL_precondition(face(border_halfedge, pmesh) == boost::graph_traits<PolygonMesh>::null_face());

    return internal::triangulate_hole_polygon_mesh(pmesh, border_halfedge, out, use_dt3).first;
  }

  template<typename PolygonMesh, typename OutputIterator>
  OutputIterator
    triangulate_hole(PolygonMesh& pmesh,
      typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
      OutputIterator out)
  {
    return triangulate_hole(pmesh, border_halfedge, out, CGAL::parameters::all_default());
  }

  /*!
  \ingroup PkgPolygonMeshProcessing
  @brief triangulates and refines a hole in a polygon mesh.

  @tparam PolygonMesh must be model of `MutableFaceGraph`
  @tparam FacetOutputIterator model of `OutputIterator`
     holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces.
  @tparam VertexOutputIterator model of `OutputIterator`
     holding `boost::graph_traits<PolygonMesh>::%vertex_descriptor` for patch vertices.

  @param pmesh polygon mesh which has the hole
  @param border_halfedge a border halfedge incident to the hole
  @param face_out output iterator over patch faces
  @param vertex_out output iterator over patch vertices without including the boundary
  @param density_control_factor factor for density where larger values cause denser refinements
  @param use_delaunay_triangulation if `true`, use the Delaunay triangulation face search space

  @return pair of `face_out` and `vertex_out`

  \sa CGAL::Polygon_mesh_processing::triangulate_hole()
  \sa CGAL::Polygon_mesh_processing::refine()

  \todo SUBMISSION: VertexPointMap
  \todo SUBMISSION: better document density_control_factor (ideally we should refer to the doc of refine)
  \todo handle islands
  */
  template<class PolygonMesh,
           class FaceOutputIterator,
           class VertexOutputIterator,
           class P, class T, class R>
  std::pair<FaceOutputIterator, VertexOutputIterator>
    triangulate_and_refine_hole(PolygonMesh& pmesh,
      typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
      FaceOutputIterator face_out,
      VertexOutputIterator vertex_out,
      const pmp_bgl_named_params<P, T, R>& p)
  {
    using boost::choose_param;
    using boost::get_param;

    std::vector<typename boost::graph_traits<PolygonMesh>::face_descriptor> patch;
    triangulate_hole(pmesh, border_halfedge, std::back_inserter(patch), p);
    face_out = std::copy(patch.begin(), patch.end(), face_out);

    return refine(pmesh, patch, face_out, vertex_out, p);
  }

  template<class PolygonMesh,
           class FaceOutputIterator,
           class VertexOutputIterator>
  std::pair<FaceOutputIterator, VertexOutputIterator>
    triangulate_and_refine_hole(PolygonMesh& pmesh,
       typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
       FaceOutputIterator face_out,
       VertexOutputIterator vertex_out)
  {
    return triangulate_and_refine_hole(pmesh, border_halfedge,
      face_out, vertex_out,
      CGAL::parameters::all_default());
  }


  /*!
  \ingroup PkgPolygonMeshProcessing
  @brief triangulates, refines and fairs a hole in a polygon mesh.

  If \ref thirdpartyEigen "Eigen" 3.2 (or greater) is available
  and `CGAL_EIGEN3_ENABLED` is defined, an overload of this function is available
  with `SparseLinearSolver` being:
  \code
  CGAL::Eigen_solver_traits<
  Eigen::SparseLU<
  CGAL::Eigen_sparse_matrix<double>::EigenType,
  Eigen::COLAMDOrdering<int> >  >
  \endcode


  @tparam SparseLinearSolver a model of `SparseLinearAlgebraTraitsWithFactor_d`
  @tparam PolygonMesh a model of `MutableFaceGraph`
  @tparam FaceOutputIterator model of `OutputIterator`
      holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces
  @tparam VertexOutputIterator model of `OutputIterator`
      holding `boost::graph_traits<PolygonMesh>::%vertex_descriptor` for patch vertices

  @param pmesh a polygon mesh which has the hole
  @param border_halfedge a border halfedge incident to the hole
  @param face_out iterator over patch faces
  @param vertex_out iterator over patch vertices without including the boundary
  @param density_control_factor factor for density where larger values cause denser refinements
  @param continuity tangential continuity, defaults to `FAIRING_C_1` and can be omitted
  @param use_delaunay_triangulation if `true`, the Delaunay triangulation face search space is used
  @param solver an instance of the sparse linear solver to use. It defaults to the 
         default construtor of the `SparseLinearSolver` template parameter

  @return tuple of
  - bool: `true` if fairing is successful
  - `face_out`
  - `vertex_out`

  \sa CGAL::Polygon_mesh_processing::triangulate_hole()
  \sa CGAL::Polygon_mesh_processing::refine()
  \sa CGAL::Polygon_mesh_processing::fair()

  \todo handle islands
  \todo SUBMISSION: VertexPointMap
  */
  template<typename PolygonMesh,
           typename FaceOutputIterator,
           typename VertexOutputIterator,
           class P, class T, class R>
  CGAL::cpp11::tuple<bool, FaceOutputIterator, VertexOutputIterator>
  triangulate_refine_and_fair_hole(PolygonMesh& pmesh,
    typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
    FaceOutputIterator face_out,
    VertexOutputIterator vertex_out,
    const pmp_bgl_named_params<P, T, R>& p)
  {
    std::vector<typename boost::graph_traits<PolygonMesh>::vertex_descriptor> patch;

    face_out = triangulate_and_refine_hole
      (pmesh, border_halfedge, face_out, std::back_inserter(patch), p).first;

    bool fair_success = fair(pmesh, patch, p);

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
      CGAL::parameters::all_default());
  }

  /*!
  \ingroup PkgPolygonMeshProcessing
  creates triangles to fill the hole defined by points in the range (`points`).
  Triangles are put into `out`
  using the indices of the input points in the range (`points`).
  Note that no degenerate triangles will be produced.
  If no triangulation can be found, then nothing is put in `out`.

  The optional range (`third_points`) indicates for each pair of consecutive points in the range (`points`),
  the third point of the facet this segment is incident to.

  Note that the ranges (`points`) and (`third_points`) may or may not contain duplicated first point at the end of sequence.

  @tparam OutputIteratorValueType value type of `OutputIterator`
    having a constructor `OutputIteratorValueType(int p0, int p1, int p2)` available.
    It defaults to `value_type_traits<OutputIterator>::%type`, and can be omitted when the default is fine.

  @tparam PointRange range of points, model of `SinglePassRange`
  @tparam OutputIterator model of `OutputIterator`
     holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces
  
  @param points the range of input points
  @param third_points the range of third points, can be omitted
  @param out iterator over output patch triangles
  @param use_delaunay_triangulation if `true`, use the Delaunay triangulation facet search space, defaults to true if omitted.

  \todo handle islands
  */
  template <typename OutputIteratorValueType,
            typename PointRange,
            typename OutputIterator>
  OutputIterator
    triangulate_hole_polyline(const PointRange& points,
                              const PointRange& third_points,
                              OutputIterator out,
                              bool use_delaunay_triangulation = true)
  {
    bool use_dt3 =
#ifdef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
      false;
#else
      use_delaunay_triangulation;
#endif

    typedef CGAL::internal::Weight_min_max_dihedral_and_area      Weight;
    typedef CGAL::internal::Weight_calculator<Weight,
                  CGAL::internal::Is_not_degenerate_triangle>  WC;
    typedef std::vector<std::pair<int, int> > Holes;
    typedef std::back_insert_iterator<Holes>  Holes_out;

    Holes holes;//just to check there is no holes left

    CGAL::internal::Tracer_polyline_incomplete<OutputIteratorValueType, OutputIterator, Holes_out>
      tracer(out, Holes_out(holes));
    triangulate_hole_polyline(points, third_points, tracer, WC(), use_dt3);
    CGAL_assertion(holes.empty());
    return tracer.out;
  }

  // overload for OutputIteratorValueType
  template <typename PointRange,
            typename OutputIterator>
  OutputIterator
  triangulate_hole_polyline(const PointRange& points,
                            const PointRange& third_points,
                            OutputIterator out,
                            bool use_delaunay_triangulation = true)
  {
    bool use_dt3 =
#ifdef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
      false;
#else
      use_delaunay_triangulation;
#endif
    return triangulate_hole_polyline<typename value_type_traits<OutputIterator>::type>
      (points, third_points, out, use_dt3);
  }

  // overload no (qbegin, qend)
  template <typename OutputIteratorValueType,
            typename PointRange,
            typename OutputIterator>
  OutputIterator
    triangulate_hole_polyline(const PointRange& points,
                              OutputIterator out,
                              bool use_delaunay_triangulation = true)
  {
    bool use_dt3 =
#ifdef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
      false;
#else
      use_delaunay_triangulation;
#endif

    return triangulate_hole_polyline<OutputIteratorValueType>
      (points, PointRange(), out, use_dt3);
  }

  // overload for OutputIteratorValueType
  template <typename PointRange,
            typename OutputIterator>
  OutputIterator
    triangulate_hole_polyline(const PointRange& points,
                              OutputIterator out,
                              bool use_delaunay_triangulation = true)
  {
    bool use_dt3 =
#ifdef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
      false;
#else
      use_delaunay_triangulation;
#endif
    return triangulate_hole_polyline<typename value_type_traits<OutputIterator>::type>
      (points, out, use_dt3);
  }

} //end namespace Polygon_mesh_processing

} //end namespace CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_TRIANGULATE_HOLE_H
