// Copyright (c) 2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Ilker O. Yaz

#ifndef CGAL_POLYGON_MESH_PROCESSING_TRIANGULATE_HOLE_H
#define CGAL_POLYGON_MESH_PROCESSING_TRIANGULATE_HOLE_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/bounding_box.h>
#include <CGAL/Polygon_mesh_processing/internal/Hole_filling/Triangulate_hole_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/internal/Hole_filling/Triangulate_hole_polyline.h>
#include <CGAL/Polygon_mesh_processing/refine.h>
#include <CGAL/Polygon_mesh_processing/fair.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <CGAL/boost/graph/helpers.h>

#include <boost/tuple/tuple.hpp>

#include <vector>

#ifdef DOXYGEN_RUNNING
#define CGAL_PMP_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_PMP_NP_CLASS NamedParameters
#endif

namespace CGAL {

namespace Polygon_mesh_processing {

  /*!
  \ingroup hole_filling_grp
  triangulates a hole in a polygon mesh.

  Depending on the choice of the underlying algorithm different preconditions apply.
  When using the 2D constrained Delaunay triangulation, the border edges of the hole
  must not intersect the surface. Otherwise, additionally, the boundary
  of the hole must not contain any non-manifold vertex. The patch generated does not
  introduce non-manifold edges nor degenerate triangles. If a hole cannot be triangulated,
  `pmesh` is not modified and nothing is recorded in `out`.

  @tparam PolygonMesh a model of `MutableFaceGraph`
  @tparam OutputIterator a model of `OutputIterator`
    holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces.
  @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  @param pmesh polygon mesh containing the hole
  @param border_halfedge a border halfedge incident to the hole
  @param out iterator over patch faces
  @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
                     as key type and `%Point_3` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      must be available in `PolygonMesh`.}
    \cgalParamNEnd

    \cgalParamNBegin{geom_traits}
      \cgalParamDescription{an instance of a geometric traits class}
      \cgalParamType{a class model of `Kernel`}
      \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
      \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
    \cgalParamNEnd

    \cgalParamNBegin{use_delaunay_triangulation}
      \cgalParamDescription{If `true`, use the Delaunay triangulation facet search space.}
      \cgalParamType{Boolean}
      \cgalParamDefault{`true`}
      \cgalParamExtra{If no valid triangulation can be found in this search space, the algorithm
                      falls back to the non-Delaunay triangulations search space to find a solution.}
    \cgalParamNEnd

    \cgalParamNBegin{use_2d_constrained_delaunay_triangulation}
      \cgalParamDescription{If `true`, the points of the boundary of the hole are used
                            to estimate a fitting plane and a 2D constrained Delaunay triangulation
                            is then used to fill the hole projected in the fitting plane.}
      \cgalParamType{Boolean}
      \cgalParamDefault{`true`}
      \cgalParamExtra{If the boundary of the hole is not planar (according to the
                      parameter `threshold_distance`) or if no valid 2D triangulation
                      can be found, the algorithm falls back to the method using
                      the 3D Delaunay triangulation. This parameter is a good choice for near planar holes.}
    \cgalParamNEnd

    \cgalParamNBegin{threshold_distance}
      \cgalParamDescription{The maximum distance between the vertices of
                            the hole boundary and the least squares plane fitted to this boundary.}
      \cgalParamType{double}
      \cgalParamDefault{one quarter of the height of the bounding box of the hole}
      \cgalParamExtra{This parameter is used only in conjunction with
                      the parameter `use_2d_constrained_delaunay_triangulation`.}
    \cgalParamNEnd
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
    using parameters::choose_parameter;
    using parameters::get_parameter;

    typedef typename GetGeomTraits<PolygonMesh,NamedParameters>::type         GeomTraits;

    bool use_dt3 =
#ifdef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
      false;
#else
      choose_parameter(get_parameter(np, internal_np::use_delaunay_triangulation), true);
#endif

    CGAL_precondition(face(border_halfedge, pmesh) == boost::graph_traits<PolygonMesh>::null_face());
    bool use_cdt =
    #ifdef CGAL_HOLE_FILLING_DO_NOT_USE_CDT2
        false;
#else
        choose_parameter(get_parameter(np, internal_np::use_2d_constrained_delaunay_triangulation), false);
#endif

    typename GeomTraits::FT max_squared_distance = typename GeomTraits::FT(-1);
    if (use_cdt) {

      std::vector<typename GeomTraits::Point_3> points;
      typedef Halfedge_around_face_circulator<PolygonMesh> Hedge_around_face_circulator;
      const auto vpmap = choose_parameter(get_parameter(np, internal_np::vertex_point), get_property_map(vertex_point, pmesh));
      Hedge_around_face_circulator circ(border_halfedge, pmesh), done(circ);
      do {
        points.push_back(get(vpmap, target(*circ, pmesh)));
      } while (++circ != done);

      const typename GeomTraits::Iso_cuboid_3 bbox = CGAL::bounding_box(points.begin(), points.end());
      typename GeomTraits::FT default_squared_distance = CGAL::abs(CGAL::squared_distance(bbox.vertex(0), bbox.vertex(5)));
      default_squared_distance /= typename GeomTraits::FT(16); // one quarter of the bbox height

      const typename GeomTraits::FT threshold_distance = choose_parameter(
        get_parameter(np, internal_np::threshold_distance), typename GeomTraits::FT(-1));
      max_squared_distance = default_squared_distance;
      if (threshold_distance >= typename GeomTraits::FT(0))
        max_squared_distance = threshold_distance * threshold_distance;
      CGAL_assertion(max_squared_distance >= typename GeomTraits::FT(0));
    }

    return internal::triangulate_hole_polygon_mesh(
      pmesh,
      border_halfedge,
      out,
      choose_parameter(get_parameter(np, internal_np::vertex_point), get_property_map(vertex_point, pmesh)),
      use_dt3,
      choose_parameter<GeomTraits>(get_parameter(np, internal_np::geom_traits)),
      use_cdt,
      max_squared_distance).first;
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
    for(typename boost::graph_traits<PM>::vertex_descriptor v0 : patch)
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
  \ingroup  hole_filling_grp
  @brief triangulates and refines a hole in a polygon mesh.

  @tparam PolygonMesh must be model of `MutableFaceGraph`
  @tparam FacetOutputIterator model of `OutputIterator`
     holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces.
  @tparam VertexOutputIterator model of `OutputIterator`
     holding `boost::graph_traits<PolygonMesh>::%vertex_descriptor` for patch vertices.
  @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  @param pmesh polygon mesh which has the hole
  @param border_halfedge a border halfedge incident to the hole
  @param face_out output iterator over patch faces
  @param vertex_out output iterator over patch vertices without including the boundary
  @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
                     as key type and `%Point_3` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      must be available in `PolygonMesh`.}
    \cgalParamNEnd

    \cgalParamNBegin{geom_traits}
      \cgalParamDescription{an instance of a geometric traits class}
      \cgalParamType{a class model of `Kernel`}
      \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
      \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
    \cgalParamNEnd

    \cgalParamNBegin{use_delaunay_triangulation}
      \cgalParamDescription{If `true`, use the Delaunay triangulation facet search space.}
      \cgalParamType{Boolean}
      \cgalParamDefault{`true`}
      \cgalParamExtra{If no valid triangulation can be found in this search space, the algorithm
                      falls back to the non-Delaunay triangulations search space to find a solution.}
    \cgalParamNEnd

    \cgalParamNBegin{use_2d_constrained_delaunay_triangulation}
      \cgalParamDescription{If `true`, the points of the boundary of the hole are used
                            to estimate a fitting plane and a 2D constrained Delaunay triangulation
                            is then used to fill the hole projected in the fitting plane.}
      \cgalParamType{Boolean}
      \cgalParamDefault{`true`}
      \cgalParamExtra{If the boundary of the hole is not planar (according to the
                      parameter `threshold_distance`) or if no valid 2D triangulation
                      can be found, the algorithm falls back to the method using
                      the 3D Delaunay triangulation. This parameter is a good choice for near planar holes.}
    \cgalParamNEnd

    \cgalParamNBegin{threshold_distance}
      \cgalParamDescription{The maximum distance between the vertices of
                            the hole boundary and the least squares plane fitted to this boundary.}
      \cgalParamType{double}
      \cgalParamDefault{one quarter of the height of the bounding box of the hole}
      \cgalParamExtra{This parameter is used only in conjunction with
                      the parameter `use_2d_constrained_delaunay_triangulation`.}
    \cgalParamNEnd

    \cgalParamNBegin{density_control_factor}
      \cgalParamDescription{factor to control density of the ouput mesh,
                            where larger values cause denser refinements, as in `refine()`}
      \cgalParamType{double}
      \cgalParamDefault{\f$ \sqrt{2}\f$}
    \cgalParamNEnd
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
  \ingroup  hole_filling_grp
  @brief triangulates, refines and fairs a hole in a polygon mesh.

  @tparam PolygonMesh a model of `MutableFaceGraph`
  @tparam FaceOutputIterator model of `OutputIterator`
      holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces
  @tparam VertexOutputIterator model of `OutputIterator`
      holding `boost::graph_traits<PolygonMesh>::%vertex_descriptor` for patch vertices
  @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  @param pmesh polygon mesh which has the hole
  @param border_halfedge a border halfedge incident to the hole
  @param face_out output iterator over patch faces
  @param vertex_out output iterator over patch vertices without including the boundary
  @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
                     as key type and `%Point_3` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      must be available in `PolygonMesh`.}
    \cgalParamNEnd

    \cgalParamNBegin{geom_traits}
      \cgalParamDescription{an instance of a geometric traits class}
      \cgalParamType{a class model of `Kernel`}
      \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
      \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
    \cgalParamNEnd

    \cgalParamNBegin{use_delaunay_triangulation}
      \cgalParamDescription{If `true`, use the Delaunay triangulation facet search space.}
      \cgalParamType{Boolean}
      \cgalParamDefault{`true`}
      \cgalParamExtra{If no valid triangulation can be found in this search space, the algorithm
                      falls back to the non-Delaunay triangulations search space to find a solution.}
    \cgalParamNEnd

    \cgalParamNBegin{use_2d_constrained_delaunay_triangulation}
      \cgalParamDescription{If `true`, the points of the boundary of the hole are used
                            to estimate a fitting plane and a 2D constrained Delaunay triangulation
                            is then used to fill the hole projected in the fitting plane.}
      \cgalParamType{Boolean}
      \cgalParamDefault{`true`}
      \cgalParamExtra{If the boundary of the hole is not planar (according to the
                      parameter `threshold_distance`) or if no valid 2D triangulation
                      can be found, the algorithm falls back to the method using
                      the 3D Delaunay triangulation. This parameter is a good choice for near planar holes.}
    \cgalParamNEnd

    \cgalParamNBegin{threshold_distance}
      \cgalParamDescription{The maximum distance between the vertices of
                            the hole boundary and the least squares plane fitted to this boundary.}
      \cgalParamType{double}
      \cgalParamDefault{one quarter of the height of the bounding box of the hole}
      \cgalParamExtra{This parameter is used only in conjunction with
                      the parameter `use_2d_constrained_delaunay_triangulation`.}
    \cgalParamNEnd

    \cgalParamNBegin{density_control_factor}
      \cgalParamDescription{factor to control density of the ouput mesh,
                            where larger values cause denser refinements, as in `refine()`}
      \cgalParamType{double}
      \cgalParamDefault{\f$ \sqrt{2}\f$}
    \cgalParamNEnd

    \cgalParamNBegin{fairing_continuity}
      \cgalParamDescription{A value controling the tangential continuity of the output surface patch.
                            The possible values are 0, 1 and 2, refering to the  C<sup>0</sup>, C<sup>1</sup>
                            and C<sup>2</sup> continuity.}
      \cgalParamType{unsigned int}
      \cgalParamDefault{`1`}
      \cgalParamExtra{The larger `fairing_continuity` gets, the more fixed vertices are required.}
    \cgalParamNEnd

    \cgalParamNBegin{sparse_linear_solver}
      \cgalParamDescription{an instance of the sparse linear solver used for fairing}
      \cgalParamType{a class model of `SparseLinearAlgebraWithFactorTraits_d`}
      \cgalParamDefault{If \ref thirdpartyEigen "Eigen" 3.2 (or greater) is available and
                        `CGAL_EIGEN3_ENABLED` is defined, then the following overload of `Eigen_solver_traits`
                        is provided as default value:\n
                        `CGAL::Eigen_solver_traits<Eigen::SparseLU<CGAL::Eigen_sparse_matrix<double>::%EigenType, Eigen::COLAMDOrdering<int> > >`}
    \cgalParamNEnd
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
  std::tuple<bool, FaceOutputIterator, VertexOutputIterator>
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
    return std::make_tuple(fair_success, face_out, vertex_out);
  }

  template<typename PolygonMesh,
           typename FaceOutputIterator,
           typename VertexOutputIterator>
  std::tuple<bool, FaceOutputIterator, VertexOutputIterator>
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
  \ingroup  hole_filling_grp
  creates triangles to fill the hole defined by points in the range `points`.
  Triangles are recorded into `out` using the indices of the input points in the range `points`.
  Note that no degenerate triangles will be produced.
  If no triangulation can be found, then nothing is recorded in `out`.

  If faces incident to the polyline outside the hole are known,
  it is recommended to use this function.
  The point range `third_points` indicates for each pair of consecutive points in the range `points`,
  the third point of the face this segment is incident to. It influences the choice
  of the best triangulation while avoiding overfolding.

  Note that the ranges `points` and `third_points` may or may not contain duplicated first point at the end of sequence.

  @pre `third_points.size() == points.size()`

  @tparam PointRange range of points, model of `Range`.
    Its iterator type is `InputIterator`.
  @tparam OutputIterator model of `OutputIterator`, to collect patch faces.
     A specialization for `CGAL::value_type_traits<OutputIterator>` must be available,
     and the corresponding value type `type` must have
     a constructor `type(int p0, int p1, int p2)` available.
     The indices correspond to the ones of input points in `points`.
  @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  @param points the range of input points
  @param third_points the range of third points
  @param out iterator over output patch triangles, described by indices of points
             in `points`
  @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamNBegin{geom_traits}
      \cgalParamDescription{an instance of a geometric traits class}
      \cgalParamType{a class model of `Kernel`}
      \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
      \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
    \cgalParamNEnd

    \cgalParamNBegin{use_delaunay_triangulation}
      \cgalParamDescription{If `true`, use the Delaunay triangulation facet search space.}
      \cgalParamType{Boolean}
      \cgalParamDefault{`true`}
      \cgalParamExtra{If no valid triangulation can be found in this search space, the algorithm
                      falls back to the non-Delaunay triangulations search space to find a solution.}
    \cgalParamNEnd

    \cgalParamNBegin{use_2d_constrained_delaunay_triangulation}
      \cgalParamDescription{If `true`, the points of the boundary of the hole are used
                            to estimate a fitting plane and a 2D constrained Delaunay triangulation
                            is then used to fill the hole projected in the fitting plane.}
      \cgalParamType{Boolean}
      \cgalParamDefault{`true`}
      \cgalParamExtra{If the boundary of the hole is not planar (according to the
                      parameter `threshold_distance`) or if no valid 2D triangulation
                      can be found, the algorithm falls back to the method using
                      the 3D Delaunay triangulation. This parameter is a good choice for near planar holes.}
    \cgalParamNEnd

    \cgalParamNBegin{threshold_distance}
      \cgalParamDescription{The maximum distance between the vertices of
                            the hole boundary and the least squares plane fitted to this boundary.}
      \cgalParamType{double}
      \cgalParamDefault{one quarter of the height of the bounding box of the hole}
      \cgalParamExtra{This parameter is used only in conjunction with
                      the parameter `use_2d_constrained_delaunay_triangulation`.}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \todo handle islands
  */
  template <typename PointRange1,
            typename PointRange2,
            typename OutputIterator,
            typename NamedParameters>
  OutputIterator
  triangulate_hole_polyline(const PointRange1& points,
                            const PointRange2& third_points,
                            OutputIterator out,
                            const NamedParameters& np)
  {
    if (points.empty()) return out;

    using parameters::choose_parameter;
    using parameters::get_parameter;

    bool use_cdt =
#ifdef CGAL_HOLE_FILLING_DO_NOT_USE_CDT2
      false;
#else
      choose_parameter(get_parameter(np, internal_np::use_2d_constrained_delaunay_triangulation), true);
#endif
bool use_dt3 =
#ifdef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
      false;
#else
      choose_parameter(get_parameter(np, internal_np::use_delaunay_triangulation), true);
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

    typedef typename PointRange1::iterator InIterator;
    typedef typename std::iterator_traits<InIterator>::value_type Point;
    typedef typename CGAL::Kernel_traits<Point>::Kernel Kernel;
#ifndef CGAL_HOLE_FILLING_DO_NOT_USE_CDT2
    struct Always_valid{
      bool operator()(const std::vector<Point>&, int,int,int)const
      {return true;}
    };
    Always_valid is_valid;

    const typename Kernel::Iso_cuboid_3 bbox = CGAL::bounding_box(points.begin(), points.end());
    typename Kernel::FT default_squared_distance = CGAL::abs(CGAL::squared_distance(bbox.vertex(0), bbox.vertex(5)));
    default_squared_distance /= typename Kernel::FT(16); // one quarter of the bbox height

    const typename Kernel::FT threshold_distance = choose_parameter(
      get_parameter(np, internal_np::threshold_distance), typename Kernel::FT(-1));
    typename Kernel::FT max_squared_distance = default_squared_distance;
    if (threshold_distance >= typename Kernel::FT(0))
      max_squared_distance = threshold_distance * threshold_distance;
    CGAL_assertion(max_squared_distance >= typename Kernel::FT(0));

    if(!use_cdt ||
       !triangulate_hole_polyline_with_cdt(
         points,
         tracer,
         is_valid,
         choose_parameter<Kernel>(get_parameter(np, internal_np::geom_traits)),
         max_squared_distance))
#endif
    triangulate_hole_polyline(points, third_points, tracer, WC(),
                              use_dt3,
                              choose_parameter<Kernel>(get_parameter(np, internal_np::geom_traits)));

    CGAL_assertion(holes.empty());
    return tracer.out;
  }

  template <typename PointRange1,
            typename PointRange2,
            typename OutputIterator>
  OutputIterator
  triangulate_hole_polyline(const PointRange1& points,
                            const PointRange2& third_points,
                            OutputIterator out)
  {
    return triangulate_hole_polyline(points, third_points, out,
      CGAL::Polygon_mesh_processing::parameters::all_default());
  }

  /*!
  \ingroup  hole_filling_grp
  Same as above but the range of third points is omitted. They are not
  taken into account in the cost computation that leads the hole filling.
*/
  template <typename PointRange,
            typename OutputIterator,
            typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
  OutputIterator
  triangulate_hole_polyline(const PointRange& points,
                            OutputIterator out,
                            const CGAL_PMP_NP_CLASS& np)
  {
    typedef typename std::iterator_traits<
      typename PointRange::iterator>::value_type Point;
    std::vector< Point > third_points;
    return triangulate_hole_polyline(points, third_points, out, np);
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

#include <CGAL/enable_warnings.h>

#endif //CGAL_POLYGON_MESH_PROCESSING_TRIANGULATE_HOLE_H
