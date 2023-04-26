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

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/boost/graph/helpers.h>

#include <tuple>

#include <vector>

namespace CGAL {

namespace Polygon_mesh_processing {

  namespace Hole_filling {
    /*! \ingroup PMP_hole_filling_grp
     *  %Default hole filling visitor model of `PMPHolefillingVisitor`.
     *  All of its functions have an empty body. This class can be used as a
     *  base class if only some of the functions of the concept require to be
     *  overridden.
     */
    struct Default_visitor{
    #ifndef DOXYGEN_RUNNING
      void start_planar_phase() const {}
      void end_planar_phase(bool) const {}
      void start_quadratic_phase(std::size_t /* n */) const {}
      void quadratic_step() const {}
      void end_quadratic_phase(bool) const {}
      void start_cubic_phase(int /* n */) const {}
      void cubic_step() const {}
      void end_cubic_phase() const {}
      void start_refine_phase() const {}
      void end_refine_phase() const {}
      void start_fair_phase() const {}
      void end_fair_phase() const {}
    #endif
    };
  } // namespace Hole_filling

  /*!
  \ingroup PMP_hole_filling_grp

  triangulates a hole in a polygon mesh.

  Depending on the choice of the underlying algorithm different preconditions apply.
  When using the 2D constrained Delaunay triangulation, the border edges of the hole
  must not intersect the surface. Otherwise, additionally, the boundary
  of the hole must not contain any non-manifold vertex. The patch generated does not
  introduce non-manifold edges nor degenerate triangles. If a hole cannot be triangulated,
  `pmesh` is not modified and nothing is recorded in the face output
  iterator.

  @tparam PolygonMesh a model of `MutableFaceGraph`
  @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  @param pmesh polygon mesh containing the hole
  @param border_halfedge a border halfedge incident to the hole
  @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin

    \cgalParamNBegin{face_output_iterator}
      \cgalParamDescription{iterator over patch faces}
      \cgalParamType{a model of `OutputIterator`
    holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces}
      \cgalParamDefault{`Emptyset_iterator`}
    \cgalParamNEnd

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

    \cgalParamNBegin{do_not_use_cubic_algorithm}
      \cgalParamDescription{Set this parameter to `true` if you only want to use the Delaunay based versions of the algorithm,
                            skipping the cubic search space one in case of failure.}
      \cgalParamType{Boolean}
      \cgalParamDefault{`false`}
      \cgalParamExtra{If `true`, `use_2d_constrained_delaunay_triangulation` or `use_delaunay_triangulation` must be set to `true`
                      otherwise nothing will be done.}
    \cgalParamNEnd

    \cgalParamNBegin{visitor}
      \cgalParamDescription{a visitor used to track when entering a given phase of the algorithm}
      \cgalParamType{A model of PMPHolefillingVisitor}
      \cgalParamType{Hole_filling::Default_visitor}
    \cgalParamNEnd

  \cgalNamedParamsEnd

  @return the face output iterator

  \todo handle islands
  @todo Replace border_halfedge by a range of border halfedges.
        The first one would describe the hole,
        the other ones would describe the islands.
  @todo Then, insert the holes vertices in the set of possibilities
        for connecting vertices together
  @todo handle the case where an island is reduced to a point
  */
  template<typename PolygonMesh,
           typename CGAL_NP_TEMPLATE_PARAMETERS>
  auto
  triangulate_hole(PolygonMesh& pmesh,
              typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
              const CGAL_NP_CLASS& np = parameters::default_values())
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;
    using parameters::get_parameter_reference;

    typedef typename GetGeomTraits<PolygonMesh,CGAL_NP_CLASS>::type         GeomTraits;

    typedef typename internal_np::Lookup_named_param_def<internal_np::face_output_iterator_t,
                                                         CGAL_NP_CLASS,
                                                         Emptyset_iterator>::type Face_output_iterator;

    Face_output_iterator out = choose_parameter<Emptyset_iterator>(get_parameter(np, internal_np::face_output_iterator));

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

    Hole_filling::Default_visitor default_visitor;

    return
      internal::triangulate_hole_polygon_mesh(
        pmesh,
        border_halfedge,
        out,
        choose_parameter(get_parameter(np, internal_np::vertex_point), get_property_map(vertex_point, pmesh)),
        use_dt3,
        choose_parameter<GeomTraits>(get_parameter(np, internal_np::geom_traits)),
        use_cdt,
        choose_parameter(get_parameter(np, internal_np::do_not_use_cubic_algorithm), false),
        choose_parameter(get_parameter_reference(np, internal_np::visitor), default_visitor),
        max_squared_distance).first;
  }

#ifndef CGAL_NO_DEPRECATED_CODE
  /*!
  \ingroup PMP_hole_filling_grp

  \deprecated This function is deprecated since \cgal 5.6 and the
  overload with the named parameter `face_output_iterator` should be
  used instead.

  \brief triangulates a hole in a polygon mesh.


  @tparam PolygonMesh a model of `MutableFaceGraph`
  @tparam OutputIterator a model of `OutputIterator`
    holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces.
  @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  */
  template<typename PolygonMesh,
           typename OutputIterator,
           typename CGAL_NP_TEMPLATE_PARAMETERS>
  CGAL_DEPRECATED
  OutputIterator
  triangulate_hole(PolygonMesh& pmesh,
              typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
              OutputIterator out,
              const CGAL_NP_CLASS& np = parameters::default_values())
  {
    return triangulate_hole(pmesh, border_halfedge,np.face_output_iterator(out));
  }
#endif // CGAL_NO_DEPRECATED_CODE

  /*!
  \ingroup PMP_hole_filling_grp
  @brief triangulates and refines a hole in a polygon mesh.

  @tparam PolygonMesh must be model of `MutableFaceGraph`
  @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  @param pmesh polygon mesh which has the hole
  @param border_halfedge a border halfedge incident to the hole
  @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin

    \cgalParamNBegin{face_output_iterator}
      \cgalParamDescription{iterator over patch faces}
      \cgalParamType{a model of `OutputIterator`
    holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces}
      \cgalParamDefault{`Emptyset_iterator`}
    \cgalParamNEnd

    \cgalParamNBegin{vertex_output_iterator}
      \cgalParamDescription{iterator over patch vertices}
      \cgalParamType{a model of `OutputIterator`
    holding `boost::graph_traits<PolygonMesh>::%vertex_descriptor` for patch vertices}
      \cgalParamDefault{`Emptyset_iterator`}
    \cgalParamNEnd

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

    \cgalParamNBegin{do_not_use_cubic_algorithm}
      \cgalParamDescription{Set this parameter to `true` if you only want to use the Delaunay based versions of the algorithm,
                            skipping the cubic search space one in case of failure.}
      \cgalParamType{Boolean}
      \cgalParamDefault{`false`}
      \cgalParamExtra{If `true`, `use_2d_constrained_delaunay_triangulation` or `use_delaunay_triangulation` must be set to `true`
                      otherwise nothing will be done.}
    \cgalParamNEnd

    \cgalParamNBegin{density_control_factor}
      \cgalParamDescription{factor to control density of the output mesh,
                            where larger values cause denser refinements, as in `refine()`}
      \cgalParamType{double}
      \cgalParamDefault{\f$ \sqrt{2}\f$}
    \cgalParamNEnd

    \cgalParamNBegin{visitor}
      \cgalParamDescription{a visitor used to track when entering a given phase of the algorithm}
      \cgalParamType{A model of PMPHolefillingVisitor}
      \cgalParamType{Hole_filling::Default_visitor}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  @return pair of face and vertex output iterator

  \sa CGAL::Polygon_mesh_processing::triangulate_hole()
  \sa CGAL::Polygon_mesh_processing::refine()

  \todo handle islands
  */
  template<typename PolygonMesh,
           typename CGAL_NP_TEMPLATE_PARAMETERS>
  auto
  triangulate_and_refine_hole(PolygonMesh& pmesh,
      typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
      const CGAL_NP_CLASS& np = parameters::default_values())
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;
    using parameters::get_parameter_reference;

    typedef typename internal_np::Lookup_named_param_def<internal_np::face_output_iterator_t,
                                                         CGAL_NP_CLASS,
                                                         Emptyset_iterator>::type Face_output_iterator;

    Face_output_iterator face_out = choose_parameter<Emptyset_iterator>(get_parameter(np, internal_np::face_output_iterator));

    typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_output_iterator_t,
                                                         CGAL_NP_CLASS,
                                                         Emptyset_iterator>::type Vertex_output_iterator;

    Vertex_output_iterator vertex_out = choose_parameter<Emptyset_iterator>(get_parameter(np, internal_np::vertex_output_iterator));

    std::vector<typename boost::graph_traits<PolygonMesh>::face_descriptor> patch;
    triangulate_hole(pmesh, border_halfedge, np.face_output_iterator(std::back_inserter(patch)));
    face_out = std::copy(patch.begin(), patch.end(), face_out);

    Hole_filling::Default_visitor default_visitor;
    typedef typename internal_np::Lookup_named_param_def<internal_np::visitor_t,
                                                         CGAL_NP_CLASS,
                                                         Hole_filling::Default_visitor>::reference Visitor;

    Visitor visitor = choose_parameter(get_parameter_reference(np, internal_np::visitor), default_visitor);
    visitor.start_refine_phase();
    std::pair<Face_output_iterator, Vertex_output_iterator> res = refine(pmesh, patch, face_out, vertex_out, np);
    visitor.end_refine_phase();
    return res;
  }


#ifndef CGAL_NO_DEPRECATED_CODE
 /*!
  \ingroup PMP_hole_filling_grp

  \deprecated This function is deprecated since \cgal 5.6 and the
  overload with the named parameters `face_output_iterator` and
  `vertex_output_iterator` should be used instead.

 @brief triangulates and refines a hole in a polygon mesh.

  @tparam PolygonMesh must be model of `MutableFaceGraph`
  @tparam FaceOutputIterator model of `OutputIterator`
     holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces.
  @tparam VertexOutputIterator model of `OutputIterator`
     holding `boost::graph_traits<PolygonMesh>::%vertex_descriptor` for patch vertices.
  @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 */

  template<typename PolygonMesh,
           typename FaceOutputIterator,
           typename VertexOutputIterator,
           typename CGAL_NP_TEMPLATE_PARAMETERS>
  CGAL_DEPRECATED
  std::pair<FaceOutputIterator, VertexOutputIterator>
    triangulate_and_refine_hole(PolygonMesh& pmesh,
      typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
      FaceOutputIterator face_out,
      VertexOutputIterator vertex_out,
      const CGAL_NP_CLASS& np = parameters::default_values())
  {
    return triangulate_and_refine_hole(pmesh, border_halfedge,
                                       np.face_output_iterator(face_out).vertex_output_iterator(vertex_out));
  }
#endif // CGAL_NO_DEPRECATED_CODE

  /*!
  \ingroup PMP_hole_filling_grp
  @brief triangulates, refines and fairs a hole in a polygon mesh.

  @tparam PolygonMesh a model of `MutableFaceGraph`
  @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  @param pmesh polygon mesh which has the hole
  @param border_halfedge a border halfedge incident to the hole

  @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin

    \cgalParamNBegin{face_output_iterator}
      \cgalParamDescription{iterator over patch faces}
      \cgalParamType{a model of `OutputIterator`
    holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces}
      \cgalParamDefault{`Emptyset_iterator`}
    \cgalParamNEnd

    \cgalParamNBegin{vertex_output_iterator}
      \cgalParamDescription{iterator over patch vertices}
      \cgalParamType{a model of `OutputIterator`
    holding `boost::graph_traits<PolygonMesh>::%vertex_descriptor` for patch vertices}
      \cgalParamDefault{`Emptyset_iterator`}
    \cgalParamNEnd

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
      \cgalParamDescription{factor to control density of the output mesh,
                            where larger values cause denser refinements, as in `refine()`}
      \cgalParamType{double}
      \cgalParamDefault{\f$ \sqrt{2}\f$}
    \cgalParamNEnd

    \cgalParamNBegin{fairing_continuity}
      \cgalParamDescription{A value controlling the tangential continuity of the output surface patch.
                            The possible values are 0, 1 and 2, referring to the  C<sup>0</sup>, C<sup>1</sup>
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

    \cgalParamNBegin{visitor}
      \cgalParamDescription{a visitor used to track when entering a given phase of the algorithm}
      \cgalParamType{A model of PMPHolefillingVisitor}
      \cgalParamType{Hole_filling::Default_visitor}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  @return tuple of `bool` with `true` if fairing is successful, and
  the face and vertex output iterator

  \sa CGAL::Polygon_mesh_processing::triangulate_hole()
  \sa CGAL::Polygon_mesh_processing::refine()
  \sa CGAL::Polygon_mesh_processing::fair()

  \todo handle islands
  */
  template<typename PolygonMesh,
           typename CGAL_NP_TEMPLATE_PARAMETERS>
  auto
  triangulate_refine_and_fair_hole(PolygonMesh& pmesh,
    typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
    const CGAL_NP_CLASS& np = parameters::default_values())
  {
    CGAL_precondition(CGAL::is_triangle_mesh(pmesh));

    using parameters::choose_parameter;
    using parameters::get_parameter;
    using parameters::get_parameter_reference;

    CGAL_precondition(is_valid_halfedge_descriptor(border_halfedge, pmesh));

    typedef typename internal_np::Lookup_named_param_def<internal_np::face_output_iterator_t,
                                                         CGAL_NP_CLASS,
                                                         Emptyset_iterator>::type Face_output_iterator;

    Face_output_iterator face_out = choose_parameter<Emptyset_iterator>(get_parameter(np, internal_np::face_output_iterator));

    typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_output_iterator_t,
                                                         CGAL_NP_CLASS,
                                                         Emptyset_iterator>::type Vertex_output_iterator;

    Vertex_output_iterator vertex_out = choose_parameter<Emptyset_iterator>(get_parameter(np, internal_np::vertex_output_iterator));

    std::vector<typename boost::graph_traits<PolygonMesh>::vertex_descriptor> patch;
    face_out = triangulate_and_refine_hole
      (pmesh, border_halfedge, np.face_output_iterator(face_out).vertex_output_iterator(std::back_inserter(patch))).first;

    CGAL_postcondition(CGAL::is_triangle_mesh(pmesh));

    Hole_filling::Default_visitor default_visitor;
    typedef typename internal_np::Lookup_named_param_def<internal_np::visitor_t,
                                                         CGAL_NP_CLASS,
                                                         Hole_filling::Default_visitor>::reference Visitor;

    Visitor visitor = choose_parameter(get_parameter_reference(np, internal_np::visitor), default_visitor);
    visitor.start_fair_phase();
    bool fair_success = fair(pmesh, patch, np);
    visitor.end_fair_phase();

    vertex_out = std::copy(patch.begin(), patch.end(), vertex_out);
    return std::make_tuple(fair_success, face_out, vertex_out);
  }

  #ifndef CGAL_NO_DEPRECATED_CODE
  /*!
  \ingroup PMP_hole_filling_grp

  \deprecated This function is deprecated since \cgal 5.6 and the
  overload with the named parameters `face_output_iterator` and
  `vertex_output_iterator` should be used instead.

  \brief triangulates, refines, and fairs a hole in a polygon mesh.

  @tparam PolygonMesh a model of `MutableFaceGraph`
  @tparam FaceOutputIterator model of `OutputIterator`
     holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces.
  @tparam VertexOutputIterator model of `OutputIterator`
     holding `boost::graph_traits<PolygonMesh>::%vertex_descriptor` for patch vertices.
  @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  */
  template<typename PolygonMesh,
           typename FaceOutputIterator,
           typename VertexOutputIterator,
           typename CGAL_NP_TEMPLATE_PARAMETERS>
  CGAL_DEPRECATED
  std::tuple<bool, FaceOutputIterator, VertexOutputIterator>
  triangulate_refine_and_fair_hole(PolygonMesh& pmesh,
    typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
    FaceOutputIterator face_out,
    VertexOutputIterator vertex_out,
    const CGAL_NP_CLASS& np = parameters::default_values())
  {
    return triangulate_refine_and_fair_hole(pmesh, border_halfedge, np.face_output_iterator(face_out).vertex_output_iterator(vertex_out));
  }
#endif // CGAL_NO_DEPRECATED_CODE

  /*!
  \ingroup PMP_hole_filling_grp
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

    \cgalParamNBegin{visitor}
      \cgalParamDescription{a visitor used to track when entering a given phase of the algorithm}
      \cgalParamType{A model of PMPHolefillingVisitor}
      \cgalParamType{Hole_filling::Default_visitor}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \todo handle islands
  */
  template <typename PointRange1,
            typename PointRange2,
            typename OutputIterator,
            typename NamedParameters = parameters::Default_named_parameters>
  OutputIterator
  triangulate_hole_polyline(const PointRange1& points,
                            const PointRange2& third_points,
                            OutputIterator out,
                            const NamedParameters& np = parameters::default_values())
  {
    if (points.empty()) return out;

    using parameters::choose_parameter;
    using parameters::get_parameter;
    using parameters::get_parameter_reference;

#ifndef CGAL_HOLE_FILLING_DO_NOT_USE_CDT2
    bool use_cdt = choose_parameter(get_parameter(np, internal_np::use_2d_constrained_delaunay_triangulation), false);
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

    Hole_filling::Default_visitor default_visitor;

#ifndef CGAL_HOLE_FILLING_DO_NOT_USE_CDT2
    if (use_cdt)
    {
      struct Always_valid
      {
        bool operator()(const std::vector<Point>&, int,int,int) const { return true; }
      };
      Always_valid is_valid;

      const typename Kernel::Iso_cuboid_3 bbox = CGAL::bounding_box(points.begin(), points.end());
      typename Kernel::FT default_squared_distance = CGAL::abs(CGAL::squared_distance(bbox.vertex(0), bbox.vertex(5)));
      default_squared_distance /= typename Kernel::FT(16); // one quarter of the bbox height

      const typename Kernel::FT threshold_distance = choose_parameter(
        get_parameter(np, internal_np::threshold_distance), typename Kernel::FT(-1));
      typename Kernel::FT max_squared_distance = default_squared_distance;
      if(threshold_distance >= typename Kernel::FT(0))
        max_squared_distance = threshold_distance * threshold_distance;

      CGAL_assertion(max_squared_distance >= typename Kernel::FT(0));
      if (triangulate_hole_polyline_with_cdt(
           points,
           tracer,
           choose_parameter(get_parameter_reference(np, internal_np::visitor), default_visitor),
           is_valid,
           choose_parameter<Kernel>(get_parameter(np, internal_np::geom_traits)),
           max_squared_distance))
      {
        CGAL_assertion(holes.empty());
        return tracer.out;
      }
    }
#endif
    triangulate_hole_polyline(points, third_points, tracer, WC(),
                              choose_parameter(get_parameter_reference(np, internal_np::visitor), default_visitor),
                              use_dt3,
                              choose_parameter(get_parameter(np, internal_np::do_not_use_cubic_algorithm), false),
                              choose_parameter<Kernel>(get_parameter(np, internal_np::geom_traits)));

    CGAL_assertion(holes.empty());
    return tracer.out;
  }

  /*!
  \ingroup PMP_hole_filling_grp
  Same as above but the range of third points is omitted. They are not
  taken into account in the cost computation that leads the hole filling.
*/
  template <typename PointRange,
            typename OutputIterator,
            typename CGAL_NP_TEMPLATE_PARAMETERS>
  OutputIterator
  triangulate_hole_polyline(const PointRange& points,
                            OutputIterator out,
                            const CGAL_NP_CLASS& np = parameters::default_values())
  {
    typedef typename std::iterator_traits<
      typename PointRange::iterator>::value_type Point;
    std::vector< Point > third_points;
    return triangulate_hole_polyline(points, third_points, out, np);
  }

} //end namespace Polygon_mesh_processing

} //end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif //CGAL_POLYGON_MESH_PROCESSING_TRIANGULATE_HOLE_H
