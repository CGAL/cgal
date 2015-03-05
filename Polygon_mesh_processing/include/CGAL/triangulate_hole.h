#ifndef CGAL_TRIANGULATE_HOLE_H
#define CGAL_TRIANGULATE_HOLE_H

#include <CGAL/internal/Hole_filling/Triangulate_hole_polygon_mesh.h>
#include <CGAL/internal/Hole_filling/Triangulate_hole_polyline.h>
#include <CGAL/Polygon_mesh_processing/refine.h>
#include <CGAL/Polygon_mesh_processing/fair.h>
#include <CGAL/Default.h>

#include <boost/tuple/tuple.hpp>

namespace CGAL {

namespace Polygon_mesh_processing {

  /*!
  \ingroup PkgPolygonMeshProcessing
  Function triangulating a hole in a polygon mesh.
  The hole should contain no non-manifold vertex. Generated patch is guaranteed to not break edge manifoldness and contain no degenerate triangle.
  If no possible patch is found, @a pmesh is not altered in any way, and no face descriptor is put into @a out.

  @tparam PolygonMesh must be a model of `MutableFaceGraph`
  @tparam OutputIterator iterator holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces.

  @param pmesh polygon mesh containing the hole
  @param border_halfedge a border halfedge incident to the hole
  @param out iterator over patch faces.
  @param use_delaunay_triangulation if `true`, use the Delaunay triangulation facet search space

  @return @a out

  \todo handle islands
  @todo Replace border_halfedge by a range of border halfedges.
        The first one would describe the hole,
        the other ones would describe the islands.
  @todo Then, insert the holes vertices in the set of possibilities
        for connecting vertices together
  @todo Handle the case where an island is reduced to a point
  */
  template<class PolygonMesh, class OutputIterator>
  OutputIterator
    triangulate_hole(PolygonMesh& pmesh,
      typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
      OutputIterator out,
      bool use_delaunay_triangulation = true)
  {
    bool use_dt3 =
#ifdef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
      false;
#else
      use_delaunay_triangulation;
#endif
    CGAL_precondition(face(border_halfedge, pmesh) == boost::graph_traits<PolygonMesh>::null_face());
    return internal::triangulate_hole_polygon_mesh
      (pmesh, border_halfedge, out, use_dt3).first;
  }

  /*!
  \ingroup PkgPolygonMeshProcessing
  @brief Function triangulating and refining a hole in a polygon mesh.

  @tparam PolygonMesh must be model of `MutableFaceGraph`
  @tparam FacetOutputIterator iterator holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces.
  @tparam VertexOutputIterator iterator holding `boost::graph_traits<PolygonMesh>::%vertex_descriptor` for patch vertices.

  @param pmesh polygon mesh which has the hole
  @param border_halfedge a border halfedge incident to the hole
  @param face_out iterator over patch facess
  @param vertex_out iterator over patch vertices without including the boundary
  @param density_control_factor factor for density where larger values cause denser refinements
  @param use_delaunay_triangulation if `true`, use the Delaunay triangulation face search space

  @return pair of @a face_out and @a vertex_out

  \todo handle islands
  */
  template<class PolygonMesh,
           class FaceOutputIterator,
           class VertexOutputIterator>
  std::pair<FaceOutputIterator, VertexOutputIterator>
    triangulate_and_refine_hole(PolygonMesh& pmesh,
      typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
      FaceOutputIterator face_out,
      VertexOutputIterator vertex_out,
      double density_control_factor = std::sqrt(2.0),
      bool use_delaunay_triangulation = true)
  {
    bool use_dt3 =
#ifdef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
      false;
#else
      use_delaunay_triangulation;
#endif
    std::vector<typename boost::graph_traits<PolygonMesh>::face_descriptor> patch;
    triangulate_hole(pmesh, border_halfedge, std::back_inserter(patch), use_dt3);
    face_out = std::copy(patch.begin(), patch.end(), face_out);
    return refine(pmesh, patch.begin(), patch.end(), face_out, vertex_out, density_control_factor);
  }


  /*!
  \ingroup PkgPolygonMeshProcessing
  @brief Function triangulating, refining and fairing a hole in a polygon mesh.

  If \ref thirdpartyEigen "Eigen" 3.2 (or greater) is available
  and `CGAL_EIGEN3_ENABLED` is defined, an overload of this function is available
  with `SparseLinearSolver` being:
  \code
  CGAL::Eigen_solver_traits<
  Eigen::SparseLU<
  CGAL::Eigen_sparse_matrix<double>::EigenType,
  Eigen::COLAMDOrdering<int> >  >
  \endcode
  .
  For using an alternative model of `FairWeightCalculator` with the default solver,
  one can pass `CGAL::Default()` as `solver`.

  @tparam SparseLinearSolver a model of `SparseLinearAlgebraTraitsWithPreFactor_d`
  @tparam PolygonMesh must be model of  `MutableFaceGraph`
  @tparam FaceOutputIterator iterator holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces.
  @tparam VertexOutputIterator iterator holding `boost::graph_traits<PolygonMesh>::%vertex_descriptor` for patch vertices.

  @param pmesh a polygon mesh which has the hole
  @param border_halfedge a border halfedge incident to the hole
  @param face_out iterator over patch faces
  @param vertex_out iterator over patch vertices without including the boundary
  @param density_control_factor factor for density where larger values cause denser refinements
  @param continuity tangential continuity, default to `FAIRING_C_1` and can be omitted
  @param use_delaunay_triangulation if `true`, use the Delaunay triangulation face search space
  @param solver An instance of the sparse linear solver to use. Note that the current implementation is
  not using the value passed but the default constructed one.

  @return tuple of
  - bool: `true` if fairing is successful
  - @a face_out
  - @a vertex_out

  \todo handle islands
  */
  template<typename SparseLinearSolver,
           typename PolygonMesh,
           typename FaceOutputIterator,
           typename VertexOutputIterator>
  boost::tuple<bool, FaceOutputIterator, VertexOutputIterator>
    triangulate_refine_and_fair_hole(PolygonMesh& pmesh,
    typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
    FaceOutputIterator face_out,
    VertexOutputIterator vertex_out,
    SparseLinearSolver
#ifdef DOXYGEN_RUNNING
    solver,
#else
    /* solver */,
#endif
    double density_control_factor = std::sqrt(2.0),
    Fairing_continuity continuity = FAIRING_C_1,
    bool use_delaunay_triangulation = true)
  {
    bool use_dt3 =
#ifdef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
      false;
#else
      use_delaunay_triangulation;
#endif
    std::vector<typename boost::graph_traits<PolygonMesh>::vertex_descriptor> patch;

    face_out = triangulate_and_refine_hole
      (pmesh, border_halfedge, face_out, std::back_inserter(patch),
      density_control_factor, use_dt3).first;

    typedef internal::Fair_default_sparse_linear_solver::Solver Default_solver;
    typedef typename Default::Get<SparseLinearSolver, Default_solver>::type Solver;

    bool fair_success = fair<Solver>(pmesh, patch.begin(), patch.end(), continuity);

    vertex_out = std::copy(patch.begin(), patch.end(), vertex_out);
    return boost::make_tuple(fair_success, face_out, vertex_out);
  }

  // use non-default weight calculator
  // WeightCalculator a model of `FairWeightCalculator`
  // weight_calculator function object to calculate weights, default to Cotangent weights and can be omitted
  template<class WeightCalculator,
           class SparseLinearSolver,
           class PolygonMesh,
           class FaceOutputIterator,
           class VertexOutputIterator>
  boost::tuple<bool, FaceOutputIterator, VertexOutputIterator>
    triangulate_refine_and_fair_hole(PolygonMesh& pmesh,
      typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
      FaceOutputIterator face_out,
      VertexOutputIterator vertex_out,
      WeightCalculator weight_calculator,
      SparseLinearSolver
      #ifdef DOXYGEN_RUNNING
        solver,
      #else
        /* solver */,
      #endif
      double density_control_factor = std::sqrt(2.0),
      Fairing_continuity continuity = FAIRING_C_1,
      bool use_delaunay_triangulation = true)
  {
    bool use_dt3 =
#ifdef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
      false;
#else
      use_delaunay_triangulation;
#endif
    std::vector<typename boost::graph_traits<PolygonMesh>::vertex_descriptor> patch;

    face_out = triangulate_and_refine_hole
      (pmesh, border_halfedge, face_out, std::back_inserter(patch),
      density_control_factor, use_dt3).first;

    typedef internal::Fair_default_sparse_linear_solver::Solver Default_solver;
    typedef typename Default::Get<SparseLinearSolver, Default_solver>::type Solver;

    bool fair_success = fair<Solver>(pmesh, patch.begin(), patch.end(), weight_calculator, continuity);

    vertex_out = std::copy(patch.begin(), patch.end(), vertex_out);
    return boost::make_tuple(fair_success, face_out, vertex_out);
  }

  //use default SparseLinearSolver and WeightCalculator
  template<class PolygonMesh,
           class FaceOutputIterator,
           class VertexOutputIterator>
  boost::tuple<bool, FaceOutputIterator, VertexOutputIterator>
    triangulate_refine_and_fair_hole(PolygonMesh& pmesh,
      typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
      FaceOutputIterator face_out,
      VertexOutputIterator vertex_out,
      double density_control_factor = std::sqrt(2.0),
      Fairing_continuity continuity = FAIRING_C_1,
      bool use_delaunay_triangulation = true)
  {
    bool use_dt3 =
#ifdef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
      false;
#else
      use_delaunay_triangulation;
#endif

    return triangulate_refine_and_fair_hole
      (pmesh, border_halfedge, face_out, vertex_out, Default(),
        density_control_factor, continuity, use_dt3);
  }

  /*!
  \ingroup PkgPolygonMeshProcessing
  Creates triangles to fill the hole defined by points in the range (@a pbegin, @a pend). Triangles are put into @a out
  using the indices of the input points in the range (@a pbegin, @a pend).
  Note that no degenerate triangle is allowed during filling. If no possible patch is found, then no triangle is put into @a out.

  The optional range (@a qbegin, @a qend) indicate for each pair of consecutive points in the range (@a pbegin, @a pend),
  the third point of the facet this segment is incident to.

  Note that the range (@a pbegin, @a pend) and (@a qbegin, @a qend) may or may not contain duplicated first point at the end of sequence.

  @tparam OutputIteratorValueType value type of OutputIterator having a constructor `OutputIteratorValueType(int p0, int p1, int p2)` available.
  It is default to value_type_traits<OutputIterator>::type, and can be omitted when the default is fine
  @tparam InputIterator iterator over input points
  @tparam OutputIterator iterator over patch triangles
  @param pbegin first iterator of the range of points
  @param pend past-the-end iterator of the range of points
  @param qbegin first iterator of the range of third points, can be omitted
  @param qend past-the-end iterator of the range of third points, can be omitted
  @param out iterator over output patch triangles
  @param use_delaunay_triangulation if `true`, use the Delaunay triangulation facet search space

  \todo handle islands
  */
  template <typename OutputIteratorValueType,
            typename InputIterator,
            typename OutputIterator>
  OutputIterator
    triangulate_hole_polyline(InputIterator pbegin, InputIterator pend,
                              InputIterator qbegin, InputIterator qend,
                              OutputIterator out,
                              bool use_delaunay_triangulation = true)
  {
    bool use_dt3 =
#ifdef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
      false;
#else
      use_delaunay_triangulation;
#endif
    typedef typename std::iterator_traits<InputIterator>::value_type Point_3;
    typedef CGAL::internal::Weight_min_max_dihedral_and_area      Weight;
    typedef CGAL::internal::Weight_calculator<Weight,
                  CGAL::internal::Is_not_degenerate_triangle>  WC;
    typedef std::vector<std::pair<int, int> > Holes;
    typedef std::back_insert_iterator<Holes>  Holes_out;

    std::vector<Point_3> P(pbegin, pend);
    std::vector<Point_3> Q(qbegin, qend);
    Holes holes;//just to check there is no holes

    CGAL::internal::Tracer_polyline_incomplete<OutputIteratorValueType, OutputIterator, Holes_out>
      tracer(out, Holes_out(holes));
    triangulate_hole_polyline(P, Q, tracer, WC(), use_dt3);
    CGAL_assertion(holes.empty());
    return tracer.out;
  }

  // overload for OutputIteratorValueType
  template <typename InputIterator,
            typename OutputIterator>
  OutputIterator
    triangulate_hole_polyline(InputIterator pbegin, InputIterator pend,
                              InputIterator qbegin, InputIterator qend,
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
      (pbegin, pend, qbegin, qend, out, use_dt3);
  }

  // overload no (qbegin, qend)
  template <typename OutputIteratorValueType,
            typename InputIterator,
            typename OutputIterator>
  OutputIterator
    triangulate_hole_polyline(InputIterator pbegin, InputIterator pend,
                              OutputIterator out,
                              bool use_delaunay_triangulation = true)
  {
    bool use_dt3 =
#ifdef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
      false;
#else
      use_delaunay_triangulation;
#endif
    typedef typename std::iterator_traits<InputIterator>::value_type Point;
    typedef std::vector<Point> Polyline_3;
    Polyline_3 Q;
    return triangulate_hole_polyline<OutputIteratorValueType>
      (pbegin, pend, Q.begin(), Q.end(), out, use_dt3);
  }

  // overload for OutputIteratorValueType
  template <typename InputIterator,
            typename OutputIterator>
  OutputIterator
    triangulate_hole_polyline(InputIterator pbegin, InputIterator pend,
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
      (pbegin, pend, out, use_dt3);
  }

} //end namespace Polygon_mesh_processing

} //end namespace CGAL

#endif //CGAL_TRIANGULATE_HOLE_H
