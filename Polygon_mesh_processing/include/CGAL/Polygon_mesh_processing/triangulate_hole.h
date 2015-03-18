#ifndef CGAL_TRIANGULATE_HOLE_H
#define CGAL_TRIANGULATE_HOLE_H

#include <CGAL/Polygon_mesh_processing/internal/Hole_filling/Triangulate_hole_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/internal/Hole_filling/Triangulate_hole_polyline.h>
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

  @tparam PolygonMesh a model of `MutableFaceGraph`
  @tparam OutputIterator a model of `OutputIterator`
    holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces.

  @param pmesh polygon mesh containing the hole
  @param border_halfedge a border halfedge incident to the hole
  @param out iterator over patch faces
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
    return refine(pmesh, patch, face_out, vertex_out, density_control_factor);
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
  - @a face_out
  - @a vertex_out

  \todo handle islands
  */
  template<typename PolygonMesh,
           typename SparseLinearSolver,
           typename FaceOutputIterator,
           typename VertexOutputIterator>
  CGAL::cpp11::tuple<bool, FaceOutputIterator, VertexOutputIterator>
  triangulate_refine_and_fair_hole(PolygonMesh& pmesh,
    typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
    FaceOutputIterator face_out,
    VertexOutputIterator vertex_out,
    SparseLinearSolver solver = CGAL::Default(),
    double density_control_factor = std::sqrt(2.0),
    unsigned int continuity = 1,
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

    bool fair_success = fair(pmesh, patch, solver, continuity);

    vertex_out = std::copy(patch.begin(), patch.end(), vertex_out);
    return CGAL::cpp11::make_tuple(fair_success, face_out, vertex_out);
  }

  // use non-default weight calculator
  // WeightCalculator a model of `FairWeightCalculator`
  // weight_calculator function object to calculate weights, default to Cotangent weights and can be omitted
  template<class WeightCalculator,
           class SparseLinearSolver,
           class PolygonMesh,
           class FaceOutputIterator,
           class VertexOutputIterator>
  CGAL::cpp11::tuple<bool, FaceOutputIterator, VertexOutputIterator>
    triangulate_refine_and_fair_hole(PolygonMesh& pmesh,
      typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
      FaceOutputIterator face_out,
      VertexOutputIterator vertex_out,
      WeightCalculator weight_calculator,
      SparseLinearSolver solver = CGAL::Default(),
      double density_control_factor = std::sqrt(2.0),
      unsigned int continuity = 1,
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

    bool fair_success = internal::fair(pmesh, patch, solver, weight_calculator, continuity);

    vertex_out = std::copy(patch.begin(), patch.end(), vertex_out);
    return CGAL::cpp11::make_tuple(fair_success, face_out, vertex_out);
  }

  //use default SparseLinearSolver and WeightCalculator
  template<class PolygonMesh,
           class FaceOutputIterator,
           class VertexOutputIterator>
  CGAL::cpp11::tuple<bool, FaceOutputIterator, VertexOutputIterator>
    triangulate_refine_and_fair_hole(PolygonMesh& pmesh,
      typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
      FaceOutputIterator face_out,
      VertexOutputIterator vertex_out,
      double density_control_factor = std::sqrt(2.0),
      unsigned int continuity = 1,
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
  Creates triangles to fill the hole defined by points in the range (@a points).
  Triangles are put into @a out
  using the indices of the input points in the range (@a points).
  Note that no degenerate triangles are allowed during filling. If no possible patch is found, then no triangle is put into @a out.

  The optional range (@a third_points) indicates for each pair of consecutive points in the range (@a points),
  the third point of the facet this segment is incident to.

  Note that the ranges (@a points) and (@a third_points) may or may not contain duplicated first point at the end of sequence.

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

#endif //CGAL_TRIANGULATE_HOLE_H
