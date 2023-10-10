// Copyright (c) 2023 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Oesau, Florent Lafarge, Dmitry Anisimov, Simon Giraudot

#ifndef CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H
#define CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H

#include <CGAL/Kinetic_shape_partition_3.h>
#include <CGAL/KSR_3/GraphCut.h>

#include <CGAL/IO/PLY.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Point_set.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/KSR/debug.h>
#include <CGAL/Shape_regularization/regularize_planes.h>

#include <boost/filesystem.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/helpers.h>

namespace CGAL
{
#ifndef DOXYGEN_RUNNING
/*!
* \ingroup PkgKineticShapePartition
  \brief Piece-wise linear reconstruction via inside/outside labeling of a kinetic partition using graph cut.

  \tparam GeomTraits
    must be a model of `KineticPartitionTraits`.

  \tparam NormalMap
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`. It must map the elements in `KineticShapePartitionTraits_3::Input_range` to `Vector_3`.
*/
template<typename GeomTraits, typename PointSet, typename PointMap, typename NormalMap, typename IntersectionKernel = CGAL::Exact_predicates_exact_constructions_kernel>
class Kinetic_shape_reconstruction_3 {
public:
  using Kernel = typename GeomTraits;
  using Intersection_kernel = typename IntersectionKernel;

  using FT = typename Kernel::FT;

  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;
  using Vector_3 = typename Kernel::Vector_3;
  using Plane_3 = typename Kernel::Plane_3;
  using Triangle_2 = typename Kernel::Triangle_2;

  using Point_set = PointSet;

  using Indices = std::vector<std::size_t>;
  using Polygon_3 = std::vector<Point_3>;

  using KSP = Kinetic_shape_partition_3<Kernel, Intersection_kernel>;

  using Point_map = typename PointMap;
  using Normal_map = typename NormalMap;

  using Region_type = CGAL::Shape_detection::Point_set::Least_squares_plane_fit_region_for_point_set<Point_set>;
  using Neighbor_query = CGAL::Shape_detection::Point_set::K_neighbor_query_for_point_set<Point_set>;
  using Sorting = CGAL::Shape_detection::Point_set::Least_squares_plane_fit_sorting_for_point_set<Point_set, Neighbor_query>;
  using Region_growing = CGAL::Shape_detection::Region_growing<Neighbor_query, Region_type>;
  using From_exact = typename CGAL::Cartesian_converter<Intersection_kernel, Kernel>;

  /*!
  \brief Creates a Kinetic_shape_reconstruction_3 object.

  \param ps
  an instance of `InputRange` with 3D points and corresponding 3D normal vectors.

  */
  template<typename NamedParameters = parameters::Default_named_parameters>
  Kinetic_shape_reconstruction_3(PointSet& ps,
    const NamedParameters& np = CGAL::parameters::default_values()) : m_points(ps), m_kinetic_partition(np), m_point_map(ps.point_map()), m_normal_map(ps.normal_map()) {}

  /*!
    \brief Detects shapes in the provided point cloud

    \tparam NamedParameters
    a sequence of \ref bgl_namedparameters "Named Parameters"

    \param np
    an instance of `NamedParameters`.

  \cgalNamedParamsBegin
    \cgalParamNBegin{point_map}
      \cgalParamDescription{a property map associating points to the elements of `input_range`}
      \cgalParamDefault{`PointMap()`}
    \cgalParamNEnd
    \cgalParamNBegin{normal_map}
      \cgalParamDescription{a property map associating normals to the elements of `input_range`}
      \cgalParamDefault{`NormalMap()`}
    \cgalParamNEnd
    \cgalParamNBegin{k_neighbors}
      \cgalParamDescription{the number of returned neighbors per each query point}
      \cgalParamType{`std::size_t`}
      \cgalParamDefault{12}
    \cgalParamNEnd
    \cgalParamNBegin{distance_tolerance}
      \cgalParamDescription{maximum distance from a point to a planar shape}
      \cgalParamType{`GeomTraits::FT`}
      \cgalParamDefault{1}
    \cgalParamNEnd
    \cgalParamNBegin{angle_tolerance}
      \cgalParamDescription{maximum angle in degrees between the normal of a point and the plane normal}
      \cgalParamType{`GeomTraits::FT`}
      \cgalParamDefault{25 degrees}
    \cgalParamNEnd
    \cgalParamNBegin{minimum_region_size}
      \cgalParamDescription{minimum number of 3D points a region must have}
      \cgalParamType{`std::size_t`}
      \cgalParamDefault{5}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  */
  template<
    typename CGAL_NP_TEMPLATE_PARAMETERS>
  std::size_t detect_planar_shapes(
    const CGAL_NP_CLASS& np = parameters::default_values()) {

    if (m_verbose)
      std::cout << std::endl << "--- DETECTING PLANAR SHAPES: " << std::endl;

    m_planes.clear();
    m_polygons.clear();
    m_region_map.clear();

    m_point_map = Point_set_processing_3_np_helper<Point_set, CGAL_NP_CLASS, Point_map>::get_point_map(m_points, np);
    m_normal_map = Point_set_processing_3_np_helper<Point_set, CGAL_NP_CLASS, Normal_map>::get_normal_map(m_points, np);

    create_planar_shapes(np);

    CGAL_assertion(m_planes.size() == m_polygons.size());
    CGAL_assertion(m_polygons.size() == m_region_map.size());

    return m_polygons.size();
  }

  /*!
  \brief Retrieves the detected shapes.

  @return
  vector with a plane equation for each detected planar shape.

  \pre `successful shape detection`
  */
  const std::vector<Plane_3>& detected_planar_shapes() {
    return m_planes;
  }

  /*!
  \brief Retrieves the indices of detected shapes.

    @return
    indices into `input_range` for each detected planar shape in vectors.

    \pre `successful shape detection`
  */
  const std::vector<std::vector<std::size_t> >& detected_planar_shape_indices() {
    return m_planar_regions;
  }

  /*!
  \brief initializes the kinetic partition.

  \param np
  a sequence of \ref bgl_namedparameters "Named Parameters"
  among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamNBegin{reorient_bbox}
      \cgalParamDescription{Use the oriented bounding box instead of the axis-aligned bounding box.}
      \cgalParamType{bool}
      \cgalParamDefault{false}
    \cgalParamNEnd
    \cgalParamNBegin{bbox_dilation_ratio}
      \cgalParamDescription{Factor for extension of the bounding box of the input data to be used for the partition.}
      \cgalParamType{FT}
      \cgalParamDefault{1.1}
    \cgalParamNEnd
    \cgalParamNBegin{angle_tolerance}
      \cgalParamDescription{The tolerance angle to snap the planes of two input polygons into one plane.}
      \cgalParamType{FT}
      \cgalParamDefault{5}
    \cgalParamNEnd
    \cgalParamNBegin{distance_tolerance}
      \cgalParamDescription{The tolerance distance to snap the planes of two input polygons into one plane.}
      \cgalParamType{FT}
      \cgalParamDefault{5}
    \cgalParamNEnd
  \cgalNamedParamsEnd

    \pre `successful shape detection`
  */
  template<typename CGAL_NP_TEMPLATE_PARAMETERS>
  void initialize_partition(const CGAL_NP_CLASS& np = parameters::default_values()) {
    m_kinetic_partition.insert(m_polygon_pts, m_polygon_indices, np);

    m_kinetic_partition.initialize(np);
  }

  /*!
  \brief Propagates the kinetic polygons in the initialized partition.

  \param k
  maximum number of allowed intersections for each input polygon before its expansion stops.

  @return
  success of kinetic partition.

  \pre `successful initialization`
  */
  void partition(std::size_t k, FT& partition_time, FT& finalization_time, FT& conformal_time) {
    m_kinetic_partition.partition(k, partition_time, finalization_time, conformal_time);
    std::cout << "Bounding box partitioned into " << m_kinetic_partition.number_of_volumes() << " volumes" << std::endl;
  }

  /*!
  \brief Access to the kinetic partition.

  @return
  created kinetic partition data structure

  \pre `successful partition`
  */
  const Kinetic_shape_partition_3<Kernel, Intersection_kernel>& partition() const {
    return m_kinetic_partition;
  }

  /*!
  \brief Creates the visibility (data-) and regularity energy terms from the input point cloud and the kinetic partition.

  @return
  success.

  \pre `successful initialization`
  */
  bool setup_energyterms() {
    CGAL::Timer timer;
    timer.start();
    if (m_kinetic_partition.number_of_volumes() == 0) {
      if (m_verbose) std::cout << "Kinetic partition is not constructed or does not have volumes" << std::endl;
      return false;
    }

    m_faces.clear();
    m_face2index.clear();

    m_face_area.clear();
    m_face_inliers.clear();
    m_face_neighbors.clear();

    m_kinetic_partition.unique_faces(std::back_inserter(m_faces));
    std::cout << "Found " << m_faces.size() << " faces between volumes" << std::endl;
    timer.reset();

    for (std::size_t i = 0; i < m_faces.size(); i++)
      m_face2index[m_faces[i]] = i;

    // Create value arrays for graph cut
    m_face_inliers.resize(m_faces.size());
    m_face_area.resize(m_faces.size());
    m_face_neighbors.resize(m_faces.size(), std::pair<int, int>(-1, -1));

    m_cost_matrix.resize(2);
    m_cost_matrix[0].resize(m_kinetic_partition.number_of_volumes() + 6);
    m_cost_matrix[1].resize(m_kinetic_partition.number_of_volumes() + 6);

    std::cout << "* computing visibility ... ";

    for (std::size_t i = 0; i < m_faces.size(); i++) {
      // Shift by 6 for accommodate outside volumes
      // Map negative numbers -1..-6 to 0..5
      std::pair<int, int> p = m_kinetic_partition.neighbors(m_faces[i]);
      assert(p.second >= -6);
      if (p.second < 0)
        m_face_neighbors[i] = std::make_pair(p.first + 6, std::size_t(-p.second - 1));
      else
        m_face_neighbors[i] = std::make_pair(p.first + 6, p.second + 6);
    }

    timer.reset();
    collect_points_for_faces3();
    timer.reset();

    std::cout << "* computing data term ... ";

    count_volume_votes();
    timer.reset();

    std::size_t max_inside = 0, max_outside = 0;
    for (std::size_t i = 0; i < m_volumes.size(); i++) {
      max_inside = (std::max<double>)(m_cost_matrix[0][i + 6], max_inside);
      max_outside = (std::max<double>)(m_cost_matrix[1][i + 6], max_outside);
    }

    // Dump volumes colored by votes
    if (false) {
      namespace fs = boost::filesystem;
      for (fs::directory_iterator end_dir_it, it("gc/i"); it != end_dir_it; ++it) {
        fs::remove_all(it->path());
      }
      for (fs::directory_iterator end_dir_it, it("gc/o"); it != end_dir_it; ++it) {
        fs::remove_all(it->path());
      }
      for (fs::directory_iterator end_dir_it, it("gc/n"); it != end_dir_it; ++it) {
        fs::remove_all(it->path());
      }
      for (fs::directory_iterator end_dir_it, it("gc/all"); it != end_dir_it; ++it) {
        fs::remove_all(it->path());
      }
      for (std::size_t i = 0; i < m_volumes.size(); i++) {
        // skip 0/0 volumes? Maybe safe them a few seconds later to be able to separate them?
        CGAL::Color c;

        int diff = int(m_cost_matrix[0][i + 6]) - int(m_cost_matrix[1][i + 6]);

        if (diff > 0) {
          std::size_t m = (std::max<int>)(50, (diff * 255) / max_inside);
          c = CGAL::Color(0, m, 0);
        }
        else {
          std::size_t m = (std::max<int>)(50, (-diff * 255) / max_outside);
          c = CGAL::Color(0, 0, m);
        }

        if (diff < 0) {
          dump_volume(i, "gc/o/" + std::to_string(i) + "-vol-" + std::to_string(m_cost_matrix[0][i + 6]) + "-" + std::to_string(m_cost_matrix[1][i + 6]), c);
          dump_volume(i, "gc/all/" + std::to_string(i) + "-vol-" + std::to_string(m_cost_matrix[0][i + 6]) + "-" + std::to_string(m_cost_matrix[1][i + 6]), c);
        }
        else if (diff > 0) {
          dump_volume(i, "gc/i/" + std::to_string(i) + "-vol-" + std::to_string(m_cost_matrix[0][i + 6]) + "-" + std::to_string(m_cost_matrix[1][i + 6]), c);
          dump_volume(i, "gc/all/" + std::to_string(i) + "-vol-" + std::to_string(m_cost_matrix[0][i + 6]) + "-" + std::to_string(m_cost_matrix[1][i + 6]), c);
        }
        else {
          dump_volume(i, "gc/n/" + std::to_string(i) + "-vol-0-0", CGAL::Color(255, 255, 255));
          dump_volume(i, "gc/all/" + std::to_string(i) + "-vol-0-0", CGAL::Color(255, 255, 255));
        }
      }
    }

    return true;
  }

  /*!
  \brief Provides the data and regularity energy terms for reconstruction via graph-cut.

  \param edges
  contains a vector of pairs of volume indices. Indicates which volumes should be connected in the graph cut formulation.

  \param edge_costs
  contains the cost for each edge specified in `edges` for two labels with different labels. For equal labels, the cost is 0. Needs to be index compatible to the `edges` parameter.

  \param cost_matrix
  provides the cost of a label for each volume cell. The first index corresponds to the label and the second index corresponds to the volume index.

  @return
  fails if the dimensions of parameters does not match the kinetic partition.

  \pre `successful initialization`
  */
  template<typename NamedParameters>
  bool setup_energyterms(
    const std::vector< std::pair<std::size_t, std::size_t> >& edges,
    const std::vector<double>& edge_costs,
    const std::vector< std::vector<double> >& cost_matrix);

  /*!
  \brief Uses graph-cut to solve an solid/empty labeling of the volumes of the kinetic partition.

  \param lambda
  trades the impact of the data term for impact of the regularization term. Should be in the range [0, 1).

  @return
  success of reconstruction.

  \pre `successful initialization`
  */
  bool reconstruct(FT lambda) {
    KSR_3::Graphcut<Kernel> gc(lambda);
    gc.solve(m_face_neighbors, m_face_area, m_cost_matrix, m_labels);

    return true;
  }

  /*!
  \brief Provides the reconstructed surface as a list of polygons.

  \param it
  an output iterator taking std::vector<Point_3>.

  \pre `successful reconstruction`
  */
  template<class OutputIterator>
  void reconstructed_model(OutputIterator it) {
    if (m_labels.empty())
      return;

    for (std::size_t i = 0; i < m_faces.size(); i++) {
      const auto& n = m_face_neighbors[i];
      // Do not extract boundary faces.
      if (n.second < 6)
        continue;

      if (m_labels[n.first] != m_labels[n.second]) {
        std::vector<Point_3> face;
        m_kinetic_partition.vertices(m_faces[i], std::back_inserter(face));
        *it++ = std::move(face);
      }
    }
  }

  /*!
  \brief Provides the reconstructed surface as a list of indexed polygons.

  \param pit
  an output iterator taking Point_3.

  \param triit
  an output iterator taking std::vector<std::size_t>.

  \pre `successful reconstruction`
  */
  template<class OutputPointIterator, class OutputPolygonIterator>
  void reconstructed_model_polylist(OutputPointIterator pit, OutputPolygonIterator polyit) {
    if (m_labels.empty())
      return;

    std::map<Point_3, std::size_t> pt2idx;

    for (std::size_t i = 0; i < m_faces.size(); i++) {
      const auto& n = m_face_neighbors[i];
      // Do not extract boundary faces.
       if (n.second < 6)
         continue;
      if (m_labels[n.first] != m_labels[n.second]) {
        std::vector<Point_3> face;
        m_kinetic_partition.vertices(m_faces[i], std::back_inserter(face));

        std::vector<std::size_t> indices(face.size());

        for (std::size_t i = 0; i < face.size(); i++) {
          auto p = pt2idx.emplace(face[i], pt2idx.size());
          if (p.second)
            *pit++ = face[i];
          indices[i] = p.first->second;
        }

        *polyit++ = std::move(indices);
      }
    }
  }

  /*!
  \brief Provides the reconstructed surface as a list of indexed triangles.

  \param pit
  an output iterator taking Point_3.

  \param triit
  an output iterator taking std::size_t.

  \pre `successful reconstruction`
  */
  template<class OutputPointIterator, class OutputTriangleIterator>
  void reconstructed_model_trilist(OutputPointIterator pit, OutputTriangleIterator triit) {
    if (m_labels.empty())
      return;

    std::map<Point_3, std::size_t> pt2idx;

    for (std::size_t i = 0; i < m_faces.size(); i++) {
      const auto& n = m_face_neighbors[i];
      // Do not extract boundary faces.
      if (n.second < 6)
        continue;
      if (m_labels[n.first] != m_labels[n.second]) {
        std::vector<Point_3> face;
        m_kinetic_partition.vertices(m_faces[i], std::back_inserter(face));

        std::vector<std::size_t> indices(face.size());

        for (std::size_t i = 0; i < face.size(); i++) {
          auto p = pt2idx.emplace(face[i], pt2idx.size());
          if (p.second)
            *pit++ = face[i];
          indices[i] = p.first->second;
        }

        for (std::size_t i = 2; i < face.size(); i++) {
          *triit++ = indices[0];
          *triit++ = indices[i - 1];
          *triit++ = indices[i];
        }
      }
    }
  }

private:

  struct Vertex_info { FT z = FT(0); };
  struct Face_info { };

  using Fbi = CGAL::Triangulation_face_base_with_info_2<Face_info, Kernel>;
  //using Fb = CGAL::Alpha_shape_face_base_2<Kernel, Fbi>;

  using Vbi = CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, Kernel>;
  //using Vb = CGAL::Alpha_shape_vertex_base_2<Kernel, Vbi>;

  using Tds = CGAL::Triangulation_data_structure_2<Vbi, Fbi>;
  using Delaunay_2 = CGAL::Delaunay_triangulation_2<Kernel, Tds>;

  using Delaunay_3 = CGAL::Delaunay_triangulation_3<Kernel>;

  struct VI
  {
    VI()
      : idx(-1)
    {}

    void set_index(std::size_t i) {
      idx = i;
    }
    std::vector<std::size_t> idx;
  };

  struct ID {
    ID()
      : id(-1)
    {}

    std::size_t id;
  };

  typedef CGAL::Triangulation_vertex_base_with_info_2<VI, Intersection_kernel> Vbi2;
  typedef CGAL::Triangulation_face_base_with_info_2<ID, Intersection_kernel> Fbi2;
  typedef CGAL::Constrained_triangulation_face_base_2<Intersection_kernel, Fbi2>    Fb;
  typedef CGAL::Triangulation_data_structure_2<Vbi2, Fb>       Tds2;
  typedef CGAL::Exact_intersections_tag                     Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<Intersection_kernel, Tds2, Itag> CDT;
  typedef CGAL::Constrained_triangulation_plus_2<CDT>       CDTplus;
  typedef typename CDTplus::Vertex_handle            Vertex_handle;
  typedef typename CDTplus::Face_handle              Face_handle;
  typedef typename CDTplus::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename CDTplus::Finite_faces_iterator    Finite_faces_iterator;

  //using Visibility = KSR_3::Visibility<Kernel, Intersection_kernel, Point_map, Normal_map>;
  using Index = typename KSP::Index;

  bool m_verbose;
  bool m_debug;

  Point_set &m_points;
  Point_map m_point_map;
  Normal_map m_normal_map;

  std::vector<std::vector<std::size_t> > m_planar_regions;
  std::vector<typename Region_growing::Primitive_and_region> m_regions;
  std::map<std::size_t, Indices> m_region_map;
  double m_detection_distance_tolerance;

  std::vector<Plane_3> m_planes;
  std::vector<Point_3> m_polygon_pts;
  std::vector<std::vector<std::size_t> > m_polygon_indices;
  std::vector<Polygon_3> m_polygons;
  KSP m_kinetic_partition;

  // Face indices are now of type Indices and are not in a range 0 to n
  std::vector<typename KSP::Index> m_faces;
  std::map<typename KSP::Index, std::size_t> m_face2index;
  std::vector<Indices> m_face_inliers;
  std::vector<FT> m_face_area;
  std::vector<std::pair<std::size_t, std::size_t> > m_face_neighbors;

  std::vector<std::pair<std::size_t, std::size_t> > m_volume_votes; // pair<inside, outside> votes
  std::vector<std::vector<double> > m_cost_matrix;
  std::vector<FT> m_volumes; // normalized volume of each kinetic volume
  std::vector<std::size_t> m_labels;

  std::size_t m_total_inliers;

  std::size_t add_convex_hull_shape(
    const std::vector<std::size_t>& region, const Plane_3& plane) {

    std::vector<Point_2> points;
    points.reserve(region.size());
    for (const std::size_t idx : region) {
      CGAL_assertion(idx < m_points.size());
      const auto& p = get(m_point_map, idx);
      const auto q = plane.projection(p);
      const auto point = plane.to_2d(q);
      points.push_back(point);
    }
    CGAL_assertion(points.size() == region.size());

    std::vector<Point_2> ch;
    CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(ch));

    std::vector<Point_3> polygon;
    for (const auto& p : ch) {
      const auto point = plane.to_3d(p);
      polygon.push_back(point);
    }

    const std::size_t shape_idx = m_polygons.size();
    m_polygons.push_back(polygon);
    m_planes.push_back(plane);

    m_polygon_indices.push_back(std::vector<std::size_t>());
    m_polygon_indices.back().resize(polygon.size());
    std::iota(std::begin(m_polygon_indices.back()), std::end(m_polygon_indices.back()), m_polygon_pts.size());
    std::copy(polygon.begin(), polygon.end(), std::back_inserter(m_polygon_pts));

    return shape_idx;
  }

  void store_convex_hull_shape(const std::string &filename,
    const std::vector<std::size_t>& region, const Plane_3& plane) {

    std::vector<Point_2> points;
    points.reserve(region.size());
    for (const std::size_t idx : region) {
      CGAL_assertion(idx < m_points.size());
      const auto& p = get(m_point_map, idx);
      const auto q = plane.projection(p);
      const auto point = plane.to_2d(q);
      points.push_back(point);
    }
    CGAL_assertion(points.size() == region.size());

    std::vector<Point_2> ch;
    CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(ch));

    std::vector<Point_3> polygon;
    for (const auto& p : ch) {
      const auto point = plane.to_3d(p);
      polygon.push_back(point);
    }

    KSR_3::dump_polygon(polygon, filename);
  }

  std::pair<int, int> make_canonical_pair(int i, int j)
  {
    if (i > j) return std::make_pair(j, i);
    return std::make_pair(i, j);
  }

  double build_cdt(CDTplus& cdt, std::vector<std::vector<std::size_t> >& faces, const std::vector<typename Intersection_kernel::Point_2> &pts, const std::vector<std::vector<std::size_t> > &indices) {
    double area = 0;

    cdt.clear();

    //check orientation of faces so that they are ccw oriented
    std::vector<std::vector<Index> > pts_idx(faces.size());
    //std::vector<std::vector<typename Intersection_kernel::Point_3> > pts(faces.size());
    for (std::size_t i = 0; i < faces.size(); ++i) {

      CGAL::Orientation res = CGAL::COLLINEAR;
      bool pos = false;
      bool neg = false;

      for (std::size_t j = 0; j < faces[i].size(); j++) {
        std::size_t k = (j + 1) % faces[i].size();
        std::size_t l = (k + 1) % faces[i].size();

        res = orientation(pts[faces[i][j]], pts[faces[i][k]], pts[faces[i][l]]);
        if (res == CGAL::LEFT_TURN)
          pos = true;
        if (res == CGAL::RIGHT_TURN)
          neg = true;
      }

      if (pos && neg)
        std::cout << "face is not convex" << std::endl;

      if (!pos && !neg)
        std::cout << "face is degenerated" << std::endl;

      if (neg)
        std::reverse(faces[i].begin(), faces[i].end());
    }

    std::vector<Vertex_handle> vertices;

    for (std::size_t v = 0; v < pts.size(); v++) {
      vertices.push_back(cdt.insert(pts[v]));
      vertices.back()->info().idx = indices[v];
    }

    typedef std::set<std::pair<int, int> > Edges;
    Edges edges;

    for (std::size_t f = 0; f < faces.size(); ++f) {
      for (std::size_t j = 0; j < faces[f].size(); ++j) {
        int vj = faces[f][j];
        int vjj = faces[f][(j + 1) % faces[f].size()];
        std::pair<Edges::iterator, bool> res = edges.insert(make_canonical_pair(vj, vjj));
#ifdef OVERLAY_2_DEBUG
        int vjjj = face2vtx[v[(j + 2) % v.size()]];
        if (orientation(vertices[vj]->point(), vertices[vjj]->point(), vertices[vjjj]->point()) != CGAL::LEFT_TURN) {
          std::cerr << "orientation( " << vertices[vj]->point() << ", " << vertices[vjj]->point() << ", " << vertices[vjjj]->point() << std::endl;
          std::cerr << orientation(vertices[vj]->point(), vertices[vjj]->point(), vertices[vjjj]->point()) << std::endl;
        }
#endif
        if (res.second) {
          cdt.insert_constraint(vertices[vj], vertices[vjj]);
        }
      }
    }

    for (CDTplus::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
      std::set<std::size_t> a, b, c;
      std::copy(fit->vertex(0)->info().idx.begin(), fit->vertex(0)->info().idx.end(), std::inserter(a, a.begin()));
      std::copy(fit->vertex(1)->info().idx.begin(), fit->vertex(0)->info().idx.end(), std::inserter(b, b.begin()));
      std::copy(fit->vertex(2)->info().idx.begin(), fit->vertex(0)->info().idx.end(), std::inserter(c, c.begin()));

      std::set<std::size_t> res, res2;
      Index common(std::size_t(-1), std::size_t(-1));
      std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::inserter(res, res.begin()));
      std::set_intersection(res.begin(), res.end(), c.begin(), c.end(), std::inserter(res2, res2.begin()));

      if (res2.size() != 1) {
        std::cout << "KSR::build_cdt: face assignment not unique!" << std::endl;
      }
      else fit->info().id = *res2.begin();
    }

    return area;
  }
  /*

  void collect_points_for_faces() {
    FT total_area = 0;
    m_total_inliers = 0;
    for (std::size_t i = 0; i < m_polygons.size(); i++) {
      std::vector<KSP::Index> faces;
      m_kinetic_partition.faces_of_input_polygon(i, std::back_inserter(faces));

      for (const auto& f : faces)
        m_face_inliers[m_face2index[f]] = std::vector<std::size_t>();

      for (std::size_t j = 0; j < faces.size(); j++) {
        std::size_t idx = m_face2index[faces[j]];
        std::vector<Point_3> face;
        m_kinetic_partition.vertices(faces[j], std::back_inserter(face));

        //multiple regions per input polygon

        Delaunay_2 tri;
        std::vector<Point_2> f2d;
        for (const Point_3& p : face) {
          f2d.push_back(m_regions[i].first.to_2d(p));
          tri.insert(m_regions[i].first.to_2d(p));
        }

        // Get area
        for (auto fit = tri.finite_faces_begin(); fit != tri.finite_faces_end(); ++fit) {
          const Triangle_2 triangle(
            fit->vertex(0)->point(),
            fit->vertex(1)->point(),
            fit->vertex(2)->point());
          m_face_area[idx] += triangle.area();
        }

        total_area += m_face_area[idx];

        for (const std::size_t index : m_regions[i].second) {
          Point_2 p = m_regions[i].first.to_2d(get(m_point_map, index));
          const auto fh = tri.locate(p);
          if (fh != nullptr && !tri.is_infinite(fh)) {
            m_face_inliers[idx].push_back(index);
            m_total_inliers++;
          }
        }
      }
    }

    // Handling face generated by the octree partition. They are not associated with an input polygon.
    for (std::size_t i = 0; i < m_faces.size(); i++) {
      if (m_face_area[i] == 0 && m_face_neighbors[i].second > 6) {
        std::vector<Point_3> face;
        m_kinetic_partition.vertices(m_faces[i], std::back_inserter(face));

        Plane_3 pl;
        CGAL::linear_least_squares_fitting_3(face.begin(), face.end(), pl, CGAL::Dimension_tag<0>());

        Delaunay_2 tri;
        for (const Point_3& p : face)
          tri.insert(pl.to_2d(p));

        // Get area
        for (auto fit = tri.finite_faces_begin(); fit != tri.finite_faces_end(); ++fit) {
          const Triangle_2 triangle(
            fit->vertex(0)->point(),
            fit->vertex(1)->point(),
            fit->vertex(2)->point());
          m_face_area[i] += triangle.area();
        }

        total_area += m_face_area[i];
      }
    }

    for (std::size_t i = 0; i < m_faces.size(); i++) {
      // If the area is 0 it is a boundary face.
      if (m_face_area[i] == 0)
        m_face_area[i] = 2.0 * m_total_inliers;
      else
        m_face_area[i] = m_face_area[i] * 2.0 * m_total_inliers / total_area;
    }
  }

  void collect_points_for_faces2() {
    FT total_area = 0;
    m_total_inliers = 0;
    From_exact from_exact;
    auto& reg2input = m_kinetic_partition.regularized_input_mapping();
    std::cout << reg2input.size() << std::endl;
    std::size_t next = 0, step = 1;
    for (std::size_t i = 0; i < reg2input.size(); i++) {

      std::vector<KSP::Index> faces;
      m_kinetic_partition.faces_of_regularized_polygon(i, std::back_inserter(faces));

      for (const auto& f : faces)
        m_face_inliers[m_face2index[f]] = std::vector<std::size_t>();

      for (std::size_t j = 0; j < faces.size(); j++) {

        std::size_t idx = m_face2index[faces[j]];
        std::vector<Point_3> face;
        m_kinetic_partition.vertices(faces[j], std::back_inserter(face));

        //multiple regions per input polygon

        Plane_3 pl = from_exact(m_kinetic_partition.regularized_plane(i));

        FT max_dist1 = 0, max_dist2 = 0;

        Delaunay_2 tri;
        std::vector<Point_2> f2d;
        for (const Point_3& p : face) {
          max_dist1 = (std::max<double>)(sqrt((pl.to_3d(pl.to_2d(p)) - p).squared_length()), max_dist1);
          f2d.push_back(pl.to_2d(p));
          tri.insert(pl.to_2d(p));
        }

        // Get area
        for (auto fit = tri.finite_faces_begin(); fit != tri.finite_faces_end(); ++fit) {
          const Triangle_2 triangle(
            fit->vertex(0)->point(),
            fit->vertex(1)->point(),
            fit->vertex(2)->point());
          m_face_area[idx] += triangle.area();
        }

        total_area += m_face_area[idx];
        for (const std::size_t& p : reg2input[i]) {
          for (const std::size_t index : m_regions[p].second) {
            Point_2 p = pl.to_2d(get(m_point_map, index));

            max_dist2 = (std::max<double>)(sqrt((pl.to_3d(p) - get(m_point_map, index)).squared_length()), max_dist2);
            const auto fh = tri.locate(p);
            if (fh != nullptr && !tri.is_infinite(fh)) {
              m_face_inliers[idx].push_back(index);
              m_total_inliers++;
            }
          }
        }
      }
    }

    set_outside_volumes(m_cost_matrix);

    // Handling face generated by the octree partition. They are not associated with an input polygon.
    for (std::size_t i = 0; i < m_faces.size(); i++) {
      if (m_face_area[i] == 0) {//}&& m_face_neighbors[i].second > 6) {
        std::vector<Point_3> face;
        m_kinetic_partition.vertices(m_faces[i], std::back_inserter(face));

        Plane_3 pl;
        CGAL::linear_least_squares_fitting_3(face.begin(), face.end(), pl, CGAL::Dimension_tag<0>());

        Delaunay_2 tri;
        for (const Point_3& p : face)
          tri.insert(pl.to_2d(p));

        // Get area
        for (auto fit = tri.finite_faces_begin(); fit != tri.finite_faces_end(); ++fit) {
          const Triangle_2 triangle(
            fit->vertex(0)->point(),
            fit->vertex(1)->point(),
            fit->vertex(2)->point());
          m_face_area[i] += triangle.area();
        }

        total_area += m_face_area[i];
      }
    }

    for (std::size_t i = 0; i < m_faces.size(); i++) {
      // Check boundary faces and if the outside node has a defined value, if not, set area to 0.

      if (m_face_neighbors[i].second < 6 && m_cost_matrix[0][m_face_neighbors[i].second] == m_cost_matrix[1][m_face_neighbors[i].second]) {
        m_face_area[i] = 0;
      }
      else
        m_face_area[i] = m_face_area[i] * 2.0 * m_total_inliers / total_area;
    }
  }
*/

  void collect_points_for_faces3() {
    FT total_area = 0;
    m_total_inliers = 0;
    From_exact from_exact;
    auto& reg2input = m_kinetic_partition.regularized_input_mapping();
    std::cout << reg2input.size() << std::endl;
    std::size_t next = 0, step = 1;
    for (std::size_t i = 0; i < reg2input.size(); i++) {

      std::vector<std::pair<KSP::Index, std::vector<std::size_t> > > mapping;
      std::size_t fusioned_input_regions = 0;
      for (const auto& p : reg2input[i])
        fusioned_input_regions += m_regions[p].second.size();

      std::vector<Point_3> pts;
      std::vector<std::size_t> pts_idx;
      pts.reserve(fusioned_input_regions);
      for (const auto& p : reg2input[i]) {
        for (std::size_t j = 0; j < m_regions[p].second.size(); j++) {
          pts.emplace_back(get(m_point_map, m_regions[p].second[j]));
          pts_idx.push_back(m_regions[p].second[j]);
        }
      }
      m_kinetic_partition.map_points_to_regularized_polygons(i, pts, mapping);

      // Still need to calculate the area
      // Remap from mapping to m_face_inliers
      for (auto p : mapping) {
        assert(m_face_inliers[m_face2index[p.first]].size() == 0);
        m_face_inliers[m_face2index[p.first]].resize(p.second.size());
        for (std::size_t i = 0; i < p.second.size(); i++)
          m_face_inliers[m_face2index[p.first]][i] = pts_idx[p.second[i]];

        m_total_inliers += p.second.size();
      }

      std::vector<KSP::Index> faces;
      m_kinetic_partition.faces_of_regularized_polygon(i, std::back_inserter(faces));

      std::vector<std::vector<std::size_t> > faces2d(faces.size());
      std::vector<std::vector<std::size_t> > indices; // Adjacent faces for each vertex
      std::map<KSP::Index, std::size_t> idx2pts; // Mapping of vertices to pts vector

      Plane_3 pl = from_exact(m_kinetic_partition.regularized_plane(i));

      for (std::size_t j = 0; j < faces.size(); j++) {
        std::size_t idx = m_face2index[faces[j]];
        std::vector<Point_3> face;
        m_kinetic_partition.vertices(faces[j], std::back_inserter(face));

        //multiple regions per input polygon

        FT max_dist1 = 0, max_dist2 = 0;

        Delaunay_2 tri;
        std::vector<Point_2> f2d;

        for (const Point_3& p : face) {
          max_dist1 = (std::max<double>)(sqrt((pl.to_3d(pl.to_2d(p)) - p).squared_length()), max_dist1);
          f2d.push_back(pl.to_2d(p));
          tri.insert(pl.to_2d(p));
        }

        // Get area
        for (auto fit = tri.finite_faces_begin(); fit != tri.finite_faces_end(); ++fit) {
          const Triangle_2 triangle(
            fit->vertex(0)->point(),
            fit->vertex(1)->point(),
            fit->vertex(2)->point());
          m_face_area[idx] += triangle.area();
        }
        total_area += m_face_area[idx];
      }
    }

    set_outside_volumes(m_cost_matrix);

    // Handling face generated by the octree partition. They are not associated with an input polygon.
    for (std::size_t i = 0; i < m_faces.size(); i++) {
      if (m_face_area[i] == 0) {//}&& m_face_neighbors[i].second > 6) {
        std::vector<Point_3> face;
        m_kinetic_partition.vertices(m_faces[i], std::back_inserter(face));

        Plane_3 pl;
        CGAL::linear_least_squares_fitting_3(face.begin(), face.end(), pl, CGAL::Dimension_tag<0>());

        Delaunay_2 tri;
        for (const Point_3& p : face)
          tri.insert(pl.to_2d(p));

        // Get area
        for (auto fit = tri.finite_faces_begin(); fit != tri.finite_faces_end(); ++fit) {
          const Triangle_2 triangle(
            fit->vertex(0)->point(),
            fit->vertex(1)->point(),
            fit->vertex(2)->point());
          m_face_area[i] += triangle.area();
        }

        total_area += m_face_area[i];
      }
    }

    for (std::size_t i = 0; i < m_faces.size(); i++) {
      // Check boundary faces and if the outside node has a defined value, if not, set area to 0.

      if (m_face_neighbors[i].second < 6 && m_cost_matrix[0][m_face_neighbors[i].second] == m_cost_matrix[1][m_face_neighbors[i].second]) {
        m_face_area[i] = 0;
      }
      else
        m_face_area[i] = m_face_area[i] * 2.0 * m_total_inliers / total_area;
    }
  }

  void count_volume_votes() {
    const int debug_volume = -1;
    FT total_volume = 0;
    std::size_t num_volumes = m_kinetic_partition.number_of_volumes();
    m_volume_votes.resize(num_volumes, std::make_pair(0, 0));

    std::size_t count_faces = 0;
    std::size_t count_points = 0;

    m_volumes.resize(num_volumes, 0);
    std::size_t next = 0, step = num_volumes / 100;
    for (std::size_t i = 0; i < num_volumes; i++) {

      std::vector<Index> faces;
      m_kinetic_partition.faces(i, std::back_inserter(faces));
      const Point_3& centroid = m_kinetic_partition.volume_centroid(i);
      std::size_t in_count = 0;
      std::size_t out_count = 0;

      std::vector<Point_3> inside, outside;
      std::vector<Vector_3> insideN, outsideN;

      // Count votes based on points on faces and their normals
      for (const Index& f : faces) {
        auto it = m_face2index.find(f);

        // If the face is not contained, it belongs to the bounding box and can be skipped.
        if (it == m_face2index.end())
          continue;

        count_faces++;

        std::size_t idx = it->second;

        for (std::size_t p : m_face_inliers[idx]) {
          const auto& point = get(m_point_map, p);
          const auto& normal = get(m_normal_map, p);

          count_points++;

          const Vector_3 vec(point, centroid);
          const FT dot_product = vec * normal;
          if (dot_product < FT(0)) {
            inside.push_back(point);
            insideN.push_back(normal);
            in_count++;
          }
          else {
            outside.push_back(point);
            outsideN.push_back(normal);
            out_count++;
          }
        }
      }

      // Calculate volume
      std::vector<Point_3> volume_vertices;

      for (const Index& f : faces)
        m_kinetic_partition.vertices(f, std::back_inserter(volume_vertices));

      Delaunay_3 tri;
      for (const Point_3& p : volume_vertices)
        tri.insert(p);

      m_volumes[i] = FT(0);
      for (auto cit = tri.finite_cells_begin(); cit != tri.finite_cells_end(); ++cit) {
        const auto& tet = tri.tetrahedron(cit);
        m_volumes[i] += tet.volume();
      }

      total_volume += m_volumes[i];

      m_volume_votes[i] = std::make_pair(in_count, out_count);
      m_cost_matrix[0][i + 6] = static_cast<double>(in_count);
      m_cost_matrix[1][i + 6] = static_cast<double>(out_count);

      if (i == debug_volume) {
        std::vector<CGAL::Color> colors;
        colors.resize(inside.size(), CGAL::Color(0, 255, 0));
        colors.resize(outside.size() + inside.size(), CGAL::Color(0, 0, 255));
        inside.reserve(inside.size() + outside.size());
        std::copy(outside.begin(), outside.end(), std::back_inserter(inside));
        insideN.reserve(inside.size() + outside.size());
        std::copy(outsideN.begin(), outsideN.end(), std::back_inserter(insideN));
        CGAL::KSR_3::dump_points(inside, insideN, colors, std::to_string(i) + "-votes");
      }
    }

    // Normalize volumes
    for (FT& v : m_volumes)
      v /= total_volume;
  }

  FT calculate_volume(std::size_t volume_index) const {
    return 0;
  }

  template<typename NamedParameters>
  void create_planar_shapes(const NamedParameters& np) {

    if (m_points.size() < 3) {
      if (m_verbose) std::cout << "* no points found, skipping" << std::endl;
      return;
    }
    if (m_verbose) std::cout << "* getting planar shapes using region growing" << std::endl;

    // Parameters.
    const std::size_t k = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::k_neighbors), 12);
    const FT max_distance_to_plane = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::maximum_distance), FT(1));
    const FT max_accepted_angle = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::maximum_angle), FT(15));
    const std::size_t min_region_size = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::minimum_region_size), 50);

    m_detection_distance_tolerance = max_distance_to_plane;

    // Region growing.
    Neighbor_query neighbor_query = CGAL::Shape_detection::Point_set::make_k_neighbor_query(
      m_points, CGAL::parameters::k_neighbors(k));

    Region_type region_type = CGAL::Shape_detection::Point_set::make_least_squares_plane_fit_region(
      m_points,
      CGAL::parameters::
      maximum_distance(max_distance_to_plane).
      maximum_angle(max_accepted_angle).
      minimum_region_size(min_region_size));

    Sorting sorting = CGAL::Shape_detection::Point_set::make_least_squares_plane_fit_sorting(m_points, neighbor_query);
    sorting.sort();

    Region_growing region_growing(
      m_points, sorting.ordered(), neighbor_query, region_type);
    region_growing.detect(std::back_inserter(m_regions));

    // Convert indices.
    m_planar_regions.clear();
    m_planar_regions.reserve(m_regions.size());

    // Copy planes for regularization.
    std::vector<Plane_3> planes(m_regions.size());
    for (std::size_t i = 0; i < m_regions.size(); i++)
      planes[i] = m_regions[i].first;

    auto range = m_regions | boost::adaptors::transformed([](typename Region_growing::Primitive_and_region& pr)->Plane_3& {return pr.first; });

    const bool regularize_axis_symmetry = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::regularize_axis_symmetry), false);
    const bool regularize_coplanarity = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::regularize_coplanarity), false);
    const bool regularize_orthogonality = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::regularize_orthogonality), false);
    const bool regularize_parallelism = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::regularize_parallelism), false);
    const FT angle_tolerance = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::angle_tolerance), 25);
    const FT maximum_offset = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::maximum_offset), 0.01);

    // Regularize detected planes.
/*
    CGAL::Shape_regularization::Planes::regularize_planes(range, m_points,
      CGAL::parameters::plane_index_map(region_growing.region_map())
      .point_map(m_point_map)
      .regularize_axis_symmetry(regularize_axis_symmetry)
      .regularize_orthogonality(regularize_orthogonality)
      .regularize_parallelism(regularize_parallelism)
      .regularize_coplanarity(regularize_coplanarity)
      .maximum_angle(angle_tolerance)
      .maximum_offset(maximum_offset));*/

    // Merge coplanar regions
    for (std::size_t i = 0; i < m_regions.size() - 1; i++) {
      for (std::size_t j = i + 1; j < m_regions.size(); j++) {
        if (m_regions[i].first == m_regions[j].first || m_regions[i].first.opposite() == m_regions[j].first) {
          std::move(m_regions[j].second.begin(), m_regions[j].second.begin(), std::back_inserter(m_regions[i].second));
          m_regions.remove(m_regions.begin() + j);
          j--
        }
      }
    }

    std::vector<Plane_3> pl;

    std::size_t idx = 0;
    for (const auto& p : range) {
      bool exists = false;
      for (std::size_t i = 0; i < pl.size(); i++)
        if (pl[i] == p || pl[i].opposite() == p) {
          //merged[i].push_back(idx);
          exists = true;
        }

      if (!exists) {
        pl.push_back(p);
      }
      idx++;
    }

    for (const auto& pair : m_regions) {
      Indices region;
      for (auto& i : pair.second)
        region.push_back(i);
      m_planar_regions.push_back(region);
      //const auto plane = fit_plane(region);
      const std::size_t shape_idx = add_convex_hull_shape(region, pair.first);
      CGAL_assertion(shape_idx != std::size_t(-1));
      m_region_map[shape_idx] = region;
    }
    CGAL_assertion(m_planar_regions.size() == m_regions.size());

    std::size_t unassigned = 0;

    region_growing.unassigned_items(m_points, boost::make_function_output_iterator([&](const auto&) { ++unassigned; }));

    std::cout << "found " << m_polygons.size() << " planar shapes regularized into " << pl.size() << std::endl;
    std::cout << "from " << m_points.size() << " input points " << unassigned << " remain unassigned" << std::endl;
  }

  const Plane_3 fit_plane(const std::vector<std::size_t>& region) const {

    std::vector<Point_3> points;
    points.reserve(region.size());
    for (const std::size_t idx : region) {
      CGAL_assertion(idx < m_points.size());
      points.push_back(get(m_point_map, idx));
    }
    CGAL_assertion(points.size() == region.size());

    Plane_3 fitted_plane;
    Point_3 fitted_centroid;
    CGAL::linear_least_squares_fitting_3(
      points.begin(), points.end(),
      fitted_plane, fitted_centroid,
      CGAL::Dimension_tag<0>());

    const Plane_3 plane(
      static_cast<FT>(fitted_plane.a()),
      static_cast<FT>(fitted_plane.b()),
      static_cast<FT>(fitted_plane.c()),
      static_cast<FT>(fitted_plane.d()));
    return plane;
  }

  void set_outside_volumes(std::vector<std::vector<double> >& cost_matrix) const {
    // Setting preferred outside label for bbox plane nodes
    // Order:
    // 0 zmin
    // 1 ymin
    // 2 xmax
    // 3 ymax
    // 4 xmin
    // 5 zmax
    const std::size_t force = m_total_inliers * 3;
    cost_matrix[0][0] = 0;
    cost_matrix[0][1] = 0;
    cost_matrix[0][2] = 0;
    cost_matrix[0][3] = 0;
    cost_matrix[0][4] = 0;
    cost_matrix[0][5] = 0;
    cost_matrix[1][0] = 0;
    cost_matrix[1][1] = 0;
    cost_matrix[1][2] = 0;
    cost_matrix[1][3] = 0;
    cost_matrix[1][4] = 0;
    cost_matrix[1][5] = 0;
  }

  void dump_volume(std::size_t i, const std::string& filename, const CGAL::Color &color) const {
    std::vector<KSP::Index> faces;
    m_kinetic_partition.faces(i, std::back_inserter(faces));

    std::vector<std::vector<Point_3> > pts(faces.size());
    std::vector<CGAL::Color> col(faces.size(), color);
    for (std::size_t j = 0; j < faces.size(); j++) {
      m_kinetic_partition.vertices(faces[j], std::back_inserter(pts[j]));
    }

    CGAL::KSR_3::dump_polygons(pts, col, filename);
  }

  void dump_face(std::size_t i, const std::string& filename, const CGAL::Color& color) const {
    std::vector<Point_3> face;
    m_kinetic_partition.vertices(m_faces[i], std::back_inserter(face));
  }
};

#endif

} // namespace CGAL


#endif // CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H
