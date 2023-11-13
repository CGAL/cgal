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
#include <CGAL/bounding_box.h>

#include <boost/filesystem.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/helpers.h>

#define WITH_LCC

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
    const NamedParameters& np = CGAL::parameters::default_values()) : m_points(ps), m_ground_polygon_index(-1), m_kinetic_partition(np) {}

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
  std::size_t detect_planar_shapes(bool estimate_ground = false,
    const CGAL_NP_CLASS& np = parameters::default_values()) {

    if (m_verbose)
      std::cout << std::endl << "--- DETECTING PLANAR SHAPES: " << std::endl;

    m_planes.clear();
    m_polygons.clear();
    m_region_map.clear();

    m_point_map = Point_set_processing_3_np_helper<Point_set, CGAL_NP_CLASS, Point_map>::get_point_map(m_points, np);
    m_normal_map = Point_set_processing_3_np_helper<Point_set, CGAL_NP_CLASS, Normal_map>::get_normal_map(m_points, np);

    create_planar_shapes(estimate_ground, np);

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

    m_kinetic_partition.get_linear_cell_complex(m_lcc);
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
    if (m_lcc.one_dart_per_cell<3>().size() == 0) {
      std::cout << "Kinetic partition is not constructed or does not have volumes" << std::endl;
      return false;
    }

    m_face_area.clear();
    m_face_inliers.clear();

    auto face_range = m_lcc.one_dart_per_cell<2>();
    m_faces_lcc.reserve(face_range.size());

    for (auto& d : face_range) {
      m_faces_lcc.push_back(m_lcc.dart_descriptor(d));

      std::size_t id = m_lcc.attribute<2>(m_faces_lcc.back());

      auto p = m_attrib2index_lcc.emplace(std::make_pair(m_lcc.attribute<2>(m_faces_lcc.back()), m_faces_lcc.size() - 1));
      CGAL_assertion(p.second);
    }

    // Create value arrays for graph cut
    m_face_inliers.resize(m_faces_lcc.size());
    m_face_area.resize(m_faces_lcc.size());
    m_face_area_lcc.resize(m_faces_lcc.size());
    m_face_neighbors_lcc.resize(m_faces_lcc.size(), std::pair<int, int>(-1, -1));

    m_cost_matrix.resize(2);
    m_cost_matrix[0].resize(m_kinetic_partition.number_of_volumes() + 6);
    m_cost_matrix[1].resize(m_kinetic_partition.number_of_volumes() + 6);

    for (std::size_t i = 0; i < m_faces_lcc.size(); i++) {
      auto n = m_lcc.one_dart_per_incident_cell<3, 2>(m_faces_lcc[i]);

      assert(n.size() == 1 || n.size() == 2);
      auto it = n.begin();

      auto& finf = m_lcc.info<2>(m_faces_lcc[i]);

      int first = m_lcc.info<3>(m_lcc.dart_descriptor(*it)).volume_index;
      auto& inf1 = m_lcc.info<3>(m_lcc.dart_descriptor(*it++));

      auto inf2 = inf1;
      if (n.size() == 2)
        inf2 = m_lcc.info<3>(m_lcc.dart_descriptor(*it));

      int second;
      if (n.size() == 2)
        second = m_lcc.info<3>(m_lcc.dart_descriptor(*it)).volume_index;

      if (n.size() == 2)
        m_face_neighbors_lcc[i] = std::make_pair(first + 6, m_lcc.info<3>(m_lcc.dart_descriptor(*it)).volume_index + 6);
      else
        m_face_neighbors_lcc[i] = std::make_pair(first + 6, -m_lcc.info<2>(m_faces_lcc[i]).input_polygon_index - 1);

      if (m_face_neighbors_lcc[i].first > m_face_neighbors_lcc[i].second)
        m_face_neighbors_lcc[i] = std::make_pair(m_face_neighbors_lcc[i].second, m_face_neighbors_lcc[i].first);

      if (m_face_neighbors_lcc[i].first < m_face_neighbors_lcc[i].second) {
        auto it = m_neighbors2index_lcc.emplace(std::make_pair(m_face_neighbors_lcc[i], i));
        assert(it.second);
      }
    }

    check_ground();

    m_face_inliers.clear();
    m_face_inliers.resize(m_faces_lcc.size());
    collect_points_for_faces_lcc();
    count_volume_votes_lcc();

    std::cout << "* computing data term ... ";

    std::size_t max_inside = 0, max_outside = 0;
    for (std::size_t i = 0; i < m_volumes.size(); i++) {
      max_inside = (std::max<double>)(m_cost_matrix[0][i + 6], max_inside);
      max_outside = (std::max<double>)(m_cost_matrix[1][i + 6], max_outside);
    }

    // Dump volumes colored by votes
/*
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
*/

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
#ifdef WITH_LCC
    gc.solve(m_face_neighbors_lcc, m_face_area_lcc, m_cost_matrix, m_labels);
#else
    gc.solve(m_face_neighbors, m_face_area, m_cost_matrix, m_labels);
#endif

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

    From_exact from_exact;

    std::map<Point_3, std::size_t> pt2idx;

    for (std::size_t i = 0; i < m_faces_lcc.size(); i++) {
      const auto& n = m_face_neighbors_lcc[i];
      // Do not extract boundary faces.
      if (n.second < 6)
        continue;
      if (m_labels[n.first] != m_labels[n.second]) {
        std::vector<Point_3> face;

        for (const auto& vd : m_lcc.one_dart_per_incident_cell<0, 2>(m_faces_lcc[i]))
          face.push_back(from_exact(m_lcc.point(m_lcc.dart_descriptor(vd))));

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
  \brief Provides the reconstructed surface as a list of indexed polygons.

  \param pit
  an output iterator taking Point_3.

  \param triit
  an output iterator taking std::vector<std::size_t>.

  \pre `successful reconstruction`
  */
  template<class OutputPointIterator, class OutputPolygonIterator>
  void reconstructed_model_polylist_lcc(OutputPointIterator pit, OutputPolygonIterator polyit) {
    if (m_labels.empty())
      return;

    From_exact from_exact;

    std::map<Point_3, std::size_t> pt2idx;

    std::vector<int> region_index(m_faces_lcc.size(), -1);
    std::size_t region = 0;

    for (std::size_t i = 0; i < m_faces_lcc.size(); i++) {
      const auto& n = m_face_neighbors_lcc[i];
      if (n.second < 6)
        continue;
      if (m_labels[n.first] != m_labels[n.second]) {
        Face_attribute fa = m_lcc.attribute<2>(m_faces_lcc[i]);

        if (region_index[fa] == -1) {
          collect_connected_component(m_faces_lcc[i], region_index, region++);
        }
      }
    }


    /*

    */
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

  typedef CGAL::Linear_cell_complex_traits<3, CGAL::Exact_predicates_exact_constructions_kernel> Traits;
  using LCC = CGAL::Linear_cell_complex_for_combinatorial_map<3, 3, Traits, typename KSP::LCC_Base_Properties>;
  using Dart_descriptor = typename LCC::Dart_descriptor;
  using Dart = typename LCC::Dart;

  //using Visibility = KSR_3::Visibility<Kernel, Intersection_kernel, Point_map, Normal_map>;
  using Index = typename KSP::Index;
  using Face_attribute = typename LCC::Base::template Attribute_descriptor<2>::type;

  bool m_verbose;
  bool m_debug;

  Point_set &m_points;
  Point_map m_point_map;
  Normal_map m_normal_map;

  std::vector<std::vector<std::size_t> > m_planar_regions;
  std::vector<typename Region_growing::Primitive_and_region> m_regions;
  std::map<std::size_t, Indices> m_region_map;
  double m_detection_distance_tolerance;

  std::size_t m_ground_polygon_index;
  Plane_3 m_ground_plane;

  std::vector<Plane_3> m_planes;
  std::vector<Point_3> m_polygon_pts;
  std::vector<std::vector<std::size_t> > m_polygon_indices;
  std::vector<Polygon_3> m_polygons;
  KSP m_kinetic_partition;

  LCC m_lcc;
  std::vector<typename LCC::Dart_const_descriptor> m_faces_lcc;
  std::map<Face_attribute, std::size_t> m_attrib2index_lcc;
  std::vector<std::size_t> lcc2index;
  std::vector<std::size_t> index2lcc;

  // Face indices are now of type Indices and are not in a range 0 to n
  std::vector<Indices> m_face_inliers;
  std::vector<FT> m_face_area, m_face_area_lcc;
  std::vector<std::pair<std::size_t, std::size_t> > m_face_neighbors_lcc;
  std::map<std::pair<std::size_t, std::size_t>, std::size_t> m_neighbors2index_lcc;

  std::vector<std::pair<std::size_t, std::size_t> > m_volume_votes; // pair<inside, outside> votes
  std::vector<bool> m_volume_below_ground;
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

  void store_convex_hull_shape(const std::vector<std::size_t>& region, const Plane_3& plane, std::vector<std::vector<Point_3> > &polys) {

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

    polys.push_back(polygon);
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

  void check_ground() {
    std::size_t num_volumes = m_kinetic_partition.number_of_volumes();
    // Set all volumes to not be below the ground, this leads to the standard 6 outside node connection.
    m_volume_below_ground.resize(num_volumes, false);
    From_exact from_exact;

    if (m_ground_polygon_index != -1)
      for (const auto &vd : m_lcc.one_dart_per_cell<3>()) {
        const auto& info = m_lcc.info<3>(m_lcc.dart_descriptor(vd));

        m_volume_below_ground[info.volume_index] = (from_exact(info.barycenter) - m_regions[m_ground_polygon_index].first.projection(from_exact(info.barycenter))).z() < 0;
      }
  }

  void collect_connected_component(typename LCC::Dart_descriptor face, std::vector<int> &region_index, std::size_t region) {
    // What does the caller have to manage? a map from face_attrib to bool to collect treated faces?

    std::queue<std::size_t> face_queue;
    face_queue.push(face);

    From_exact from_exact;

    std::vector<Dart_descriptor> border_edges;

    std::vector<Point_3> pts;

    int ip = m_lcc.info<2>(face).input_polygon_index;
    typename Intersection_kernel::Plane_3 pl = m_lcc.info<2>(face).plane;

    while (!face_queue.empty()) {
      Dart_descriptor cur_fdh(face_queue.front());
      Face_attribute cur_fa = m_lcc.attribute<2>(cur_fdh);
      face_queue.pop();

      if (region_index[cur_fa] == region)
        continue;

      region_index[cur_fa] = region;

      // Iterate over edges of face.
      for (auto& ed : m_lcc.one_dart_per_incident_cell<1, 2>(cur_fdh)) {
        Dart_descriptor edh = m_lcc.dart_descriptor(ed);

        for (auto &fd : m_lcc.one_dart_per_incident_cell<2, 1, 3>(edh)) {
          Dart_descriptor fdh = m_lcc.dart_descriptor(fd);
          Face_attribute fa = m_lcc.attribute<2>(fdh);
          auto& inf = m_lcc.info<2>(fdh);
          bool added = false;

          const auto& n = m_face_neighbors_lcc[m_attrib2index_lcc[fa]];

          // Do not segment bbox surface
          if (n.second < 6)
            continue;

          // Belongs to reconstruction?
          if (m_labels[n.first] == m_labels[n.second])
            continue;

          // Already segmented?
          if (region_index[fa] != -1)
            continue;

          if (ip != -7) {
            if (m_lcc.info<2>(fdh).input_polygon_index == ip) {
              added = true;
              face_queue.push(fdh);
            }
          }
          else
            if (m_lcc.info<2>(fdh).plane == pl) {
              added = true;
              face_queue.push(fdh);
            }

          if (!added)
            border_edges.push_back(edh);
        }
      }
    }

    pts.clear();
    for (Dart_descriptor &edh : border_edges)
      for (auto& vd : m_lcc.one_dart_per_incident_cell<0, 1>(edh)) {
        pts.push_back(from_exact(m_lcc.point(m_lcc.dart_descriptor(vd))));
      }

    std::ofstream vout(std::to_string(region) + ".xyz");
    for (Point_3& p : pts) {
      vout << p.x() << " " << p.y() << " " << p.z() << std::endl;
    }
    vout << std::endl;
    vout.close();

    // Iterate over all edges and collect all face descriptors and border edges.
    // For each edge, find neighbor face m_lcc.one_dart_of_incident_cell<2, 1>
     // To find neighbor face, m_lcc.one_dart_of_incident_cell<1, 2>?
      // Exclude faces by region_index (including own face)
    // Identifying other face by comparing input_polygon (special case -7)
    // Two sets -> faces set<attribute_descriptor> and edges std::set?

    // After collecting faces, I can collect border edges and just start from one to get the loop.
  }

  void collect_points_for_faces_lcc() {
    FT total_area = 0;
    m_total_inliers = 0;
    From_exact from_exact;

    std::vector<std::vector<Dart_descriptor> > poly2faces(m_kinetic_partition.input_planes().size());
    std::vector<Dart_descriptor> other_faces;
    for (auto& d : m_lcc.one_dart_per_cell<2>()) {
      Dart_descriptor dh = m_lcc.dart_descriptor(d);
      if (m_lcc.info<2>(dh).input_polygon_index >= 0)
        poly2faces[m_lcc.info<2>(dh).input_polygon_index].push_back(dh);
      else
        other_faces.push_back(dh); // Contains faces originating from the octree decomposition as well as bbox faces
    }

    assert(m_kinetic_partition.input_planes().size() == m_regions.size());

    std::size_t next = 0, step = 1;
    for (std::size_t i = 0; i < m_kinetic_partition.input_planes().size(); i++) {

      std::vector<std::pair<Dart_descriptor, std::vector<std::size_t> > > mapping;

      std::vector<Point_3> pts;
      pts.reserve(m_regions[i].second.size());

      for (std::size_t j = 0; j < m_regions[i].second.size(); j++)
        pts.emplace_back(get(m_point_map, m_regions[i].second[j]));

      map_points_to_faces(i, pts, mapping);

      // Remap from mapping to m_face_inliers
      for (auto p : mapping) {
        //m_face_inliers[m_attrib2index_lcc[m_lcc.attribute<2>(p.first)]].clear();
        Face_attribute& f = m_lcc.attribute<2>(p.first);
        std::size_t id = m_attrib2index_lcc[f];
        assert(m_face_inliers[id].size() == 0);

        m_face_inliers[m_attrib2index_lcc[m_lcc.attribute<2>(p.first)]].resize(p.second.size());
        for (std::size_t k = 0; k < p.second.size(); k++)
          m_face_inliers[m_attrib2index_lcc[m_lcc.attribute<2>(p.first)]][k] = m_regions[i].second[p.second[k]];

        m_total_inliers += p.second.size();
      }

      Plane_3 pl = from_exact(m_kinetic_partition.input_planes()[i]);

      for (std::size_t j = 0; j < poly2faces[i].size(); j++) {
        std::size_t idx = m_attrib2index_lcc[m_lcc.attribute<2>(poly2faces[i][j])];
        m_face_area_lcc[idx] = 0;

        //multiple regions per input polygon

        Delaunay_2 tri;

        Dart_descriptor n = poly2faces[i][j];
        do {
          tri.insert(pl.to_2d(from_exact(m_lcc.point(n))));
          n = m_lcc.beta(n, 0);
        } while (n != poly2faces[i][j]);

        // Get area
        for (auto fit = tri.finite_faces_begin(); fit != tri.finite_faces_end(); ++fit) {
          const Triangle_2 triangle(
            fit->vertex(0)->point(),
            fit->vertex(1)->point(),
            fit->vertex(2)->point());
          m_face_area_lcc[idx] += triangle.area();
        }

        total_area += m_face_area_lcc[idx];
      }
    }

    set_outside_volumes(m_cost_matrix);

    // Handling face generated by the octree partition. They are not associated with an input polygon.
    for (std::size_t i = 0; i < other_faces.size(); i++) {
      std::vector<Point_3> face;
      std::size_t idx = m_attrib2index_lcc[m_lcc.attribute<2>(other_faces[i])];

      Dart_descriptor n = other_faces[i];
      do {
        face.push_back(from_exact(m_lcc.point(n)));
        n = m_lcc.beta(n, 0);
      } while (n != other_faces[i]);

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
        m_face_area_lcc[idx] += triangle.area();
      }

      total_area += m_face_area_lcc[idx];
    }

    m_face_area_lcc.resize(m_faces_lcc.size(), 0);

    for (std::size_t i = 0; i < m_faces_lcc.size(); i++)
      m_face_area_lcc[i] = m_face_area_lcc[i] * 2.0 * m_total_inliers / total_area;

    std::size_t redirected = 0;
    for (std::size_t i = 0; i < m_face_neighbors_lcc.size(); i++) {
      // Check if the face is on a bbox face besides the top face.
      // If so and the connected volume is below the ground, redirect the face to the bottom face node.
      if (m_face_neighbors_lcc[i].second < 5 && m_volume_below_ground[m_face_neighbors_lcc[i].first - 6]) {
        redirected++;
        m_face_neighbors_lcc[i].second = 0;
      }
    }
    std::cout << redirected << " faces redirected to below ground" << std::endl;
  }

  void count_volume_votes_lcc() {
    const int debug_volume = -1;
    FT total_volume = 0;
    std::size_t num_volumes = m_kinetic_partition.number_of_volumes();
    m_volume_votes.clear();
    m_volume_votes.resize(num_volumes, std::make_pair(0, 0));

    m_volumes.resize(num_volumes, 0);

    for (std::size_t i = 0; i < num_volumes; i++) {
      m_cost_matrix[0][i] = m_cost_matrix[1][i] = 0;
      m_volumes[i] = 0;
    }

    std::size_t count_faces = 0;
    std::size_t count_points = 0;

    From_exact from_exact;

    std::size_t idx = 0;

    for (std::size_t i = 0; i < m_faces_lcc.size(); i++) {
      std::size_t v[] = { -1, -1 };
      Point_3 c[2];
      std::size_t in[] = {0, 0}, out[] = {0, 0};

      std::size_t idx = 0;
      for (auto& vd : m_lcc.one_dart_per_incident_cell<3, 2>(m_faces_lcc[i])) {
        typename LCC::Dart_descriptor vdh = m_lcc.dart_descriptor(vd);
        v[idx] = m_lcc.info<3>(vdh).volume_index;
        c[idx] = from_exact(m_lcc.info<3>(vdh).barycenter);
        idx++;
      }

      for (std::size_t p : m_face_inliers[i]) {
        const auto& point = get(m_point_map, p);
        const auto& normal = get(m_normal_map, p);

        count_points++;

        for (std::size_t j = 0; j < idx; j++) {
          const Vector_3 vec(point, c[j]);
          const FT dot_product = vec * normal;
          if (dot_product < FT(0))
            in[j]++;
          else
            out[j]++;
        }
      }

      for (std::size_t j = 0; j < idx; j++) {
        m_volume_votes[v[j]].first += in[j];
        m_volume_votes[v[j]].second += out[j];
        m_cost_matrix[0][v[j] + 6] += static_cast<double>(in[j]);
        m_cost_matrix[1][v[j] + 6] += static_cast<double>(out[j]);
      }
    }

    for (auto &d : m_lcc.one_dart_per_cell<3>()) {
      typename LCC::Dart_descriptor dh = m_lcc.dart_descriptor(d);

      std::vector<Point_3> volume_vertices;

      std::size_t volume_index = m_lcc.info<3>(dh).volume_index;

      // Collect all vertices of volume to calculate volume
      for (auto &fd : m_lcc.one_dart_per_incident_cell<2, 3>(dh)) {
        typename LCC::Dart_descriptor fdh = m_lcc.dart_descriptor(fd);

        for (const auto &vd : m_lcc.one_dart_per_incident_cell<0, 2>(fdh))
          volume_vertices.push_back(from_exact(m_lcc.point(m_lcc.dart_descriptor(vd))));
      }

      Delaunay_3 tri;
      for (const Point_3& p : volume_vertices)
        tri.insert(p);

      m_volumes[volume_index] = FT(0);
      for (auto cit = tri.finite_cells_begin(); cit != tri.finite_cells_end(); ++cit) {
        const auto& tet = tri.tetrahedron(cit);
        m_volumes[volume_index] += tet.volume();
      }

      total_volume += m_volumes[volume_index];
    }

    // Normalize volumes
    for (FT& v : m_volumes)
      v /= total_volume;

//     for (std::size_t i = 0; i < m_volumes.size(); i++)
//       std::cout << i << ": " << m_cost_matrix[0][i] << " o: " << m_cost_matrix[1][i] << " v: " << m_volumes[i] << std::endl;
  }

  template<typename NamedParameters>
  void create_planar_shapes(bool estimate_ground, const NamedParameters& np) {

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

    std::size_t unassigned = 0;
    region_growing.unassigned_items(m_points, boost::make_function_output_iterator([&](const auto&) { ++unassigned; }));

    std::vector<std::vector<Point_3> > polys_debug;

    for (std::size_t i = 0; i < m_regions.size(); i++) {

      Indices region;
      for (auto& j : m_regions[i].second)
        region.push_back(j);

      store_convex_hull_shape(region, m_regions[i].first, polys_debug);
      //KSR_3::dump_polygon(polys_debug[i], std::to_string(i) + "-detected-region.ply");
    }

    //KSR_3::dump_polygons(polys_debug, "detected-" + std::to_string(m_regions.size()) + "-polygons.ply");

    // Convert indices.
    m_planar_regions.clear();
    m_planar_regions.reserve(m_regions.size());

    // Copy planes for regularization.
    std::vector<Plane_3> planes(m_regions.size());
    for (std::size_t i = 0; i < m_regions.size(); i++)
      planes[i] = m_regions[i].first;

    auto range = m_regions | boost::adaptors::transformed([](typename Region_growing::Primitive_and_region& pr)->Plane_3& {return pr.first; });

    std::size_t num_shapes = m_regions.size();

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

    CGAL::Shape_regularization::Planes::regularize_planes(range, m_points,
      CGAL::parameters::plane_index_map(region_growing.region_map())
      .point_map(m_point_map)
      .regularize_axis_symmetry(regularize_axis_symmetry)
      .regularize_orthogonality(regularize_orthogonality)
      .regularize_parallelism(regularize_parallelism)
      .regularize_coplanarity(regularize_coplanarity)
      .maximum_angle(angle_tolerance)
      .maximum_offset(maximum_offset));

    polys_debug.clear();

    for (std::size_t i = 0; i < m_regions.size(); i++) {

      Indices region;
      for (auto& j : m_regions[i].second)
        region.push_back(j);

      store_convex_hull_shape(region, m_regions[i].first, polys_debug);
      //KSR_3::dump_polygon(polys_debug[i], std::to_string(i) + "-detected-region.ply");
    }

    //KSR_3::dump_polygons(polys_debug, "regularized-" + std::to_string(m_regions.size()) + "-polygons.ply");

    // Merge coplanar regions
    for (std::size_t i = 0; i < m_regions.size() - 1; i++) {
      for (std::size_t j = i + 1; j < m_regions.size(); j++) {
        if (m_regions[i].first == m_regions[j].first || m_regions[i].first.opposite() == m_regions[j].first) {
          std::move(m_regions[j].second.begin(), m_regions[j].second.end(), std::back_inserter(m_regions[i].second));
          m_regions.erase(m_regions.begin() + j);
          j--;
        }
      }
    }

    if (estimate_ground) {
      // How to estimate ground plane? Find low z peak
      std::vector<std::size_t> candidates;
      FT low_z_peak = (std::numeric_limits<FT>::max)();
      FT bbox_min[] = { (std::numeric_limits<FT>::max)(), (std::numeric_limits<FT>::max)(), (std::numeric_limits<FT>::max)() };
      FT bbox_max[] = { -(std::numeric_limits<FT>::max)(), -(std::numeric_limits<FT>::max)(), -(std::numeric_limits<FT>::max)() };
      for (const auto& p : m_points) {
        const auto& point = get(m_point_map, p);
        for (std::size_t i = 0; i < 3; i++) {
          bbox_min[i] = (std::min)(point[i], bbox_min[i]);
          bbox_max[i] = (std::max)(point[i], bbox_max[i]);
        }
      }

      FT bbox_center[] = { 0.5 * (bbox_min[0] + bbox_max[0]), 0.5 * (bbox_min[1] + bbox_max[1]), 0.5 * (bbox_min[2] + bbox_max[2]) };

      for (std::size_t i = 0; i < m_regions.size(); i++) {
        Vector_3 d = m_regions[i].first.orthogonal_vector();
        if (abs(d.z()) > 0.98) {
          candidates.push_back(i);
          FT z = m_regions[i].first.projection(Point_3(bbox_center[0], bbox_center[1], bbox_center[2])).z();
          low_z_peak = (std::min<FT>)(z, low_z_peak);
        }
      }

      m_ground_polygon_index = -1;
      std::vector<std::size_t> other_ground;
      for (std::size_t i = 0; i < candidates.size(); i++) {
        Vector_3 d = m_regions[candidates[i]].first.orthogonal_vector();
        FT z = m_regions[candidates[i]].first.projection(Point_3(bbox_center[0], bbox_center[1], bbox_center[2])).z();
        if (z - low_z_peak < max_distance_to_plane) {
          if (m_ground_polygon_index == -1)
            m_ground_polygon_index = candidates[i];
          else
            other_ground.push_back(candidates[i]);
        }
      }

      for (std::size_t i = 0; i < other_ground.size(); i++)
        std::move(m_regions[other_ground[i]].second.begin(), m_regions[other_ground[i]].second.end(), std::back_inserter(m_regions[m_ground_polygon_index].second));

      std::cout << "ground polygon " << m_ground_polygon_index << ", merging other polygons";
      while (other_ground.size() != 0) {
        m_regions.erase(m_regions.begin() + other_ground.back());
        std::cout << " " << other_ground.back();
        other_ground.pop_back();
      }
      std::cout << std::endl;

      std::vector<Point_3> ground_plane;
      ground_plane.reserve(m_regions[m_ground_polygon_index].second.size());
      for (std::size_t i = 0; i < m_regions[m_ground_polygon_index].second.size(); i++) {
        ground_plane.push_back(get(m_point_map, m_regions[m_ground_polygon_index].second[i]));
      }

      CGAL::linear_least_squares_fitting_3(ground_plane.begin(), ground_plane.end(), m_regions[m_ground_polygon_index].first, CGAL::Dimension_tag<0>());

      if (m_regions[m_ground_polygon_index].first.orthogonal_vector().z() < 0)
        m_regions[m_ground_polygon_index].first = m_regions[m_ground_polygon_index].first.opposite();
    }

    polys_debug.clear();

    for (std::size_t i = 0; i < m_regions.size(); i++) {

      Indices region;
      for (auto& j : m_regions[i].second)
        region.push_back(j);

      store_convex_hull_shape(region, m_regions[i].first, polys_debug);
      //KSR_3::dump_polygon(polys_debug[i], std::to_string(i) + "-detected-region.ply");
    }

    //KSR_3::dump_polygons(polys_debug, "merged-" + std::to_string(m_regions.size()) + "-polygons.ply");

    std::vector<Plane_3> pl;

    std::size_t idx = 0;
    for (const auto& p : m_regions) {
      bool exists = false;
      for (std::size_t i = 0; i < pl.size(); i++)
        if (pl[i] == p.first || pl[i].opposite() == p.first) {
          //merged[i].push_back(idx);
          exists = true;
        }

      if (!exists) {
        pl.push_back(p.first);
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

    std::cout << "found " << num_shapes << " planar shapes regularized into " << m_planar_regions.size() << std::endl;
    std::cout << "from " << m_points.size() << " input points " << unassigned << " remain unassigned" << std::endl;
  }

  void map_points_to_faces(const std::size_t polygon_index, const std::vector<Point_3>& pts, std::vector<std::pair<typename LCC::Dart_descriptor, std::vector<std::size_t> > >& face_to_points) {
    std::vector<Index> faces;

    if (polygon_index >= m_kinetic_partition.input_planes().size())
      assert(false);

    From_exact from_exact;

    const typename Intersection_kernel::Plane_3& pl = m_kinetic_partition.input_planes()[polygon_index];
    const Plane_3 inexact_pl = from_exact(pl);
    std::vector<Point_2> pts2d;
    pts2d.reserve(pts.size());

    for (const Point_3& p : pts)
      pts2d.push_back(inexact_pl.to_2d(p));

    //std::cout << "switch to Data_structure::m_face2sp" << std::endl;
    //ToDo I need to check whether the current way provides all faces as some faces may have been added during the make_conformal step

    // Iterate over all faces of the lcc
    for (Dart& d : m_lcc.one_dart_per_cell<2>()) {
      Dart_descriptor dd = m_lcc.dart_descriptor(d);
      if (m_lcc.info<2>(m_lcc.dart_descriptor(d)).input_polygon_index != polygon_index || !m_lcc.info<2>(m_lcc.dart_descriptor(d)).part_of_initial_polygon)
        continue;

      // No filtering of points per partition

      face_to_points.push_back(std::make_pair(m_lcc.dart_descriptor(d), std::vector<std::size_t>()));

      auto& info = m_lcc.info<2>(m_lcc.dart_descriptor(d));

      std::vector<Point_2> vts2d;
      vts2d.reserve(m_lcc.one_dart_per_incident_cell<0, 2>(m_lcc.dart_descriptor(d)).size());

      typename LCC::Dart_descriptor n = dd;
      do {
        vts2d.push_back(inexact_pl.to_2d(from_exact(m_lcc.point(n))));
        n = m_lcc.beta(n, 0);
      } while (n != dd);

      Polygon_2<Kernel> poly(vts2d.begin(), vts2d.end());
      if (poly.is_clockwise_oriented())
        std::reverse(vts2d.begin(), vts2d.end());

      for (std::size_t i = 0; i < pts2d.size(); i++) {
        const auto& pt = pts2d[i];
        bool outside = false;

        // poly, vertices and edges in IFace are oriented ccw
        std::size_t idx = 0;
        for (std::size_t i = 0; i < vts2d.size(); i++) {
          Vector_2 ts = (vts2d[(i + vts2d.size() - 1) % vts2d.size()]) - pt;
          Vector_2 tt = (vts2d[i]) - pt;

          bool ccw = (tt.x() * ts.y() - tt.y() * ts.x()) <= 0;
          if (!ccw) {
            outside = true;
            break;
          }
        }

        if (outside)
          continue;

        face_to_points.back().second.push_back(i);
      }
    }
  }

/*
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
*/

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
    // 0 - cost for labelled as outside
    cost_matrix[0][0] = 0;
    cost_matrix[0][1] = 0;
    cost_matrix[0][2] = 0;
    cost_matrix[0][3] = 0;
    cost_matrix[0][4] = 0;
    cost_matrix[0][5] = 0;
    // 1 - cost for labelled as inside
    cost_matrix[1][0] = 0;
    cost_matrix[1][1] = 0;
    cost_matrix[1][2] = 0;
    cost_matrix[1][3] = 0;
    cost_matrix[1][4] = 0;
    cost_matrix[1][5] = 0;

    if (m_ground_polygon_index != -1) {
      std::cout << "using estimated ground plane for reconstruction" << std::endl;
      cost_matrix[0][1] = 0;
      cost_matrix[0][2] = 0;
      cost_matrix[0][3] = 0;
      cost_matrix[0][4] = 0;
      cost_matrix[1][1] = force;
      cost_matrix[1][2] = force;
      cost_matrix[1][3] = force;
      cost_matrix[1][4] = force;
    }
  }

/*
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
  }*/
};

#endif

} // namespace CGAL


#endif // CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H
