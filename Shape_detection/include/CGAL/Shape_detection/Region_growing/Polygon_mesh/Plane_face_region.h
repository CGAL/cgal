// Copyright (c) 2024-2025 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sébastien Loriot
//

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_PLANE_FACE_REGION_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_PLANE_FACE_REGION_H

#include <CGAL/license/Shape_detection.h>

#include <CGAL/Shape_detection/Region_growing/internal/property_map.h>
#include <CGAL/Shape_detection/Region_growing/internal/utils.h>
#ifdef CGAL_SD_RG_USE_PMP
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#endif

namespace CGAL {
namespace Shape_detection {
namespace Polygon_mesh {

/*!
  \ingroup PkgShapeDetectionRGOnMesh

  \brief Region type based on the plane of the first face selected.

  This class uses the supporting plane of the first face picked for the region
  and expands it for all faces with a normal close to that of the first face
  ("close" being defined by the `maximum_distance` parameter), and such that vertices are
  not far from that supporting plane (far being defined by the `maximum_angle` or `cosine_of_maximum_angle` parameter).

  \tparam GeomTraits
  a model of `Kernel`

  \tparam PolygonMesh
  a model of `FaceListGraph`

  \tparam VertexToPointMap
  a model of `ReadablePropertyMap` whose key type is the vertex type of a polygon mesh
  (`boost::graph_traits<PolygonMesh>::vertex_descriptor`) and value type is `GeomTraits::Point_3`

  \cgalModels{RegionType}
*/
template<typename GeomTraits,
         typename PolygonMesh,
         typename VertexToPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::const_type>
class Plane_face_region {
public:
  /// \name Types
  /// @{

  /// \cond SKIP_IN_MANUAL
  using Face_graph = PolygonMesh;
  using Vertex_to_point_map = VertexToPointMap;

  using face_descriptor = typename boost::graph_traits<PolygonMesh>::face_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
  using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
  /// \endcond

  /// Number type.
  using FT = typename GeomTraits::FT;

  /// Item type.
  using Item = face_descriptor;
  using Region = std::vector<Item>;

  /// Primitive
  using Primitive = typename GeomTraits::Plane_3;

  /// Region map
  using Region_index_map = typename boost::property_map<Face_graph, CGAL::dynamic_face_property_t<std::size_t> >::const_type;

  /// @}

private:
  using Point_3 = typename GeomTraits::Point_3;
  using Vector_3 = typename GeomTraits::Vector_3;
  using Plane_3 = typename GeomTraits::Plane_3;
  using Triangle_3 = typename GeomTraits::Triangle_3;

  using Squared_length_3 = typename GeomTraits::Compute_squared_length_3;
  using Squared_distance_3 = typename GeomTraits::Compute_squared_distance_3;
  using Scalar_product_3 = typename GeomTraits::Compute_scalar_product_3;
  using Cross_product_3 = typename GeomTraits::Construct_cross_product_vector_3;

public:
  /// \name Initialization
  /// @{

  /*!
    \brief initializes all internal data structures.

    \tparam NamedParameters
    a sequence of \ref bgl_namedparameters "Named Parameters"

    \param pmesh
    an instance of `PolygonMesh` that represents a polygon mesh

    \param np
    a sequence of \ref bgl_namedparameters "Named Parameters"
    among the ones listed below

    \cgalNamedParamsBegin
      \cgalParamNBegin{maximum_distance}
        \cgalParamDescription{the maximum distance from a point to the plane of the primitive}
        \cgalParamType{`GeomTraits::FT`}
        \cgalParamDefault{1}
      \cgalParamNEnd
      \cgalParamNBegin{maximum_angle}
        \cgalParamDescription{the maximum angle in degrees between
        the normal of a face and the normal of the plane of the primitive}
        \cgalParamType{`GeomTraits::FT`}
        \cgalParamDefault{25 degrees}
      \cgalParamNEnd
      \cgalParamNBegin{cosine_of_maximum_angle}
        \cgalParamDescription{the cosine value `cos(maximum_angle * PI / 180)`, to be used instead of the parameter `maximum_angle()`}
        \cgalParamType{`GeomTraits::FT`}
        \cgalParamDefault{`cos(25 * PI / 180)`}
      \cgalParamNEnd
      \cgalParamNBegin{minimum_region_size}
        \cgalParamDescription{the minimum number of faces a region must have}
        \cgalParamType{`std::size_t`}
        \cgalParamDefault{1}
      \cgalParamNEnd
      \cgalParamNBegin{vertex_point_map}
        \cgalParamDescription{an instance of `VertexToPointMap` that maps a polygon mesh
        vertex (`boost::graph_traits<PolygonMesh>::vertex_descriptor`) to `GeomTraits::Point_3`}
        \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
      \cgalParamNEnd
      \cgalParamNBegin{face_normal_map}
        \cgalParamDescription{a property map associating normal vectors to the faces of `pmesh`.}
        \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%face_descriptor`
                       as key type and `GeomTraits::Vector_3` as value type.}
        \cgalParamDefault{If this parameter is omitted, face normals will be estimated using crossproducts of vectors created
                          from consecutive vertices of the face.}
      \cgalParamNEnd
      \cgalParamNBegin{geom_traits}
        \cgalParamDescription{an instance of `GeomTraits`}
        \cgalParamDefault{`GeomTraits()`}
      \cgalParamNEnd
    \cgalNamedParamsEnd

    \pre `faces(tmesh).size() > 0`
    \pre `maximum_distance >= 0`
    \pre `maximum_angle >= 0 && maximum_angle <= 90`
    \pre `cosine_of_maximum_angle >= 0 && cosine_of_maximum_angle <= 1`
    \pre `minimum_region_size > 0`
  */
  template<typename CGAL_NP_TEMPLATE_PARAMETERS>
  Plane_face_region(
    const PolygonMesh& pmesh,
    const CGAL_NP_CLASS& np = parameters::default_values()) :
  m_face_graph(pmesh),
  m_vpm(parameters::choose_parameter(parameters::get_parameter(
    np, internal_np::vertex_point), get_const_property_map(CGAL::vertex_point, pmesh))),
  m_traits(parameters::choose_parameter<GeomTraits>(parameters::get_parameter(np, internal_np::geom_traits))),
  m_squared_length_3(m_traits.compute_squared_length_3_object()),
  m_squared_distance_3(m_traits.compute_squared_distance_3_object()),
  m_scalar_product_3(m_traits.compute_scalar_product_3_object()),
  m_cross_product_3(m_traits.construct_cross_product_vector_3_object()),
  m_face_normals(get(CGAL::dynamic_face_property_t<Vector_3>(), pmesh)),
  m_face_triangulations( get(CGAL::dynamic_face_property_t<std::vector<Triangle_3>>(), pmesh) ) {
    static constexpr bool use_input_face_normal =
      !parameters::is_default_parameter<CGAL_NP_CLASS, internal_np::face_normal_t>::value;

#ifdef CGAL_SD_RG_USE_PMP
    auto get_face_normal = [this](Item face, const PolygonMesh& pmesh)
    {
      return Polygon_mesh_processing::compute_face_normal(face, pmesh, parameters::vertex_point_map(m_vpm));
    };
#else
    auto get_face_normal = [this](const Item &face, const PolygonMesh& pmesh) -> Vector_3
    {
      const auto hedge = halfedge(face, pmesh);
      const auto vertices = vertices_around_face(hedge, pmesh);
      CGAL_precondition(vertices.size() >= 3);

      auto vertex = vertices.begin();
      const boost::property_traits<Vertex_to_point_map>::reference p1 = get(m_vpm, *vertex); ++vertex;
      const boost::property_traits<Vertex_to_point_map>::reference p2 = get(m_vpm, *vertex); ++vertex;
      Point_3 p3 = get(m_vpm, *vertex);
      while(collinear(p1, p2, p3)) {
        if (++vertex == vertices.end()) return NULL_VECTOR;
        p3 = get(m_vpm, *vertex);
      }

      const Vector_3 u = p2 - p1;
      const Vector_3 v = p3 - p1;
      return m_cross_product_3(u, v);
    };
#endif

    if constexpr (!use_input_face_normal) {
      for (const Item &i : faces(pmesh))
        put(m_face_normals, i, get_face_normal(i, pmesh));
    }
    else {
      auto fnm = parameters::get_parameter(np, internal_np::face_normal);
      for (const Item &i : faces(pmesh))
        put(m_face_normals, i, get(fnm, i));
    }

    for (const Item &i : faces(pmesh)) {
      std::vector<Point_3> pts;
      auto h = halfedge(i, pmesh);
      auto s = h;

      do {
        pts.push_back(get(m_vpm, target(h, pmesh)));
        h = next(h, pmesh);
      } while (h != s);

      std::vector<Triangle_3> face_triangulation;
      internal::triangulate_face<GeomTraits>(pts, face_triangulation);
      put(m_face_triangulations, i, face_triangulation);
    }

    CGAL_precondition(faces(m_face_graph).size() > 0);
    const FT max_distance = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::maximum_distance), FT(1));
    CGAL_precondition(max_distance >= FT(0));
    m_distance_threshold = max_distance;

    const FT max_angle = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::maximum_angle), FT(25));
    CGAL_precondition(max_angle >= FT(0) && max_angle <= FT(90));

    m_min_region_size = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::minimum_region_size), 1);
    CGAL_precondition(m_min_region_size > 0);

    const FT default_cos_value = static_cast<FT>(std::cos(CGAL::to_double(
      (max_angle * static_cast<FT>(CGAL_PI)) / FT(180))));
    const FT cos_value = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::cosine_of_maximum_angle), default_cos_value);
    CGAL_precondition(cos_value >= FT(0) && cos_value <= FT(1));
    m_cos_value_threshold = cos_value;
  }

  /// @}

  /// \name Access
  /// @{

  /*!
    \brief implements `RegionType::region_index_map()`.

    This function creates an empty property map that maps each face to a `std::size_t`.
  */
  Region_index_map region_index_map() {
    return get(CGAL::dynamic_face_property_t<std::size_t>(), m_face_graph);
  }

  /*!
    \brief implements `RegionType::primitive()`.

    This function provides the support plane of the seed face.

    \return Primitive parameters that fits the region.

    \pre is_valid_region(region)
    \pre update(region)
  */
  Primitive primitive() const {
    return m_plane;
  }

  /*!
    \brief implements `RegionType::is_part_of_region()`.

    This function determines if the face `query` is within
    the `maximum_distance` from the corresponding plane and if the angle between
    its normal and the plane's normal is below the `maximum_angle`. If both conditions
    are satisfied, it returns `true`, otherwise `false`.

    \param query item of the query face

    \pre `query` is a valid face_descriptor (`boost::graph_traits<PolygonMesh>::%face_descriptor`) of `pmesh`
  */
  bool is_part_of_region(
    const Item query,
    const Region&) const
  {
    if (m_cos_value_threshold==1 || m_distance_threshold == 0)
    {
      halfedge_descriptor h = halfedge(query, m_face_graph);
      for (vertex_descriptor v : vertices_around_face(h, m_face_graph))
      {
        if (!coplanar(m_p, m_q, m_r, get(m_vpm, v)))
          return false;
      }
      return true;
    }
    else
    {
      // test on distance of points to the plane of the seed face
      const FT squared_distance_threshold = m_distance_threshold * m_distance_threshold;
      halfedge_descriptor h = halfedge(query, m_face_graph);
      for (vertex_descriptor v : vertices_around_face(h, m_face_graph))
      {
        //TODO: that's a bit dummy that we retest points that are already in the region...
        //      not sure caching in a vpm does worth it (need reset for each region)
        if (typename GeomTraits::Compare_squared_distance_3()(m_p, m_q, m_r,get(m_vpm, v), squared_distance_threshold) != SMALLER)
          return false;
      }

      if (m_cos_value_threshold == 1)
        return true;

      const typename GeomTraits::Point_3& p2=get(m_vpm,source(h, m_face_graph));
      const typename GeomTraits::Point_3& q2=get(m_vpm,target(h, m_face_graph));
      typename GeomTraits::Point_3 r2;

      halfedge_descriptor guard = prev(h, m_face_graph);
      do{
        h=next(h, m_face_graph);
        if (h == guard) return true;
        r2=get(m_vpm,target(h, m_face_graph));
      }
      while(collinear(p2,q2,r2));

      // test on the normal of the query face to the normal of the seed face
      return typename GeomTraits::Compare_angle_3()(m_p,m_q,m_r,
                                                    p2,q2,r2,
                                                    m_cos_value_threshold) == SMALLER;
    }
  }

  /*!
    \brief implements `RegionType::is_valid_region()`.

    This function controls if the `region` contains at least `minimum_region_size` faces.

    \param region
    Faces of the region represented as `Items`.
  */
  inline bool is_valid_region(const Region& region) const {
    return (region.size() >= m_min_region_size);
  }

  /*!
    \brief implements `RegionType::update()`.

    This function uses the support plane of the seed face as primitive.

    \param region
    Faces of the region represented as `Items`.

    \return Boolean `true` if the seed face is not degenerated and `false` otherwise

    \pre `region.size() > 0`
  */
  bool update(const Region& region) {

    CGAL_precondition(region.size() > 0);
    if (region.size() == 1) { // init reference plane and normal
      m_seed_face = region[0];

      halfedge_descriptor h = halfedge(m_seed_face, m_face_graph);

      //safety check for degenerate faces
      halfedge_descriptor guard = prev(h, m_face_graph);
      m_p = get(m_vpm, source(h, m_face_graph));
      m_q = get(m_vpm, target(h, m_face_graph));

      do {
        h = next(h, m_face_graph);
        // If all vertices are collinear, the face is degenerate and not suitable as seed
        if (h == guard) return false;
        m_r = get(m_vpm, target(h, m_face_graph));
      } while (collinear(m_p, m_q, m_r));

      const Vector_3 face_normal = get(m_face_normals, m_seed_face);
      if (face_normal == CGAL::NULL_VECTOR) return false;

      CGAL_precondition(face_normal != CGAL::NULL_VECTOR);
      m_plane = Plane_3(m_p, face_normal);
      m_normal = face_normal;
    }

    return true;
  }

  /// @}

private:
  const Face_graph& m_face_graph;
  const Vertex_to_point_map m_vpm;
  GeomTraits m_traits;

  FT m_distance_threshold;
  FT m_cos_value_threshold;
  std::size_t m_min_region_size;

  const Squared_length_3 m_squared_length_3;
  const Squared_distance_3 m_squared_distance_3;
  const Scalar_product_3 m_scalar_product_3;
  const Cross_product_3 m_cross_product_3;

  typename boost::property_map<Face_graph, CGAL::dynamic_face_property_t<Vector_3> >::const_type m_face_normals;
  typename boost::property_map<Face_graph, CGAL::dynamic_face_property_t<std::vector<Triangle_3>> >::const_type m_face_triangulations;

  Plane_3 m_plane;
  Vector_3 m_normal;
  face_descriptor m_seed_face;
  Point_3 m_p, m_q, m_r; // Three non-collinear points of the m_seed_face

  // Compute centroid of the face.
  template<typename Face>
  Point_3 get_face_centroid(const Face& face) const {

    const auto hedge = halfedge(face, m_face_graph);
    const auto vertices = vertices_around_face(hedge, m_face_graph);
    CGAL_precondition(vertices.size() > 0);

    FT sum = FT(0), x = FT(0), y = FT(0), z = FT(0);
    for (const auto vertex : vertices) {
      const Point_3& point = get(m_vpm, vertex);
      x += point.x();
      y += point.y();
      z += point.z();
      sum += FT(1);
    }
    CGAL_precondition(sum > FT(0));
    x /= sum;
    y /= sum;
    z /= sum;
    return Point_3(x, y, z);
  }
};

} // namespace Polygon_mesh
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_PLANE_FACE_REGION_H
