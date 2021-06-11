#ifndef CGAL_BARYCENTRIC_INTERNAL_UTILS_3_H
#define CGAL_BARYCENTRIC_INTERNAL_UTILS_3_H

// STL includes 
#include <tuple>

// Internal includes
#include <CGAL/Barycentric_coordinates_3/barycentric_enum_3.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>

namespace CGAL{
namespace Barycentric_coordinates{
namespace internal{

// Get default values.
  template<typename OutputIterator>
  void get_default(
    const std::size_t n, OutputIterator output) {

    for (std::size_t i = 0; i < n; ++i) {
      *(output++) = 0;
    }
  }

// Compute barycentric coordinates in the space.
  template<
  typename OutputIterator,
  typename GeomTraits>
  OutputIterator planar_coordinates_3(
    const typename GeomTraits::Point_3& p0,
    const typename GeomTraits::Point_3& p1,
    const typename GeomTraits::Point_3& p2,
    const typename GeomTraits::Point_3& p3,
    const typename GeomTraits::Point_3& query,
    OutputIterator coordinates,
    const GeomTraits& traits) {

    // Number type.
    using FT = typename GeomTraits::FT;

    // Functions.
    const auto volume_3 = traits.compute_volume_3_object();
    const FT total_volume = volume_3(p0, p1, p2, p3);

    CGAL_precondition(total_volume != FT(0));
    if (total_volume == FT(0)) {
      get_default(4, coordinates);
      return coordinates;
    }

    // Compute some related sub-volumes.
    const FT V1 = volume_3(p1, p3, p2, query);
    const FT V2 = volume_3(p2, p3, p0, query);
    const FT V3 = volume_3(p3, p1, p0, query);

    // Compute the inverted total volume of the tetrahedron.
    CGAL_assertion(total_volume != FT(0));
    const FT inverted_total_volume = FT(1) / total_volume;

    // Compute coordinates.
    const FT b0 = V1 * inverted_total_volume;
    const FT b1 = V2 * inverted_total_volume;
    const FT b2 = V3 * inverted_total_volume;
    const FT b3 = FT(1) - b0 - b1 - b2;

    // Return coordinates.
    *(coordinates++) = b0;
    *(coordinates++) = b1;
    *(coordinates++) = b2;
    *(coordinates++) = b3;

    return coordinates;
  }

  // Compute wp coordinates for a given vertex v and a query q
  template<
  typename GeomTraits, 
  typename PolygonMesh,
  typename VertexToPointMap, 
  typename FaceCirculator>
  typename GeomTraits::FT calculate_wp_vertex_query(
  const typename GeomTraits::Point_3& vertex,
  const typename GeomTraits::Point_3& query,
  const PolygonMesh& mesh,
  const VertexToPointMap& vertex_to_point_map,
  FaceCirculator& face_circulator,
  const GeomTraits& traits){

    using FT = typename GeomTraits::FT;
    using Point_3 = typename GeomTraits::Point_3;
    using Polygon_mesh = PolygonMesh;

    const auto det3 = traits.compute_determinant_3_object();
    const auto dot3 = traits.compute_scalar_product_3_object();

    // Vertices of triangle
    Point_3 p0;
    Point_3 p1;
    Point_3 p2;

    // Use face circulator 
    CGAL::Face_around_target_circulator<Polygon_mesh> done(face_circulator);
    do{

      //Figure out how to pick every vertex ofr one face;
      auto vertices_of_face = mesh.vertex(mesh.edge(mesh.halfedge(*face_circulator)), 0);

      std::cout << p0 << std::endl;
    
      face_circulator++;
    }while(face_circulator!=done);
    

    // Vector between query and vertex
    typename GeomTraits::Vector_3 query_vertex(query, vertex);

    return FT(0);
  }

}
}
}

#endif
