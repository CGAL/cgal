#ifndef CGAL_IO_READ_GLTF_H
#define CGAL_IO_READ_GLTF_H

#include <CGAL/IO/GLTF/read_gltf.h>
#include <string>

namespace CGAL {
namespace IO {

// Public entry point for polygon soup (points and polygons)
template <typename PointRange, typename PolygonRange>
bool read_GLTF(const std::string& filename,
               PointRange& points,
               PolygonRange& polygons) {
    return internal::read_GLTF(filename, points, polygons);
}

// Public entry point that wraps your internal implementation
template <typename PolygonMesh>
bool read_GLTF(const std::string& filename, PolygonMesh& mesh) {
    return internal::read_GLTF(filename, mesh);
}

} // namespace IO
} // namespace CGAL

#endif
