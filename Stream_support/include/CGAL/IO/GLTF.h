#ifndef CGAL_IO_GLTF_H
#define CGAL_IO_GLTF_H

#include <CGAL/IO/GLTF/GLTF_reader.h>

namespace CGAL {
namespace IO {

// Public entry point that wraps your internal implementation
template <typename PolygonMesh>
bool read_GLTF(const std::string& filename, PolygonMesh& mesh) {
    return internal::read_GLTF(filename, mesh);
}

} // namespace IO
} // namespace CGAL

#endif
