#ifndef CGAL_IO_GLTF_H
#define CGAL_IO_GLTF_H

#include <CGAL/IO/io.h>
#include <CGAL/IO/helpers.h>
#include <CGAL/Container_helper.h>
#include <boost/range/value_type.hpp>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <iostream>
#include <string>
#include <vector>

// You must have tinygltf.h in your include path
#include "tiny_gltf.h"

namespace CGAL {
namespace IO {
namespace internal {

// The core implementation function
template <typename PointRange, typename PolygonRange>
bool read_GLTF(const std::string& filename,
               PointRange& points,
               PolygonRange& polygons,
               const bool verbose = false)
{
    // Use the Point type from the PointRange container
    typedef typename boost::range_value<PointRange>::type Point;

    tinygltf::Model model;
    tinygltf::TinyGLTF loader;
    std::string err;
    std::string warn;

    bool ret = loader.LoadASCIIFromFile(&model, &err, &warn, filename);

    if (!warn.empty() && verbose) {
        std::cout << "GLTF warning: " << warn << std::endl;
    }
    if (!err.empty()) {
        if (verbose) std::cerr << "GLTF error: " << err << std::endl;
        // try loading as binary glTF (.glb) as a fallback
        ret = loader.LoadBinaryFromFile(&model, &err, &warn, filename);
        if (!ret) {
          if (verbose) std::cerr << "Failed to load GLTF/GLB file." << std::endl;
          return false;
        }
    }
    if (!ret) {
        if (verbose) std::cerr << "Failed to parse GLTF file: " << filename << std::endl;
        return false;
    }

    // Iterate over all meshes in the model
    for (const auto& mesh : model.meshes) {
        // Iterate over all primitives in the mesh (a primitive is a single draw call with one material)
        for (const auto& primitive : mesh.primitives) {
            // We are only interested in triangles for creating simple meshes
            if (primitive.mode != TINYGLTF_MODE_TRIANGLES) {
                continue;
            }

            // --- 1. Get VERTICES ---
            const float* position_buffer = nullptr;
            size_t vertex_count = 0;
            size_t position_stride = 0;

            // Find the "POSITION" attribute
            if (primitive.attributes.find("POSITION") != primitive.attributes.end()) {
                const tinygltf::Accessor& accessor = model.accessors[primitive.attributes.find("POSITION")->second];
                const tinygltf::BufferView& bufferView = model.bufferViews[accessor.bufferView];
                const tinygltf::Buffer& buffer = model.buffers[bufferView.buffer];

                vertex_count = accessor.count;
                position_buffer = reinterpret_cast<const float*>(&buffer.data[bufferView.byteOffset + accessor.byteOffset]);

                // Get the stride (bytes between consecutive vertices)
                position_stride = accessor.ByteStride(bufferView);
            } else {
                if (verbose) std::cerr << "Primitive has no POSITION attribute." << std::endl;
                continue; // Skip this primitive if it has no vertex positions
            }

            // Remember the number of points we had before adding new ones.
            // This is crucial for correctly indexing the faces.
            size_t vertex_offset = points.size();

            // Add the vertices from this primitive to our global points container
            for (size_t i = 0; i < vertex_count; ++i) {
                // Assuming points are tightly packed VEC3 of float
                if (position_stride == 0 || position_stride == sizeof(float) * 3) {
                     points.push_back(Point(position_buffer[i * 3 + 0],
                                           position_buffer[i * 3 + 1],
                                           position_buffer[i * 3 + 2]));
                } else { // Handle strided data
                    const unsigned char* byte_ptr = reinterpret_cast<const unsigned char*>(position_buffer) + i * position_stride;
                    const float* float_ptr = reinterpret_cast<const float*>(byte_ptr);
                     points.push_back(Point(float_ptr[0], float_ptr[1], float_ptr[2]));
                }
            }

            // --- 2. Get INDICES (Faces) ---
            if (primitive.indices >= 0) {
                const tinygltf::Accessor& accessor = model.accessors[primitive.indices];
                const tinygltf::BufferView& bufferView = model.bufferViews[accessor.bufferView];
                const tinygltf::Buffer& buffer = model.buffers[bufferView.buffer];

                const unsigned char* base_ptr = &buffer.data[bufferView.byteOffset + accessor.byteOffset];
                size_t index_count = accessor.count;

                // glTF indices can be unsigned short or unsigned int. We need to handle both.
                switch (accessor.componentType) {
                    case TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT: {
                        const unsigned short* index_buffer = reinterpret_cast<const unsigned short*>(base_ptr);
                        for (size_t i = 0; i < index_count; i += 3) {
                            polygons.emplace_back();
                            auto& current_polygon = polygons.back();
                            // Add the vertex_offset to get the correct global index
                            current_polygon.push_back(index_buffer[i + 0] + vertex_offset);
                            current_polygon.push_back(index_buffer[i + 1] + vertex_offset);
                            current_polygon.push_back(index_buffer[i + 2] + vertex_offset);
                        }
                        break;
                    }
                    case TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT: {
                        const unsigned int* index_buffer = reinterpret_cast<const unsigned int*>(base_ptr);
                        for (size_t i = 0; i < index_count; i += 3) {
                            polygons.emplace_back();
                            auto& current_polygon = polygons.back();
                            current_polygon.push_back(index_buffer[i + 0] + vertex_offset);
                            current_polygon.push_back(index_buffer[i + 1] + vertex_offset);
                            current_polygon.push_back(index_buffer[i + 2] + vertex_offset);
                        }
                        break;
                    }
                    default:
                        if (verbose) std::cerr << "Unsupported index component type." << std::endl;
                        break;
                }
            }
        }
    }

    return true;
}

} // namespace internal


// Public-facing function, similar to read_OBJ
template <typename PointRange, typename PolygonRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_GLTF(const std::string& fname,
               PointRange& points,
               PolygonRange& polygons,
               const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
               , std::enable_if_t<internal::is_Range<PolygonRange>::value>* = nullptr
#endif
              )
{
    const bool verbose = parameters::choose_parameter(parameters::get_parameter(np, internal_np::verbose), false);
    // Note: We don't clear the containers, we append, to match read_OBJ behavior.
    return internal::read_GLTF(fname, points, polygons, verbose);
}

} // namespace IO
} // namespace CGAL

#endif // CGAL_IO_GLTF_H
