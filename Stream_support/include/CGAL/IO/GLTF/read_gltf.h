// Copyright (c) 2025
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Dema Nasrawt and Aidan Pawlak

#ifndef CGAL_IO_GLTF_READ_GLTF_H
#define CGAL_IO_GLTF_READ_GLTF_H

#include <CGAL/IO/io.h>
#include <CGAL/IO/helpers.h>
#include <CGAL/Container_helper.h>
#include <boost/range/value_type.hpp>
#include <CGAL/Named_function_parameters.h>

#include <iostream>
#include <string>
#include <vector>
#include <limits>

// Include the tinygltf single-header implementation exactly once.
// The implementation macros are undefined afterwards so that a second include
// of tiny_gltf.h (e.g. through write_gltf.h) skips the implementation section
// and avoids redefinition errors with the stb libraries.
#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "tiny_gltf.h"
#undef TINYGLTF_IMPLEMENTATION
#undef STB_IMAGE_IMPLEMENTATION
#undef STB_IMAGE_WRITE_IMPLEMENTATION

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
        // Iterate over all primitives in the mesh
        for (const auto& primitive : mesh.primitives) {
            // We are only interested in triangles for creating simple meshes
            if (primitive.mode != TINYGLTF_MODE_TRIANGLES) {
                if (verbose) {
                    std::cout << "Skipping non-triangle primitive (mode: "
                              << primitive.mode << ")" << std::endl;
                }
                continue;
            }

            // --- 1. Get VERTICES ---
            size_t vertex_count = 0;
            size_t position_stride = 0;

            // Find the "POSITION" attribute
            auto pos_it = primitive.attributes.find("POSITION");
            if (pos_it == primitive.attributes.end()) {
                if (verbose) std::cerr << "Primitive has no POSITION attribute." << std::endl;
                continue;
            }

            const tinygltf::Accessor& accessor = model.accessors[pos_it->second];

            // Validate accessor
            if (accessor.type != TINYGLTF_TYPE_VEC3) {
                if (verbose) std::cerr << "POSITION accessor must be VEC3." << std::endl;
                continue;
            }

            if (accessor.bufferView < 0 || size_t(accessor.bufferView) >= model.bufferViews.size()) {
                if (verbose) std::cerr << "Invalid or missing bufferView for POSITION." << std::endl;
                continue;
            }

            const tinygltf::BufferView& bufferView = model.bufferViews[accessor.bufferView];
            if (size_t(bufferView.buffer) >= model.buffers.size()) {
                if (verbose) std::cerr << "Invalid buffer index." << std::endl;
                continue;
            }

            const tinygltf::Buffer& buffer = model.buffers[bufferView.buffer];
            vertex_count = accessor.count;
            size_t total_offset = bufferView.byteOffset + accessor.byteOffset;

            // Get the stride (bytes between consecutive vertices)
            position_stride = accessor.ByteStride(bufferView);
            
            // Validate and default stride based on component type
            if (position_stride == 0) {
                if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_DOUBLE) {
                    position_stride = sizeof(double) * 3;
                } else {
                    position_stride = sizeof(float) * 3;
                }
            }

            // Pointer to the start of the data
            const unsigned char* base_ptr = &buffer.data[total_offset];
            size_t vertex_offset = points.size();

            // Add the vertices from this primitive
            for (size_t i = 0; i < vertex_count; ++i) {
                const unsigned char* element_ptr = base_ptr + (i * position_stride);
                
                // Bounds check
                if (element_ptr + (accessor.componentType == TINYGLTF_COMPONENT_TYPE_DOUBLE ? sizeof(double)*3 : sizeof(float)*3) > buffer.data.data() + buffer.data.size()) {
                    if (verbose) std::cerr << "Vertex data out of bounds." << std::endl;
                    return false;
                }

                if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT) {
                    const float* f_ptr = reinterpret_cast<const float*>(element_ptr);
                    points.push_back(Point(static_cast<double>(f_ptr[0]), 
                                          static_cast<double>(f_ptr[1]), 
                                          static_cast<double>(f_ptr[2])));
                } 
                else if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_DOUBLE) {
                    // This path ensures bit-perfect identity for high-precision data
                    const double* d_ptr = reinterpret_cast<const double*>(element_ptr);
                    points.push_back(Point(d_ptr[0], d_ptr[1], d_ptr[2]));
                }
                else {
                    if (verbose) std::cerr << "Unsupported component type for POSITION." << std::endl;
                    return false;
                }
            }

            // --- 2. Get INDICES (Faces) ---
            if (primitive.indices >= 0) {
                if (size_t(primitive.indices) >= model.accessors.size()) {
                    if (verbose) std::cerr << "Invalid indices accessor index." << std::endl;
                    continue;
                }

                const tinygltf::Accessor& idx_accessor = model.accessors[primitive.indices];

                if (idx_accessor.bufferView < 0) {
                    if (verbose) std::cerr << "Indices accessor has no bufferView." << std::endl;
                    continue;
                }

                if (size_t(idx_accessor.bufferView) >= model.bufferViews.size()) {
                    if (verbose) std::cerr << "Invalid indices bufferView." << std::endl;
                    continue;
                }

                const tinygltf::BufferView& idx_bufferView = model.bufferViews[idx_accessor.bufferView];

                if (size_t(idx_bufferView.buffer) >= model.buffers.size()) {
                    if (verbose) std::cerr << "Invalid indices buffer." << std::endl;
                    continue;
                }

                const tinygltf::Buffer& idx_buffer = model.buffers[idx_bufferView.buffer];

                size_t idx_total_offset = idx_bufferView.byteOffset + idx_accessor.byteOffset;
                const unsigned char* base_ptr = &idx_buffer.data[idx_total_offset];
                size_t index_count = idx_accessor.count;

                // Validate index count is multiple of 3 for triangles
                if (index_count % 3 != 0) {
                    if (verbose) std::cerr << "Index count not multiple of 3 for triangles." << std::endl;
                    continue;
                }

                // glTF indices can be unsigned byte, short, or int
                switch (idx_accessor.componentType) {
                    case TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE: {
                        const unsigned char* index_buffer = base_ptr;
                        for (size_t i = 0; i < index_count; i += 3) {
                            polygons.emplace_back();
                            auto& current_polygon = polygons.back();
                            current_polygon.push_back(index_buffer[i + 0] + vertex_offset);
                            current_polygon.push_back(index_buffer[i + 1] + vertex_offset);
                            current_polygon.push_back(index_buffer[i + 2] + vertex_offset);
                        }
                        break;
                    }
                    case TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT: {
                        const unsigned short* index_buffer = reinterpret_cast<const unsigned short*>(base_ptr);
                        for (size_t i = 0; i < index_count; i += 3) {
                            polygons.emplace_back();
                            auto& current_polygon = polygons.back();
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
                            // Check for overflow when adding offset
                            if (index_buffer[i + 0] > std::numeric_limits<size_t>::max() - vertex_offset ||
                                index_buffer[i + 1] > std::numeric_limits<size_t>::max() - vertex_offset ||
                                index_buffer[i + 2] > std::numeric_limits<size_t>::max() - vertex_offset) {
                                if (verbose) std::cerr << "Index overflow detected." << std::endl;
                                return false;
                            }
                            current_polygon.push_back(index_buffer[i + 0] + vertex_offset);
                            current_polygon.push_back(index_buffer[i + 1] + vertex_offset);
                            current_polygon.push_back(index_buffer[i + 2] + vertex_offset);
                        }
                        break;
                    }
                    default:
                        if (verbose) std::cerr << "Unsupported index component type: "
                                               << idx_accessor.componentType << std::endl;
                        continue;
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

#endif // CGAL_IO_GLTF_READ_GLTF_H