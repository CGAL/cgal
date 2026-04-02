// Author(s)     : Dema Nasrawt and Aidan Pawlak

#ifndef CGAL_IO_GLTF_WRITE_GLTF_H
#define CGAL_IO_GLTF_WRITE_GLTF_H

#include "tiny_gltf.h"

#include <CGAL/IO/helpers.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/number_utils.h>

#include <algorithm>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace CGAL {
namespace IO {
namespace internal {

template <typename PointRange, typename PolygonRange>
bool write_GLTF(const std::string& filename,
                const PointRange&  points,
                const PolygonRange& polygons,
                const bool verbose = false)
{
  const std::size_t vertex_count = points.size();

  if (vertex_count == 0 && verbose)
    std::cerr << "Warning: writing GLTF with no vertices." << std::endl;

  std::vector<unsigned char> buffer_data;

  // --- Vertex positions: VEC3 of DOUBLE ---
  const std::size_t vertex_byte_offset = 0;

  double min_x =  std::numeric_limits<double>::max();
  double min_y =  std::numeric_limits<double>::max();
  double min_z =  std::numeric_limits<double>::max();
  double max_x = -std::numeric_limits<double>::max();
  double max_y = -std::numeric_limits<double>::max();
  double max_z = -std::numeric_limits<double>::max();

  for (const auto& p : points)
  {
    double x = CGAL::to_double(p.x());
    double y = CGAL::to_double(p.y());
    double z = CGAL::to_double(p.z());

    min_x = std::min(min_x, x);  max_x = std::max(max_x, x);
    min_y = std::min(min_y, y);  max_y = std::max(max_y, y);
    min_z = std::min(min_z, z);  max_z = std::max(max_z, z);

    unsigned char bytes[sizeof(double)];
    std::memcpy(bytes, &x, sizeof(double));
    buffer_data.insert(buffer_data.end(), bytes, bytes + sizeof(double));
    std::memcpy(bytes, &y, sizeof(double));
    buffer_data.insert(buffer_data.end(), bytes, bytes + sizeof(double));
    std::memcpy(bytes, &z, sizeof(double));
    buffer_data.insert(buffer_data.end(), bytes, bytes + sizeof(double));
  }

  const std::size_t vertex_byte_length = buffer_data.size();

  // Pad to 4-byte boundary before index data (glTF requirement).
  while (buffer_data.size() % 4 != 0)
    buffer_data.push_back(0x00);

  // --- Triangle indices: SCALAR of UNSIGNED_INT ---
  const std::size_t index_byte_offset = buffer_data.size();
  std::size_t index_count = 0;

  for (const auto& poly : polygons)
  {
    if (poly.size() != 3)
    {
      if (verbose)
        std::cerr << "Warning: skipping non-triangle face (size="
                  << poly.size() << ")." << std::endl;
      continue;
    }

    for (const auto& idx : poly)
    {
      unsigned int index = static_cast<unsigned int>(idx);
      unsigned char bytes[sizeof(unsigned int)];
      std::memcpy(bytes, &index, sizeof(unsigned int));
      buffer_data.insert(buffer_data.end(), bytes, bytes + sizeof(unsigned int));
    }
    index_count += 3;
  }

  const std::size_t index_byte_length = buffer_data.size() - index_byte_offset;

  // --- Assemble the glTF model ---
  tinygltf::Model model;
  tinygltf::TinyGLTF writer;

  tinygltf::Buffer buf;
  buf.data = buffer_data;
  model.buffers.push_back(buf);

  tinygltf::BufferView pos_bv;
  pos_bv.buffer     = 0;
  pos_bv.byteOffset = vertex_byte_offset;
  pos_bv.byteLength = vertex_byte_length;
  pos_bv.target     = TINYGLTF_TARGET_ARRAY_BUFFER;
  model.bufferViews.push_back(pos_bv);   // index 0

  tinygltf::BufferView idx_bv;
  idx_bv.buffer     = 0;
  idx_bv.byteOffset = index_byte_offset;
  idx_bv.byteLength = index_byte_length;
  idx_bv.target     = TINYGLTF_TARGET_ELEMENT_ARRAY_BUFFER;
  model.bufferViews.push_back(idx_bv);   // index 1

  tinygltf::Accessor pos_acc;
  pos_acc.bufferView    = 0;
  pos_acc.byteOffset    = 0;
  pos_acc.componentType = TINYGLTF_COMPONENT_TYPE_DOUBLE;
  pos_acc.count         = vertex_count;
  pos_acc.type          = TINYGLTF_TYPE_VEC3;
  if (vertex_count > 0)
  {
    pos_acc.minValues = { min_x, min_y, min_z };
    pos_acc.maxValues = { max_x, max_y, max_z };
  }
  model.accessors.push_back(pos_acc);    // index 0

  tinygltf::Accessor idx_acc;
  idx_acc.bufferView    = 1;
  idx_acc.byteOffset    = 0;
  idx_acc.componentType = TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT;
  idx_acc.count         = index_count;
  idx_acc.type          = TINYGLTF_TYPE_SCALAR;
  model.accessors.push_back(idx_acc);    // index 1

  tinygltf::Primitive prim;
  prim.attributes["POSITION"] = 0;
  prim.indices                = 1;
  prim.mode                   = TINYGLTF_MODE_TRIANGLES;

  tinygltf::Mesh mesh;
  mesh.name = "mesh";
  mesh.primitives.push_back(prim);
  model.meshes.push_back(mesh);

  tinygltf::Node node;
  node.mesh = 0;
  model.nodes.push_back(node);

  tinygltf::Scene scene;
  scene.nodes.push_back(0);
  model.scenes.push_back(scene);
  model.defaultScene = 0;

  model.asset.version   = "2.0";
  model.asset.generator = "CGAL";

  const bool ok = writer.WriteGltfSceneToFile(
      &model, filename,
      /*embed_images=*/  true,
      /*embed_buffers=*/ true,
      /*pretty_print=*/  true,
      /*write_binary=*/  false);

  if (!ok && verbose)
    std::cerr << "Failed to write GLTF file: " << filename << std::endl;

  return ok;
}

} // namespace internal


/// \ingroup PkgStreamSupportIoFuncsGLTF
/// \brief writes `points` and `polygons` to `fname` as an ASCII .gltf file.
///
/// Only triangle faces are written; non-triangles are skipped with an optional
/// warning.  Positions are stored as 64-bit doubles for lossless round-tripping
/// with `read_GLTF`.
///
/// \returns `true` on success.
template <typename PointRange, typename PolygonRange,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_GLTF(const std::string&  fname,
                const PointRange&   points,
                const PolygonRange& polygons,
                const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
                , std::enable_if_t<internal::is_Range<PolygonRange>::value>* = nullptr
#endif
               )
{
  const bool verbose = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::verbose), false);
  return internal::write_GLTF(fname, points, polygons, verbose);
}

} // namespace IO
} // namespace CGAL

#endif // CGAL_IO_GLTF_WRITE_GLTF_H
