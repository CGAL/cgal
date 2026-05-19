// Copyright (c) 2025 Universitaet Bremen
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Rene Weller, Sven Oesau
//

#ifndef CGAL_SPHERE_PACKING_H_
#define CGAL_SPHERE_PACKING_H_

#include <CGAL/license/Sphere_packing.h>

#include <CGAL/disable_warnings.h>
#include <CGAL/Named_function_parameters.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#pragma warning( push )
#pragma warning( disable : 4251 )
#include <ProtoSphereExtended.h>
#pragma warning( pop )

namespace CGAL {

namespace internal {

template<typename TriangleMesh, typename VertexPointMap>
class FaceGraphWrapper : public ProtoSphere::ITriMesh {
public:
  using uint4 = ProtoSphere::uint4;
  using float3 = ProtoSphere::float3;
  using float4 = ProtoSphere::float4;

  FaceGraphWrapper(const TriangleMesh& tm, const VertexPointMap vpm) : mesh(tm), vpm(vpm) {
    aabb_min = { (std::numeric_limits<float>::max)(), (std::numeric_limits<float>::max)(), (std::numeric_limits<float>::max)() };
    aabb_max = { -(std::numeric_limits<float>::max)(), -(std::numeric_limits<float>::max)(), -(std::numeric_limits<float>::max)() };

    for (const auto& v : vertices(mesh)) {
      auto p = get(vpm, v);
      aabb_min.x = (std::min<float>)(aabb_min.x, float(p.x()));
      aabb_min.y = (std::min<float>)(aabb_min.y, float(p.y()));
      aabb_min.z = (std::min<float>)(aabb_min.z, float(p.z()));
      aabb_max.x = (std::max<float>)(aabb_max.x, float(p.x()));
      aabb_max.y = (std::max<float>)(aabb_max.y, float(p.y()));
      aabb_max.z = (std::max<float>)(aabb_max.z, float(p.z()));
    }

    float3 center = { (aabb_min.x + aabb_max.x) / 2.0f, (aabb_min.y + aabb_max.y) / 2.0f, (aabb_min.z + aabb_max.z) / 2.0f };

    points.reserve(vertices(mesh).size());
    for (const auto& v : vertices(mesh)) {
      auto p = get(vpm, v);
      points.push_back(float4{float(p.x()), float(p.y()), float(p.z()), 1.0f});
    }
  }
  ~FaceGraphWrapper() {}

	virtual float getVolume()	const {
    return static_cast<float>(Polygon_mesh_processing::volume(mesh));
  }

	virtual float3 const getAABBMin() const { return aabb_min; }
	virtual float3 const getAABBMax() const { return aabb_max; }

	virtual unsigned int getNumVertices() const {
    return static_cast<unsigned int>(vertices(mesh).size());
  }

	virtual unsigned int getNumTriangles() const {
    return static_cast<unsigned int>(faces(mesh).size());
  }

	virtual std::vector<uint4> const& getIndices() const {
#pragma message("ToDo: Does this work with garbage collection?")
    if (indices.empty()) {
      indices.reserve(faces(mesh).size() * 3);
      for (const auto& f : faces(mesh)) {
        auto h = halfedge(f, mesh);
        indices.push_back(uint4{target(h, mesh), target(next(h, mesh), mesh),
          target(next(next(h, mesh), mesh), mesh), 0});
      }
    }
    return indices;
  }

	virtual std::vector<float4>	const& getVertices() const { return points; }

	virtual std::vector<uint4>	const& getNormalIndices()	const {
    return getIndices();
  }

  virtual std::vector<float4>	const& getNormals()	const {
    if (normals.empty()) {
      normals.reserve(faces(mesh).size());
      for (const auto& f : faces(mesh)) {
        auto n = Polygon_mesh_processing::compute_face_normal(f, mesh, CGAL::parameters::vertex_point_map(vpm));
        normals.push_back(float4{float(n.x()), float(n.y()), float(n.z()), 0.0f});
      }
    }

    return normals;
  }

  virtual void moveToCenter(const float cellSize) {
    std::cout << "moveToCenter(" << cellSize << ") called" << std::endl;
    if (!points.empty()) {
      for (float4 &v : points) {
        v.x += cellSize;
        v.y += cellSize;
        v.z += cellSize;
      }
    }
  }

private:
  float3 aabb_min, aabb_max;
  const TriangleMesh &mesh;
  const VertexPointMap vpm;
  mutable std::vector<uint4> indices;
  std::vector<float4> points;
  mutable std::vector<float4> normals;
};
} //namespace internal

/**
 * \ingroup PkgSpherePackingAlgorithms
 * \brief Packs spheres into a closed mesh.
 * Computes a packing of spheres into a closed and self-intersection free triangle mesh. The possible range of radii of spheres can be chosen; large spheres are preferred during the packing.
 * The optimal solution is not guaranteed as it has NP-hard complexity.
 *
 * Note that the precision of the method is limited by the precision of the GPU, i.e., non-exact arithmetic with 32-bit floating point numbers is used.
 *
 * \tparam TriangleMesh a model of `HalfedgeListGraph`, `FaceListGraph`, and `MutableFaceGraph`
 *
 * \tparam OutputIterator must be an output iterator accepting variables of type `geom_traits::Sphere_3`.
 *
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param tm the input tm to pack spheres into. It has to be a closed triangle mesh and may not have self-intersections.
 * \param out output iterator into which the packed spheres are written.
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *
 *   \cgalParamNBegin{number_of_spheres}
 *     \cgalParamDescription{target number of spheres to be packed into `tm`.}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{1,000}
 *     \cgalParamExtra{May be exceeded as one packing of iteration may insert many spheres.
 *         May also not be reached if `target_coverage` is attained first or no further spheres satisfying `minimum_radius` can be inserted.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{target_coverage}
 *     \cgalParamDescription{ratio of the volume of `tm` to be covered by spheres.}
 *     \cgalParamType{float}
 *     \cgalParamDefault{0.92}
 *     \cgalParamExtra{May be exceeded as one packing of iteration may insert many spheres.
 *         May also not be reached if `number_of_target_spheres` is attained first or no further spheres satisfying `minimum_radius` can be inserted.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{minimum_radius}
 *     \cgalParamDescription{minimum radius of packed spheres}
 *     \cgalParamType{float}
 *     \cgalParamDefault{0}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{initial_grid_resolution}
 *     \cgalParamDescription{number of grid cells for the longest side of the bounding box. Grid cells are cubic.}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{2}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{maximum_splits}
 *     \cgalParamDescription{how often the grid is subdivided at most.}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{3}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{iterations_between_splits}
 *     \cgalParamDescription{iterations of placing spheres between splitting the grid into a higher resolution.}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{30}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `TriangleMesh`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `Kernel`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type of `TriangleMesh`, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
 *
 * \cgalNamedParamsEnd
 *
 * \return `out`, the output iterator
 */


template<typename TriangleMesh, typename OutputIterator, typename NamedParameters = parameters::Default_named_parameters>
OutputIterator pack_spheres(const TriangleMesh& tm, OutputIterator out, const NamedParameters& np = parameters::default_values()) {
  using parameters::choose_parameter;
  using parameters::get_parameter;

  if (faces(tm).size() == 0)
    return out;

  using VPM = typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type;
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point), get_const_property_map(CGAL::vertex_point, tm));
  using Geom_traits = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;
  using Point_3 = typename Geom_traits::Point_3;
  using Sphere_3 = typename Geom_traits::Sphere_3;

  // named parameters
  const unsigned int number_of_spheres = parameters::choose_parameter(parameters::get_parameter(np, internal_np::number_of_spheres), 1000);
  const unsigned int initial_grid_resolution = parameters::choose_parameter(parameters::get_parameter(np, internal_np::initial_grid_resolution), 2);
  const unsigned int maximum_splits = parameters::choose_parameter(parameters::get_parameter(np, internal_np::maximum_splits), 3);
  const unsigned int iterations_between_splits = parameters::choose_parameter(parameters::get_parameter(np, internal_np::iterations_between_splits), 30);
  const float minimum_radius = parameters::choose_parameter(parameters::get_parameter(np, internal_np::minimum_radius), 0.0f);
  const float target_coverage = parameters::choose_parameter(parameters::get_parameter(np, internal_np::target_coverage), 0.92f);
  const bool verbose = parameters::choose_parameter(parameters::get_parameter(np, internal_np::verbose), false);

  ProtoSphere::Profiler::getProfilerPtr()->disable(true);

  internal::FaceGraphWrapper<TriangleMesh, VPM> wrapper(tm, vpm);

  if (!verbose) {
    ProtoSphere::Log::getLog("ProtoSphere", false).setSilentMode(true);
    ProtoSphere::Log::getLog("Profiler", false).setSilentMode(true);
  }

  if (!ProtoSphere::GPUWrapper::init()) {
    std::cout << "GPUWrapper initialization failed!\n";
    return out;
  }

#ifdef USE_OPENCL
  if (ProtoSphere::GPUWrapper::getSingletonPtr()->getType() == GPOpenCL)
    ProtoSphere::GPUWrapper::getSingletonPtr()->getRawOpenCL()->setKernelFolder("../oclKernels/");
#endif // USE_OPENCL

  if (verbose)
    ProtoSphere::Profiler::getProfilerPtr()->startProfiling();

  // Allocate memory for the parameters. The memory is freed by the GridManager.
  ProtoSphere::HybridGridParameters *params = new ProtoSphere::HybridGridParameters;
  params->setNumCellsOnLongestSide(initial_grid_resolution);
  params->setMaxSplits(maximum_splits);
  params->setMaxIterations(iterations_between_splits);
  params->setExpImpSplitGap(2);
  params->setSplitGapGrowthFactor(2.0f);

  // create the gridmanager that will take care of the grid
  ProtoSphere::GridManager pGridManager(ProtoSphere::GTHybridGrid, params, number_of_spheres);
  pGridManager.setMinVolumeDiffPerSphere(1.26f);
  pGridManager.setMaxVolumeDiffPerSphere(8.0f);
  pGridManager.setMinVolumeDiffPerStep(0.02f);
  pGridManager.setMaxVolumeDiffPerStep(0.25f);
  pGridManager.setTargetVolumeCovered(target_coverage * 100.0f);

  // initialize the grid with the data from the loaded object
  if (!pGridManager.init(&wrapper))
    return out;

  pGridManager.setIterationFileExport(false);

  // start the stepping of the grid
  pGridManager.run();

  ProtoSphere::float4* pSpheres = static_cast<ProtoSphere::float4*>(pGridManager.getSpheres());

  for (std::size_t i = 0;i < pGridManager.getCurrentNumberOfSpheres(); i++)
    out++ = Sphere_3(Point_3(pSpheres[i].x, pSpheres[i].y, pSpheres[i].z), pSpheres[i].w * pSpheres[i].w);

  return out;
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SPHERE_PACKING_H_
