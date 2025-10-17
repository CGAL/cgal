namespace CGAL {
namespace IO {

/** \file VTK.h
 * Functions to import/export 3D Linear_cell_complex from/to VTK legacy ASCII
 * format.
 *
 * Only supports:
 * - Linear_cell_complex_for_combinatorial_map<3,3>
 * - VTK legacy ASCII format (.vtk files)
 * - Optional scalar fields for vertices and volumes
 *
 * Supported VTK cell types:
 * - VTK_TETRA (10): Tetrahedron
 * - VTK_VOXEL (11): Voxel (special hexahedron ordering)
 * - VTK_HEXAHEDRON (12): Hexahedron
 * - VTK_WEDGE (13): Prism/Wedge
 * - VTK_PYRAMID (14): Pyramid
 * - VTK_PENTAGONAL_PRISM (15): Pentagonal prism
 * - VTK_HEXAGONAL_PRISM (16): Hexagonal prism
 * - VTK_POLYHEDRON (42): Generic polyhedron
 */

/**
 * \brief Read a VTK legacy ASCII file and load it into a 3D
 *        Linear_cell_complex.
 * \ingroup PkgLinearCellComplexRefIOVTK
 *
 * \tparam LCC must be a Linear_cell_complex_for_combinatorial_map<3,3>
 * \tparam VertexScalarType Type for vertex scalar data (default: float)
 * \tparam VolumeScalarType Type for volume scalar data (default: float)
 * \param filename Path to the VTK file
 * \param alcc The Linear_cell_complex to populate (will be cleared first)
 * \param vertex_scalars Optional output vector to store per-vertex scalar values.
 *                      If provided, will be resized to match number of vertices.
 * \param volume_scalars Optional output vector to store per-volume scalar values.
 *                      If provided, will be resized to match number of volumes.
 * \return `true` if loading was successful, `false` otherwise
 */
template <typename LCC, typename VertexScalarType, typename VolumeScalarType>
bool read_VTK(const char* filename,
              LCC& alcc,
              std::vector<VertexScalarType>* vertex_scalars,
              std::vector<VolumeScalarType>* volume_scalars);

/**
 * \brief Write a 3D Linear_cell_complex to a VTK legacy ASCII file.
 * \ingroup PkgLinearCellComplexRefIOVTK
 *
 * \tparam LCC must be a Linear_cell_complex_for_combinatorial_map<3,3>
 * \tparam VertexScalarType Type for vertex scalar data (default: float)
 * \tparam VolumeScalarType Type for volume scalar data (default: float)
 * \param filename Path to the output VTK file
 * \param alcc The Linear_cell_complex to export
 * \param vertex_scalars Optional per-vertex scalar data. If provided, must have
 *                      same size as number of vertex attributes in the LCC.
 * \param volume_scalars Optional per-volume scalar data. If provided, must have
 *                      same size as number of 3-cells in the LCC.
 * \return `true` if writing was successful, `false` otherwise
 */
template <typename LCC, typename VertexScalarType, typename VolumeScalarType>
bool write_VTK(const char* filename,
               const LCC& alcc,
               const std::vector<VertexScalarType>* vertex_scalars,
               const std::vector<VolumeScalarType>* volume_scalars);

} // namespace IO
} // namespace CGAL
