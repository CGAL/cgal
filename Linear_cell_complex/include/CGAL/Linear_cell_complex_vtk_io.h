// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_LINEAR_CELL_COMPLEX_VTK_IO_H
#define CGAL_LINEAR_CELL_COMPLEX_VTK_IO_H 1

#include <vector>

namespace CGAL {

/** \file Linear_cell_complex_vtk_io.h
 * Functions to import/export 3D Linear_cell_complex from/to VTK legacy ASCII format.
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
 * \brief Read a VTK legacy ASCII file and load it into a 3D Linear_cell_complex.
 * \ingroup PkgLinearCellComplexExamples
 *
 * \tparam LCC must be a Linear_cell_complex_for_combinatorial_map<3,3>
 * \param alcc The Linear_cell_complex to populate (will be cleared first)
 * \param filename Path to the VTK file
 * \param vertex_scalars Optional output vector to store per-vertex scalar values.
 *                      If provided, will be resized to match number of vertices.
 * \param volume_scalars Optional output vector to store per-volume scalar values.
 *                      If provided, will be resized to match number of volumes.
 * \return `true` if loading was successful, `false` otherwise
 */
template <typename LCC>
bool read_vtk(LCC& alcc,
              const char* filename,
              std::vector<float>* vertex_scalars = nullptr,
              std::vector<float>* volume_scalars = nullptr);

/**
 * \brief Write a 3D Linear_cell_complex to a VTK legacy ASCII file.
 * \ingroup PkgLinearCellComplexExamples
 *
 * \tparam LCC must be a Linear_cell_complex_for_combinatorial_map<3,3>
 * \param alcc The Linear_cell_complex to export
 * \param filename Path to the output VTK file
 * \param vertex_scalars Optional per-vertex scalar data. If provided, must have
 *                      same size as number of vertex attributes in the LCC.
 * \param volume_scalars Optional per-volume scalar data. If provided, must have
 *                      same size as number of 3-cells in the LCC.
 * \return `true` if writing was successful, `false` otherwise
 */
template <typename LCC>
bool write_vtk(const LCC& alcc,
               const char* filename,
               const std::vector<float>* vertex_scalars = nullptr,
               const std::vector<float>* volume_scalars = nullptr);

} // namespace CGAL

#include <CGAL/Linear_cell_complex_vtk_io_impl.h>

#endif // CGAL_LINEAR_CELL_COMPLEX_VTK_IO_H