//=============================================================================
// Copyright (C) 2001-2005 by Computer Graphics Group, RWTH Aachen
// Copyright (C) 2011 by Graphics & Geometry Group, Bielefeld University
//
// SPDX-License-Identifier: GPL-2.0-only
//
//=============================================================================


#ifndef SURFACE_MESH_IO_H
#define SURFACE_MESH_IO_H


//== INCLUDES =================================================================

#include <string>
#include "Surface_mesh.h"

//=============================================================================


bool read_mesh(Surface_mesh& mesh, const std::string& filename);
template<typename NamedParameters>
bool read_off(Surface_mesh& mesh, const std::string& filename, NamedParameters& np);
bool read_off(Surface_mesh& mesh, const std::string& filename);

bool write_mesh(const Surface_mesh& mesh, const std::string& filename);
template<typename NamedParameters>
bool write_off(const Surface_mesh& mesh, const std::string& filename, const NamedParameters& np);
bool write_off(const Surface_mesh& mesh, const std::string& filename);


//=============================================================================
#endif // SURFACE_MESH_IO_H
//=============================================================================
