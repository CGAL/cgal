// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef ALGO_3D_OUTWARDMESHOFFSET_H
#define ALGO_3D_OUTWARDMESHOFFSET_H

#include "data/3d/ptrs.h"
#include "cgal_kernel.h"

#include <CGAL/Surface_mesh.h>

#include <algorithm>
#include <filesystem>
#include <list>

namespace algo { namespace _3d {

using namespace data::_3d;

class OutwardMeshOffset
{
    using Mesh = CGAL::Surface_mesh<Point3>;

    using vertex_descriptor = typename boost::graph_traits<Mesh>::vertex_descriptor;
    using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;
    using edge_descriptor = boost::graph_traits<Mesh>::edge_descriptor;
    using face_descriptor = typename boost::graph_traits<Mesh>::face_descriptor;
    using vertex_iterator = typename boost::graph_traits<Mesh>::vertex_iterator;

    virtual ~OutwardMeshOffset();

public:
    // @todo customer specific, belongs into a .cpp
    // here currently because it's shared between `offset_mesh.cpp` and `convert_to_weighted_PLY.cpp`
    static bool assign_weights(Mesh& sm,
                               const char* weights_filename);

    static bool invert_and_add_bbox(Mesh& sm);
    static bool remove_bbox_and_invert(Mesh& sm);
    static PolyhedronSPtr convert(Mesh& sm);
    static PolyhedronSPtr preprocess(Mesh& sm);

public:
    // @todo add missing parameters:
    // - config file path
    static bool run(const char* mesh_filename,
                    const char* weights_filename,
                    const std::list<CGAL::FT>& save_offsets,
                    const std::filesystem::path save_path = std::filesystem::current_path());
};

} }

#endif /* ALGO_3D_OUTWARDMESHOFFSET_H */
