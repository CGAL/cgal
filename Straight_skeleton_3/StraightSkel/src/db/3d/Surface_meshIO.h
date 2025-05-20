// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef DB_3D_SURFACE_MESH_H
#define DB_3D_SURFACE_MESH_H

#include "db/3d/AbstractFile.h"
#include "cgal_kernel.h"

#include "data/3d/skel/ptrs.h"

#include <CGAL/Surface_mesh.h>

#include <map>
#include <string>

namespace db { namespace _3d {

using namespace data::_3d;

/**
 * PLY file format
 */
class Surface_meshIO : public AbstractFile {

    using Mesh = CGAL::Surface_mesh<Point3>;
    using vertex_descriptor = typename boost::graph_traits<Mesh>::vertex_descriptor;
    using halfedge_descriptor = typename boost::graph_traits<Mesh>::halfedge_descriptor;
    using edge_descriptor = typename boost::graph_traits<Mesh>::edge_descriptor;
    using face_descriptor = typename boost::graph_traits<Mesh>::face_descriptor;

public:
    virtual ~Surface_meshIO();

    static PolyhedronSPtr load(const Mesh& sm,
                               const std::string& description,
                               std::map<edge_descriptor, EdgeWPtr>& e2e);
    static PolyhedronSPtr load(const Mesh& sm,
                               const std::string& description = {});

    static bool save(const PolyhedronSPtr& polyhedron,
                     CGAL::Surface_mesh<Point3>& sm,
                     bool do_triangulate = true,
                     bool convert_to_double = true);

protected:
    Surface_meshIO();
};

} }

#endif /* DB_3D_SURFACE_MESH_H */
