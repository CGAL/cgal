// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef DB_3D_PLYFILE_H
#define DB_3D_PLYFILE_H

#include "db/3d/AbstractFile.h"
#include "cgal_kernel.h"

#include <CGAL/Surface_mesh.h>

#include <string>

namespace db { namespace _3d {

using namespace data::_3d;

/**
 * PLY file format
 */
class PLYFile : public AbstractFile {

    using Mesh = CGAL::Surface_mesh<Point3>;
    using vertex_descriptor = typename boost::graph_traits<Mesh>::vertex_descriptor;
    using halfedge_descriptor = typename boost::graph_traits<Mesh>::halfedge_descriptor;
    using edge_descriptor = typename boost::graph_traits<Mesh>::edge_descriptor;
    using face_descriptor = typename boost::graph_traits<Mesh>::face_descriptor;

public:
    virtual ~PLYFile();

    static PolyhedronSPtr load(const std::string& filename);

protected:
    PLYFile();
};

} }

#endif /* DB_3D_PLYFILE_H */
