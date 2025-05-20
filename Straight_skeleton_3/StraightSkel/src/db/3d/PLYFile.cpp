// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#include "db/3d/PLYFile.h"

#include "debug.h"
#include "data/3d/Polyhedron.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"
#include "data/3d/Triangle.h"
#include "data/3d/KernelFactory.h"
#include "util/StringFuncs.h"
#include "util/Configuration.h"

#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/SkelFacetData.h"

#include "db/3d/Surface_meshIO.h"

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/IO/PLY.h>

#include <cmath>
#include <exception>
#include <fstream>
#include <sstream>
#include <vector>

namespace db { namespace _3d {

PLYFile::PLYFile() {
    // intentionally does nothing
}

PLYFile::~PLYFile() {
    // intentionally does nothing
}

// This particular reader can read weights if they are stored as a property in the PLY file
PolyhedronSPtr PLYFile::load(const std::string& filename) {
    PolyhedronSPtr result = PolyhedronSPtr();

    std::ifstream ifs(filename.c_str());
    if (ifs.is_open()) {
        // @todo avoid building the intermediate mesh (but still use CGAL readers)
        CGAL::Surface_mesh<Point3> sm;
        bool success = CGAL::IO::read_PLY(ifs, sm);
        if (!success) {
          return {};
        }

        result = Surface_meshIO::load(sm, "filename='"+filename+"'; ");
    }

    return result;
}

} }
