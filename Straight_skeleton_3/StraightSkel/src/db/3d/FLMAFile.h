// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * @file:   db/3d/FLMAFile.h
 * @author: Gernot Walzl
 * @date:   2013-12-16
 */

#ifndef DB_3D_FLMAFILE_H
#define DB_3D_FLMAFILE_H

#include "db/3d/AbstractFile.h"
#include <string>

namespace db { namespace _3d {

using namespace data::_3d;

/**
 * AVL flma file format
 */
class FLMAFile : public AbstractFile {
public:
    virtual ~FLMAFile();
    static PolyhedronSPtr load(const std::string& filename);
    static bool save(const std::string& filename, PolyhedronSPtr polyhedron);
protected:
    FLMAFile();
};

} }

#endif /* DB_3D_FLMAFILE_H */
