// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef DB_3D_TRIANGLEDAO_H
#define DB_3D_TRIANGLEDAO_H

#include "data/3d/ptrs.h"
#include "data/3d/Triangle.h"
#include "db/ptrs.h"
#include "db/3d/DAOFactory.h"
#include <string>

namespace db { namespace _3d {

using data::_3d::Triangle;
using data::_3d::TriangleSPtr;

class TriangleDAO {
friend class DAOFactory;
public:
    virtual ~TriangleDAO();
    std::string getTableSchema() const;
    int insert(TriangleSPtr triangle);
    bool del(TriangleSPtr triangle);
    TriangleSPtr find(int polyhedronid, int fid, int tid);
    bool update(TriangleSPtr triangle);
protected:
    TriangleDAO();
    int nextTID(int polyhedronid, int fid);
};

} }

#endif /* DB_3D_TRIANGLEDAO_H */
