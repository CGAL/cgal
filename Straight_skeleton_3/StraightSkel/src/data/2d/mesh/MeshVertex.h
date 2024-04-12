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
 * @file   data/2d/mesh/MeshVertex.h
 * @author Gernot Walzl
 * @date   2014-01-24
 */

#ifndef DATA_2D_MESH_MESHVERTEX_H
#define DATA_2D_MESH_MESHVERTEX_H

#include "data/2d/ptrs.h"
#include "data/2d/mesh/ptrs.h"
#include <list>
#include <string>

namespace data { namespace _2d { namespace mesh {

class MeshVertex {
public:
    virtual ~MeshVertex();
    static MeshVertexSPtr create(Point2SPtr point);

    Point2SPtr getPoint() const;
    void setPoint(Point2SPtr point);
    MeshSPtr getMesh() const;
    void setMesh(MeshSPtr mesh);

    void addCell(MeshCellSPtr cell);
    bool removeCell(MeshCellSPtr cell);
    bool containsCell(MeshCellSPtr cell) const;

    MeshCellSPtr firstCell() const;

    MeshVertexSPtr next(MeshCellSPtr cell) const;
    MeshVertexSPtr prev(MeshCellSPtr cell) const;

    std::list<MeshCellWPtr>& cells();

    unsigned int countCells() const;

#ifdef USE_CGAL
    CGAL::FT getX() const;
    CGAL::FT getY() const;
#else
    double getX() const;
    double getY() const;
#endif

    std::string toString() const;

protected:
    MeshVertex(Point2SPtr point);
    Point2SPtr point_;
    MeshWPtr mesh_;
    std::list<MeshCellWPtr> cells_;
};

} } }

#endif /* DATA_2D_MESH_MESHVERTEX_H */
