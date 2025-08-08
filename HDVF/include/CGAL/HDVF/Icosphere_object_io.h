// Copyright (c) 2025 LIS Marseille (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alexandra Bac <alexandra.bac@univ-amu.fr>

#ifndef CGAL_HDVF_ICOSPHERE_OBJECT_H
#define CGAL_HDVF_ICOSPHERE_OBJECT_H

#include <CGAL/license/HDVF.h>

#include <iostream>
#include <cassert>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <CGAL/HDVF/Mesh_object_io.h>

namespace CGAL {
namespace HDVF {

using namespace std ;

/* Class used to build an icosphere embedding a triangular mesh (for Alexander Duality) */

class Icosphere_object_io : public Mesh_object_io
{
public:
    using Index=size_t ;
    using Lookup=std::map<std::pair<Index, Index>, Index>;

    Icosphere_object_io(size_t subdivisions, const IONodeType &c = IONodeType({0, 0, 0}), double r=1.) : Mesh_object_io(2, vertices_ico, triangles_ico)
    {
        for (size_t i=0; i<subdivisions; ++i)
        {
            subdivide();
        }
        rigid_transformation(c, r) ;
    }


    // Icosahedron
    inline static const float X=.525731112119133606f;
    inline static const float Z=.850650808352039932f;
    inline static const float N=0.f;

    inline static const std::vector<IONodeType> vertices_ico =
    {
        IONodeType({-X,N,Z}), IONodeType({X,N,Z}), IONodeType({-X,N,-Z}), IONodeType({X,N,-Z}),
        IONodeType({N,Z,X}), IONodeType({N,Z,-X}), IONodeType({N,-Z,X}), IONodeType({N,-Z,-X}),
        IONodeType({Z,X,N}), IONodeType({-Z,X, N}), IONodeType({Z,-X,N}), IONodeType({-Z,-X, N})
    };

    inline static const std::vector<IOCellType> triangles_ico =
    {
        {0,4,1},{0,9,4},{9,5,4},{4,5,8},{4,8,1},
        {8,10,1},{8,3,10},{5,3,8},{5,2,3},{2,7,3},
        {7,10,3},{7,6,10},{7,11,6},{11,0,6},{0,1,6},
        {6,1,10},{9,0,11},{9,11,2},{9,2,5},{7,2,11}
    };

    // Methods
    Index vertex_for_edge(Lookup& lookup, Index first, Index second)
    {
        Lookup::key_type key(first, second);
        if (key.first>key.second)
            std::swap(key.first, key.second);

        auto inserted=lookup.insert({key, nodes.size()});
        if (inserted.second)
        {
            IONodeType& edge0(nodes[first]);
            IONodeType& edge1(nodes[second]);
            IONodeType& point(edge0+edge1);
            normalize(point) ;
            add_node(point);
        }

        return inserted.first->second;
    }

    void subdivide()
    {
        Lookup lookup;
        std::vector<IOCellType> result ;

        for (IOCellType & each:cells)
        {
            std::array<Index, 3> mid;
            std::vector<size_t> each_vertex ;
            for (size_t c : each)
                each_vertex.push_back(c) ;
            for (size_t edge=0; edge<3; ++edge)
            {
                mid[edge]=vertex_for_edge(lookup,
                                          each_vertex[edge], each_vertex[(edge+1)%3]);
            }
            result.push_back({each_vertex[0], mid[0], mid[2]});
            result.push_back({each_vertex[1], mid[1], mid[0]});
            result.push_back({each_vertex[2], mid[2], mid[1]});
            result.push_back({mid[0], mid[1], mid[2]});
        }
        clear_cells() ;
        for (IOCellType c : result)
            add_cell(c) ;
    }

    void rigid_transformation(const IONodeType &c, double r)
    {
        for (size_t i=0; i<nvertices; ++i)
        {
            nodes[i] *= r ;
            nodes[i] += c ;
        }
    }
};

} /* end namespace HDVF */
} /* end namespace CGAL */


#endif // CGAL_HDVF_ICOSPHERE_OBJECT_H
