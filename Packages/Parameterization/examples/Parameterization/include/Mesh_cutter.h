/***************************************************************************

    begin                : jan 02
    copyright            : (C) 2002 by Pierre Alliez
    email                : pierre.alliez@sophia.inria.fr
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef MESH_CUTTER_H
#define MESH_CUTTER_H

#include "Polyhedron_ex.h"
#include "Feature_skeleton.h"
#include <list>


class Mesh_cutter
{
    // Private types
    typedef My_kernel::Vector_3 Vector_3;
    typedef My_kernel::Point_3  Point_3;
    typedef Feature_backbone<Polyhedron_ex::Vertex_handle,
                             Polyhedron_ex::Halfedge_handle> backbone;
    enum {FREE,DONE,FIXED};

    // Fields
    Polyhedron_ex *m_pPolyhedron;  // the model to cut
    backbone *m_pBackbone;
    Polyhedron_ex::Facet_handle m_pSeedFacet;
    Polyhedron_ex::Vertex_handle m_pSeedVertex;

public:

    // life cycle
    Mesh_cutter(Polyhedron_ex *pPolyhedron)
    {
        CGAL_assertion(pPolyhedron != NULL);
        m_pPolyhedron = pPolyhedron;
        m_pBackbone = NULL;
    }
    ~Mesh_cutter() {}

    void cut(backbone *pBackbone);
    void cut_genus(backbone *pBackbone);

private:

    // genus > 0
    bool init();
    bool simplify();
    bool extend();
    void precompute_distances();
    Polyhedron_ex::Halfedge_handle pick_best_halfedge(
                    std::list<Polyhedron_ex::Halfedge_handle>::iterator &pos);
    void recursive_tag(Polyhedron_ex::Facet_handle pFacet,int index);
};


#endif // MESH_CUTTER_H
