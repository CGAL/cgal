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

#include <list>


class Mesh_cutter
{
// Public types
public:

    typedef std::list<Polyhedron_ex::Halfedge_handle> 
                                            Backbone;

// Private types
private:

    typedef My_kernel::Vector_3 Vector_3;
    typedef My_kernel::Point_3  Point_3;
    enum {FREE,DONE,FIXED};

// Public operations
public:

    // life cycle
    Mesh_cutter(Polyhedron_ex& polyhedron)
    {
        m_pPolyhedron = &polyhedron;
        m_pBackbone = NULL;
    }
    ~Mesh_cutter() {}

    void cut(Backbone& backbone);
    void cut_genus(Backbone& backbone);

// Private operations
private:

    // genus > 0
    bool init();
    bool simplify();
    bool extend();
    void precompute_distances();
    Polyhedron_ex::Halfedge_handle pick_best_halfedge(
                    std::list<Polyhedron_ex::Halfedge_handle>::iterator &pos);
    void recursive_tag(Polyhedron_ex::Facet_handle pFacet,int index);

// Fields
private:

    Polyhedron_ex*              m_pPolyhedron;  // the model to cut
    Backbone*                   m_pBackbone;    // the backbone to fill
    Polyhedron_ex::Facet_handle m_pSeedFacet;
    Polyhedron_ex::Vertex_handle m_pSeedVertex;
};


#endif // MESH_CUTTER_H
