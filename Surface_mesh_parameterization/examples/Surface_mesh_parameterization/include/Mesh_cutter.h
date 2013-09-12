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
#include <cassert>


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


//***************************************************
// simple cut for genus 0 mesh
//***************************************************
inline void Mesh_cutter::cut(Backbone& backbone)
{
    m_pBackbone = &backbone;

    // special init -> tag all vertices, but two
    m_pPolyhedron->tag_vertices(FREE);
    Polyhedron_ex::Vertex_handle pVertexMin,pVertexMax;
    m_pPolyhedron->farthest_point_aligned(pVertexMin,pVertexMax);
    pVertexMin->tag(FIXED);
    pVertexMax->tag(FIXED);
    init();

    // cutting
    while(extend()) {}
}

/////////////////////////////////////////////////////
// GENUS > 0
/////////////////////////////////////////////////////

//***************************************************
// cut for genus>0 mesh
//***************************************************
inline void Mesh_cutter::cut_genus(Backbone& backbone)
{
    m_pBackbone = &backbone;

    // init
    m_pPolyhedron->tag_vertices(FREE); // all free
    init();

    // cutting
    while(extend()) {}
}

//***************************************************
// init
//***************************************************
inline bool Mesh_cutter::init()
{
    // tag facets
    m_pPolyhedron->tag_facets(FREE);

    // compute bounding box and center
    double xmin = m_pPolyhedron->minimum(0);
    double ymin = m_pPolyhedron->minimum(1);
    double zmin = m_pPolyhedron->minimum(2);
    double xmax = m_pPolyhedron->maximum(0);
    double ymax = m_pPolyhedron->maximum(1);
    double zmax = m_pPolyhedron->maximum(2);
    double xcenter = 0.5*(xmin+xmax);
    double ycenter = 0.5*(ymin+ymax);
    double zcenter = 0.5*(zmin+zmax);
    Point_3 center(xcenter,ycenter,zcenter);

    // get closest facet
    m_pSeedFacet = m_pPolyhedron->get_closest_inner_facet(center);
    assert(m_pSeedFacet != NULL);

    Polyhedron_ex::Halfedge_handle he = m_pSeedFacet->halfedge();
    assert(he != NULL);
    assert(m_pBackbone != NULL);
    m_pBackbone->push_back(he);
    m_pBackbone->push_back(he->next());
    m_pBackbone->push_back(he->next()->next());

    precompute_distances();
    m_pSeedFacet->tag(DONE);

    return true;
}

//***************************************************
// extend
//***************************************************
inline bool Mesh_cutter::extend()
{
    std::list<Polyhedron_ex::Halfedge_handle>::iterator pos;
    Polyhedron_ex::Halfedge_handle pHalfedge = pick_best_halfedge(pos);
    if(pHalfedge == NULL)
    return false;

    // flag facet
    pHalfedge->opposite()->facet()->tag(DONE);

    // insert halfedge
    std::list<Polyhedron_ex::Halfedge_handle>::iterator tmp =
    m_pBackbone->insert(pos,pHalfedge->opposite()->next()->next());
    m_pBackbone->insert(tmp,pHalfedge->opposite()->next());

    // remove this one
    m_pBackbone->remove(pHalfedge);

    // simplify current backbone
    while(simplify()) {}
    return true;
}

//***************************************************
// simplify
//***************************************************
inline bool Mesh_cutter::simplify()
{
    // cleanup
    std::list<Polyhedron_ex::Halfedge_handle>::iterator iter;
    for(iter  = m_pBackbone->begin();
        iter != m_pBackbone->end();
        iter++)
    {
        Polyhedron_ex::Halfedge_handle pHalfedge = (*iter);
        Polyhedron_ex::Halfedge_handle opposite = pHalfedge->opposite();

        // get next halfedge in the list
        iter++;
        Polyhedron_ex::Halfedge_handle pNext = NULL;
        if(iter == m_pBackbone->end()) // loop
            pNext = (*m_pBackbone->begin());
        else
            pNext = (*iter);

        if(pNext == opposite &&
            pHalfedge->vertex()->tag() == FREE)
        {
            m_pBackbone->remove(pHalfedge);
            m_pBackbone->remove(opposite);
            return true;
        }

        iter--; // restore
    }
    return false;
}

//***************************************************
// precompute_distances
//***************************************************
inline void Mesh_cutter::precompute_distances()
{
    Polyhedron_ex::Halfedge_iterator pHalfedge;
    for(pHalfedge = m_pPolyhedron->halfedges_begin();
        pHalfedge != m_pPolyhedron->halfedges_end();
        pHalfedge++)
    pHalfedge->distance(m_pPolyhedron->distance(m_pSeedFacet,pHalfedge));
}

//***************************************************
// pick_best_halfedge
//***************************************************
inline Polyhedron_ex::Halfedge_handle Mesh_cutter::pick_best_halfedge(
    std::list<Polyhedron_ex::Halfedge_handle>::iterator &pos)
{
    Polyhedron_ex::Halfedge_handle pBest = NULL;
    double min_distance = 1e308; //

    // cleanup
    std::list<Polyhedron_ex::Halfedge_handle>::iterator iter;
    for(iter  = m_pBackbone->begin();
        iter != m_pBackbone->end();
        iter++)
    {
        Polyhedron_ex::Halfedge_handle pHalfedge = (*iter);
        Polyhedron_ex::Halfedge_handle opposite = pHalfedge->opposite();
        Polyhedron_ex::Facet_handle pFacet = opposite->facet();

        // check
        if(pHalfedge->is_border() ||
            pFacet == NULL)
            continue;

        if(pFacet->tag() == DONE)
            continue;

        // no border vertex
        Polyhedron_ex::Vertex_handle pVertex = opposite->next()->vertex();
        if(m_pPolyhedron->is_border(pVertex))
            continue;

        // precomputed distance
        double distance = pHalfedge->distance();
        if(distance < min_distance)
        {
            pos = iter;
            pBest = pHalfedge;
            min_distance = distance;
        }
    }
    return pBest;
}


#endif // MESH_CUTTER_H
