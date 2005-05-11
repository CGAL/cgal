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

#ifndef FEATURE_SKELETON_H
#define FEATURE_SKELETON_H

#include <CGAL/basic.h>

#include <list>
#include <vector>


template<class VERTEX_HANDLE,class HALFEDGE_HANDLE>
class Feature_backbone
{
// Public types
public:
    typedef VERTEX_HANDLE                           Vertex_handle;
    typedef HALFEDGE_HANDLE                         Halfedge_handle;
    typedef typename std::list<HALFEDGE_HANDLE>     Halfedge_list;
    typedef typename Halfedge_list::iterator        Halfedge_list_iterator;
    typedef typename Halfedge_list::const_iterator  Halfedge_list_const_iterator;

// Public operations
public:
    Feature_backbone()
    {
        m_begin = NULL;
        m_end = NULL;
    }
    ~Feature_backbone()
    {
        clear();
    }

    void clear()
    {
        m_halfedges.clear();
        m_begin = NULL;
        m_end = NULL;
    }

    // copy from
    void copy_from(Feature_backbone *pBackbone)
    {
        m_halfedges.clear();
        Halfedge_list_iterator iter;
        for(iter  = pBackbone->halfedges()->begin();
            iter != pBackbone->halfedges()->end();
            iter++)
            m_halfedges.push_back(*iter);
    }

    // data access
    Halfedge_list*       halfedges()        { return &m_halfedges; }
    const Halfedge_list* halfedges() const  { return &m_halfedges; }
    Vertex_handle begin()                   { return m_begin; }
    Vertex_handle end()                     { return m_end;   }
    void begin(Vertex_handle hVertex)       { m_begin = hVertex; }
    void end(Vertex_handle hVertex)         { m_end = hVertex;   }

// Fields
private:
    Vertex_handle  m_begin;
    Vertex_handle  m_end;
    Halfedge_list m_halfedges;

}; // Feature_backbone


template<class VERTEX_HANDLE,class HALFEDGE_HANDLE>
class Feature_skeleton
{
private:
    typedef Feature_backbone<VERTEX_HANDLE,HALFEDGE_HANDLE> backbone;
    std::vector<backbone*> m_backbones;

public:

    // data access
    std::vector<backbone*> *backbones() { return &m_backbones; }

    // life cycle
    Feature_skeleton() {}
    ~Feature_skeleton()
    {
        free();
    }

    void free()
    {
        unsigned int i;
        for(i=0;i<m_backbones.size();i++)
        {
            delete m_backbones[i];
            m_backbones[i] = NULL;
        }
        m_backbones.clear();
    }

    bool erase(backbone *pBackboneToErase)
    {
        typedef typename std::vector<backbone*>::iterator backbone_iterator;
        backbone_iterator ppBackbone;
        for(ppBackbone  = m_backbones.begin();
            ppBackbone != m_backbones.end();
            ppBackbone++)
        {
            backbone *pBackbone = *ppBackbone;
            if(pBackbone == pBackboneToErase)
            {
            delete pBackbone;
            m_backbones.erase(ppBackbone);
            return true;
            }
        }
        return false;
    }

}; // Feature_skeleton


#endif // FEATURE_SKELETON_H



