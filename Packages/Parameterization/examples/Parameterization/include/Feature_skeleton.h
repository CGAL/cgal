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
private:
  VERTEX_HANDLE m_begin;
  VERTEX_HANDLE m_end;
  std::list<HALFEDGE_HANDLE> m_halfedges;

public:

  // type of backbone
  enum backbone_type {SEAMING,
                      SHARP,
                      BOUNDARY};
private:
  backbone_type m_type;

public:
  Feature_backbone(backbone_type t = SEAMING)
  {
    m_begin = NULL;
    m_end = NULL;
    m_type = t;
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
    typedef typename std::list<HALFEDGE_HANDLE> list_halfedges;
    typedef typename list_halfedges::iterator iter_list_halfedges;
    iter_list_halfedges iter;
    for(iter  = pBackbone->halfedges()->begin();
        iter != pBackbone->halfedges()->end();
        iter++)
      m_halfedges.push_back(*iter);
  }

  // type
  backbone_type type() { return m_type; }
  void type(backbone_type t) { m_type = t; }

  // data access
  std::list<HALFEDGE_HANDLE>* halfedges() { return &m_halfedges; }
  const std::list<HALFEDGE_HANDLE>* halfedges() const { return &m_halfedges; }

  // misc
  bool is_closed() const { return (m_begin == m_end); }
  bool has_no_corner() const { return (m_begin == NULL && m_end == NULL); }
  VERTEX_HANDLE begin()  { return m_begin; }
  VERTEX_HANDLE end()    { return m_end;   }
  void begin(VERTEX_HANDLE hVertex)  { m_begin = hVertex; }
  void end(VERTEX_HANDLE hVertex)    { m_end = hVertex;   }
};


template<class VERTEX_HANDLE,class HALFEDGE_HANDLE>
class Feature_skeleton
{
private:
  typedef Feature_backbone<VERTEX_HANDLE,HALFEDGE_HANDLE> backbone;
  typedef typename backbone::backbone_type backbone_type;
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

  bool erase_one_type(backbone_type type)
  {
    typedef typename std::vector<backbone*>::iterator backbone_iterator;
    backbone_iterator ppBackbone;
    for(ppBackbone  = m_backbones.begin();
        ppBackbone != m_backbones.end();
        ppBackbone++)
    {
      backbone *pBackbone = *ppBackbone;
      if(pBackbone->type() == type)
      {
        delete pBackbone;
        m_backbones.erase(ppBackbone);
        return true;
      }
    }
    return false;
  }

  backbone *get_backbone_type(backbone_type type)
  {
    typedef typename std::vector<backbone*>::iterator backbone_iterator;
    backbone_iterator ppBackbone;
    for(ppBackbone  = m_backbones.begin();
        ppBackbone != m_backbones.end();
        ppBackbone++)
    {
      backbone *pBackbone = *ppBackbone;
      if(pBackbone->type() == type)
        return pBackbone;
    }
    return NULL;
  }

  void erase_type(backbone_type type)
  {
    while(erase_one_type(type));
  }

};

#endif // FEATURE_SKELETON_H



