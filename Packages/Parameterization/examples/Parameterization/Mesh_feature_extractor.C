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

#include "Mesh_feature_extractor.h"

//*********************************************
// life cycle
//*********************************************
Mesh_feature_extractor::Mesh_feature_extractor(Polyhedron_ex *pPolyhedron)
{
  m_pPolyhedron = pPolyhedron;
  CGAL_assertion(m_pPolyhedron != NULL);
  m_pSkeleton = m_pPolyhedron->get_skeleton();
}
Mesh_feature_extractor::~Mesh_feature_extractor()
{
}

//*************************************************************
// extract boundaries
// return the number of boundary backbones
//*************************************************************
int Mesh_feature_extractor::extract_boundaries(bool sort)
{
  const int tag_free = 0;
  const int tag_processed = 1;
  m_pPolyhedron->tag_halfedges(tag_free); // unprocessed

  int nb = 0;

  // add isolated boundary backbones
  while(add_boundary_backbone(tag_free,tag_processed))
    nb++;

  // sort them by decreasing length if required
  if(nb>1 && sort)
  {
    int index = m_pPolyhedron->get_index_longest_backbone();
    skeleton *pSkeleton = m_pPolyhedron->get_skeleton();
    // put longest backbone first
    backbone *tmp = (*pSkeleton->backbones())[0];
    (*pSkeleton->backbones())[0] = (*pSkeleton->backbones())[index];
    (*pSkeleton->backbones())[index] = tmp;
  }

  std::cerr << "  " << nb << " boundary backbones added" << std::endl;
  return nb;
}


//*************************************************************
// add closed boundary backbone
//*************************************************************
bool Mesh_feature_extractor::add_boundary_backbone(int tag_free,
                                                   int tag_processed)
{
  Halfedge_handle seed_halfedge = m_pPolyhedron->get_border_halfedge_tag(tag_free);
  if(seed_halfedge == NULL)
    return false;

  // add one backbone
  std::cerr << "  add one closed boundary backbone...";
  backbone *pNewBackbone = new backbone;
  m_pSkeleton->backbones()->push_back(pNewBackbone);

  // add the seed halfedge
  pNewBackbone->halfedges()->push_back(seed_halfedge);
  seed_halfedge->tag(tag_processed); // processed

  // set the begin/end vertex NULL)
  Vertex_handle seed_vertex = seed_halfedge->prev()->vertex();
  pNewBackbone->begin(NULL);
  pNewBackbone->end(NULL);

  // fill it
  int size = 1;
  Halfedge_handle current_halfedge = seed_halfedge;
  do
  {
    Halfedge_handle next_halfedge = current_halfedge->next();
    if(next_halfedge->prev()->vertex() == seed_vertex)
      break;
    pNewBackbone->halfedges()->push_back(next_halfedge);
    next_halfedge->tag(tag_processed); // processed
    current_halfedge = next_halfedge;
    size++;
  }
  while(1);

  std::cerr << size << " halfedges" << std::endl;
  return true;
}


