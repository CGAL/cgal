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

//*********************************************
// build sharp skeleton
// return the number of sharp crease backbones
//*********************************************
int Mesh_feature_extractor::extract_features()
{
  const int tag_free = 0;
  const int tag_processed = 1;
  m_pPolyhedron->tag_halfedges(tag_free); // unprocessed

  // from corners first
  int nb = 0;
  while(add_feature_backbone_from_corner(tag_free,tag_processed))
    nb++;

  // isolated feature backbones
  while(add_feature_backbone(tag_free,tag_processed))
    nb++;

  std::cerr << "  overall " << nb << " feature backbone(s)" << std::endl;
  return nb;
}

//*************************************************************
// extract boundaries
// return the number of boundary backbones
//*************************************************************
int Mesh_feature_extractor::extract_boundaries(bool open_backbones,
                                               bool sort)
{
  const int tag_free = 0;
  const int tag_processed = 1;
  m_pPolyhedron->tag_halfedges(tag_free); // unprocessed

  // add boundary backbones from corners
  int nb = 0;
  if(open_backbones)
    while(add_boundary_backbone_from_corner(tag_free,tag_processed))
      nb++;

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

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

//*************************************************************
// add feature backbone issued from corner
//*************************************************************
bool Mesh_feature_extractor::add_feature_backbone_from_corner(int tag_free,
                                                              int tag_processed)
{
  // find sharp backbone issued from a corner
  Polyhedron_ex::Vertex_handle seed_corner =
    m_pPolyhedron->get_corner_sharp_or_seaming_halfedge_tag(tag_free);
  if(seed_corner == NULL)
    return false;

  Polyhedron_ex::Halfedge_handle halfedge =
    m_pPolyhedron->get_sharp_or_seaming_halfedge_tag(seed_corner,tag_free);
  CGAL_assertion(halfedge != NULL);

  // reverse it (start from seed_corner)
  // note that both are tagged as proceed in order
  // to avoid redundant sharp backbones
  halfedge->tag(tag_processed); // processed
  halfedge = halfedge->opposite();
  CGAL_assertion(halfedge != NULL);
  halfedge->tag(tag_processed); // processed opposite also

  // add one backbone
  std::cerr << "  add one feature backbone attached to a corner...";
  backbone::backbone_type type = backbone::SHARP;
  if(halfedge->is_seaming())
    type = backbone::SEAMING;

  backbone *pNewBackbone = new backbone(type);
  CGAL_assertion(pNewBackbone != NULL);
  m_pPolyhedron->get_skeleton()->backbones()->push_back(pNewBackbone);

  // add halfedge
  pNewBackbone->halfedges()->push_back(halfedge);
  Polyhedron_ex::Vertex_handle seed_vertex = halfedge->prev()->vertex();
  Polyhedron_ex::Vertex_handle current_vertex = halfedge->vertex();
  pNewBackbone->begin(seed_vertex);
  int size = 1;
  while(!current_vertex->is_corner())
  {
    // pick the next halfedge
    Polyhedron_ex::Halfedge_handle next =
      m_pPolyhedron->get_sharp_or_seaming_halfedge_tag(current_vertex,tag_free);
    CGAL_assertion(next != NULL);
    next->tag(tag_processed); // tag this one
    next = next->opposite(); // reverse it
    CGAL_assertion(next != NULL);
    next->tag(tag_processed); // tag also the opposite

    // add it to the backbone and update the end
    pNewBackbone->halfedges()->push_back(next);
    current_vertex = next->vertex();
    CGAL_assertion(current_vertex != NULL);
    size++;
  }

  // set the end (might be the same)
  pNewBackbone->end(current_vertex);
  std::cerr << size << " halfedges" << std::endl;

  return true;
}

//*************************************************************
// add isolated feature backbone
//*************************************************************
bool Mesh_feature_extractor::add_feature_backbone(int tag_free,
                                                  int tag_processed)
{
  // some (closed) sharp backbones may not be linked
  // to any corner, we have to grab them also
  Halfedge_handle seed_halfedge =
    m_pPolyhedron->get_sharp_or_seaming_halfedge_tag(tag_free);
  if(seed_halfedge == NULL)
    return false;

  // add one backbone
  std::cerr << "  add one closed sharp backbone...";
  backbone *pNewBackbone = new backbone(backbone::SHARP);
  m_pSkeleton->backbones()->push_back(pNewBackbone);

  // add the seed halfedge
  pNewBackbone->halfedges()->push_back(seed_halfedge);
  seed_halfedge->tag(tag_processed); // processed
  seed_halfedge = seed_halfedge->opposite();
  CGAL_assertion(seed_halfedge != NULL);
  seed_halfedge->tag(tag_processed); // processed

  // set the begin/end vertex the same)
  Vertex_handle seed_vertex = seed_halfedge->prev()->vertex();
  CGAL_assertion(seed_vertex != NULL);

  int size = 1;
  Vertex_handle current_vertex = seed_halfedge->vertex();
  //CGAL_assertion(current_vertex->is_crease());
  Halfedge_handle next = NULL;
  // pick the next halfedge
  while((next =
    m_pPolyhedron->get_sharp_or_seaming_halfedge_tag(current_vertex,0)) != NULL)
  {
    next->tag(tag_processed); // tag this one
    next = next->opposite(); // reverse it
    CGAL_assertion(next != NULL);
    next->tag(tag_processed); // tag also the opposite

    // add it to the backbone and update the end
    pNewBackbone->halfedges()->push_back(next);
    current_vertex = next->vertex();
    CGAL_assertion(current_vertex != NULL);
    size++;
  }
  std::cerr << size << " halfedges" << std::endl;

  // set no begin/end
  pNewBackbone->begin(NULL);
  pNewBackbone->end(NULL);

  return true;
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
  backbone *pNewBackbone = new backbone(backbone::BOUNDARY);
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

//*************************************************************
// add boundary backbone from corner
//*************************************************************
bool Mesh_feature_extractor::add_boundary_backbone_from_corner(int tag_free,
                                                               int tag_processed)
{
  // find one border backbone issued from a seed corner
  Vertex_handle seed_corner = m_pPolyhedron->get_corner_border_halfedge_tag(tag_free);
  if(seed_corner == NULL)
    return false;

  // careful, the halfedge is issued from seed_corner
  Halfedge_handle current_halfedge = m_pPolyhedron->get_border_halfedge_tag(seed_corner,tag_free);
  current_halfedge->tag(tag_processed); // processed

  // add one backbone
  std::cerr << "  add one boundary backbone attached to a corner...";
  backbone *pNewBackbone = new backbone(backbone::BOUNDARY);
  CGAL_assertion(pNewBackbone != NULL);
  m_pSkeleton->backbones()->push_back(pNewBackbone);

  // add halfedge
  pNewBackbone->halfedges()->push_back(current_halfedge);
  Vertex_handle seed_vertex = current_halfedge->prev()->vertex();
  CGAL_assertion(seed_vertex == seed_corner);
  pNewBackbone->begin(seed_vertex);

  // fill it
  int size = 1;
  Vertex_handle current_vertex = NULL;
  do
  {
    // detect corner or loop
    current_vertex = current_halfedge->vertex();
    Halfedge_handle next_halfedge = current_halfedge->next();
    CGAL_assertion(next_halfedge != NULL);

    if(current_vertex->is_corner() ||
       current_vertex == seed_vertex)
       // || next_halfedge->tag() == tag_processed
      break;

    pNewBackbone->halfedges()->push_back(next_halfedge);
    next_halfedge->tag(tag_processed); // processed
    current_halfedge = next_halfedge;
    size++;
  }
  while(1);

  // set the end (might be the same or/and another corner)
  pNewBackbone->end(current_vertex);
  std::cerr << size << " halfedges, ";
  if(pNewBackbone->is_closed())
    std::cerr << "closed backbone";
  else
    std::cerr << "opened backbone";
  std::cerr << std::endl;

  return true;
}

