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

#ifndef MESH_FEATURE_EXTRACTOR_H
#define MESH_FEATURE_EXTRACTOR_H

#include "cgal_types.h"


class Mesh_feature_extractor
{
public:
	Mesh_feature_extractor(Polyhedron_ex *pPolyhedron);
	~Mesh_feature_extractor();
	
	typedef Feature_skeleton<Polyhedron_ex::Vertex_handle,
	                         Polyhedron_ex::Halfedge_handle> skeleton;
	typedef Feature_backbone<Polyhedron_ex::Vertex_handle,
	                         Polyhedron_ex::Halfedge_handle> backbone;
	typedef Polyhedron_ex::Vertex_handle Vertex_handle;
	typedef Polyhedron_ex::Halfedge_handle Halfedge_handle;

private:
  Polyhedron_ex *m_pPolyhedron;
  skeleton *m_pSkeleton;

public: 

	// main functions	
	int extract_boundaries(bool sort = true);
	
private:
  bool add_boundary_backbone(int tag_free,int tag_processed);
};


#endif // MESH_FEATURE_EXTRACTOR_H
