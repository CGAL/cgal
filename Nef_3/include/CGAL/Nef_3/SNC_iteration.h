// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
#ifndef CGAL_SNC_ITERATION_H
#define CGAL_SNC_ITERATION_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/Nef_S2/SM_iteration.h>

#define CGAL_forall_vertices(x,SNC)\
for(x = (SNC).vertices_begin(); x != (SNC).vertices_end(); ++x) 

#define CGAL_forall_halfedges(x,SNC)\
for(x = (SNC).halfedges_begin(); x != (SNC).halfedges_end(); ++x)

#define CGAL_forall_edges(x,SNC)\
for(x = (SNC).halfedges_begin(); x != (SNC).halfedges_end(); ++x) \
if ( x->is_twin() ) continue; else

#define CGAL_forall_halffacets(x,SNC)\
for(x = (SNC).halffacets_begin(); x != (SNC).halffacets_end(); ++x) 
 
#define CGAL_forall_facets(x,SNC)\
for(x = (SNC).halffacets_begin(); x != (SNC).halffacets_end(); ++x) \
if ( x->is_twin() ) continue; else

#define CGAL_forall_volumes(x,SNC)\
for(x = (SNC).volumes_begin(); x != (SNC).volumes_end(); ++x) 

#define CGAL_forall_facet_cycles_of(x,F)\
for(x = (F)->facet_cycles_begin(); x != (F)->facet_cycles_end(); ++x) 

#define CGAL_forall_shells_of(x,C)\
for(x = (C)->shells_begin(); x != (C)->shells_end(); ++x) 

#define CGAL_forall_svertices_of(x,V)\
for(x = (V)->svertices_begin(); x != (V)->svertices_end(); ++x)

#define CGAL_forall_sedges_of(x,V)\
for(x = (V)->shalfedges_begin(); x != (V)->shalfedges_end(); ++(++x))

#define CGAL_forall_shalfedges_of(x,V)\
for(x = (V)->shalfedges_begin(); x != (V)->shalfedges_end(); ++x)

#define CGAL_forall_sfaces_of(x,V)\
for(x = (V)->sfaces_begin(); x != (V)->sfaces_end(); ++x)

#endif //CGAL_SNC_ITERATION_H
