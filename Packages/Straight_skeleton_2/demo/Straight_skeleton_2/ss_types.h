// Copyright (c) 2002  Max Planck Institut fuer Informatik (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Radu Ursu

#include <CGAL/Straight_skeleton_builder_2.h>

typedef CGAL::Straight_skeleton_2<K>                          Ssds;
typedef CGAL::Straight_skeleton_builder_traits_2<K>           BuilderTraits;
typedef CGAL::Straight_skeleton_builder_2<BuilderTraits,Ssds> Builder;

typedef Ssds::Halfedge_iterator     Halfedge_iterator;
typedef Ssds::Vertex_handle         Vertex_handle;
typedef Ssds::Face_const_iterator   Face_const_iterator;
typedef Ssds::Halfedge_const_handle Halfedge_const_handle ;
typedef Ssds::Vertex_const_handle   Vertex_const_handle ;
