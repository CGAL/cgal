// Copyright (c) 2001-2004  ENS of Paris (France).
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
// Author(s)     : Pierre Angelier, Michel Pocchiola

#ifndef CGAL_VISIBILITY_COMPLEX_ITEMS_H
#define CGAL_VISIBILITY_COMPLEX_ITEMS_H

#include <CGAL/Visibility_complex_vertex_base.h>
#include <CGAL/Visibility_complex_edge_base.h>
#include <CGAL/Visibility_complex_face_base.h>

CGAL_BEGIN_NAMESPACE

class Visibility_complex_items {
public:
    template <class VCA>
    struct Vertex_wrapper {
	typedef Visibility_complex_vertex_base<VCA> Vertex;
    };
    template <class VCA>
    struct Edge_wrapper {
	typedef Visibility_complex_edge_base<VCA>   Edge;
    };
    template <class VCA>
    struct Face_wrapper {
	typedef Visibility_complex_face_base<VCA>   Face;
    };
};

CGAL_END_NAMESPACE

#endif
