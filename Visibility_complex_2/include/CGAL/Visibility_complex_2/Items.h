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
// $URL$
// $Id$
//
// Author(s)     : Pierre Angelier, Michel Pocchiola

#ifndef CGAL_VISIBILITY_COMPLEX_2_ITEMS_H
#define CGAL_VISIBILITY_COMPLEX_2_ITEMS_H

#include <CGAL/Visibility_complex_2/Vertex_base.h>
#include <CGAL/Visibility_complex_2/Edge_base.h>
#include <CGAL/Visibility_complex_2/Face_base.h>

CGAL_BEGIN_NAMESPACE
namespace Visibility_complex_2_details {

class Items {
public:
    template <class VCA>
    struct Vertex_wrapper {
	typedef Vertex_base<VCA> Vertex;
    };
    template <class VCA>
    struct Edge_wrapper {
	typedef Edge_base<VCA>   Edge;
    };
    template <class VCA>
    struct Face_wrapper {
	typedef Face_base<VCA>   Face;
    };
};
}

typedef Visibility_complex_2_details::Items Visibility_complex_2_items;

CGAL_END_NAMESPACE


#endif // CGAL_VISIBILITY_COMPLEX_2_ITEMS_H
