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

#ifndef CGAL_CIRCLE_2_BITANGENT_2_INTERSECTION_H
#define CGAL_CIRCLE_2_BITANGENT_2_INTERSECTION_H

#include <CGAL/basic.h>
#include <CGAL/Point_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Visibility_complex/Bitangent_2.h>

CGAL_BEGIN_NAMESPACE

// -----------------------------------------------------------------------------

template < class R_ , class C_ >
bool do_intersect( const Circle_2<R_>& c1, const Bitangent_2<C_>& c2 )
{
    cerr << "Not implemented Circle_2 - Bitangent_2 intersection !" << endl;
    return false;
}

// -----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // CGAL_CIRCLE_2_BITANGENT_2_INTERSECTION_H
