// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_REGULAR_TRIANGULATION_3_WATCHER_BASE_H
#define CGAL_KINETIC_REGULAR_TRIANGULATION_3_WATCHER_BASE_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/Delaunay_triangulation_visitor_base_3.h>

namespace CGAL { namespace Kinetic {

struct Regular_triangulation_visitor_base_3:
public Delaunay_triangulation_visitor_base_3
{
//typedef Tr Triangulation;
    Regular_triangulation_visitor_base_3(){}

    template <class Key, class Cell>
        void pre_move(Key, Cell){}

    template <class Key, class Cell>
        void post_move(Key, Cell){}
};

} } //namespace CGAL::Kinetic
#endif
