// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$ $Date$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KDS_DELAUNAY_TRIANGULATION_2_WATCHER_BASE_H
#define CGAL_KDS_DELAUNAY_TRIANGULATION_2_WATCHER_BASE_H
#include <CGAL/KDS/basic.h>

CGAL_KDS_BEGIN_NAMESPACE

struct Delaunay_triangulation_visitor_base_2
{
    Delaunay_triangulation_visitor_base_2(){}
    template <class Ok>
    void remove_vertex(Ok) {
    }

    template <class Ok>
    void create_vertex(Ok) {
    }

    template <class Ok>
    void modify_vertex(Ok) {
    }

    template <class It>
    void create_faces(It, It) {
    }

    template <class It>
    void remove_faces(It, It) {
    }

    template <class E>
    void before_flip(E) {

    }
    template <class E>
    void after_flip(E) {

    }
};

CGAL_KDS_END_NAMESPACE
#endif
