// Copyright (c) 2016-2017 INRIA Nancy Grand-Est (France).
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
// Author(s)     : Iordan Iordanov <iordan.iordanov@loria.fr>

#include <iostream>
#include <CGAL/CORE_Expr.h>
#include <CGAL/exact_complex.h>
#include <CGAL/Point_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Hyperbolic_octagon_translation_matrix.h>
#include <CGAL/Hyperbolic_octagon_word_4.h>
#include <vector>

#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>

using namespace CGAL;
using namespace std;

int main(void) {

    typedef CORE::Expr                                                          NT;
    typedef CGAL::Cartesian<NT>                                                 Kernel;
    typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel> Traits;
    typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>        Triangulation;

    Triangulation tr;
    tr.insert_dummy_points(false);
    CGAL_assertion(tr.is_valid());

    cout << "triangulation works!" << 
    cout << "nb of vertices: " << tr.number_of_vertices() << endl;
    cout << "nb of faces: " << tr.number_of_faces() << endl;

    return 0;
}