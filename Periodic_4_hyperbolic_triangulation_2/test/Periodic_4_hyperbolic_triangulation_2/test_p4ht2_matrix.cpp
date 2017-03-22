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
#include <vector>

using namespace CGAL;
using namespace std;

int main(void) {

    typedef CORE::Expr                                  NT;
    typedef Hyperbolic_octagon_translation_matrix<NT>   Matrix;
    typedef Point_2< Cartesian<NT> >                    Point;
    // typedef Circle_2< Cartesian<NT> >   Circle;

    Matrix m;
    cout << "Identity matrix: " << m << endl;

    vector<Matrix> gens;
    get_generators(gens);
    for (int i = 0; i < gens.size(); i++) {
        cout << "g[" << i << "] = " << gens[i] << endl;
    }

    CGAL_assertion(gens[0]*gens[4] == m);
    CGAL_assertion(gens[1]*gens[5] == m);
    CGAL_assertion(gens[2]*gens[6] == m);
    CGAL_assertion(gens[3]*gens[7] == m);

    Point o(NT(0), NT(0));
    vector<Point> imp;
    for (int i = 0; i < gens.size(); i++) {
        imp.push_back(gens[i].apply(o));
        cout << "imp[" << i << "] = " << imp[i] << endl;
    }

    CGAL_assertion(imp[0] == Point(CGAL::sqrt(NT(2))/(CGAL::sqrt(NT(1)+CGAL::sqrt(NT(2)))), NT(0)));
    CGAL_assertion(imp[1] == Point(NT(1)/(CGAL::sqrt(NT(1)+CGAL::sqrt(NT(2)))), NT(1)/(CGAL::sqrt(NT(1)+CGAL::sqrt(NT(2))))));

    cout << "test concluded successfully!" << endl;
    return 0;
}