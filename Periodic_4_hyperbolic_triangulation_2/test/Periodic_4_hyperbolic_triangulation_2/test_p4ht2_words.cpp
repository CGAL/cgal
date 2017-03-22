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

using namespace CGAL;
using namespace std;

int main(void) {

    typedef CORE::Expr                                          NT;
    typedef Hyperbolic_octagon_translation_matrix<NT>           Matrix;
    typedef Point_2< Cartesian<NT> >                            Point;
    typedef Hyperbolic_octagon_word_4<unsigned short int, NT>   Word;

    Word w;
    cout << "empty word: " << w << ", matrix: " << w.get_matrix() << endl;

    Word a(0), ab(0, 5), abc(0, 5, 2), abcd(0, 5, 2, 7), dcb(7, 2, 5), dc(7, 2), d(7);
    cout << "a    = " << a << ",    matrix: " << a.get_matrix()     << endl;
    cout << "ab   = " << ab << ",   matrix: " << ab.get_matrix()    << endl;
    cout << "abc  = " << abc << ",  matrix: " << abc.get_matrix()   << endl;
    cout << "abcd = " << abcd << ", matrix: " << abcd.get_matrix()  << endl;
    cout << "dcb  = " << dcb << ",  matrix: " << dcb.get_matrix()   << endl;
    cout << "dc   = " << dc << ",   matrix: " << dc.get_matrix()    << endl;
    cout << "d    = " << d << ",    matrix: " << d.get_matrix()     << endl;    

    return 0;
}