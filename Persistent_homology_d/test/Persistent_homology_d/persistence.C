// Copyright (c) 2003,2004  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Lutz Kettner, Afra Zomorodian

#include <CGAL/basic.h>
#include <iostream>
#include <gmpxx.h>
#include <CGAL/Persistent_homology_d/persistence.h>

typedef int Coeff;

// Define example filtration: simplices 0 to 10; we store their boundary

int filtration[11][4] = {
    {-1},
    {-1},
    { 0, 1, -1},
    { -1},
    { -1},
    { 0, 3, -1},
    { 3, 4, -1},
    { 1, 4, -1},
    { 1, 3, -1},
    { 6, 7, 8, -1},
    { 2, 5, 8, -1}
};

class Simple_boundary_operator {
public:
    typedef CGAL::Persistent_homology_d::Boundary_index_tag  rep_tag;
    typedef int                                 (*Simplex_handle)[4];
    typedef std::pair<int, Coeff>                         value_type;
    template <class OutputIterator>
    OutputIterator operator()(Simplex_handle h, OutputIterator out) {
        int* sh = *h;
        int sign = 1;
        while(*sh != -1) {
            value_type p(*sh++, Coeff(sign));
            *out++ = p;
            sign = -sign;
        }
        return out;
    }
};

int main() {
    CGAL::Persistent_homology_d::persistence( filtration, filtration + 11, 
	       CGAL::Persistent_homology_d::Z2_coefficients(),
	       CGAL::Persistent_homology_d::Our_storage(), 
	       CGAL::Persistent_homology_d::Pairs_only(),
	       Simple_boundary_operator());
    return 0;
}

