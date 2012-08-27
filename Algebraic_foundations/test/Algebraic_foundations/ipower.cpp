// Copyright (c) 2006-2007 Max-Planck-Institute Saarbruecken (Germany).
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
// $URL: svn+ssh://hemmer@scm.gforge.inria.fr/svn/cgal/trunk/Algebraic_foundations/test/Algebraic_foundations/extended_euclidean_algorithm.cpp $
// $Id: extended_euclidean_algorithm.cpp 47265 2008-12-08 06:26:27Z hemmer $
//
// Author(s)     : Michael Hemmer
//
// ============================================================================

#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/ipower.h>

template<class Integer>
void test_ipower() {
  {
    Integer i = 3;
    assert(CGAL::ipower(i,0)==Integer(1));
    assert(CGAL::ipower(i,1)==Integer(3));
    assert(CGAL::ipower(i,2)==Integer(9));
    assert(CGAL::ipower(i,3)==Integer(27));
    assert(CGAL::ipower(i,4)==Integer(81));
    assert(CGAL::ipower(i,5)==Integer(243));
  }
  {
    Integer i = -3;
    assert(CGAL::ipower(i,0)==Integer(1));
    assert(CGAL::ipower(i,1)==Integer(-3));
    assert(CGAL::ipower(i,2)==Integer(9));
    assert(CGAL::ipower(i,3)==Integer(-27));
    assert(CGAL::ipower(i,4)==Integer(81));
    assert(CGAL::ipower(i,5)==Integer(-243));
  }
}

int main(){
    test_ipower<long>();     
    return 0;
}

// EOF 
