// Copyright (c) 2006-2007 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Michael Hemmer
//
// ============================================================================

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
