// Copyright (c) 2009 Inria Lorraine (France). All rights reserved.
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
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_TEST_IO_H
#define CGAL_TEST_IO_H

#include <CGAL/basic.h>
#ifdef NDEBUG
#  undef NDEBUG
#  include <cassert>
#  define NDEBUG 1
#endif
#include <sstream>

namespace CGAL{

// construct a _NT from a value of type _CT, output it to a stream, read
// it again and check that the number is the same
template <class _NT,class _CT>
void test_io(_CT x){
        typedef _NT     NT;

        NT a(x);
        std::stringstream ss;
        ss<<CGAL::oformat(a);
        ss>>CGAL::iformat(a);
        assert(a==NT(x));
}

// construct by default, write it to stdout and read it again
// (use with caution: when constructing by default, the value of the
// number may be undefined and the test may fail, what does not mean
// that the i/o functions are incorrect)
template <class _NT>
void test_io(){
        typedef _NT     NT;

        NT a;
        std::stringstream ss;
        ss<<CGAL::oformat(a);
        ss>>CGAL::iformat(a);
        assert(a==NT());
}

// construct an interval of type _NT from a value of type _CT, output it to
// a stream, read it again and check that the bounds of the number
// (obtained with member functions inf() and sup()) are the same
template <class _NT,class _CT>
void test_interval_io(_CT x){
        typedef _NT     NT;

        NT a(x),b(x);
        std::stringstream ss;
        ss<<CGAL::oformat(a);
        ss>>CGAL::iformat(a);
        assert(a.inf()==b.inf());
        assert(a.sup()==b.sup());
}

} // namespace CGAL

#endif

// vim: tabstop=8: softtabstop=8: smarttab: shiftwidth=8: expandtab
