// Copyright (c) 2009 Inria Lorraine (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_TEST_IO_H
#define CGAL_TEST_IO_H

#include <CGAL/config.h>
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
        ss<<CGAL::IO::oformat(a);
        ss>>CGAL::IO::iformat(a);
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
        ss<<CGAL::IO::oformat(a);
        ss>>CGAL::IO::iformat(a);
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
        ss<<CGAL::IO::oformat(a);
        ss>>CGAL::IO::iformat(a);
        assert(a.inf()==b.inf());
        assert(a.sup()==b.sup());
}

} // namespace CGAL

#endif

// vim: tabstop=8: softtabstop=8: smarttab: shiftwidth=8: expandtab
