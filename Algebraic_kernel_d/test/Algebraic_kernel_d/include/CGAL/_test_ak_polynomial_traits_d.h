// Copyright (c) 2006-2013 INRIA Nancy-Grand Est (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-2.1-only
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#include <CGAL/Test/_test_polynomial_traits_d.h>

namespace CGAL{

// This structure is empty because it will be partially specialized.
template <class Is_exact,class PT>
struct test_ak_polynomial_traits_d{};

template <class PT>
struct test_ak_polynomial_traits_d<CGAL::Tag_true,PT>{
        void operator()(const PT& traits){
                return test_polynomial_traits_d(traits);
        }
};

template <class PT>
struct test_ak_polynomial_traits_d<CGAL::Tag_false,PT>{
        void operator()(const PT& /* traits */){
                std::cout<<
                        "\nATTENTION: not testing inexact polynomial traits"<<
                        std::endl;
                return;
        }
};

} // namespace CGAL
