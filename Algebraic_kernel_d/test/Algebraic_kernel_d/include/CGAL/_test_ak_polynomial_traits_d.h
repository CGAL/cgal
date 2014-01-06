// Copyright (c) 2006-2013 INRIA Nancy-Grand Est (France). All rights reserved.
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
// $Id$
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
