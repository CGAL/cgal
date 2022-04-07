// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Marc Glisse

#ifndef CGAL_VECTOR_VECTOR_H
#define CGAL_VECTOR_VECTOR_H
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <CGAL/Dimension.h>
#include <CGAL/NewKernel_d/utils.h>
#include <vector>
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
namespace CGAL {

//Derive from a class that doesn't depend on Dim, or still use Dim for checking?
template<class NT_,class Dim_,class Max_dim_=Dim_> struct Vector_vector {
        typedef NT_ NT;
        typedef Dim_ Dimension;
        typedef Max_dim_ Max_dimension;
        typedef std::vector<NT> Vector;
        template< class D2, class D3=D2 >
        struct Rebind_dimension {
          typedef Vector_vector< NT, D2, D3 > Other;
        };
        template<class> struct Property : boost::false_type {};

        struct Construct_vector {
                struct Dimension {
                        Vector operator()(int d) const {
                                return Vector(d);
                        }
                };

                struct Iterator {
                        template<typename Iter>
                        Vector operator()(int CGAL_assertion_code(d),Iter const& f,Iter const& e) const {
                                        CGAL_assertion(d==std::distance(f,e));
                                        return Vector(f,e);
                                }
                };

                // unneeded thanks to Iterator_and_last?
#if 0
                struct Iterator_add_one {
                        template<typename Iter>
                        Vector operator()(int  CGAL_assertion_code(d),Iter const& f,Iter const& e) const {
                                        CGAL_assertion(d==std::distance(f,e)+1);
                                        Vector a;
                                        a.reserve(d+1);
                                        a.insert(a.end(),f,e);
                                        a.push_back(1);
                                        return a;
                                }
                };
#endif

                struct Iterator_and_last {
                        template<typename Iter,typename T>
                                Vector operator()(int d,Iter const& f,Iter const& e,T&& t) const {
                                        CGAL_assertion(d==std::distance(f,e)+1);
                                        Vector a;
                                        a.reserve(d+1);
                                        a.insert(a.end(),f,e);
                                        a.push_back(std::forward<T>(t));
                                        return a;
                                }
                };

                // useless, use a transform_iterator?
#if 0
                struct Iterator_and_last_divide {
                        template<typename Iter,typename T>
                                Vector operator()(int d,Iter f,Iter const& e,T const&t) const {
                                        CGAL_assertion(d==std::distance(f,e)+1);
                                        Vector a;
                                        a.reserve(d+1);
                                        for(;f!=e;++f){
                                                a.push_back(*f/t);
                                        }
                                        return a;
                                }
                };
#endif

                struct Values {
                        template<class...U>
                                Vector operator()(U&&...u) const {
                                        //TODO: check the right number of {}, g++ accepts one and two
                                        Vector a={forward_safe<NT,U>(u)...};
                                        return a;
                                }
                };

                struct Values_divide {
                        template<class H,class...U>
                                Vector operator()(H const&h,U&&...u) const {
                                        //TODO: do we want to cast at some point?
                                        //e.g. to avoid 1/2 in integers
                                        // ==> use Rational_traits<NT>().make_rational(x,y) ?
                                        Vector a={Rational_traits<NT>().make_rational(std::forward<U>(u),h)...};
                                        return a;
                                }
                };
        };
        typedef typename Vector::const_iterator Vector_const_iterator;
        static Vector_const_iterator vector_begin(Vector const&a){
                return a.begin();
        }
        static Vector_const_iterator vector_end(Vector const&a){
                return a.end();
        }
        static int size_of_vector(Vector const&a){
                return (int)a.size();
        }
};


}
#endif

