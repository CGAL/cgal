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

#ifndef CGAL_VECTOR_ARRAY_H
#define CGAL_VECTOR_ARRAY_H
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <CGAL/Dimension.h>
#include <CGAL/NewKernel_d/utils.h>
#include <CGAL/array.h>
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/repetition/enum.hpp>

#include <CGAL/NewKernel_d/Vector/determinant_of_points_from_vectors.h>
#include <CGAL/NewKernel_d/Vector/determinant_of_vectors_small_dim.h>
#include <CGAL/NewKernel_d/Vector/determinant_of_iterator_to_points_from_iterator_to_vectors.h>
#include <CGAL/NewKernel_d/Vector/determinant_of_iterator_to_points_from_points.h>
#include <CGAL/NewKernel_d/Vector/determinant_of_iterator_to_vectors_from_vectors.h>



namespace CGAL {

// May not be safe to use with dim!=max_dim.
// In that case, we should store the real dim next to the array.
template<class NT_,class Dim_,class Max_dim_=Dim_> struct Array_vector {
        typedef NT_ NT;
        typedef Dim_ Dimension;
        typedef Max_dim_ Max_dimension;
        template< class D2, class D3=D2 >
        struct Rebind_dimension {
          typedef Array_vector< NT, D2, D3 > Other;
        };
        template<class> struct Property : boost::false_type {};

        static const unsigned d_=Max_dim_::value;
        CGAL_static_assertion(d_ != (unsigned)UNKNOWN_DIMENSION);

        typedef std::array<NT,d_> Vector;
        struct Construct_vector {
                struct Dimension {
                        // Initialize with NaN if possible?
                  Vector operator()(unsigned CGAL_assertion_code(d)) const {
                                CGAL_assertion(d<=d_);
                                return Vector();
                        }
                };

                struct Iterator {
                        template<typename Iter>
                                Vector operator()(unsigned CGAL_assertion_code(d),Iter const& f,Iter const& e) const {
                                        CGAL_assertion(d==(unsigned) std::distance(f,e));
                                        CGAL_assertion(d<=d_);
                                        //TODO: optimize for forward iterators
                                        Vector a;
                                        CGAL_assume(f!=e);
                                        std::copy(f,e,a.begin());
                                        return a;
                                }
                };

#if 0
                struct Iterator_add_one {
                        template<typename Iter>
                                Vector operator()(unsigned d,Iter const& f,Iter const& e) const {
                                        CGAL_assertion(d==std::distance(f,e)+1);
                                        CGAL_assertion(d<=d_);
                                        //TODO: optimize
                                        Vector a;
                                        std::copy(f,e,a.begin());
                                        a.back()=1;
                                        return a;
                                }
                };
#endif

                struct Iterator_and_last {
                        template<typename Iter,typename T>
                        Vector operator()(unsigned  CGAL_assertion_code(d),Iter const& f,Iter const& e,T&& t) const {
                                        CGAL_assertion(d==std::distance(f,e)+1);
                                        CGAL_assertion(d<=d_);
                                        //TODO: optimize for forward iterators
                                        Vector a;
                                        std::copy(f,e,a.begin());
                                        a.back()=std::forward<T>(t);
                                        return a;
                                }
                };

                struct Values {
                        template<class...U>
                                Vector operator()(U&&...u) const {
                                        static_assert(sizeof...(U)<=d_,"too many arguments");
                                        Vector a={{forward_safe<NT,U>(u)...}};
                                        return a;
                                }
                };

                struct Values_divide {
                        template<class H,class...U>
                                Vector operator()(H const& h,U&&...u) const {
                                        static_assert(sizeof...(U)<=d_,"too many arguments");
                                        Vector a={{Rational_traits<NT>().make_rational(std::forward<U>(u),h)...}};
                                        return a;
                                }
                };
        };

        typedef NT const* Vector_const_iterator;
        static Vector_const_iterator vector_begin(Vector const&a){
                return &a[0];
        }
        static Vector_const_iterator vector_end(Vector const&a){
                return &a[0]+d_; // Don't know the real size
        }
        static unsigned size_of_vector(Vector const&){
                return d_; // Don't know the real size
        }

};
// Do not try to instantiate the above
template<class NT_,class Max_dim_> struct Array_vector<NT_,Dynamic_dimension_tag,Max_dim_> {};

}
#endif
