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

#ifndef CGAL_LA_EIGEN_CONSTRUCTORS_H
#define CGAL_LA_EIGEN_CONSTRUCTORS_H
#include <CGAL/config.h>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4003) // not enough actual parameters for macro 'BOOST_PP_EXPAND_I'
                                // https://lists.boost.org/boost-users/2014/11/83291.php
#endif

#ifndef CGAL_EIGEN3_ENABLED
#error Requires Eigen
#endif
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <CGAL/Dimension.h>
#include <Eigen/Dense>
#include <CGAL/iterator_from_indices.h>
#include <CGAL/NewKernel_d/utils.h>
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>

namespace CGAL {
  template <class Vector_> struct Construct_eigen {
    typedef Vector_ result_type;
    typedef typename Vector_::Scalar NT;

    private:
    static void check_dim(int CGAL_assertion_code(d)){
      CGAL_assertion_code(int m = result_type::MaxSizeAtCompileTime;)
      CGAL_assertion((m == Eigen::Dynamic) || (d <= m));
    }
    public:

    struct Dimension {
      // Initialize with NaN if possible?
      result_type operator()(int d) const {
        check_dim(d);
        return result_type(d);
      }
    };

    struct Iterator {
      template<typename Iter>
        result_type operator()(int d,Iter const& f,Iter const& e) const {
          check_dim(d);
          CGAL_assertion(d==std::distance(f,e));
          result_type a(d);
          // TODO: check the right way to do this
          std::copy(f,e,&a[0]);
          return a;
        }
    };

#if 0
    struct Iterator_add_one {
      template<typename Iter>
        result_type operator()(int d,Iter const& f,Iter const& e) const {
          check_dim(d);
          CGAL_assertion(d==std::distance(f,e)+1);
          result_type a(d);
          std::copy(f,e,&a[0]);
          a[d-1]=1;
          return a;
        }
    };
#endif

    struct Iterator_and_last {
      template<typename Iter,typename T>
        result_type operator()(int d,Iter const& f,Iter const& e,T&& t) const {
          check_dim(d);
          CGAL_assertion(d==std::distance(f,e)+1);
          result_type a(d);
          std::copy(f,e,&a[0]);
          a[d-1]=std::forward<T>(t);
          return a;
        }
    };

    struct Initializer_list {
      // Fix T==NT?
      template<class T>
        result_type operator()(std::initializer_list<T> l) const {
          return Iterator()(l.size(),l.begin(),l.end());
        }
    };

    struct Values {
      // TODO avoid going through Initializer_list which may cause extra copies. Possibly use forward_as_tuple.
      template<class...U>
        result_type operator()(U&&...u) const {
          check_dim(sizeof...(U)); // TODO: use static_assert
          return Initializer_list()({forward_safe<NT,U>(u)...});
        }
    };

    struct Values_divide {
      template<class H,class...U>
        result_type operator()(H const&h,U&&...u) const {
          check_dim(sizeof...(U)); // TODO: use static_assert
          return Initializer_list()({Rational_traits<NT>().make_rational(std::forward<U>(u),h)...});
        }
    };
  };
}
#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif
