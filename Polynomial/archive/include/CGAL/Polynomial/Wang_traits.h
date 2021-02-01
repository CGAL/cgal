// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hemmer
//
// ============================================================================

#ifndef CGAL_POLYNOMIAL_WANG_TRAITS_H
#define CGAL_POLYNOMIAL_WANG_TRAITS_H 1

#include <CGAL/basic.h>
#include <CGAL/Polynomial/wang.h>

/*! \file CGAL/Polynomial/Wang_traits.h
 *  \brief definition of traits class CGAL::Wang_traits.
 */

namespace CGAL{

// fwd
template <class A > class Polynomial;

namespace internal{

/*! \nosubgrouping
 *  \brief traits class for rational reconstrcution based on wangs
 *  algorithm
 *
 *  This is experimental, and should serve as a design study, i.e.,
 *  It may be joint with Scalar_factor_traits.
 *
 *  This is the default implementation of CGAL::Wang_traits.
 *  It is valid for scalar types beeing a EuclideanRing, e.g.,  Integer
 */
template <class NT_>
class Wang_traits {
public:
    // the supported number type
    typedef NT_ NT;
    // NT is also
    typedef NT Scalar;

    struct Wang {
        bool
        operator()
            (const NT& u, const Scalar& m, NT& n, Scalar& d) const {
            n = d = NT(0);
            return CGAL::internal::wang(u,m,n,d);
        }
    };
};


template <class AS >
class Wang_traits< Polynomial<AS> >{

    typedef Wang_traits<AS> WT;
public:
    // the supported number type
    typedef Polynomial<AS> NT;
    // the scalar type (same as Scalar factor traits ?)
    typedef typename WT::Scalar Scalar;

    struct Wang {
        bool operator()
            (const NT& p, const Scalar& m, NT& result_n, Scalar& result_d) const {
            typename Algebraic_structure_traits<Scalar>::Integral_division idiv;
            typename Algebraic_structure_traits<Scalar>::Gcd          gcd;
            typename WT::Wang wang;

            result_n = NT(0);
            result_d = Scalar(0);
//            std::cout<<"Poly "<<p<<" m "<<m<<std::endl;
            const int d = p.degree();
            std::vector<AS>     nums(d+1);
            std::vector<Scalar> denoms(d+1);
            for (int i = 0; i <= d; i++) {
//                bool w = wang(p[i], m, nums[i], denoms[i]);
//                wang(p[i], m, nums[i], denoms[i]);
//                std::cout<<i<<" "<<p[i]<<" "<<w<<std::endl;
              if(!wang(p[i], m, nums[i], denoms[i])) return false;
//              if(!w) return false;    !!!!!!
            }

            // c = lcm(denoms[0], ..., denoms[d])
            result_d = denoms[0];
            for (int i = 1; i <= d; i++) {
                result_d *= idiv(denoms[i], gcd(result_d, denoms[i]));
            }

            // expand each (nums[i], denoms[i]) pair to common denominator
            for (int i = 0; i <= d; i++) {
                nums[i] *= AS(idiv(result_d, denoms[i]));
            }
            result_n = NT(nums.begin(),nums.end());
            return true;
        }
    };
};

} // namespace internal
} // namespace CGAL

#endif // CGAL_POLYNOMIAL_WANG_TRAITS_H
// EOF
