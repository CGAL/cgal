// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : CGAL
// File          : include/CGAL/Wang_traits.h
// CGAL_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Michael Hemmer   <mhemmer@uni-mainz.de>
//
// ============================================================================

#ifndef CGAL_WANG_TRAITS_H
#define CGAL_WANG_TRAITS_H 1

#include <CGAL/basic.h>

/*! \file CGAL/Wang_traits.h
 *  \brief Definition of traits class CGAL::Wang_traits.
 */

namespace CGAL{
// fwd 
template <class A, class B> class Sqrt_extension;
} //namespace CGAL
namespace CGAL {
template <class A > class Polynomial;

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
            return CGAL::wang(u,m,n,d);
        }
    };
};

template <class AS, class ROOT>
class Wang_traits< CGAL::Sqrt_extension<AS,ROOT> >{
    typedef Wang_traits<AS> WT; 
public:
    // the supported number type
    typedef  CGAL::Sqrt_extension<AS,ROOT> NT;
    // the scalar type (same as Scalar factor traits ?) 
    typedef typename WT::Scalar Scalar;

    struct Wang { 
        bool 
        operator()
            (const NT& ext, const Scalar& m, NT& n, Scalar& d) const {
            typename Algebraic_structure_traits<Scalar>::Integral_division idiv;
            typename WT::Wang wang; 
            
            AS     a0,a1;
            Scalar d0,d1; 
            ROOT root;
            n = NT(0);
            d = Scalar(0);
            
            if(!wang(ext.a0(),m,a0,d0)) return false;
            
            if(ext.is_extended()){
                if(!wang(ext.a1(),m,a1,d1)) return false;
                d  = d0 * idiv(d1,CGAL::gcd(d0,d1));
                a0 = a0 * idiv(d,d0);
                a1 = a1 * idiv(d,d1);
                n  = NT(a0,a1,ext.root());
            }else{
                d = d0; 
                n = NT(a0);
            }
            return true; 
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



} // namespace CGAL

#endif // CGAL_WANG_TRAITS_H
// EOF
