// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>


#ifndef CGAL_SQRT_EXTENSION_WANG_TRAITS_H
#define CGAL_SQRT_EXTENSION_WANG_TRAITS_H


#include <CGAL/basic.h>
#include <CGAL/Sqrt_extension/Sqrt_extension_type.h>

namespace CGAL {
namespace internal{

template <class NT_> class Wang_traits;

template <class AS, class ROOT, class ACDE_TAG, class FP_TAG>
class Wang_traits< CGAL::Sqrt_extension<AS,ROOT,ACDE_TAG,FP_TAG> >{
    typedef Wang_traits<AS> WT;
public:
    // the supported number type
  typedef  CGAL::Sqrt_extension<AS,ROOT,ACDE_TAG,FP_TAG> NT;
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


} // namespace internal
} //namespace CGAL

#endif // CGAL_SQRT_EXTENSION_WANG_TRAITS_H
