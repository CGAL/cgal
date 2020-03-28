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


#ifndef CGAL_SQRT_EXTENSION_EIGEN_NUMTRAITS_H
#define CGAL_SQRT_EXTENSION_EIGEN_NUMTRAITS_H

namespace Eigen {
  template<class> struct NumTraits;
  template <class NT,class ROOT, class ACDE_TAG, class FP_TAG>
  struct NumTraits<CGAL::Sqrt_extension<NT, ROOT, ACDE_TAG, FP_TAG> >
  {
    enum {
      IsInteger = 0,
      IsSigned = 1,
      IsComplex = 0,
      RequireInitialization = 1,
      ReadCost = 2*NumTraits<NT>::ReadCost+NumTraits<ROOT>::ReadCost,
      AddCost = 2*NumTraits<NT>::AddCost+NumTraits<ROOT>::ReadCost,
      MulCost = 5*NumTraits<NT>::MulCost+2*NumTraits<NT>::AddCost
    };

    typedef CGAL::Sqrt_extension<NT, ROOT, ACDE_TAG, FP_TAG> Real;
    typedef Real NonInteger;
    typedef Real Nested;
    typedef Real Literal;

    static inline Real epsilon() { return NumTraits<NT>::epsilon(); }
    static inline int digits10() { return NumTraits<NT>::digits10(); }
    static inline Real dummy_precision() { return NumTraits<NT>::dummy_precision(); }
    static inline Real highest() { return NumTraits<NT>::highest(); }
    static inline Real lowest() { return NumTraits<NT>::lowest(); }
    static inline Real infinity() { return NumTraits<NT>::infinity(); }
    static inline Real quiet_NaN() { return NumTraits<NT>::quiet_NaN(); }
  };
}

#endif
