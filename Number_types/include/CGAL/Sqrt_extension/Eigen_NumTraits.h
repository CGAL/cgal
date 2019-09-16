// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
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
    typedef CGAL::Sqrt_extension<NT, ROOT, ACDE_TAG, FP_TAG> Real;
    typedef Real NonInteger;
    typedef Real Nested;
    typedef Real Literal;

    static inline Real epsilon() { return NumTraits<NT>::epsilon(); }

    enum {
      IsInteger = 0,
      IsSigned = 1,
      IsComplex = 0,
      RequireInitialization = 1,
      ReadCost = 2*NumTraits<NT>::ReadCost+NumTraits<ROOT>::ReadCost,
      AddCost = 2*NumTraits<NT>::AddCost+NumTraits<ROOT>::ReadCost,
      MulCost = 5*NumTraits<NT>::MulCost+2*NumTraits<NT>::AddCost
    };
  };
}

#endif
