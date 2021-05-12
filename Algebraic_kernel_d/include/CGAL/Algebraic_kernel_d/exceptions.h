// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================



#ifndef CGAL_ALGEBRAIC_KERNEL_EXCEPTIONS_H
#define CGAL_ALGEBRAIC_KERNEL_EXCEPTIONS_H


namespace CGAL {

  namespace internal {

    /*!
     * \brief Exception class for not sufficiently generic positions.
     *
     * Must be thrown whenever a curve cannot be analysed because its position
     * is not "good enough".
     */
    class Non_generic_position_exception {

    public:

      //! Default constructible
      Non_generic_position_exception() {}

    };

    /*!
     * \brief Exception class for not sufficiently generic positions.
     *
     * Must be thrown whenever a curve cannot be analysed because its position
     * is not "good enough".
     */
    template<typename Polynomial>
    class Zero_resultant_exception {

      Polynomial curve1,curve2;
      bool one_curve_failure;

    public:

      Zero_resultant_exception(Polynomial c)
        : curve1(c), curve2(c),one_curve_failure(true)
        {}

      Zero_resultant_exception(Polynomial c1,Polynomial c2)
        : curve1(c1),curve2(c2),one_curve_failure(false)
        {}

    };

  } // namespace internal

} //namespace CGAL


#endif
