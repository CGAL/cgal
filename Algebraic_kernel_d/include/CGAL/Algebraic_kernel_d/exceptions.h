// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
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
