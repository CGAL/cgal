// Copyright (c) 1997-2001  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>

#ifndef CGAL_COMPUTE_QUOTIENT_H
#define CGAL_COMPUTE_QUOTIENT_H

#include <CGAL/Quotient.h>
#include <functional>

CGAL_BEGIN_NAMESPACE

template < class NT >
struct Compute_quotient {
    typedef  NT                         argument1_type;
    typedef  NT                         argument2_type;
    typedef  CGAL::Quotient<NT>         result_type;

    result_type  operator () ( const NT& numer, const NT& denom) const
        { return result_type( numer, denom); }
};

CGAL_END_NAMESPACE
  
#endif

// ===== EOF ==================================================================
