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

#ifndef CGAL_IDENTITY_H
#define CGAL_IDENTITY_H

#include <CGAL/basic.h>
#include <functional>

CGAL_BEGIN_NAMESPACE

template < class T >
class identity
    : public CGAL_STD::unary_function<T,T> {
  public:
    T  operator () ( T t) const { return t; }
};

CGAL_END_NAMESPACE
  
#endif // CGAL_IDENTITY_H

// ===== EOF ==================================================================
