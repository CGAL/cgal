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

#ifndef CGAL_ACCESS_BY_INDEX_H
#define CGAL_ACCESS_BY_INDEX_H

#include <CGAL/basic.h>
#include <functional>
#include <iterator>

CGAL_BEGIN_NAMESPACE

template < class RndAccIt, bool check_lower = false, bool check_upper = false>
class Access_by_index {
  public:
    typedef  int                        argument_type;
    typedef  typename std::iterator_traits<RndAccIt>::value_type
                                        result_type;

    Access_by_index( RndAccIt it = RndAccIt(),
		     const result_type& default_result = result_type(),
		     int lower = 0, int upper = 0)
        : a( it), r( default_result)
        {
	    if ( check_lower) l = lower;
	    if ( check_upper) u = upper;
	}

    result_type  operator () ( int i) const
        {
	    if ( check_lower && i <  l) return r;
	    if ( check_upper && i >= u) return r;
	    return a[ i];
	}

  private:
    RndAccIt     a;
    int          l;
    int          u;
    result_type  r;
};

CGAL_END_NAMESPACE
  
#endif

// ===== EOF ==================================================================
