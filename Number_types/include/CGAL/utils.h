// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id:$
// 
//
// Author(s)     : Sebastian Limbach  <slimbach@mpi-sb.mpg.de>

#ifndef CGAL_UTILS_H
#define CGAL_UTILS_H

#include <CGAL/basic.h>
#include <CGAL/utils_classes.h>

CGAL_BEGIN_NAMESPACE


template< class Number_type >
inline bool is_valid( const Number_type& x ) {
  return Is_valid< Number_type >()( x );
};

CGAL_END_NAMESPACE

#endif // CGAL_UTILS_H
