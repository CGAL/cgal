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

#include <CGAL/number_type_basic.h>

CGAL_BEGIN_NAMESPACE

template< class T >
inline bool is_valid( const T& x ) {
  return Is_valid< T >()( x );
};

CGAL_END_NAMESPACE

#endif // CGAL_UTILS_H
