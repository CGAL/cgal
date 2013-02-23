// Copyright (c) 2013
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
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
//
//
// Author(s)     : Sebastien Loriot <sebastien.loriot@cgal.org>


#ifndef CGAL_RESULT_OF_H
#define CGAL_RESULT_OF_H

#include <boost/utility/result_of.hpp>
#include <boost/version.hpp>

namespace CGAL{

// Even if for now we use the tr1 result_of implementation, we use the cpp11
// namespace since in the future, that's the decltype version that will be used
namespace cpp11{

template<typename F>
struct result_of
{
  // from boost 1.44 release note http://www.boost.org/users/history/version_1_44_0.html :
  //    New template boost::tr1_result_of that implements the TR1 ResultOf protocol even if boost::result_of uses the C++0x decltype-based implementation.
  #if BOOST_VERSION < 104400
  typedef typename boost::result_of<F>::type type;
  #else
  typedef typename boost::tr1_result_of<F>::type type;
  #endif
};

}

}

#endif //CGAL_RESULT_OF_H
