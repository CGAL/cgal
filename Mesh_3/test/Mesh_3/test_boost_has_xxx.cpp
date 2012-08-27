// Copyright (c) 2010  GeometryFactory Sarl (France).
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
//
// Author(s)     : Laurent Rineau

//| If Boost MPL cannot support the introspection macro
//| BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF with that compiler, then failure!
//| See http://www.boost.org/libs/mpl/doc/refmanual/cfg-no-has-xxx.html

// One needs that introspection macro, to ensure backward compatibility
// with specification as of CGAL-3.7.

#include <boost/config.hpp>

int main()
{
#ifdef BOOST_MPL_CFG_NO_HAS_XXX
  return 1;
#endif
  return 0;
}
