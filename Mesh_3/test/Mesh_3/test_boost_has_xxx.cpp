// Copyright (c) 2010  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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
