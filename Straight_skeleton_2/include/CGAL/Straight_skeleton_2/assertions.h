// Copyright (c) 2007 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>

#ifndef CGAL_STRAIGHT_SKELETON_ASSERTIONS_H
#define CGAL_STRAIGHT_SKELETON_ASSERTIONS_H

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/assertions.h>


namespace CGAL {

namespace {

template<class Handle> inline bool handle_assigned ( Handle const& aH )
{
  Handle null ;
  return aH != null ;
}

}

} // end namespace CGAL

#endif // CGAL_STRAIGHT_SKELETON_ASSERTIONS_H //
// EOF //
