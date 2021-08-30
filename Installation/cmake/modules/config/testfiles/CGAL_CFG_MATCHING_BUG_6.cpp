// Copyright (c) 2005
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

//| VC 7.3 does not recognize when an operator in a class
//| redefines the operator with the same signature in a base class
//| It happens with the regular triangulation.
//| No minimal testcase yet

#if (defined _MSC_VER && ! defined __INTEL_COMPILER) || defined __SUNPRO_CC

#error "this should not compile and that is good so"

#endif

int main()
{
  return 0;
}
