// Copyright (c) 2012 Inria Saclay (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author: Marc Glisse

#include <CGAL/int.h>
#include <CGAL/long_long.h>

int main(){
  unsigned int  a = 42;
  (void)CGAL::compare(a,a);
  unsigned long b = 42;
  (void)CGAL::compare(b,b);
#ifdef CGAL_USE_LONG_LONG
  unsigned long long c = 42;
  (void)CGAL::compare(c,c);
#endif
}
