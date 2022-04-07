// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Manuel Caroli

#include "types.h"

#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>

typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<K> GT;
typedef CGAL::Periodic_3_Delaunay_triangulation_3<GT> Triang;

int main(int argc, char* argv[]) {
  Triang T;
  return bench_triang<Triang>(argc,argv,T,true);
}
