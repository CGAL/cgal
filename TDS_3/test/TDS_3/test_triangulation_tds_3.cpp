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
// Author(s)     : Francois Rebufat

#include <CGAL/Triangulation_data_structure_3.h>

#include <CGAL/_test_cls_tds_3.h>

typedef CGAL::Triangulation_data_structure_3<>               Tds;

// Explicit instantiation :
template class CGAL::Triangulation_data_structure_3<>;

int main()
{
  _test_cls_tds_3(Tds());
  return 0;
}
