// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mariette Yvinec

#ifndef CGAL_TRIANGULATION_DUMMY_TDS_2_H
#define CGAL_TRIANGULATION_DUMMY_TDS_2_H

#include <CGAL/license/TDS_2.h>


#include <CGAL/config.h>

namespace CGAL {

// Dummy TDS which provides all types that a vertex_base or cell_base can use.
struct Dummy_tds_2 {
  struct Vertex {};
  struct Face {};
  struct Edge {};

  struct Vertex_handle {};
  struct Face_handle {};

  struct Vertex_iterator {};
  struct Face_iterator {};
  struct Edge_iterator {};

  struct Edge_circulator {};
  struct Facet_circulator {};
  struct Face_circulator {};
  struct Vertex_circulator {};

  typedef std::size_t size_type;
};

} //namespace CGAL

#endif // CGAL_TRIANGULATION_DUMMY_TDS_2_H
