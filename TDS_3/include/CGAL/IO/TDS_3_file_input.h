#ifndef TDS_3_FILE_INPUT_H
#define TDS_3_FILE_INPUT_H
// Copyright (c) 1997-2010  INRIA Sophia-Antipolis (France).
// Copyright (c) 2011, 2020 GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau, Maxime Gimeno

// Adapted from operator>>(std::istream&, Triangulation_3&) from
// <CGAL/Triangulation_3.h>

#include <CGAL/license/TDS_3.h>


#include <CGAL/internal/Tr_or_tds_file_input.h>

namespace CGAL {

template <typename Tr_src,
          typename Tr_tgt,
          typename ConvertVertex,
          typename ConvertCell>
std::istream& file_input(std::istream& is, Tr_tgt &tr,
                         ConvertVertex convert_vertex = ConvertVertex(),
                         ConvertCell convert_cell = ConvertCell())
{
  return internal::file_input(is, tr, true, convert_vertex, convert_cell);
}

} //end CGAL
#endif // TDS_3_FILE_INPUT_H
