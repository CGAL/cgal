// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_APOLLONIUS_GRAPH_DATA_STRUCTURE_2_H
#define CGAL_APOLLONIUS_GRAPH_DATA_STRUCTURE_2_H

#include <CGAL/license/Apollonius_graph_2.h>


#include <CGAL/Triangulation_data_structure_2.h>


namespace CGAL {

// For backward compatibility

template <class Vb, class Fb>
class Apollonius_graph_data_structure_2
  : public Triangulation_data_structure_2<Vb, Fb>
{};


} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_DATA_STRUCTURE_2_H
