// Copyright (c) 2019 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Jane Tournois


#ifndef CGAL_TETRAHEDRAL_REMESHING_TRIANGULATION_H
#define CGAL_TETRAHEDRAL_REMESHING_TRIANGULATION_H

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Default.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_cell_base.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_vertex_base.h>


namespace CGAL
{
namespace Tetrahedral_remeshing
{
  template<typename K,
           typename Info,
           typename Input_tr = CGAL::Default,
           typename Cb = CGAL::Triangulation_cell_base_3<K> >
  class Remeshing_triangulation_3
    : public CGAL::Triangulation_3<K,
        CGAL::Triangulation_data_structure_3<
          Remeshing_vertex_base<K>,
          Remeshing_cell_base<K,
            Info,
            typename CGAL::Default::Get<Input_tr, CGAL::Triangulation_3<K> >::type::Cell,
            Cb
          >
        >
      >
  {
    typedef Remeshing_vertex_base<K>                       RVb;

    typedef typename CGAL::Default::Get<Input_tr,
      CGAL::Triangulation_3<K> >::type                     Input_tr;
    typedef typename Input_tr::Cell                        Input_cell;
    typedef Remeshing_cell_base<K, Info, Input_cell, Cb>   RCb;

  public:
    typedef CGAL::Triangulation_data_structure_3<RVb, RCb> Tds;
    typedef CGAL::Triangulation_3<K, Tds>                  Self;
    typedef Self                                           type;
  };


}//end namespace Tetrahedral_remeshing
}//end namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_TRIANGULATION_H
