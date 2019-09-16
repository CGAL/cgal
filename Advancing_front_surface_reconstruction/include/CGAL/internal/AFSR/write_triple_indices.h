// Copyright (c) 2015  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Frank Da, David Cohen-Steiner, Andreas Fabri

#ifndef CGAL_AFSR_WRITE_TRIPLE_INDICES_H
#define CGAL_AFSR_WRITE_TRIPLE_INDICES_H

#include <CGAL/license/Advancing_front_surface_reconstruction.h>


#include <CGAL/array.h>

namespace CGAL {

  template <class Triangulation, class Filter>
  class Advancing_front_surface_reconstruction;



  template <class OutputIterator, class Triangulation, class Filter>
  OutputIterator
  write_triple_indices(OutputIterator out, const Advancing_front_surface_reconstruction<Triangulation,Filter>& S)
  {
    typedef Advancing_front_surface_reconstruction<Triangulation,Filter> Surface;
    typedef typename Surface::TDS_2 TDS_2;
    typedef typename TDS_2::Face_iterator Face_iterator;

    if(S.triangulation_3().dimension() < 3){
      std::cerr << "not 3D\n";
      return out;
    }
    const TDS_2& tds = S.triangulation_data_structure_2();

    for(Face_iterator fit = tds.faces_begin(); fit != tds.faces_end(); ++fit){

      if(fit->is_on_surface()){
        *out++ = CGAL::make_array(fit->vertex(0)->vertex_3()->id(),
                                  fit->vertex(1)->vertex_3()->id(),
                                  fit->vertex(2)->vertex_3()->id());
      }
    }
    return out;
  }

}


#endif
