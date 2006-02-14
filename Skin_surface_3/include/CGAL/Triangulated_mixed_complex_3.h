// Copyright (c) 2005 RuG (Netherlands)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$ $Date$
// 
//
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef TRIANGULATED_MIXED_COMPLEX_3
#define TRIANGULATED_MIXED_COMPLEX_3

#include <CGAL/Skin_surface_traits_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulated_mixed_complex_cell_3.h>

CGAL_BEGIN_NAMESPACE

template <
  class SkinSurfaceTraits_3,
  class GT=typename SkinSurfaceTraits_3::Triangulated_mixed_complex_traits,
  class PolyhedronKernel_3=typename SkinSurfaceTraits_3::Polyhedron_traits,
  class Tds = Triangulation_data_structure_3 <
    Triangulation_vertex_base_3<GT>,
    Triangulated_mixed_complex_cell_3<GT, PolyhedronKernel_3> > >
class Triangulated_mixed_complex_3 : public Triangulation_3<GT, Tds> {
};

CGAL_END_NAMESPACE

#endif // TRIANGULATED_MIXED_COMPLEX_3
