// Copyright (c) 2006-2009 INRIA Sophia-Antipolis (France).
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
//
//
// Author(s)     : Laurent Rineau, Stephane Tayeb


#ifndef CGAL_MESH_TRIANGULATION_3_H
#define CGAL_MESH_TRIANGULATION_3_H

#include <CGAL/Kernel_traits.h>

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Mesh_3/Robust_weighted_circumcenter_filtered_traits_3.h>

#include <CGAL/Mesh_vertex_base_3.h>
#include <CGAL/Mesh_cell_base_3.h>

namespace CGAL {
  
  namespace details {
    
    template<typename K>
    struct Mesh_geom_traits_generator
    {
    private:
      typedef Robust_weighted_circumcenter_filtered_traits_3<K>
        Geom_traits;
      
    public:
      typedef Geom_traits type;
      typedef type Type;
    };  // end struct Mesh_geom_traits_generator
    
  }  // end namespace details
  
  
// Struct Mesh_triangulation_3
//
template< class MD,
          class K=typename Kernel_traits<MD>::Kernel,
          class GT = typename details::Mesh_geom_traits_generator<K>::type,
          class Cb = CGAL::Regular_triangulation_cell_base_3
          <
            GT, 
            CGAL::Triangulation_cell_base_with_circumcenter_3<GT> 
          >
        >
struct Mesh_triangulation_3
{
private:
  typedef GT                                                    Geom_traits;
  typedef Mesh_vertex_base_3<Geom_traits, MD>                   Vertex_base;
  typedef Mesh_cell_base_3<Geom_traits, MD, Cb, Sequential_tag> Cell_base;
  typedef Triangulation_data_structure_3<Vertex_base,Cell_base> Tds;
  typedef Regular_triangulation_3<Geom_traits, Tds>             Triangulation;

public:
  typedef Triangulation type;
  typedef type Type;
};  // end struct Mesh_triangulation_3

// Struct Parallel_mesh_triangulation_3
//
template< class MD,
          class K=typename Kernel_traits<MD>::Kernel,
          class GT = typename details::Mesh_geom_traits_generator<K>::type,
          class Cb = CGAL::Regular_triangulation_cell_base_3
          <
            GT, 
            CGAL::Triangulation_cell_base_with_circumcenter_3
            <
              GT
#ifdef LINKED_WITH_TBB
              , Triangulation_lazy_ds_cell_base_3<Parallel_tag>
#endif // LINKED_WITH_TBB
            >
          > 
        >
struct Parallel_mesh_triangulation_3
{
private:
  typedef GT                                                    Geom_traits;
  typedef Mesh_vertex_base_3<Geom_traits, MD>                   Vertex_base;

#ifdef LINKED_WITH_TBB
  typedef Mesh_cell_base_3<Geom_traits, MD, Cb, Parallel_tag>   Cell_base;
  typedef Triangulation_data_structure_3<
                            Vertex_base,Cell_base, true>        Tds;
  typedef Regular_triangulation_3<Geom_traits, Tds, true>       Triangulation;

#else // !LINKED_WITH_TBB
  typedef Mesh_cell_base_3<Geom_traits, MD, Cb, Sequential_tag> Cell_base;
  typedef Triangulation_data_structure_3<
                            Vertex_base,Cell_base>              Tds;
  typedef Regular_triangulation_3<Geom_traits, Tds>             Triangulation;
#endif // LINKED_WITH_TBB

public:
  typedef Triangulation type;
  typedef type Type;
};  // end struct Mesh_triangulation_3


}  // end namespace CGAL

#endif // CGAL_MESH_TRIANGULATION_3_H
