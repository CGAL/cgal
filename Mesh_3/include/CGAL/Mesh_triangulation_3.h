// Copyright (c) 2006-2009 INRIA Sophia-Antipolis (France).
// Copyright (c) 2011      GeometryFactory Sarl (France)
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

#ifdef CGAL_COMPACT_MESH_VERTEX_CELL
#include <CGAL/Compact_mesh_vertex_base_3.h>
#include <CGAL/Compact_mesh_cell_base_3.h>
#endif // CGAL_COMPACT_MESH_VERTEX_CELL

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
template<class MD, 
         class Concurrency_tag = Sequential_tag,
         class K=typename Kernel_traits<MD>::Kernel>
struct Mesh_triangulation_3;

// Sequential version (default)
template<class MD, class Concurrency_tag, class K>
struct Mesh_triangulation_3
{
private:
  typedef typename details::Mesh_geom_traits_generator<K>::type Geom_traits;

#ifdef CGAL_COMPACT_MESH_VERTEX_CELL
  typedef Compact_mesh_vertex_base_3<Geom_traits, MD>           Vertex_base;
  typedef Compact_mesh_cell_base_3<Geom_traits, MD>             Cell_base;
#else // NOT CGAL_COMPACT_MESH_VERTEX_CELL
  typedef Mesh_vertex_base_3<Geom_traits, MD>                   Vertex_base;
  typedef Mesh_cell_base_3<Geom_traits, MD>                     Cell_base;
#endif // NOT CGAL_COMPACT_MESH_VERTEX_CELL

#if defined(CGAL_MESH_3_USE_LAZY_SORTED_REFINEMENT_QUEUE)\
 || defined(CGAL_MESH_3_USE_LAZY_UNSORTED_REFINEMENT_QUEUE)
  typedef Triangulation_data_structure_3<
    Vertex_base, Cell_base,
    Compact_container_strategy_with_counter, 
    Compact_container_strategy_with_counter>                    Tds;
#else
  typedef Triangulation_data_structure_3<Vertex_base,Cell_base> Tds;
#endif
  typedef Regular_triangulation_3<Geom_traits, Tds>             Triangulation;

public:
  typedef Triangulation type;
  typedef type Type;
};  // end struct Mesh_triangulation_3

#ifdef CGAL_LINKED_WITH_TBB
// Parallel version (specialization)
//
template<class MD, class K>
struct Mesh_triangulation_3<MD, Parallel_tag, K>
{
private:
  typedef typename details::Mesh_geom_traits_generator<K>::type Geom_traits;

# ifdef CGAL_COMPACT_MESH_VERTEX_CELL
  typedef Compact_mesh_vertex_base_3<Geom_traits, MD>           Vertex_base;
  typedef Compact_mesh_cell_base_3<Geom_traits,MD>              Cell_base;
# else // NOT CGAL_COMPACT_MESH_VERTEX_CELL
  typedef Mesh_vertex_base_3<Geom_traits, MD>                   Vertex_base;
  typedef Mesh_cell_base_3<Geom_traits, MD>                     Cell_base;
# endif // NOT CGAL_COMPACT_MESH_VERTEX_CELL

  typedef Triangulation_data_structure_3<
    Vertex_base, Cell_base, 
    Compact_container_strategy_with_counter, 
    Compact_container_strategy_with_counter,
    Parallel_tag>                                               Tds;
  typedef Regular_triangulation_3<Geom_traits, Tds>             Triangulation;

public:
  typedef Triangulation type;
  typedef type Type;
};  // end struct Mesh_triangulation_3
#endif // CGAL_LINKED_WITH_TBB

}  // end namespace CGAL

#endif // CGAL_MESH_TRIANGULATION_3_H
