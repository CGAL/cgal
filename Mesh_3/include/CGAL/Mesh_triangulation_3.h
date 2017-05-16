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

#include <CGAL/license/Mesh_3.h>


#include <CGAL/Mesh_3/config.h>

#include <CGAL/Kernel_traits.h>

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Robust_weighted_circumcenter_filtered_traits_3.h>

#include <CGAL/Mesh_vertex_base_3.h>
#include <CGAL/Compact_mesh_cell_base_3.h>

namespace CGAL {
  
  namespace details {
    
    template<typename K>
    struct Mesh_geom_traits_generator
    {
    private:
      typedef Robust_weighted_circumcenter_filtered_traits_3<K> Geom_traits;

    public:
      typedef Geom_traits type;
      typedef type Type;
    };  // end struct Mesh_geom_traits_generator
    
  }  // end namespace details
  
  
// Struct Mesh_triangulation_3
//
template<class MD,
         class K_ = Default,
         class Concurrency_tag = Sequential_tag,
         class Vertex_base_ = Default,
         class Cell_base_   = Default>
struct Mesh_triangulation_3;

// Sequential version (default)
template<class MD, class K_, class Concurrency_tag,
         class Vertex_base_, class Cell_base_>
struct Mesh_triangulation_3
{
private:
  typedef typename Default::Get<K_, typename Kernel_traits<MD>::Kernel>::type K;

  typedef typename details::Mesh_geom_traits_generator<K>::type Geom_traits;

  typedef typename Default::Get<
    Vertex_base_, 
    Mesh_vertex_base_3<Geom_traits, MD> >::type                 Vertex_base;
  typedef typename Default::Get<
    Cell_base_, 
    Compact_mesh_cell_base_3<Geom_traits, MD> >::type           Cell_base;

  typedef Triangulation_data_structure_3<Vertex_base,Cell_base> Tds;
  typedef Regular_triangulation_3<Geom_traits, Tds>             Triangulation;

public:
  typedef Triangulation type;
  typedef type Type;
};  // end struct Mesh_triangulation_3

#ifdef CGAL_LINKED_WITH_TBB
// Parallel version (specialization)
//
template<class MD, class K_,
         class Vertex_base_, class Cell_base_>
struct Mesh_triangulation_3<MD, K_, Parallel_tag, Vertex_base_, Cell_base_>
{
private:
  typedef typename Default::Get<K_, typename Kernel_traits<MD>::Kernel>::type K;

  typedef typename details::Mesh_geom_traits_generator<K>::type Geom_traits;

  typedef typename Default::Get<
    Vertex_base_, 
    Mesh_vertex_base_3<Geom_traits, MD> >::type                 Vertex_base;
  typedef typename Default::Get<
    Cell_base_, 
    Compact_mesh_cell_base_3<Geom_traits, MD> >::type           Cell_base;

  typedef Triangulation_data_structure_3<
    Vertex_base, Cell_base, Parallel_tag>                       Tds;
  typedef Regular_triangulation_3<Geom_traits, Tds>             Triangulation;

public:
  typedef Triangulation type;
  typedef type Type;
};  // end struct Mesh_triangulation_3
#endif // CGAL_LINKED_WITH_TBB

}  // end namespace CGAL

#endif // CGAL_MESH_TRIANGULATION_3_H
