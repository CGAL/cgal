// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// $Id$
// 
//
// Author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>
#ifndef CGAL_TD_TRAITS_H
#include <CGAL/Arr_point_location/Td_traits.h>
#endif

namespace CGAL {


/*
template <class Traits,class Arrangement_on_surface_2>
typename Td_traits<Traits,Arrangement_on_surface_2>::Td_map_item
Td_traits<Traits,Arrangement_on_surface_2>::m_blank ;
*/
template <class Traits,class Arrangement_on_surface_2>
typename Td_traits<Traits,Arrangement_on_surface_2>::Vertex_const_handle
Td_traits<Traits,Arrangement_on_surface_2>::m_empty_vtx_handle = Vertex_const_handle();

template <class Traits,class Arrangement_on_surface_2>
typename Td_traits<Traits,Arrangement_on_surface_2>::Halfedge_const_handle
Td_traits<Traits,Arrangement_on_surface_2>::m_empty_he_handle = Halfedge_const_handle();

//template <class Traits,class Arrangement_on_surface_2>
//typename Td_traits<Traits,Arrangement_on_surface_2>::Td_map_item
//Td_traits<Traits,Arrangement_on_surface_2>::m_empty_map_item = Td_map_item();

} //namespace CGAL
