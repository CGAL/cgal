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

template <class Traits,class Arrangement_on_surface_2>
typename Td_traits<Traits,Arrangement_on_surface_2>::Vertex_const_handle
Td_traits<Traits,Arrangement_on_surface_2>::vtx_at_left_infinity()
{
  //  static Vertex_const_handle* m_p_vtx_at_lt_inf;
  if (!m_p_vtx_at_lt_inf)
    m_p_vtx_at_lt_inf = new Vertex_const_handle();
  return *m_p_vtx_at_lt_inf;
}

template <class Traits,class Arrangement_on_surface_2>
typename Td_traits<Traits,Arrangement_on_surface_2>::Vertex_const_handle
Td_traits<Traits,Arrangement_on_surface_2>::vtx_at_right_infinity()
{
  //  static Vertex_const_handle* m_p_vtx_at_rt_inf;
  if (!m_p_vtx_at_rt_inf)
    m_p_vtx_at_rt_inf = new Vertex_const_handle();
  return *m_p_vtx_at_rt_inf;
}

template <class Traits,class Arrangement_on_surface_2>
typename Td_traits<Traits,Arrangement_on_surface_2>::Halfedge_const_handle
Td_traits<Traits,Arrangement_on_surface_2>::he_at_bottom_infinity()
{
  //  static Halfedge_const_handle* m_p_he_at_btm_inf;
  if (!m_p_he_at_btm_inf)
    m_p_he_at_btm_inf = new Halfedge_const_handle();
  return *m_p_he_at_btm_inf;
}

template <class Traits,class Arrangement_on_surface_2>
typename Td_traits<Traits,Arrangement_on_surface_2>::Halfedge_const_handle
Td_traits<Traits,Arrangement_on_surface_2>::he_at_top_infinity()
{
  //  static Halfedge_const_handle* m_p_he_at_top_inf;
  if (!m_p_he_at_top_inf)
    m_p_he_at_top_inf = new Halfedge_const_handle();
  return *m_p_he_at_top_inf;
}

template <class Traits,class Arrangement_on_surface_2>
typename Td_traits<Traits,Arrangement_on_surface_2>::Vertex_const_handle *
Td_traits<Traits,Arrangement_on_surface_2>::m_p_vtx_at_lt_inf = 0;

template <class Traits,class Arrangement_on_surface_2>
typename Td_traits<Traits,Arrangement_on_surface_2>::Vertex_const_handle *
Td_traits<Traits,Arrangement_on_surface_2>::m_p_vtx_at_rt_inf = 0;

template <class Traits,class Arrangement_on_surface_2>
typename Td_traits<Traits,Arrangement_on_surface_2>::Halfedge_const_handle *
Td_traits<Traits,Arrangement_on_surface_2>::m_p_he_at_btm_inf = 0;

template <class Traits,class Arrangement_on_surface_2>
typename Td_traits<Traits,Arrangement_on_surface_2>::Halfedge_const_handle *
Td_traits<Traits,Arrangement_on_surface_2>::m_p_he_at_top_inf = 0;

} //namespace CGAL
