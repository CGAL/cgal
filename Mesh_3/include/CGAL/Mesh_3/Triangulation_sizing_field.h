// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : Defines a sizing field which use an internal triangulation
// to store the sizes 
//******************************************************************************

#ifndef CGAL_MESH_3_TRIANGULATION_SIZING_FIELD_H
#define CGAL_MESH_3_TRIANGULATION_SIZING_FIELD_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <boost/iterator/transform_iterator.hpp>

namespace CGAL {

namespace Mesh_3
{
  

/**
 * @class Triangulation_sizing_field
 */
template <typename Tr>
class Triangulation_sizing_field
{
  // Types
  typedef typename Tr::Geom_traits    Gt;
  typedef typename Tr::Weighted_point Weighted_point;
  typedef typename Gt::FT             FT;
  
  typedef Triangulation_vertex_base_with_info_3<FT, Gt>   Vb;
  typedef Triangulation_cell_base_3<Gt>                   Cb;
  typedef Triangulation_data_structure_3<Vb, Cb>          Tds;
  typedef Regular_triangulation_3<Gt,Tds>                 Compact_triangulation;
  typedef Compact_triangulation                           Ctr;
  
  typedef typename Tr::Vertex_handle      Vertex_handle;
  typedef typename Tr::Cell_handle        Cell_handle;
  typedef typename Ctr::Cell_handle       CCell_handle;
  typedef typename Ctr::Vertex_handle     CVertex_handle;
  
public:
  // Vertices of mesh triangulation do not need to be updated 
  static const bool is_vertex_update_needed = false;
  
public:
  /**
   * Constructor
   */
  Triangulation_sizing_field(const Tr& tr);
  
  /**
   * Fill sizing field, using size associated to point in \c value_map.
   */
  void fill(const std::map<Weighted_point, FT>& value_map);
  
  /**
   * Returns size at point \c p.
   */
  FT operator()(const Weighted_point& p) const;
  
  /**
   * Returns size at point \c p. (needed for compatibility)
   */
  template <typename Handle>
  FT operator()(const Weighted_point& p, const Handle&) const
  { return this->operator()(p); }
  
private:
  /**
   * Returns size at point \c p, by interpolation into tetrahedron.
   */
  FT interpolate_on_cell_vertices(const Weighted_point& p,
                                  const CCell_handle& cell) const;
  
  /**
   * Returns size at point \c p, by interpolation into facet (\c cell is assumed
   * to be an infinite cell).
   */
  FT interpolate_on_facet_vertices(const Weighted_point& p,
                                   const CCell_handle& cell) const;
  
  /**
   * Returns an hint for \c p location.
   */
  CCell_handle get_hint(const Weighted_point& p) const
  { return last_cell_; }

  /**
   * A functor which extract the point from a vertex handle.
   * Used by boost transform iterator
   */
  struct Extract_point :
  public std::unary_function<typename Tr::Vertex,Weighted_point>
  {
    Weighted_point operator()(const typename Tr::Vertex& v) const { return v.point(); }
  };
  
private:
  Compact_triangulation ctr_;
  mutable CCell_handle last_cell_;
};
  
  
template <typename Tr>
Triangulation_sizing_field<Tr>::
Triangulation_sizing_field(const Tr& tr)
  : ctr_(boost::make_transform_iterator(tr.finite_vertices_begin(), 
                                        Extract_point()),
         boost::make_transform_iterator(tr.finite_vertices_end(),
                                        Extract_point()) )
  , last_cell_(CCell_handle())
{
}  
  
  
  
template <typename Tr>
void
Triangulation_sizing_field<Tr>::
fill(const std::map<Weighted_point, FT>& value_map)
{
  typedef typename Ctr::Finite_vertices_iterator  Fvi;
  
  for ( Fvi vit = ctr_.finite_vertices_begin() ;
        vit != ctr_.finite_vertices_end() ;
        ++ vit )
  {
    typename std::map<Weighted_point, FT>::const_iterator find_result = 
      value_map.find(vit->point());
    
    if ( find_result != value_map.end() )
      vit->info() = find_result->second;
    else
      vit->info() = FT(0);
  }
}  
  
  
template <typename Tr>
typename Triangulation_sizing_field<Tr>::FT
Triangulation_sizing_field<Tr>::
operator()(const Weighted_point& p) const  
{  
  CCell_handle hint = get_hint(p);
  CCell_handle cell = ctr_.locate(p,hint);
  last_cell_ = cell;
  
  if ( !ctr_.is_infinite(cell) )
    return interpolate_on_cell_vertices(p,cell);
  else
    return interpolate_on_facet_vertices(p,cell);
}
  

template <typename Tr>
typename Triangulation_sizing_field<Tr>::FT
Triangulation_sizing_field<Tr>::
interpolate_on_cell_vertices(const Weighted_point& p, const CCell_handle& cell) const
{
  typename Gt::Compute_volume_3 volume =
    ctr_.geom_traits().compute_volume_3_object();
  
  // Interpolate value using tet vertices values
  const FT& va = cell->vertex(0)->info();
  const FT& vb = cell->vertex(1)->info();
  const FT& vc = cell->vertex(2)->info();
  const FT& vd = cell->vertex(3)->info();
  
  const Weighted_point& a = cell->vertex(0)->point();
  const Weighted_point& b = cell->vertex(1)->point();
  const Weighted_point& c = cell->vertex(2)->point();
  const Weighted_point& d = cell->vertex(3)->point();
  
  const FT abcp = CGAL::abs(volume(a,b,c,p));
  const FT abdp = CGAL::abs(volume(a,b,d,p));
  const FT acdp = CGAL::abs(volume(a,c,d,p));
  const FT bcdp = CGAL::abs(volume(b,c,d,p));
  
  // TODO: improve this (static filter ?)
  // If volume is 0, then compute the average value
  if ( is_zero(abcp+abdp+acdp+bcdp) )
    return (va+vb+vc+vd)/4.;
  
  return ( (abcp*vd + abdp*vc + acdp*vb + bcdp*va) / (abcp+abdp+acdp+bcdp) );
}
  
template <typename Tr>
typename Triangulation_sizing_field<Tr>::FT
Triangulation_sizing_field<Tr>::
interpolate_on_facet_vertices(const Weighted_point& p, const CCell_handle& cell) const
{
  typename Gt::Compute_area_3 area =
    ctr_.geom_traits().compute_area_3_object();
  
  // Find infinite vertex and put it in k0
  int k0 = 0;
  int k1 = 1;
  int k2 = 2;
  int k3 = 3;
  
  if ( ctr_.is_infinite(cell->vertex(1)) )
    std::swap(k0,k1);
  if ( ctr_.is_infinite(cell->vertex(2)) )
    std::swap(k0,k2);
  if ( ctr_.is_infinite(cell->vertex(3)) )
    std::swap(k0,k3);
  
  // Interpolate value using tet vertices values
  const FT& va = cell->vertex(k1)->info();
  const FT& vb = cell->vertex(k2)->info();
  const FT& vc = cell->vertex(k3)->info();
  
  const Weighted_point& a = cell->vertex(k1)->point();
  const Weighted_point& b = cell->vertex(k2)->point();
  const Weighted_point& c = cell->vertex(k3)->point();
  
  const FT abp = area(a,b,p);
  const FT acp = area(a,c,p);
  const FT bcp = area(b,c,p);
  
  // TODO: improve this (static filter ?)
  // If area is 0, then compute the average value
  if ( is_zero(abp+acp+bcp) )
    return (va+vb+vc)/3.;
  
  return ( (abp*vc + acp*vb + bcp*va ) / (abp+acp+bcp) );
}

} // end namespace Mesh_3


} //namespace CGAL

#endif // CGAL_MESH_3_TRIANGULATION_SIZING_FIELD_H
