// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : Defines a sizing field stored into an external
// mesh triangulation
//******************************************************************************

#ifndef CGAL_MESH_3_MESH_SIZING_FIELD_H
#define CGAL_MESH_3_MESH_SIZING_FIELD_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/enumerable_thread_specific.h>
#endif

namespace CGAL {

namespace Mesh_3
{

/**
 * @class Mesh_sizing_field_base
 */
// Sequential
template <typename Cell_handle, typename Concurrency_tag>
class Mesh_sizing_field_base
{
protected:
  Cell_handle get_last_cell() const
  {
    return last_cell_;
  }

  void set_last_cell(Cell_handle c) const
  {
    last_cell_ = c;
  }

private:
  /// A cell that is used to accelerate location queries
  mutable Cell_handle last_cell_;
};

#ifdef CGAL_LINKED_WITH_TBB
/**
 * @class Mesh_sizing_field_base specialization
 */
// Parallel
template <typename Cell_handle>
class Mesh_sizing_field_base<Cell_handle, Parallel_tag>
{
protected:
  Cell_handle get_last_cell() const
  {
    return last_cell_.local();
  }

  void set_last_cell(Cell_handle c) const
  {
    last_cell_.local() = c;
  }

private:
  /// A cell that is used to accelerate location queries
  mutable tbb::enumerable_thread_specific<Cell_handle> last_cell_;
};
#endif // CGAL_LINKED_WITH_TBB

/**
 * @class Mesh_sizing_field
 */
template <typename Tr, bool Need_vertex_update = true>
class Mesh_sizing_field
  : public Mesh_sizing_field_base<typename Tr::Cell_handle,
                                  typename Tr::Concurrency_tag>
{
  // Types
  typedef typename Tr::Geom_traits              Gt;
  typedef typename Tr::Bare_point               Bare_point;
  typedef typename Tr::Weighted_point           Weighted_point;
  typedef typename Gt::FT                       FT;

  typedef typename Tr::Vertex_handle            Vertex_handle;
  typedef typename Tr::Cell_handle              Cell_handle;

public:
  // update vertices of mesh triangulation ?
  static const bool is_vertex_update_needed = Need_vertex_update;

public:
  /**
   * Constructor
   */
  Mesh_sizing_field(Tr& tr);

  /**
   * Fill sizing field, using size associated to point in \c value_map
   */
  void fill(const std::map<Bare_point, FT>& value_map);

  /**
   * Returns size at point \c p.
   */
  FT operator()(const Bare_point& p) const
  { return this->operator()(p, this->get_last_cell()); }

  /**
   * Returns size at point \c p, using \c v to accelerate \c p location
   * in triangulation
   */
  FT operator()(const Bare_point& p, const Vertex_handle& v) const
  { return this->operator()(p, v->cell()); }

  /**
   * Returns size at point \c p.
   */
  FT operator()(const Bare_point& p, const Cell_handle& c) const;

  /**
   * Returns size at point \c p. Assumes that p is the centroid of c.
   */
  FT operator()(const Bare_point& p, const std::pair<Cell_handle, bool>& c) const;

private:
  /**
   * Returns size at point \c p, by interpolation into tetrahedron.
   */
  FT interpolate_on_cell_vertices(const Bare_point& p,
                                  const Cell_handle& cell) const;

  /**
   * Returns size at point \c p, by interpolation into facet (\c cell is assumed
   * to be an infinite cell).
   */
  FT interpolate_on_facet_vertices(const Bare_point& p,
                                   const Cell_handle& cell) const;

private:
  /// The triangulation
  Tr& tr_;
};


template <typename Tr, bool B>
Mesh_sizing_field<Tr,B>::
Mesh_sizing_field(Tr& tr)
  : tr_(tr)
{
}


template <typename Tr, bool B>
void
Mesh_sizing_field<Tr,B>::
fill(const std::map<Bare_point, FT>& value_map)
{
  typedef typename Tr::Finite_vertices_iterator  Fvi;

  typename Gt::Construct_point_3 cp = tr_.geom_traits().construct_point_3_object();

  for ( Fvi vit = tr_.finite_vertices_begin(); vit != tr_.finite_vertices_end(); ++ vit )
  {
    const Weighted_point& position = tr_.point(vit);

    typename std::map<Bare_point, FT>::const_iterator find_result =
      value_map.find(cp(position));

    if ( find_result != value_map.end() )
    {
      vit->set_meshing_info(find_result->second);
    }
    else
    {
      CGAL_assertion(false);
      vit->set_meshing_info(FT(0));
    }
  }
}


template <typename Tr, bool B>
typename Mesh_sizing_field<Tr,B>::FT
Mesh_sizing_field<Tr,B>::
operator()(const Bare_point& p, const Cell_handle& c) const
{
  typename Gt::Construct_weighted_point_3 cwp = tr_.geom_traits().construct_weighted_point_3_object();

#ifdef CGAL_MESH_3_SIZING_FIELD_INEXACT_LOCATE
  //use the inexact locate (much faster than locate) to get a hint
  //and then use locate to check whether p is really inside hint
  // if not, an exact locate will be performed
  Cell_handle hint = tr_.inexact_locate(cwp(p),c);
  const Cell_handle cell = tr_.locate(cwp(p), hint);
#else
  const Cell_handle cell = tr_.locate(cwp(p),c);
#endif
  this->set_last_cell(cell);

  if ( !tr_.is_infinite(cell) )
    return interpolate_on_cell_vertices(p, cell);
  else
    return interpolate_on_facet_vertices(p, cell);
}


template <typename Tr, bool B>
typename Mesh_sizing_field<Tr,B>::FT
Mesh_sizing_field<Tr,B>::
operator()(const Bare_point&, const std::pair<Cell_handle,bool>& c) const
{
  // Assumes that p is the centroid of c
  const Cell_handle& cell = c.first;

  // Interpolate value using tet vertices values
  const FT& va = cell->vertex(0)->meshing_info();
  const FT& vb = cell->vertex(1)->meshing_info();
  const FT& vc = cell->vertex(2)->meshing_info();
  const FT& vd = cell->vertex(3)->meshing_info();

  return ( (va+vb+vc+vd)/4 );
}


template <typename Tr, bool B>
typename Mesh_sizing_field<Tr,B>::FT
Mesh_sizing_field<Tr,B>::
interpolate_on_cell_vertices(const Bare_point& p, const Cell_handle& cell) const
{
  typename Gt::Construct_point_3 cp = tr_.geom_traits().construct_point_3_object();
  typename Gt::Compute_volume_3 volume = tr_.geom_traits().compute_volume_3_object();

  // Interpolate value using tet vertices values
  const FT& va = cell->vertex(0)->meshing_info();
  const FT& vb = cell->vertex(1)->meshing_info();
  const FT& vc = cell->vertex(2)->meshing_info();
  const FT& vd = cell->vertex(3)->meshing_info();

  const Weighted_point& wa = tr_.point(cell, 0);
  const Weighted_point& wb = tr_.point(cell, 1);
  const Weighted_point& wc = tr_.point(cell, 2);
  const Weighted_point& wd = tr_.point(cell, 3);
  const Bare_point& a = cp(wa);
  const Bare_point& b = cp(wb);
  const Bare_point& c = cp(wc);
  const Bare_point& d = cp(wd);

  const FT abcp = CGAL::abs(volume(a,b,c,p));
  const FT abdp = CGAL::abs(volume(a,d,b,p));
  const FT acdp = CGAL::abs(volume(a,c,d,p));
  const FT bcdp = CGAL::abs(volume(b,d,c,p));

  // If volume is 0, then compute the average value
  if ( is_zero(abcp+abdp+acdp+bcdp) )
    return (va+vb+vc+vd) / 4.;

  return ( (abcp*vd + abdp*vc + acdp*vb + bcdp*va) / (abcp+abdp+acdp+bcdp) );
}


template <typename Tr, bool B>
typename Mesh_sizing_field<Tr,B>::FT
Mesh_sizing_field<Tr,B>::
interpolate_on_facet_vertices(const Bare_point& p, const Cell_handle& cell) const
{
  typename Gt::Compute_area_3 area =  tr_.geom_traits().compute_area_3_object();

  typename Gt::Construct_point_3 cp = tr_.geom_traits().construct_point_3_object();
  // Find infinite vertex and put it in k0
  int k0 = 0;
  int k1 = 1;
  int k2 = 2;
  int k3 = 3;

  if ( tr_.is_infinite(cell->vertex(1)) )
    std::swap(k0,k1);
  if ( tr_.is_infinite(cell->vertex(2)) )
    std::swap(k0,k2);
  if ( tr_.is_infinite(cell->vertex(3)) )
    std::swap(k0,k3);

  // Interpolate value using tet vertices values
  const FT& va = cell->vertex(k1)->meshing_info();
  const FT& vb = cell->vertex(k2)->meshing_info();
  const FT& vc = cell->vertex(k3)->meshing_info();

  const Weighted_point& wa = tr_.point(cell, k1);
  const Weighted_point& wb = tr_.point(cell, k2);
  const Weighted_point& wc = tr_.point(cell, k3);
  const Bare_point& a = cp(wa);
  const Bare_point& b = cp(wb);
  const Bare_point& c = cp(wc);

  const FT abp = area(a, b, p);
  const FT acp = area(a, c, p);
  const FT bcp = area(b, c, p);

  CGAL_assertion(abp >= 0);
  CGAL_assertion(acp >= 0);
  CGAL_assertion(bcp >= 0);

  // If area is 0, then compute the average value
  if ( is_zero(abp+acp+bcp) )
    return (va+vb+vc)/3.;

  return ( (abp*vc + acp*vb + bcp*va ) / (abp+acp+bcp) );
}

} // end namespace Mesh_3


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MESH_3_MESH_SIZING_FIELD_H
