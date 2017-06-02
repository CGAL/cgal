// Copyright (c) 2004-2009  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Mikhail Bogdanov
//
#ifndef CGAL_PERIODIC_3_MESH_CELL_CRITERIA_3_H
#define CGAL_PERIODIC_3_MESH_CELL_CRITERIA_3_H

#include <CGAL/license/Periodic_3_mesh_3.h>

#include <CGAL/Periodic_3_mesh_3/config.h>

#include <CGAL/Periodic_3_mesh_3/periodic_3_mesh_standard_cell_criteria.h>

namespace CGAL {

template <typename Tr,
          typename Visitor_ = Mesh_3::Periodic_3_mesh_3::Cell_criteria_visitor_with_features<Tr> >
class Periodic_3_mesh_cell_criteria_3
{
  typedef Periodic_3_mesh_cell_criteria_3<Tr>     Self;

public:
  typedef Visitor_                                Visitor;
  typedef typename Visitor::Cell_quality          Cell_quality;
  typedef typename Visitor::Cell_badness          Cell_badness;

  typedef Mesh_3::Abstract_criterion<Tr, Visitor> Abstract_criterion;

private:
  typedef Mesh_3::Criteria<Tr,Visitor>            Criteria;

  typedef typename Tr::Cell_handle                Cell_handle;
  typedef typename Tr::Geom_traits::FT            FT;

  typedef typename Tr::Iso_cuboid                 Iso_cuboid;

public:
  /**
   * @brief Constructor
   *
   * @param radius_edge_bound the radius-edge bound
   * @param radius_bound the radius bound (tet sizing)
   */
  Periodic_3_mesh_cell_criteria_3(const Iso_cuboid& iso_cuboid,
                                  const FT& radius_edge_bound,
                                  const FT& radius_bound) :
    helper_(iso_cuboid)
  {
    if ( FT(0) != radius_bound )
      init_radius(radius_bound);

    if ( FT(0) != radius_edge_bound )
      init_radius_edge(radius_edge_bound);
  }

  /**
   * @brief Constructor
   *
   * @param radius_edge_bound the radius-edge bound
   * @param radius_bound the radius bound (tet sizing)
   */
  template<typename MD>
  Periodic_3_mesh_cell_criteria_3(const MD& periodic_domain,
                                  const FT& radius_edge_bound,
                                  const FT& radius_bound)
    : helper_(periodic_domain.periodic_bounding_box())
  {
    if ( FT(0) != radius_bound )
      init_radius(radius_bound);

    if ( FT(0) != radius_edge_bound )
      init_radius_edge(radius_edge_bound);
  }

  // Nb: SFINAE (dummy) to avoid wrong matches with built-in numerical types
  // as int.
  template <typename MD, typename Sizing_field>
  Periodic_3_mesh_cell_criteria_3(const MD& periodic_domain,
                                  const FT& radius_edge_bound,
                                  const Sizing_field& radius_bound,
                                  typename Sizing_field::FT /* dummy */ = 0)
    : helper_(periodic_domain.periodic_bounding_box())
  {
    init_radius(radius_bound);

    if ( FT(0) != radius_edge_bound )
      init_radius_edge(radius_edge_bound);
  }

  /// Destructor
  ~Periodic_3_mesh_cell_criteria_3() { }

  /**
   * @brief returns the badness of cell \c cell
   * @param cell the cell
   * @return the badness of \c cell
   */
  Cell_badness operator()(const Cell_handle& cell) const
  {
    return criteria_(cell);
  }

  void add(Abstract_criterion* criterion)
  {
    criteria_.add(criterion);
  }

private:
  void init_radius_edge(const FT& radius_edge_bound)
  {
    typedef Mesh_3::Periodic_3_mesh_3::Cell_radius_edge_criterion<Tr,Visitor> Radius_edge_criterion;
    criteria_.add(new Radius_edge_criterion(helper_, radius_edge_bound));
  }

  void init_radius(const FT& radius_bound)
  {
    typedef Mesh_3::Periodic_3_mesh_3::Cell_uniform_size_criterion<Tr,Visitor> Radius_criterion;
    criteria_.add(new Radius_criterion(helper_, radius_bound));
  }

  template < typename Sizing_field>
  void init_radius(const Sizing_field& radius_bound)
  {
    typedef Mesh_3::Periodic_3_mesh_3::Cell_variable_size_criterion<Tr,Visitor,Sizing_field>
      Radius_criterion;

    criteria_.add(new Radius_criterion(helper_, radius_bound));
  }

private:
  Criteria criteria_;

  Tr helper_;
};  // end class Periodic_3_mesh_cell_criteria_3

}  // end namespace CGAL


#endif // CGAL_PERIODIC_3_MESH_CELL_CRITERIA_3_H
