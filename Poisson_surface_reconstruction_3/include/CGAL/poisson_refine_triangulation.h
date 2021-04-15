// Copyright (c) 2007-2008  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent RINEAU, Laurent Saboret

#ifndef CGAL_POISSON_REFINE_TRIANGULATION_H
#define CGAL_POISSON_REFINE_TRIANGULATION_H

#include <CGAL/license/Poisson_surface_reconstruction_3.h>


// CGAL
#include <CGAL/IO/trace.h>
#include <CGAL/Mesher_level.h>
#include <CGAL/Mesh_3/Poisson_refine_cells_3.h>
#include <CGAL/Poisson_mesh_cell_criteria_3.h>
#include <CGAL/surface_reconstruction_points_assertions.h>
#include <CGAL/Surface_mesh_traits_generator_3.h>

namespace CGAL {

/// \cond SKIP_IN_MANUAL

/// Utility class for poisson_refine_triangulation():
/// implements Delaunay refinement in a loose bounding
/// box of point set (break bad tetrahedra, where
/// bad means badly shaped or too big).
///
/// This class must be derived to inherit from Mesher_level.
template <class Tr,
          class Criteria,
          class Surface,
          class Oracle,
          class Container = Meshes::Double_map_container<typename Tr::Cell_handle,
                                                         typename Criteria::Cell_quality>
>
class Poisson_mesher_level_impl_base :
    public Mesh_3::Poisson_refine_tets_with_oracle_base<Tr, Criteria, Surface, Oracle, Container>
{
  typedef Mesh_3::Poisson_refine_tets_with_oracle_base<Tr, Criteria, Surface, Oracle, Container> Base;

public:
  // Inherited methods and fields used below
  using Base::triangulation_ref_impl;
  using Base::oracle;
  using Base::surface;
  using Base::should_be_refined;

  typedef typename Tr::Geom_traits Geom_traits;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Point Point;
  typedef typename Base::Cell_quality Cell_quality;

public:
  /** \name CONSTRUCTORS */

  Poisson_mesher_level_impl_base(Tr& t, Criteria crit, unsigned int max_vert, Surface& surface, Oracle& oracle)
    : Base(t, crit, surface, oracle),
      max_vertices(max_vert) ///< number of vertices bound (ignored if zero)
  {}

protected:
  /* --- protected functions --- */

  bool test_if_cell_is_bad(const Cell_handle c)
  {
    Cell_quality q;
    if( is_in_domain(c) && should_be_refined(c, q) )
    {
      this->add_bad_element(c, q);
      return true;
    }
    return false;
  }

  bool is_in_domain(const Cell_handle c)
  {
    return oracle.is_in_volume(surface, triangulation_ref_impl().dual(c));
  }

public:
  /* Overriden functions of this level: */
  /* we override all methods that call test_if_cell_is_bad() */

  void scan_triangulation_impl()
  {
    for(typename Tr::Finite_cells_iterator cit = triangulation_ref_impl().finite_cells_begin(),
        eit = triangulation_ref_impl().finite_cells_end();
        cit != eit;
        ++cit)
      test_if_cell_is_bad(cit);
  }

  void after_insertion_impl(const Vertex_handle& v)
  {
    update_star(v);
  }

  void update_star(const Vertex_handle& v)
  {
    // scan refiner
    typedef std::vector<Cell_handle> Cells;
    typedef typename Cells::iterator Cell_iterator;
    Cells incident_cells;
    incident_cells.reserve(32);
    triangulation_ref_impl().incident_cells(v, std::back_inserter(incident_cells));

    for(Cell_iterator cit = incident_cells.begin();
        cit != incident_cells.end();
        ++cit)
    {
      if( ! triangulation_ref_impl().is_infinite(*cit) )
        test_if_cell_is_bad(*cit);
    }
  }

  /// Tells if the algorithm is done.
  bool no_longer_element_to_refine_impl() const
  {
    return Base::no_longer_element_to_refine_impl() ||
           (max_vertices > 0 && triangulation_ref_impl().number_of_vertices() >= max_vertices);
  }


private:
  /* --- private datas --- */
  unsigned int max_vertices; ///< number of vertices bound (ignored if zero)

}; // end Poisson_mesher_level_impl_base


/// Utility class for poisson_refine_triangulation():
/// glue class that inherits from both Mesher_level
/// and Poisson_mesher_level_impl_base.
template <typename Tr,
          typename Criteria,
          typename Surface,
          typename Oracle = typename CGAL::Surface_mesh_traits_generator_3<Surface>::type,
          typename PreviousLevel = Null_mesher_level
 >
class Poisson_mesher_level :
  public Poisson_mesher_level_impl_base<Tr, Criteria, Surface, Oracle>,
  public Mesher_level <
    Tr,
    Poisson_mesher_level<Tr, Criteria, Surface, Oracle, PreviousLevel>,
    typename Tr::Cell_handle,
    PreviousLevel,
    Triangulation_mesher_level_traits_3<Tr>
  >
{
  typedef Poisson_mesher_level_impl_base<Tr, Criteria, Surface, Oracle> Base;

public:
  typedef Mesher_level <
    Tr,
    Poisson_mesher_level<Tr, Criteria, Surface, Oracle, PreviousLevel>,
    typename Tr::Cell_handle,
    PreviousLevel,
    Triangulation_mesher_level_traits_3<Tr>
  > Mesher;

  Poisson_mesher_level(Tr& t, Criteria criteria, unsigned int max_vertices, Surface& surface, Oracle& oracle, PreviousLevel& previous_level)
    : Base(t, criteria, max_vertices, surface, oracle),
      Mesher(previous_level)
  {
  }

}; // end class Poisson_mesher_level


/// Delaunay refinement in a loose bounding box
/// of input point set (break bad tetrahedra, where
/// bad means badly shaped or too big).
/// @return the number of vertices inserted.
///
/// @commentheading Preconditions:
/// - Tr must use a geometric traits with robust circumcenter computation.
/// - convergence is guaranteed if radius_edge_ratio_bound >= 1.0.
///
/// @commentheading Template Parameters:
/// @param Tr 3D Delaunay triangulation.
/// @param Surface Sphere_3 or Iso_cuboid_3.
/// @param Sizing_field A sizing field functor type
/// @param Second_sizing_field A sizing field functor type
///
/// @commentheading Sizing fields
/// - The first sizing field is the real sizing field that is targeted by
/// the refinement process. It may be costly to use.
/// - The second sizing field is supposed to be a sizing field that is less
/// costly to use (such as a constant sizing field). If a cell has a radius
/// that is smaller than the value of the second sizing field at its
/// center, then the cell is considered as small enough and the first
/// sizing field is not evaluated for that cell.
template <typename Tr,
          typename Surface,
          typename Sizing_field,
          typename Second_sizing_field>
unsigned int poisson_refine_triangulation(
  Tr& tr,
  double radius_edge_ratio_bound, ///< radius edge ratio bound (>= 1.0)
  const Sizing_field& sizing_field, ///< sizing field for cell radius bound
  const Second_sizing_field& second_sizing_field, ///< second sizing field for cell radius bound
  unsigned int max_vertices, ///< number of vertices bound (ignored if zero)
  Surface& enlarged_bbox) ///< new bounding sphere or box
{

  // Convergence is guaranteed if radius_edge_ratio_bound >= 1.0
  CGAL_surface_reconstruction_points_precondition(radius_edge_ratio_bound >= 1.0);

  // Mesher_level types
  typedef Poisson_mesh_cell_criteria_3<
    Tr
    , Sizing_field
    , Second_sizing_field
    > Tets_criteria;
  typedef typename CGAL::Surface_mesh_traits_generator_3<Surface>::type Oracle;
  typedef Poisson_mesher_level<Tr, Tets_criteria, Surface, Oracle, Null_mesher_level> Refiner;


  std::size_t nb_vertices = tr.number_of_vertices(); // get former #vertices

  // Delaunay refinement
  Tets_criteria tets_criteria(radius_edge_ratio_bound*radius_edge_ratio_bound,
                              sizing_field,
                              second_sizing_field);
  Oracle oracle;
  Null_mesher_level null_mesher_level;
  Refiner refiner(tr, tets_criteria, max_vertices, enlarged_bbox, oracle, null_mesher_level);
  refiner.scan_triangulation(); // Push bad cells to the queue
  refiner.refine(Null_mesh_visitor()); // Refine triangulation until queue is empty

  nb_vertices = tr.number_of_vertices() - nb_vertices;


  return static_cast<unsigned int>(nb_vertices);
}

/// \endcond


} //namespace CGAL

#endif // CGAL_POISSON_REFINE_TRIANGULATION_H
