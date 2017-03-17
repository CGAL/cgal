// Copyright (c) 2010  INRIA Sophia-Antipolis (France).
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
// $URL:$
// $Id:$
//
//
// Author(s)     : Mikhail Bogdanov

#ifndef CGAL_MESH_3_PERIODIC_TRIANGULATION_3_H
#define CGAL_MESH_3_PERIODIC_TRIANGULATION_3_H

// traits class
#include <CGAL/Kernel_traits.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Robust_circumcenter_traits_3.h>
#include <CGAL/Robust_weighted_circumcenter_filtered_traits_3.h>

// periodic issues
#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/Periodic_3_regular_triangulation_3.h>

// vertex and cell bases
#include <CGAL/Mesh_vertex_base_3.h>
#include <CGAL/Mesh_cell_base_3.h>
#include <CGAL/Triangulation_cell_base_with_circumcenter_3.h>

// mesh domain
#include <CGAL/Implicit_mesh_domain_3.h>

// IO
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

#include <CGAL/Periodic_3_mesh_3/config.h>
#include <CGAL/array.h>
#include <CGAL/tags.h>

namespace CGAL {

/// This class currently provides an interface between the classe
/// `CGAL::Periodic_3_regular_triangulation_3` and the mesher `Mesh_3`.
/// As periodic triangulations are parallelized, a lot of these functions will
/// become obsolete.
template<class Gt, class Tds>
class Periodic_3_regular_triangulation_3_mesher_3
    : public Periodic_3_regular_triangulation_3<Gt, Tds> {
public:

  typedef Sequential_tag   Concurrency_tag;
  typedef void             Lock_data_structure;

  void *get_lock_data_structure() const { return 0; }
  void set_lock_data_structure(void *) const { }

  typedef Periodic_3_regular_triangulation_3<Gt, Tds>  Base;

  typedef Gt                              Geometric_traits;
  typedef Tds                             Triangulation_data_structure;

  typedef typename Base::Cell_iterator    Finite_cells_iterator;
  typedef typename Base::Facet_iterator   Finite_facets_iterator;
  typedef typename Base::Edge_iterator    Finite_edges_iterator;
  typedef typename Base::Vertex_iterator  Finite_vertices_iterator;

  typedef typename Base::Vertex_handle    Vertex_handle;
  typedef typename Base::Cell_handle      Cell_handle;
  typedef typename Base::Facet            Facet;
  typedef typename Base::Edge             Edge;

  typedef typename Base::Point            Point;
  typedef typename Base::Bare_point       Bare_point;
  typedef typename Base::Weighted_point   Weighted_point;
  typedef typename Base::Periodic_point   Periodic_point;

  typedef typename Base::Locate_type      Locate_type;

  typedef typename Base::Segment          Segment;
  typedef typename Base::Periodic_segment Periodic_segment;

  typedef typename Base::Tetrahedron      Tetrahedron;
  typedef typename Base::Periodic_tetrahedron Periodic_tetrahedron;

  typedef typename Base::Offset           Offset;
  typedef typename Base::Iso_cuboid       Iso_cuboid;
  typedef typename Base::Conflict_tester  Conflict_tester;
  typedef typename Base::Covering_sheets  Covering_sheets;

  using Base::dual;
  using Base::tds;
  using Base::get_offset;
  using Base::incident_cells;
  using Base::incident_edges;
  using Base::incident_facets;
  using Base::segment;
  using Base::set_offsets;
  using Base::point;
  using Base::tetrahedron;
  using Base::periodic_tetrahedron;
#ifndef CGAL_NO_STRUCTURAL_FILTERING
  using Base::inexact_locate;
#endif

  /// Constructor
  Periodic_3_regular_triangulation_3_mesher_3(
      const Iso_cuboid& domain = Iso_cuboid(0,0,0,1,1,1),
      const Geometric_traits& gt = Geometric_traits())
    : Base(domain, gt)
  {
    this->insert_dummy_points();
  }

  /// Concurrency related
  template <typename Cell_handle>
  bool try_lock_cell(const Cell_handle &, int = 0) const { return true; }

  bool try_lock_and_get_incident_cells(Vertex_handle /* v */,
                                       std::vector<Cell_handle>& /* cells */) const
  {
    std::cerr << "ERROR: implement try_lock_and_get_incident_cells()"<< std::endl;
    return true;
  }

  /// Basic setters-getters
  // there are no infinite elements in a periodic triangulation
  bool is_infinite(const Vertex_handle) const { return false; }
  bool is_infinite(const Edge&) const { return false; }
  bool is_infinite(const Facet&) const { return false; }
  bool is_infinite(const Cell_handle) const { return false; }
  bool is_infinite(const Cell_handle, int) const { return false; }
  bool is_infinite(const Cell_handle c, int i, int j) const;

  Cell_handle infinite_cell() const
  {
    // there are no infinite cells in a periodic triangulation
    assert(false);
    return Cell_handle();
  }

  Vertex_handle infinite_vertex() const
  {
    // there is no infinite vertex in a periodic triangulation
    assert(false);
    return Vertex_handle();
  }

  int dimension() const
  {
    return (this->number_of_vertices() == 0) ? -2 : 3;
  }

  /// transform a point (living anywhere in space) into a point living inside
  /// the canonical iso cuboid and an offset
  Point canonicalize_point(const Point& p) const
  {
    return point(this->periodic_point(p));
  }

  void set_domain(const Iso_cuboid& domain)
  {
    Base::set_domain(domain);
    this->insert_dummy_points();
  }

  Tetrahedron tetrahedron(const Cell_handle c) const
  {
    Periodic_tetrahedron ptet = periodic_tetrahedron(c);
    return tetrahedron(ptet);
  }

  /*!
  Copies all finite `Edge`s incident to `v` to the
  output iterator `edges`. Returns the resulting output iterator.

  Since there are no infinite edges in a periodic triangulation, this
  function is simply a wrapper around `incident_edges`

  \pre `t.dimension() > 0`, `v != Vertex_handle()`, `t.is_vertex(v)`.
  */
  template<class OutputIterator>
  OutputIterator
  finite_incident_edges(Vertex_handle v, OutputIterator edges) const
  {
    return incident_edges(v, edges);
  }

  /*!
  Copies the `Cell_handle`s of all finite cells incident to `v` to the output
  iterator `cells`.
  Returns the resulting output iterator.

  Since there are no infinite cells in a periodic triangulation, this
  function is simply a wrapper around `incident_cells`

\pre `t.dimension() == 3`, `v != Vertex_handle()`, `t.is_vertex(v)`.
  */
  template<class OutputIterator>
  OutputIterator
  finite_incident_cells(Vertex_handle v, OutputIterator cells) const
  {
    return incident_cells(v, cells);
  }

  /*!
  Copies all finite `Facet`s incident to `v` to the output iterator
  `facets`.
  Returns the resulting output iterator.

  Since there are no infinite facets in a periodic triangulation, this
  function is simply a wrapper around `incident_facets`

  \pre `t.dimension() > 1`, `v != Vertex_handle()`, `t.is_vertex(v)`.
  */
  template<class OutputIterator>
  OutputIterator
  finite_incident_facets(Vertex_handle v, OutputIterator facets) const
  {
    return incident_facets(v, facets);
  }

  void clear_v_offsets() const
  {
    for (typename std::vector<Vertex_handle>::iterator voit = this->v_offsets.begin();
         voit != this->v_offsets.end() ; ++voit) {
      (*voit)->clear_offset();
    }
    this->v_offsets.clear();
  }

  /// Call `CGAL::side_of_power_sphere` with a canonicalized point
  Bounded_side side_of_power_sphere(const Cell_handle& c, const Point& p,
                                    bool perturb = false) const
  {
    Point point = this->canonicalize_point(p);
    return Base::side_of_power_sphere(c, point, Offset(), perturb);
  }

  /// \name Locate functions
  ///
  /// Locate points within a periodic triangulation.
  ///
  /// These functions are temporarily here to interface between Mesh_3 and
  /// the periodic triangulations, until the latter are made parallel.
  ///
  /// \sa `CGAL::Regular_triangulation_3::locate`
  /// @{
  Cell_handle locate(const Point & p,
                     Cell_handle start = Cell_handle(),
                     bool* could_lock_zone = NULL) const
  {
    assert(could_lock_zone == NULL);
    return Base::locate(p, start);
  }

  Cell_handle locate(const Point & p,
                     Vertex_handle hint,
                     bool* could_lock_zone = NULL) const
  {
    assert(could_lock_zone == NULL);
    // compared to the non-periodic version in T3, the infinite cell cannot
    // be used; `Cell_handle()` is used instead
    return Base::locate(p, hint == Vertex_handle() ? Cell_handle() : hint->cell());
  }

  Cell_handle locate(const Point& p,
                     Locate_type& l, int& i, int& j,
                     Cell_handle start = Cell_handle(),
                     bool* could_lock_zone = NULL) const
  {
    assert(could_lock_zone == NULL);
    return Base::locate(p, l, i, j, start);
  }

  Cell_handle locate(const Point& p,
                     Locate_type& l, int& i, int& j,
                     Vertex_handle hint,
                     bool* could_lock_zone = NULL) const
  {
    assert(could_lock_zone == NULL);
    return Base::locate(p, l, i, j,
                        hint == Vertex_handle() ? Cell_handle() : hint->cell());
  }
  /// @}

  /// \name Conflict functions
  /// Returns the vertices on the interior of the conflict hole.
  ///
  /// These functions are temporarily here to interface between Mesh_3 and
  /// the periodic triangulations, until the latter are made parallel.
  ///
  /// @{
  template <class OutputIterator>
  OutputIterator
  vertices_inside_conflict_zone(const /*Weighted_point*/Point& /* p */,
                                Cell_handle /* c */,
                                OutputIterator res) const
  {
    return res;
    assert(false); // not yet supported
  }

  //template < class Gt, class Tds >
  template <class OutputIteratorBoundaryFacets, class OutputIteratorCells,
            class OutputIteratorInternalFacets>
  Triple<OutputIteratorBoundaryFacets, OutputIteratorCells, OutputIteratorInternalFacets>
  find_conflicts(const Point &p,
                 Cell_handle c,
                 OutputIteratorBoundaryFacets bfit,
                 OutputIteratorCells cit,
                 OutputIteratorInternalFacets ifit,
                 bool* could_lock_zone = NULL,
                 const Facet* /* this_facet_must_be_in_the_cz */ = NULL,
                 bool* /* the_facet_is_in_its_cz */ = NULL ) const
  {
    assert(could_lock_zone == NULL);

    clear_v_offsets();

    CGAL_triangulation_precondition( this->number_of_vertices() != 0);

    const Point canonic_p = canonicalize_point(p);

    // #warning rewrite these lines
    Locate_type lt;
    int li, lj;
    locate( p, lt, li, lj, Cell_handle());
    if(lt == 0 /*Locate_type::VERTEX*/ ) {
      return make_triple(bfit, cit, ifit);
    }

    std::vector<Facet> facets;
    facets.reserve(64);
    std::vector<Cell_handle> cells;
    cells.reserve(32);

    Conflict_tester tester(canonic_p, this);
    Triple<typename std::back_insert_iterator<std::vector<Facet> >,
           typename std::back_insert_iterator<std::vector<Cell_handle> >,
           OutputIteratorInternalFacets> tit =
             Base::find_conflicts(c, tester,
                                  make_triple(std::back_inserter(facets),
                                              std::back_inserter(cells), ifit));
    ifit = tit.third;

    // Reset the conflict flag on the boundary.
    for(typename std::vector<Facet>::iterator fit=facets.begin();
        fit != facets.end(); ++fit) {
      fit->first->neighbor(fit->second)->tds_data().clear();
      *bfit++ = *fit;
    }

    // Reset the conflict flag in the conflict cells.
    for(typename std::vector<Cell_handle>::iterator ccit=cells.begin();
        ccit != cells.end(); ++ccit) {
      (*ccit)->tds_data().clear();
      *cit++ = *ccit;
    }

    return make_triple(bfit, cit, ifit);
  }

  template <class OutputIteratorBoundaryFacets, class OutputIteratorCells>
  std::pair<OutputIteratorBoundaryFacets, OutputIteratorCells>
  find_conflicts(const Point &p, Cell_handle c,
                 OutputIteratorBoundaryFacets bfit,
                 OutputIteratorCells cit,
                 bool* could_lock_zone = NULL) const
  {
    assert(could_lock_zone == NULL);

    Triple<OutputIteratorBoundaryFacets,
           OutputIteratorCells,
           Emptyset_iterator> t = find_conflicts(p, c, bfit, cit,
                                                 Emptyset_iterator(),
                                                 could_lock_zone);
    return std::make_pair(t.first, t.second);
  }
  /// @}

  /// \name Insert functions
  ///
  /// Insert points in the triangulation.
  ///
  /// These functions are temporarily here to interface between Mesh_3 and
  /// the periodic triangulations, until the latter are made parallel.
  ///
  /// @{
  template <class CellIt>
  Vertex_handle insert_in_hole(const Point & p, CellIt cell_begin, CellIt cell_end,
                               Cell_handle begin, int i)
  {
    Vertex_handle v = tds().insert_in_hole(cell_begin, cell_end, begin, i);
    v->set_point(canonicalize_point(p));

    std::vector<Cell_handle> nbs;
    incident_cells(v, std::back_inserter(nbs));

    // For all neighbors of the newly added vertex v: fetch their offsets from
    // the tester and reset them in the triangulation data structure.
    for (typename std::vector<Cell_handle>::iterator cit = nbs.begin();
         cit != nbs.end(); cit++) {
      Offset off[4];
      for (int i=0 ; i<4 ; i++)
        off[i] = (*cit)->vertex(i)->offset();

      set_offsets(*cit, off[0], off[1], off[2], off[3]);
    }

    clear_v_offsets();

    return v;
  }

  Vertex_handle insert(const Point& p,
                       Cell_handle start = Cell_handle(),
                       bool* could_lock_zone = NULL)
  {
    assert(could_lock_zone == NULL);
    return Base::insert(canonicalize_point(p), start);
  }

  Vertex_handle insert(const Point& p,
                       Vertex_handle hint,
                       bool* could_lock_zone = NULL)
  {
    assert(could_lock_zone == NULL);
    // compared to the non-periodic version in T3, the infinite cell cannot
    // be used; `Cell_handle()` is used instead
    return Base::insert(canonicalize_point(p),
                        hint == Vertex_handle() ? Cell_handle() : hint->cell());
  }

  Vertex_handle insert(const Point& p,
                       Locate_type lt, Cell_handle loc, int li, int lj,
                       bool* could_lock_zone = NULL)
  {
    assert(could_lock_zone == NULL);
    return Base::insert(canonicalize_point(p), lt, loc, li, lj);
  }
  /// @}

  /// Remove function
  void remove(Vertex_handle v,
              bool* could_lock_zone = NULL)
  {
    assert(could_lock_zone == NULL);
    return Base::remove(v);
  }

  /// Dual computations
  Point dual(Cell_handle c) const
  {
    // return the canonical point
    //return canonicalize_point(periodic_circumcenter(c));

    // return the point with respect to the canonical cell c
    return this->geom_traits().construct_weighted_circumcenter_3_object()(
          c->vertex(0)->point(), c->vertex(1)->point(),
          c->vertex(2)->point(), c->vertex(3)->point(),
          get_offset(c,0), get_offset(c,1),
          get_offset(c,2), get_offset(c,3));
  }

  Object dual(const Facet & f) const
  {
    Segment s = segment(Base::dual(f));
    return make_object(s);
  }
};

namespace details {

template<typename K>
struct Periodic_mesh_geom_traits_generator
{
private:
  // Mesh_3 requires the Bare_point / Weighted_point typedef
  // but `Robust_weighted_circumcenter_filtered_traits_3` no longer
  // derives from `Regular_triangulation_euclidean_traits_3`.
  // Since it's just a typedef and the order does not matter, we just
  // wrap here instead.
  typedef Regular_triangulation_euclidean_traits_3<
            Periodic_3_regular_triangulation_traits_3<
              Robust_weighted_circumcenter_filtered_traits_3<K> > > Geom_traits;

public:
  typedef Geom_traits type;
  typedef type Type;
};  // end struct Periodic_mesh_geom_traits_generator

}  // end namespace details

template<class MD, class K=typename Kernel_traits<MD>::Kernel>
struct Mesh_periodic_3_triangulation_3
{
private:
  // traits
  typedef typename details::Periodic_mesh_geom_traits_generator<K>::type Geom_traits;

  // Periodic vertex and cell bases
  typedef Periodic_3_triangulation_ds_vertex_base_3<> VbDS;
  typedef Triangulation_vertex_base_3<Geom_traits, VbDS> PVb;

  typedef Periodic_3_triangulation_ds_cell_base_3<> CbDS;
  typedef Triangulation_cell_base_3<Geom_traits, CbDS> PTCb;
  typedef Triangulation_cell_base_with_circumcenter_3<Geom_traits, PTCb> PCb;

  // Mesh vertex and cell bases (inherits from periodic ones)
  typedef Mesh_vertex_base_3<Geom_traits, MD, PVb> Vertex_base;
  typedef Mesh_cell_base_3<Geom_traits, MD, PCb> Cell_base;

  // Triangulation and tds
  typedef Triangulation_data_structure_3<Vertex_base, Cell_base> Tds;
  typedef Periodic_3_regular_triangulation_3_mesher_3<Geom_traits, Tds> Triangulation;

public:
  typedef Triangulation type;
  typedef type Type;
};  // end struct Mesh_triangulation_3

// helper function moving periodic
// triangles in some canonical positions in order to get a surface with
// less "holes"
template<class Triangulation>
typename Triangulation::Periodic_triangle
canonicalize(const typename Triangulation::Periodic_triangle& pt)
{
  typedef typename Triangulation::Offset Offset;

  Offset o0 = pt[0].second;
  Offset o1 = pt[1].second;
  Offset o2 = pt[2].second;
  int diffx = std::min(o0.x(),std::min(o1.x(),o2.x()));
  int diffy = std::min(o0.y(),std::min(o1.y(),o2.y()));
  int diffz = std::min(o0.z(),std::min(o1.z(),o2.z()));
  Offset diff_off(diffx,diffy,diffz);

  return make_array(std::make_pair(pt[0].first,o0 - diff_off),
                    std::make_pair(pt[1].first,o1 - diff_off),
                    std::make_pair(pt[2].first,o2 - diff_off));
}

template<class Triangulation>
typename Triangulation::Periodic_tetrahedron
canonicalize_tetrahedron(const typename Triangulation::Periodic_tetrahedron& pt)
{
  typedef typename Triangulation::Offset Offset;

  Offset o0 = pt[0].second;
  Offset o1 = pt[1].second;
  Offset o2 = pt[2].second;
  Offset o3 = pt[3].second;

  int diffx = std::min(std::min(o0.x(),o1.x()),std::min(o2.x(),o3.x()));
  int diffy = std::min(std::min(o0.y(),o1.y()),std::min(o2.y(),o3.y()));
  int diffz = std::min(std::min(o0.z(),o1.z()),std::min(o2.z(),o3.z()));
  Offset diff_off(diffx,diffy,diffz);

  return make_array(std::make_pair(pt[0].first,o0 - diff_off),
                    std::make_pair(pt[1].first,o1 - diff_off),
                    std::make_pair(pt[2].first,o2 - diff_off),
                    std::make_pair(pt[3].first,o3 - diff_off));
}

template <class Stream, class C3t3>
Stream &vertices_medit(Stream &out, C3t3 &c3t3)
{
  out << std::setprecision(20);
  out << "MeshVersionFormatted 1" << std::endl;
  out << "Dimension 3" << std::endl;
  out << "Vertices " << c3t3.triangulation().nb_of_vertices() * 8 << std::endl;
}

// Writing a restricted Delaunay triangulation to the .mesh file format
// Writing the triangulation to 8 domains.
// Can be used with the medit (a viewer)
template <class Stream, class C3t3>
Stream &write_complex_to_medit(Stream &out, C3t3 &c3t3,
                               unsigned occurence_count = 8)
{
  typedef typename C3t3::Triangulation           Triangulation;
  typedef Triangulation                          Tr;

  typedef typename Triangulation::Iso_cuboid     Iso_cuboid;

  typedef typename Triangulation::Triangle       Triangle;
  typedef typename Triangulation::Tetrahedron    Tetrahedron;

  typedef typename C3t3::Facet_iterator          Facet_iterator;
  typedef typename C3t3::Cell_iterator           Cell_iterator;

  Triangulation& t = c3t3.triangulation();
  int number_of_facets = static_cast<int>(c3t3.number_of_facets());
  int number_of_cells = static_cast<int>(c3t3.number_of_cells());
  int number_of_vertices = 3 * number_of_facets + 4 * number_of_cells;
  out << std::setprecision(20);
  out << "MeshVersionFormatted 1\nDimension 3\nVertices"
      << "\n" << number_of_vertices * occurence_count
      << std::endl;

  Iso_cuboid cb = t.domain();

  for(unsigned j = 0; j < occurence_count; j++ ) {
    for (Facet_iterator it =c3t3.facets_begin(); it!=c3t3.facets_end(); it++) {
      Triangle tri = t.triangle(canonicalize<Tr>(t.periodic_triangle(*it)));
      for(int i = 0; i < 3; i++) {
        out << tri[i].x() + (j&1) << " "
            << tri[i].y() + ((j&2) >> 1) << " "
            << tri[i].z() + ((j&4) >> 2) << " "
            << 32*j + 1 << std::endl;
      }
    }
  }

  for(unsigned j = 0; j < occurence_count; j++ ) {
    for (Cell_iterator it = c3t3.cells_begin(); it !=c3t3.cells_end(); it++) {
      Tetrahedron tet = t.tetrahedron(canonicalize_tetrahedron<Tr>(t.periodic_tetrahedron( it )));
      for(int i = 0; i < 4; i++) {
        out << tet[i].x() + (j&1) << " "
            << tet[i].y() + ((j&2) >> 1) << " "
            << tet[i].z() + ((j&4) >> 2) << " "
            << 32*j + 1 << std::endl;
      }
    }
  }

  int first_vertex = 1;
  out << "Triangles\n"
      << number_of_facets * occurence_count
      << std::endl;
  const int number_of_vertices_on_facets = number_of_facets * 3;
  for(unsigned j = 0; j < occurence_count; j++ ) {

    for( int i = 0; i < number_of_facets; i++) {
      out << i * 3 + j * number_of_vertices_on_facets + first_vertex  << " "
          << i * 3+1 + j * number_of_vertices_on_facets + first_vertex << " "
          << i * 3+2 + j * number_of_vertices_on_facets + first_vertex << " "
          << 128 * j + 1 << std::endl;
    }
  }

  const int shift = number_of_vertices_on_facets * occurence_count;
  out << "Tetrahedra\n"
      << number_of_cells * occurence_count
      << std::endl;
  const int number_of_vertices_on_cells = number_of_cells * 4;
  for(unsigned j = 0; j < occurence_count; j++ ) {

    Cell_iterator it = c3t3.cells_begin();

    for (int i = 0; i < number_of_cells; i++) {
      out << i * 4 +     j * number_of_vertices_on_cells + first_vertex + shift << " "
          << i * 4 + 1 + j * number_of_vertices_on_cells + first_vertex + shift << " "
          << i * 4 + 2 + j * number_of_vertices_on_cells + first_vertex + shift << " "
          << i * 4 + 3 + j * number_of_vertices_on_cells + first_vertex + shift << " "
          << it->subdomain_index() << std::endl;
      //  << 128 * j + 64 + 1 << std::endl;

      it++;
    }
  }

  out << "0\nEnd";
  return out;
}

} // namespace CGAL

#endif // CGAL_MESH_3_PERIODIC_TRIANGULATION_3_H

