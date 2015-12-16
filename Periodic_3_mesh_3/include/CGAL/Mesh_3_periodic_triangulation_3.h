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
// $URL: svn+ssh://mbogdanov@scm.gforge.inria.fr/svn/cgal/trunk/Periodic_3_triangulation_3/include/CGAL/Periodic_3_Delaunay_triangulation_3.h $
// $Id: Mesh_3_periodic_3_triangulation_3.h 55473 2010-04-14 16:20:51Z mbogdanov $
// 
//
// Author(s)     : Mikhail Bogdanov

#ifndef CGAL_MESH_3_PERIODIC_TRIANGULATION_3_H
#define CGAL_MESH_3_PERIODIC_TRIANGULATION_3_H

// traits class
#include <CGAL/Kernel_traits.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Robust_circumcenter_traits_3.h>
#include <CGAL/Mesh_3/Robust_weighted_circumcenter_filtered_traits_3.h>

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

#include <CGAL/Periodic_mesh_3/config.h>
#include <CGAL/tags.h>

namespace CGAL {

template<class Gt, class Tds>
class Periodic_3_Delaunay_triangulation_3_Mesher_3 :
  public Periodic_3_regular_triangulation_3<Gt, Tds> {
public:

  typedef Sequential_tag   Concurrency_tag;
  typedef void             Lock_data_structure;

  void *get_lock_data_structure() const
  {
    return 0;
  }

  void set_lock_data_structure(void *) const
  {
  }

  typedef Periodic_3_regular_triangulation_3<Gt, Tds>  Base;
  typedef typename Base::Base Base_Base;

  typedef Gt                     Geometric_traits;
  typedef Tds                    Triangulation_data_structure;

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
  
  //typedef typename Gt::Bare_point         Bare_point;
  //typedef typename Gt::Weighted_point_3   Weighted_point;
  typedef typename Base::Locate_type Locate_type;
    

  typedef typename Base::Segment          Segment;
  typedef typename Base::Periodic_segment Periodic_segment;

  typedef typename Base::Tetrahedron      Tetrahedron;
  typedef typename Base::Periodic_tetrahedron Periodic_tetrahedron;

  typedef typename Base::Offset           Offset;
  typedef typename Base::Iso_cuboid       Iso_cuboid;
  typedef typename Base::Conflict_tester  Conflict_tester;
  typedef typename Base::Covering_sheets  Covering_sheets;

  using Base::tds;
  using Base::get_offset;
  using Base::incident_cells;
  using Base::incident_edges;
  using Base::incident_facets;
  using Base::segment;
  using Base::set_offsets;
  using Base::point;
  using Base::periodic_tetrahedron;
//  using Base::nearest_vertex;

  Periodic_3_Delaunay_triangulation_3_Mesher_3(
    const Iso_cuboid& domain = Iso_cuboid(0,0,0,1,1,1),
    const Geometric_traits& gt = Geometric_traits())
    : Base(domain, gt) {
      this->insert_dummy_points();      
  }

  template <typename Cell_handle>
  bool try_lock_cell(const Cell_handle &, int = 0) const
  { return true; }
  bool
  try_lock_and_get_incident_cells(Vertex_handle v,
                                  std::vector<Cell_handle>& cells) const
  {
    std::cerr << "ERROR: implement try_lock_and_get_incident_cells()"<< std::endl;
    return true;
  }


  bool is_infinite(const Vertex_handle v) const
  { return false; }

  bool is_infinite(const Cell_handle c) const
  { return false; }

  bool is_infinite(const Cell_handle c, int i) const
  { return false; }

  bool is_infinite(const Facet & f) const
  { return false; }

  bool is_infinite(const Cell_handle c, int i, int j) const;

  bool is_infinite(const Edge & e) const
  { return false; }
    
  int dimension() const 
  { 
    if (this->number_of_vertices() == 0) {
      return 0;
    }
    return 3; 
  }

  using Base::dual;

  Point dual(Cell_handle c) const {
    // it returns the canonical point
    //return point(periodic_circumcenter(c));

    // it returns the point with respect to the canonical cell c
    return this->geom_traits().construct_weighted_circumcenter_3_object()(
                                                    c->vertex(0)->point(), c->vertex(1)->point(),
                                                    c->vertex(2)->point(), c->vertex(3)->point(),
                                                    get_offset(c,0), get_offset(c,1),
                                                    get_offset(c,2), get_offset(c,3));
  }

  Object dual(const Facet & f) const {
    Segment s = segment(Base::dual(f));
    return make_object(s);
  }  

  Vertex_handle push_back(const Point& p) {
    return insert(p);
  }

  Vertex_handle insert(const Point & p, Cell_handle start = Cell_handle()) {
    return Base::insert(canonicalize_point(p), start);
  }

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
        for (int i=0 ; i<4 ; i++) {
          off[i] = (*cit)->vertex(i)->offset();
        }
        set_offsets(*cit, off[0], off[1], off[2], off[3]);
    }

    clear_v_offsets(); 

    return v;
  }

  void clear_v_offsets() const
  {
    for (typename std::vector<Vertex_handle>::iterator voit = this->v_offsets.begin();
      voit != this->v_offsets.end() ; ++voit) {
        (*voit)->clear_offset();
    }
    this->v_offsets.clear();
  }

  //template < class Gt, class Tds >
  template <class OutputIteratorBoundaryFacets, class OutputIteratorCells,
  class OutputIteratorInternalFacets>
    Triple<OutputIteratorBoundaryFacets, OutputIteratorCells,
    OutputIteratorInternalFacets> 
    
    find_conflicts( const Point &p,
    Cell_handle c, OutputIteratorBoundaryFacets bfit,
    OutputIteratorCells cit, OutputIteratorInternalFacets ifit
    , bool *could_lock_zone = NULL
    , const Facet *this_facet_must_be_in_the_cz = NULL
    , bool *the_facet_is_in_its_cz = NULL
    ) const /**/ {
      clear_v_offsets();

      CGAL_triangulation_precondition( this->number_of_vertices() != 0);

      const Point canonic_p = canonicalize_point(p);
      
      //#warning rewrite these lines
      Locate_type lt;
      int li, lj;
      locate( p, lt, li, lj, Cell_handle());
      if(lt == 0/*Locate_type::VERTEX*/) {
        return make_triple(bfit, cit, ifit);
      }

      std::vector<Facet> facets;
      facets.reserve(64);
      std::vector<Cell_handle> cells;
      cells.reserve(32);

      Conflict_tester tester(canonic_p, this);
      Triple<typename std::back_insert_iterator<std::vector<Facet> >,
        typename std::back_insert_iterator<std::vector<Cell_handle> >,
        OutputIteratorInternalFacets> tit = Base_Base::find_conflicts(c, tester,
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
                 OutputIteratorCells cit, bool*) const
  {
      Triple<OutputIteratorBoundaryFacets,
             OutputIteratorCells,
	     Emptyset_iterator> t = find_conflicts(p, c, bfit, cit,
		                                   Emptyset_iterator());
      return std::make_pair(t.first, t.second);
  }

  //TODO: integrate into P3DT3, it almost corresponds to periodic_circumcenter
  Point canonicalize_point(const Point& p) const {
    Iso_cuboid dom = this->domain();
    if (   !(p.x() < dom.xmin()) && p.x()<dom.xmax()
      && !(p.y() < dom.ymin()) && p.y()<dom.ymax()
      && !(p.z() < dom.zmin()) && p.z()<dom.zmax() )
      return p;

    int ox=-1, oy=-1, oz=-1;
    if (p.x() < dom.xmin()) ox = 1;
    else if (p.x() < dom.xmax()) ox = 0;
    if (p.y() < dom.ymin()) oy = 1;
    else if (p.y() < dom.ymax()) oy = 0;
    if (p.z() < dom.zmin()) oz = 1;
    else if (p.z() < dom.zmax()) oz = 0;
    Offset transl_offx(0,0,0);
    Offset transl_offy(0,0,0);
    Offset transl_offz(0,0,0);
    Point dp(p);

    // Find the right offset such that the translation will yield a
    // point inside the original domain.
    while ( dp.x() < dom.xmin() || !(dp.x() < dom.xmax()) ) {
      transl_offx.x() = transl_offx.x() + ox;
      dp = point(std::make_pair(p,transl_offx));
    }
    while ( dp.y() < dom.ymin() || !(dp.y() < dom.ymax()) ) {
      transl_offy.y() = transl_offy.y() + oy;
      dp = point(std::make_pair(p,transl_offy));
    }
    while ( dp.z() < dom.zmin() || !(dp.z() < dom.zmax()) ) {
      transl_offz.z() = transl_offz.z() + oz;
      dp = point(std::make_pair(p,transl_offz));
    }

    Offset transl_off(transl_offx.x(),transl_offy.y(),transl_offz.z());
    Periodic_point ppp(std::make_pair(p,transl_off));
    CGAL_triangulation_assertion_code(Point rv(point(ppp));)
      CGAL_triangulation_assertion( !(rv.x()<dom.xmin()) && rv.x()<dom.xmax() );
    CGAL_triangulation_assertion( !(rv.y()<dom.ymin()) && rv.y()<dom.ymax() );
    CGAL_triangulation_assertion( !(rv.z()<dom.zmin()) && rv.z()<dom.zmax() );
    return point(ppp);
  }

  void set_domain(const Iso_cuboid & domain)
  {
    Base::set_domain(domain);
    this->insert_dummy_points();
  }

  using Base::tetrahedron;

  Tetrahedron tetrahedron(const Cell_handle c) const
  {
    Periodic_tetrahedron pt =  periodic_tetrahedron(c);
    return tetrahedron(pt);
  }

  template<class OutputIterator>
  OutputIterator
    finite_incident_edges(Vertex_handle v, OutputIterator edges) const
  {
    return incident_edges(v, edges);
  }

  template<class OutputIterator>
  OutputIterator
    finite_incident_cells(Vertex_handle v, OutputIterator cells) const
  {
    return incident_cells(v, cells);
  }

  template<class OutputIterator>
  OutputIterator
    finite_incident_facets(Vertex_handle v, OutputIterator facets) const
  {
    return incident_facets(v, facets);
  }
 
  Bounded_side side_of_power_sphere(const Cell_handle& c, const Point& p, bool perturb = false) const
  {
    Point point = this->canonicalize_point(p);
    return Base::side_of_power_sphere(c, point, Offset(), perturb);
  }
    
  Vertex_handle nearest_power_vertex(const Bare_point& p, Cell_handle start) const
  {
    return Base::nearest_power_vertex(p, start);
    
    // not yet implemented
    assert(false);
  }
  
  using Base::locate;

  Cell_handle locate(const Point& p, Locate_type& l, int& i, int& j, Cell_handle start = Cell_handle(), bool* could_lock_zone = NULL)
  {
    assert(could_lock_zone == NULL);
    return Base::locate(p,l,i,j,start);
  }
    
  Cell_handle locate(const Point & p, Vertex_handle hint) const
  {
    return Base::locate(p, hint == Vertex_handle() ? infinite_cell() : hint->cell());

    assert(false); // not yet supported
  }
   
  // Returns the vertices on the interior of the conflict hole.
  template <class OutputIterator>
  OutputIterator
  vertices_inside_conflict_zone(const /*Weighted_point*/Point& p, Cell_handle c,
                                OutputIterator res) const
  {
    return res;
    
    assert(false); // not yet supported
  }
    
  Cell_handle infinite_cell() const
  {
    assert(false);
    return Cell_handle();
  }
  
  Vertex_handle infinite_vertex() const
  {
    assert(false);
    return Vertex_handle();
  }

};

namespace details {

  template<typename K>
  struct Periodic_mesh_geom_traits_generator
  {
  private:
    //TODO: avoid this trick
    //class K2 : public K {};

    typedef Robust_weighted_circumcenter_filtered_traits_3<K> K3;
    //typedef Regular_triangulation_euclidean_traits_3<K> K3;
    typedef Periodic_3_regular_triangulation_traits_3<K3> K4;
    //typedef Robust_weighted_circumcenter_filtered_traits_3<K3> Geom_traits;
    class K5 : public K4
    {
    public:
      typedef K5 Self;
      typedef typename  K4::Weighted_point_3 Weighted_point_3;
      typedef Regular_traits_with_offsets_adaptor<Self, typename K3::In_smallest_orthogonal_sphere_3> In_smallest_orthogonal_sphere_3;
      typedef Regular_traits_with_offsets_adaptor<Self, typename K3::Side_of_bounded_orthogonal_sphere_3> Side_of_bounded_orthogonal_sphere_3;
      typedef Regular_traits_with_offsets_adaptor<Self, typename K3::Does_simplex_intersect_dual_support_3> Does_simplex_intersect_dual_support_3;
      typedef Regular_traits_with_offsets_adaptor<Self, typename K3::Construct_weighted_circumcenter_3> Construct_weighted_circumcenter_3;
      typedef Regular_traits_with_offsets_adaptor<Self, typename K3::Compute_squared_radius_smallest_orthogonal_sphere_3> Compute_squared_radius_smallest_orthogonal_sphere_3;
      typedef Regular_traits_with_offsets_adaptor<Self, typename K3::Compute_power_product_3> Compute_power_product_3;
      typedef Regular_traits_with_offsets_adaptor<Self, typename K3::Compute_critical_squared_radius_3> Compute_critical_squared_radius_3;
      typedef Regular_traits_with_offsets_adaptor<Self, typename K3::Compare_weighted_squared_radius_3> Compare_weighted_squared_radius_3;
//      typedef Regular_traits_with_offsets_adaptor<Self, typename K3::Construct_circumcenter_3> Construct_circumcenter_3;

//      Construct_circumcenter_3 construct_circumcenter_3_object () const
//      {
//        return Construct_circumcenter_3(&_domain);
//      }

      In_smallest_orthogonal_sphere_3 in_smallest_orthogonal_sphere_3_object () const
      {
        return In_smallest_orthogonal_sphere_3(&this->_domain);
      }

      Side_of_bounded_orthogonal_sphere_3 side_of_bounded_orthogonal_sphere_3_object () const
      {
        return Side_of_bounded_orthogonal_sphere_3(&this->_domain);
      }

      Does_simplex_intersect_dual_support_3 does_simplex_intersect_dual_support_3_object () const
      {
        return Does_simplex_intersect_dual_support_3(&this->_domain);
      }

      Construct_weighted_circumcenter_3 construct_weighted_circumcenter_3_object () const
      {
        return Construct_weighted_circumcenter_3(&this->_domain);
      }

      Compute_power_product_3 compute_power_product_3_object () const
      {
        return Compute_power_product_3(&this->_domain);
      }

      Compute_squared_radius_smallest_orthogonal_sphere_3 compute_squared_radius_smallest_orthogonal_sphere_3_object () const
      {
        return Compute_squared_radius_smallest_orthogonal_sphere_3(&this->_domain);
      }

      Compute_critical_squared_radius_3 compute_critical_squared_radius_3_object () const
      {
        return Compute_critical_squared_radius_3(&this->_domain);
      }

      Compare_weighted_squared_radius_3 compare_weighted_squared_radius_3_object () const
      {
        return Compare_weighted_squared_radius_3(&this->_domain);
      }
    };
    typedef K5 Geom_traits;
    
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
  typedef Periodic_3_Delaunay_triangulation_3_Mesher_3<Geom_traits, Tds> Triangulation;

public:
  typedef Triangulation type;
  typedef type Type;
};  // end struct Mesh_triangulation_3

// helper function moving periodic
// triangles in some canonical positions in order to get a surface with
// less "holes"
template<class Triangulation>
typename Triangulation::Periodic_triangle canonicalize(const typename Triangulation::Periodic_triangle& pt) {
  typedef typename Triangulation::Offset Offset;

  Offset o0 = pt[0].second;
  Offset o1 = pt[1].second;
  Offset o2 = pt[2].second;
  int diffx = std::min(o0.x(),std::min(o1.x(),o2.x()));
  int diffy = std::min(o0.y(),std::min(o1.y(),o2.y()));
  int diffz = std::min(o0.z(),std::min(o1.z(),o2.z()));
  Offset diff_off(diffx,diffy,diffz);

  return make_array(
    std::make_pair(pt[0].first,o0 - diff_off),
    std::make_pair(pt[1].first,o1 - diff_off),
    std::make_pair(pt[2].first,o2 - diff_off));
}

template<class Triangulation>
typename Triangulation::Periodic_tetrahedron canonicalize_tetrahedron(const typename Triangulation::Periodic_tetrahedron& pt) {
  typedef typename Triangulation::Offset Offset;

  Offset o0 = pt[0].second;
  Offset o1 = pt[1].second;
  Offset o2 = pt[2].second;
  Offset o3 = pt[3].second;

  int diffx = std::min(std::min(o0.x(),o1.x()),std::min(o2.x(),o3.x()));
  int diffy = std::min(std::min(o0.y(),o1.y()),std::min(o2.y(),o3.y()));
  int diffz = std::min(std::min(o0.z(),o1.z()),std::min(o2.z(),o3.z()));
  Offset diff_off(diffx,diffy,diffz);

  return make_array(
    std::make_pair(pt[0].first,o0 - diff_off),
    std::make_pair(pt[1].first,o1 - diff_off),
    std::make_pair(pt[2].first,o2 - diff_off),
    std::make_pair(pt[3].first,o3 - diff_off));
}
  
template <class Stream, class C3t3>
Stream &vertices_medit(Stream &out, C3t3 &c3t3) {

  out << std::setprecision(20);
  out << "MeshVersionFormatted 1" << std::endl;
  out << "Dimension 3" << std::endl;
  out << "Vertices " << c3t3.triangulation().nb_of_vertices() * 8 << std::endl;
}

// Writing a restricted Delaunay triangulation to the .mesh file format
// Writing the triangulation to 8 domains.
// Can be used with the medit (a viewer)
template <class Stream, class C3t3>
Stream &write_complex_to_medit(Stream &out, C3t3 &c3t3, unsigned occurence_count = 8) {
  typedef typename C3t3::Triangulation Triangulation;
  typedef Triangulation Tr;
  
  typedef typename Triangulation::Geom_traits Gt;
  
  typedef typename Triangulation::Iso_cuboid Iso_cuboid;
  
  typedef typename Triangulation::Triangle Triangle;
  typedef typename Triangulation::Tetrahedron Tetrahedron;
  
  typedef typename C3t3::Facet_iterator Facet_iterator;
  typedef typename C3t3::Cell_iterator Cell_iterator;
  
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
      //<< 128 * j + 64 + 1 << std::endl;
      
      it++;
    }
  }
  
  out << "0\nEnd";     
  return out;
}
  
} //namespace CGAL

#endif // CGAL_MESH_3_PERIODIC_TRIANGULATION_3_H

