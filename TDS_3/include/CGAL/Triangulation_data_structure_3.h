// Copyright (c) 1999-2005  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Sylvain Pion

// combinatorial triangulation of the boundary of a polytope
// of dimension d in dimension d+1
// for -1 <= d <= 3

#ifndef CGAL_TRIANGULATION_DATA_STRUCTURE_3_H
#define CGAL_TRIANGULATION_DATA_STRUCTURE_3_H

#include <CGAL/license/TDS_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/basic.h>

#include <utility>
#include <map>
#include <set>
#include <vector>
#include <stack>

#include <boost/unordered_set.hpp>
#include <CGAL/utility.h>
#include <CGAL/iterator.h>
#include <CGAL/internal/Has_member_visited.h>

#include <CGAL/Unique_hash_map.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_utils_3.h>

#include <CGAL/Concurrent_compact_container.h>
#include <CGAL/Compact_container.h>

#include <CGAL/Triangulation_ds_cell_base_3.h>
#include <CGAL/Triangulation_ds_vertex_base_3.h>
#include <CGAL/Triangulation_simplex_3.h>

#include <CGAL/internal/Triangulation_ds_iterators_3.h>
#include <CGAL/internal/Triangulation_ds_circulators_3.h>

#ifdef CGAL_LINKED_WITH_TBB
#  include <tbb/scalable_allocator.h>
#endif

#include <boost/type_traits/is_convertible.hpp>

namespace CGAL {

// TODO : noms : Vb != Vertex_base : clarifier.

template < class Vb = Triangulation_ds_vertex_base_3<>,
           class Cb = Triangulation_ds_cell_base_3<>,
           class Concurrency_tag_ = Sequential_tag
>
class Triangulation_data_structure_3
  : public Triangulation_utils_3
{
  typedef Triangulation_data_structure_3<Vb, Cb, Concurrency_tag_> Tds;

public:
  typedef Concurrency_tag_            Concurrency_tag;

  // This tag is used in the parallel operations of RT_3 to access some functions
  // of the TDS (tds.vertices().is_used(Vertex_handle)) that are much more efficient
  // than what is exposed by the TDS concept (tds.is_vertex(Vertex_handle)).
  typedef CGAL::Tag_true              Is_CGAL_TDS_3;

  // Tools to change the Vertex and Cell types of the TDS.
  template < typename Vb2 >
  struct Rebind_vertex {
    typedef Triangulation_data_structure_3<Vb2, Cb, Concurrency_tag> Other;
  };

  template < typename Cb2 >
  struct Rebind_cell {
    typedef Triangulation_data_structure_3<Vb, Cb2, Concurrency_tag> Other;
  };

  // Put this TDS inside the Vertex and Cell types.
  typedef typename Vb::template Rebind_TDS<Tds>::Other  Vertex;
  typedef typename Cb::template Rebind_TDS<Tds>::Other  Cell;

  class Cell_data {
    unsigned char conflict_state;
  public:
    Cell_data() : conflict_state(0) {}

    void clear()            { conflict_state = 0; }
    void mark_in_conflict() { conflict_state = 1; }
    void mark_on_boundary() { conflict_state = 2; }
    void mark_processed()   { conflict_state = 1; }

    bool is_clear()       const { return conflict_state == 0; }
    bool is_in_conflict() const { return conflict_state == 1; }
    bool is_on_boundary() const { return conflict_state == 2; }
    bool processed() const { return conflict_state == 1; }
  };

private:

  friend class internal::Triangulation_ds_facet_iterator_3<Tds>;
  friend class internal::Triangulation_ds_edge_iterator_3<Tds>;

  friend class internal::Triangulation_ds_cell_circulator_3<Tds>;
  friend class internal::Triangulation_ds_facet_circulator_3<Tds>;

public:
    
  // Cells
  // N.B.: Concurrent_compact_container requires TBB
#ifdef CGAL_LINKED_WITH_TBB
  typedef typename boost::mpl::if_c
  <
    boost::is_convertible<Concurrency_tag, Parallel_tag>::value,
    Concurrent_compact_container<Cell, tbb::scalable_allocator<Cell> >,
    Compact_container<Cell>
  >::type                                                Cell_range;

# else
  typedef Compact_container<Cell>                        Cell_range;
#endif

  // Vertices
  // N.B.: Concurrent_compact_container requires TBB
#ifdef CGAL_LINKED_WITH_TBB
  typedef typename boost::mpl::if_c
  <
    boost::is_convertible<Concurrency_tag, Parallel_tag>::value,
    Concurrent_compact_container<Vertex, tbb::scalable_allocator<Vertex> >,
    Compact_container<Vertex>
  >::type                                                Vertex_range;

# else
  typedef Compact_container<Vertex>                      Vertex_range;
#endif

  
  typedef typename Cell_range::size_type       size_type;
  typedef typename Cell_range::difference_type difference_type;

  typedef typename Cell_range::iterator        Cell_iterator;
  typedef typename Vertex_range::iterator      Vertex_iterator;

  typedef internal::Triangulation_ds_facet_iterator_3<Tds>   Facet_iterator;
  typedef internal::Triangulation_ds_edge_iterator_3<Tds>    Edge_iterator;

  typedef internal::Triangulation_ds_cell_circulator_3<Tds>  Cell_circulator;
  typedef internal::Triangulation_ds_facet_circulator_3<Tds> Facet_circulator;

  typedef Iterator_range<Prevent_deref<Cell_iterator> >    Cell_handles;
  typedef Iterator_range<Prevent_deref<Vertex_iterator> >  Vertex_handles;

  typedef Iterator_range<Facet_iterator> Facets;
  typedef Iterator_range<Edge_iterator> Edges;
  
//private: // In 2D only :
  typedef internal::Triangulation_ds_face_circulator_3<Tds>  Face_circulator;

  typedef Vertex_iterator                          Vertex_handle;
  typedef Cell_iterator                            Cell_handle;

  typedef std::pair<Cell_handle, int>              Facet;
  typedef Triple<Cell_handle, int, int>            Edge;

  typedef Triangulation_simplex_3<Tds>             Simplex;
//#ifndef CGAL_TDS_USE_RECURSIVE_CREATE_STAR_3
  //internally used for create_star_3 (faster than a tuple)
  struct iAdjacency_info{
    int v1;
    Cell_handle v2;
    int v3;
    Cell_handle v4;
    int v5;
    int v6;
    iAdjacency_info(){}
    iAdjacency_info(int a1,Cell_handle a2,int a3,Cell_handle a4,int a5 ,int a6):
      v1(a1),v2(a2),v3(a3),v4(a4),v5(a5),v6(a6) {}
    void update_variables(int& a1,Cell_handle& a2,int& a3,Cell_handle& a4,int& a5 ,int& a6)
    {
      a1=v1;
      a2=v2;
      a3=v3;
      a4=v4;
      a5=v5;
      a6=v6;
    }
  };
//#endif  
 

public:
  Triangulation_data_structure_3()
    : _dimension(-2)
  {}

  Triangulation_data_structure_3(const Tds & tds)
  {
    copy_tds(tds);
  }

  Tds & operator= (const Tds & tds)
  {
    if (&tds != this) {
      Tds tmp(tds);
      swap(tmp);
    }
    return *this;
  }

  size_type number_of_vertices() const { return vertices().size(); }

  int dimension() const {return _dimension;}

  size_type number_of_cells() const
    {
      if ( dimension() < 3 ) return 0;
      return cells().size();
    }

  size_type number_of_facets() const
    {
      if ( dimension() < 2 ) return 0;
      return std::distance(facets_begin(), facets_end());
    }

  size_type number_of_edges() const
    {
      if ( dimension() < 1 ) return 0;
      return std::distance(edges_begin(), edges_end());
    }

  // USEFUL CONSTANT TIME FUNCTIONS

  // SETTING

  void set_dimension(int n) { _dimension = n; }

  Vertex_handle create_vertex(const Vertex &v)
  {
      return vertices().insert(v);
  }

  Vertex_handle create_vertex()
  {
      return vertices().emplace();
  }

  Vertex_handle create_vertex(Vertex_handle v)
  {
      return create_vertex(*v);
  }

  Cell_handle create_cell(const Cell &c)
    {
      return cells().insert(c);
    }

  Cell_handle create_cell()
    {
      return cells().emplace();
    }

  Cell_handle create_cell(Cell_handle c)
    {
      return create_cell(*c);
    }

  Cell_handle create_cell(Vertex_handle v0, Vertex_handle v1,
                          Vertex_handle v2, Vertex_handle v3)
    {
      return cells().emplace(v0, v1, v2, v3);
    }

  Cell_handle create_cell(Vertex_handle v0, Vertex_handle v1,
                          Vertex_handle v2, Vertex_handle v3,
                          Cell_handle n0, Cell_handle n1,
                          Cell_handle n2, Cell_handle n3)
    {
      return cells().emplace(v0, v1, v2, v3, n0, n1, n2, n3);
    }

  Cell_handle create_face()
    {
      CGAL_triangulation_precondition(dimension()<3);
      return create_cell();
    }

  Cell_handle create_face(Vertex_handle v0, Vertex_handle v1,
                          Vertex_handle v2)
    {
      CGAL_triangulation_precondition(dimension()<3);
      return cells().emplace(v0, v1, v2, Vertex_handle());
    }

  // The following functions come from TDS_2.
  Cell_handle create_face(Cell_handle f0, int i0,
                          Cell_handle f1, int i1,
                          Cell_handle f2, int i2)
    {
      CGAL_triangulation_precondition(dimension() <= 2);
      Cell_handle newf = create_face(vertex(f0, cw(i0)),
                                     vertex(f1, cw(i1)),
                                     vertex(f2, cw(i2)));
      set_adjacency(newf, 2, f0, i0);
      set_adjacency(newf, 0, f1, i1);
      set_adjacency(newf, 1, f2, i2);
      return newf;
    }

  Cell_handle create_face(Cell_handle f0, int i0,
                          Cell_handle f1, int i1)
    {
      CGAL_triangulation_precondition(dimension() <= 2);
      Cell_handle newf = create_face(vertex(f0, cw(i0)),
                                     vertex(f1, cw(i1)),
                                     vertex(f1, ccw(i1)));
      set_adjacency(newf, 2, f0, i0);
      set_adjacency(newf, 0, f1, i1);
      return newf;
    }

  Cell_handle create_face(Cell_handle f, int i, Vertex_handle v)
    {
      CGAL_triangulation_precondition(dimension() <= 2);
      Cell_handle newf = create_face(vertex(f, cw(i)),
                                     vertex(f, ccw(i)),
                                     v);
      set_adjacency(newf, 2, f, i);
      return newf;
    }

  // not documented
  void read_cells(std::istream& is, const std::vector< Vertex_handle > &V,
                  std::size_t & m, std::vector< Cell_handle > &C);
  // not documented
  void print_cells(std::ostream& os,
                   const Unique_hash_map<Vertex_handle, std::size_t> &V ) const;

  // ACCESS FUNCTIONS

  void delete_vertex( Vertex_handle v )
  {
      CGAL_triangulation_expensive_precondition( is_vertex(v) );
      vertices().erase(v);
  }

  void delete_cell( Cell_handle c )
  {
      CGAL_triangulation_expensive_precondition( is_simplex(c) );
      cells().erase(c);
  }

  template <class InputIterator>
  void delete_vertices(InputIterator begin, InputIterator end)
  {
      for(; begin != end; ++begin)
          delete_vertex(*begin);
  }

  template <class InputIterator>
  void delete_cells(InputIterator begin, InputIterator end)
  {
      for(; begin != end; ++begin)
          delete_cell(*begin);
  }

  // QUERIES

  bool is_simplex(Cell_handle c) const; // undocumented for now
  bool is_vertex(Vertex_handle v) const;
  bool is_edge(Cell_handle c, int i, int j) const;
  bool is_edge(Vertex_handle u, Vertex_handle v, Cell_handle & c,
               int & i, int & j) const;
  bool is_edge(Vertex_handle u, Vertex_handle v) const;
  bool is_facet(Cell_handle c, int i) const;
  bool is_facet(Vertex_handle u, Vertex_handle v,
                Vertex_handle w,
                Cell_handle & c, int & i, int & j, int & k) const;
  bool is_cell(Cell_handle c) const;
  bool is_cell(Vertex_handle u, Vertex_handle v,
               Vertex_handle w, Vertex_handle t,
               Cell_handle & c, int & i, int & j, int & k, int & l) const;
  bool is_cell(Vertex_handle u, Vertex_handle v,
               Vertex_handle w, Vertex_handle t) const;

  bool has_vertex(const Facet & f, Vertex_handle v, int & j) const;
  bool has_vertex(Cell_handle c, int i,
                  Vertex_handle v, int & j) const;
  bool has_vertex(const Facet & f, Vertex_handle v) const;
  bool has_vertex(Cell_handle c, int i, Vertex_handle v) const;

  bool are_equal(Cell_handle c, int i,
                 Cell_handle n, int j) const;
  bool are_equal(const Facet & f, const Facet & g) const;
  bool are_equal(const Facet & f, Cell_handle n, int j) const;


  // Functions that were member functions of the vertex and cell class
  
  Cell_data& tds_data(Cell_handle c)
  {
    return c->itds_data();
  }
  
  const Cell_data& tds_data(Cell_handle c) const
  {
    return c->itds_data(); }


  bool& visited_for_vertex_extractor(Vertex_handle v)
  {
    return v->visited_for_vertex_extractor;
  }
  
  Vertex_handle vertex(Cell_handle c, int i) const
  {
    CGAL_triangulation_precondition( i >= 0 && i <= 3 );
    CGAL_assume( i >= 0 && i <= 3 );
    //return V[i];
    return c->ivertex(i);
  }

  
  int index(Cell_handle c, Vertex_handle v) const
  {
    /*
    if (v == V[0]) { return 0; }
    if (v == V[1]) { return 1; }
    if (v == V[2]) { return 2; }
    CGAL_triangulation_assertion( v == V[3] );
    return 3;
    */
    return c->iindex(v);
  }


  int index(Cell_handle c, Cell_handle n) const
  {
    /*
    if (n == N[0]) return 0;
    if (n == N[1]) return 1;
    if (n == N[2]) return 2;
    CGAL_triangulation_assertion( n == N[3] );
    return 3;
    */
    return c->iindex(n);
  }

  Cell_handle neighbor(Cell_handle c, int i) const
  {
    /*
    CGAL_triangulation_precondition( i >= 0 && i <= 3);
    return N[i];
    */
    return c->ineighbor(i);
  }

  
  bool has_vertex(Cell_handle c, Vertex_handle v) const
  {
    return c->ihas_vertex(v);
  }

  bool has_vertex(Cell_handle c, Vertex_handle v, int& i) const
  {
    return c->ihas_vertex(v,i);
  }

  
  void set_vertex(Cell_handle c, int i, Vertex_handle v)
  {
    /*
    CGAL_triangulation_precondition( i >= 0 && i <= 3);
    V[i] = v;
    */
    c->iset_vertex(i,v);
  }
  
  void set_vertices(Cell_handle c)
  {
    c->iset_vertices();
  }
  
  void set_vertices(Cell_handle c,
                    Vertex_handle v0, Vertex_handle v1,
                    Vertex_handle v2, Vertex_handle v3)
  {
    c->iset_vertices(v0,v1,v2,v3);
  }
  
  void set_neighbor(Cell_handle c, int i, Cell_handle n)
  {
    /*
    CGAL_triangulation_precondition( i >= 0 && i <= 3);
    CGAL_triangulation_precondition( this != n.operator->() );
    N[i] = n;
    */
    c->iset_neighbor(i,n);
  }

  void set_neighbors(Cell_handle c)
  {
    c->iset_neighbors();
  }
  
  void set_neighbors(Cell_handle c,
                      Cell_handle n0, Cell_handle n1,
                      Cell_handle n2, Cell_handle n3)
  {
    c->iset_neighbors(n0,n1,n2,n3);
  }
  
  bool has_neighbor(Cell_handle c, Cell_handle n) const
  {
    return c->ihas_neighbor(n);
  }

  bool has_neighbor(Cell_handle c, Cell_handle n, int& i) const
  {
    return c->ihas_neighbor(n,i);
  }

  
  Cell_handle cell(Vertex_handle v) const 
  { return v->icell(); }  

  
  void set_cell(Vertex_handle v, Cell_handle c)
  {
    v->iset_cell(c);
  }
  
  // MODIFY

  bool flip(Cell_handle c, int i);
  bool flip(const Facet &f)
  { return flip( f.first, f.second); }

  void flip_flippable(Cell_handle c, int i);
  void flip_flippable(const Facet &f)
  { flip_flippable( f.first, f.second ); }

  bool flip(Cell_handle c, int i, int j);
  bool flip(const Edge &e)
  { return flip( e.first, e.second, e.third ); }

  void flip_flippable(Cell_handle c, int i, int j);
  void flip_flippable(const Edge &e)
  { flip_flippable( e.first, e.second, e.third ); }

private:
  // common to flip and flip_flippable
  void flip_really(Cell_handle c, int i, Cell_handle n, int in);
  void flip_really(Cell_handle c, int i, int j,
                   Cell_handle c1, Vertex_handle v1,
                   int i1, int j1, int next1,
                   Cell_handle c2, Vertex_handle v2,
                   int i2, int j2, int next2,
                   Vertex_handle v3);

  #ifdef CGAL_TDS_USE_RECURSIVE_CREATE_STAR_3
  Cell_handle create_star_3(Vertex_handle v, Cell_handle c,
                            int li, int prev_ind2 = -1);
  #else
  Cell_handle recursive_create_star_3(Vertex_handle v, Cell_handle c, int li, int prev_ind2,int depth);
  Cell_handle non_recursive_create_star_3(Vertex_handle v, Cell_handle c, int li, int prev_ind2);

  Cell_handle create_star_3(Vertex_handle v, Cell_handle c,
                            int li, int prev_ind2 = -1)
  {
    return recursive_create_star_3(v,c,li,prev_ind2,0);
  }
  #endif

  Cell_handle create_star_2(Vertex_handle v,
                            Cell_handle c, int li);

public:

  // Internal function : assumes the conflict cells are marked.
  template <class CellIt>
  Vertex_handle _insert_in_hole(CellIt cell_begin, CellIt cell_end,
                                Cell_handle begin, int i,
                                Vertex_handle newv)
  {
      CGAL_triangulation_precondition(begin != Cell_handle());
      // if begin == nullptr (default arg), we could compute one by walking in
      // CellIt.  At the moment, the functionality is not available, you have
      // to specify a starting facet.

      Cell_handle cnew;
      if (dimension() == 3)
          cnew = create_star_3(newv, begin, i);
      else
          cnew = create_star_2(newv, begin, i);

      set_cell(newv,cnew);
      delete_cells(cell_begin, cell_end);
      return newv;
  }

  // Internal function : assumes the conflict cells are marked.
  template <class CellIt>
  Vertex_handle _insert_in_hole(CellIt cell_begin, CellIt cell_end,
                                Cell_handle begin, int i)
  {
      return _insert_in_hole(cell_begin, cell_end, begin, i, create_vertex());
  }

  // Mark the cells in conflict, then calls the internal function.
  template <class CellIt>
  Vertex_handle insert_in_hole(CellIt cell_begin, CellIt cell_end,
                               Cell_handle begin, int i,
                               Vertex_handle newv)
  {
      for (CellIt cit = cell_begin; cit != cell_end; ++cit)
          tds_data(*cit).mark_in_conflict();

      return _insert_in_hole(cell_begin, cell_end, begin, i, newv);
  }

  // Mark the cells in conflict, then calls the internal function.
  template <class CellIt>
  Vertex_handle insert_in_hole(CellIt cell_begin, CellIt cell_end,
                               Cell_handle begin, int i)
  {
      return insert_in_hole(cell_begin, cell_end, begin, i, create_vertex());
  }

  //INSERTION
  
  // Create a finite cell with v1, v2, v3 and v4
  // Precondition: v1, v2, v3 and v4 MUST BE positively oriented
  Vertex_handle insert_first_finite_cell(
    Vertex_handle &v1, Vertex_handle &v2, Vertex_handle &v3, Vertex_handle &v4,
    Vertex_handle v_infinite = Vertex_handle());

  Vertex_handle insert_in_cell(Cell_handle c);

  Vertex_handle insert_in_facet(const Facet & f)
    { return insert_in_facet(f.first, f.second); }

  Vertex_handle insert_in_facet(Cell_handle c, int i);

  Vertex_handle insert_in_edge(const Edge & e)
    { return insert_in_edge(e.first, e.second, e.third); }

  Vertex_handle insert_in_edge(Cell_handle c, int i, int j);

  Vertex_handle insert_increase_dimension(Vertex_handle star =Vertex_handle());

  // REMOVAL

private:
  Cell_handle remove_degree_4(Vertex_handle v);
  Cell_handle remove_degree_3(Vertex_handle v);
  Cell_handle remove_degree_2(Vertex_handle v);
public:
  Cell_handle remove_from_maximal_dimension_simplex(Vertex_handle v);
  void remove_decrease_dimension(Vertex_handle v)
  {
      remove_decrease_dimension (v, v);
  }
  void remove_decrease_dimension(Vertex_handle v, Vertex_handle w);
  void decrease_dimension(Cell_handle f, int i);

  // Change orientation of the whole TDS.
  void reorient()
  {
      CGAL_triangulation_precondition(dimension() >= 1);
      for (Cell_iterator i = cells().begin();
              i != cells().end(); ++i)
          change_orientation(i);
  }

  // ITERATOR METHODS

  Cell_iterator cells_begin() const
  {
    if ( dimension() < 3 )
        return cells_end();
    return cells().begin();
  }

  Cell_iterator cells_end() const
  {
    return cells().end();
  }

  Cell_handles cell_handles() const
  {
    return make_prevent_deref_range(cells_begin(), cells_end()); 
  }
  
  Cell_iterator raw_cells_begin() const
  {
    return cells().begin();
  }

  Cell_iterator raw_cells_end() const
  {
    return cells().end();
  }

  Facet_iterator facets_begin() const
  {
    if ( dimension() < 2 )
        return facets_end();
    return Facet_iterator(this);
  }

  Facet_iterator facets_end() const
  {
    return Facet_iterator(this, 1);
  }

  Facets facets() const
  {
    return Facets(facets_begin(), facets_end());
  }
  
  Edge_iterator edges_begin() const
  {
    if ( dimension() < 1 )
        return edges_end();
    return Edge_iterator(this);
  }

  Edge_iterator edges_end() const
  {
    return Edge_iterator(this,1);
  }

  Edges edges() const
  {
    return Edges(edges_begin(), edges_end());
  }
  
  Vertex_iterator vertices_begin() const
  {
    return vertices().begin();
  }

  Vertex_iterator vertices_end() const
  {
    return vertices().end();
  }

  Vertex_handles vertex_handles() const
  {
    return make_prevent_deref_range(vertices_begin(), vertices_end()); 
  }
  
  // CIRCULATOR METHODS

  // cells around an edge
  Cell_circulator incident_cells(const Edge & e) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    return Cell_circulator(this, e);
  }
  Cell_circulator incident_cells(Cell_handle ce, int i, int j) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    return Cell_circulator(this, ce, i, j);
  }

  Cell_circulator incident_cells(const Edge &e, Cell_handle start) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    return Cell_circulator(this, e, start);
  }
  Cell_circulator incident_cells(Cell_handle ce, int i, int j,
                                 Cell_handle start) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    return Cell_circulator(this, ce, i, j, start);
  }

  //facets around an edge
  Facet_circulator incident_facets(const Edge & e) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    return Facet_circulator(this, e);
  }
  Facet_circulator incident_facets(Cell_handle ce, int i, int j) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    return Facet_circulator(this, ce, i, j);
  }
  Facet_circulator incident_facets(const Edge & e, const Facet & start) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    return Facet_circulator(this, e, start);
  }
  Facet_circulator incident_facets(Cell_handle ce, int i, int j,
                                   const Facet & start) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    return Facet_circulator(this, ce, i, j, start);
  }
  Facet_circulator incident_facets(const Edge & e,
                                   Cell_handle start, int f) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    return Facet_circulator(this, e, start, f);
  }
  Facet_circulator incident_facets(Cell_handle ce, int i, int j,
                                   Cell_handle start, int f) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    return Facet_circulator(this, ce, i, j, start, f);
  }

  // 2D : circulates on the faces adjacent to a vertex.
  Face_circulator incident_faces(Vertex_handle v) const
  {
    CGAL_triangulation_precondition( dimension() == 2 );
    return Face_circulator(this, v, cell(v));
  }

  // around a vertex
private:

  template <class IncidentCellIterator, class IncidentFacetIterator>
  std::pair<IncidentCellIterator, IncidentFacetIterator>
  incident_cells_3(Vertex_handle v, Cell_handle d,
                                 std::pair<IncidentCellIterator,
                                 IncidentFacetIterator> it) const
  {
        CGAL_triangulation_precondition(dimension() == 3);
        
        std::stack<Cell_handle> cell_stack;
        cell_stack.push(d);
        const_cast<Tds*>(this)->tds_data(d).mark_in_conflict();
        *it.first++ = d;
        
        do {
                Cell_handle c = cell_stack.top();
                cell_stack.pop();
                
                for (int i=0; i<4; ++i) {
                  if (vertex(c, i) == v)
                                continue;
                        Cell_handle next = neighbor(c,i);
                        if (c < next)
                                *it.second++ = Facet(c, i); // Incident facet.
                        if (! tds_data(next).is_clear())
                                continue;
                        cell_stack.push(next);
                        const_cast<Tds*>(this)->tds_data(next).mark_in_conflict();
                        *it.first++ = next;
                }
        } while(!cell_stack.empty());

        return it;
  }

  template <class IncidentFacetIterator>
  void
  incident_cells_3_threadsafe(Vertex_handle v, Cell_handle d,
                              std::vector<Cell_handle> &cells,
                              IncidentFacetIterator facet_it) const
  {
    boost::unordered_set<Cell_handle, Handle_hash_function> found_cells;

    cells.push_back(d);
    found_cells.insert(d);
    int head=0;
    int tail=1;
    do {
      Cell_handle c = cells[head];

      for (int i=0; i<4; ++i) {
        if (vertex(c, i) == v)
          continue;
        Cell_handle next = neighbor(c, i);
        if (c < next)
          *facet_it++ = Facet(c, i); // Incident facet
        if (! found_cells.insert(next).second )
          continue;
        cells.push_back(next);
        ++tail;
      }
      ++head;
    } while(head != tail);
  }
  
  void just_incident_cells_3(Vertex_handle v,
                             std::vector<Cell_handle>& cells) const
  {
    CGAL_triangulation_precondition(dimension() == 3);

    Cell_handle d = cell(v);
    cells.push_back(d);
    const_cast<Tds*>(this)->tds_data(d).mark_in_conflict(); // AF: why did this work before
    int head=0;
    int tail=1;
    do {
      Cell_handle c = cells[head];

      for (int i=0; i<4; ++i) {
        if (vertex(c, i) == v)
          continue;
        Cell_handle next = neighbor(c, i);
        if (! tds_data(next).is_clear())
          continue;
        cells.push_back(next);
        ++tail;
        const_cast<Tds*>(this)->tds_data(next).mark_in_conflict();
      }
      ++head;
    } while(head != tail);
  }

  template <class OutputIterator>
  void
  incident_cells_2(Vertex_handle v, Cell_handle,
                   OutputIterator cells) const
  {
      CGAL_triangulation_precondition(dimension() == 2);

      Face_circulator fc = incident_faces(v);
      Face_circulator done(fc);
      do {
        *cells++ = fc;
        ++fc;
      } while (fc != done);
  }

public:

  class False_filter {
    public:
    False_filter() {}
    template<class T>
    bool operator() (T) {
      return false;
    }
  };

        // Visitor for visit_incident_cells:
        // outputs the facets
  template <class OutputIterator, class Filter>
  class Facet_extractor {
    OutputIterator output;
    Filter filter;
    public:
    Facet_extractor(Vertex_handle, OutputIterator _output, const Tds*, Filter _filter):
    output(_output), filter(_filter){}
    void operator() (Cell_handle) {}

    OutputIterator result() {
      return output;
    }

    class Facet_it {
      OutputIterator& output;
      Filter& filter;
      public:
      Facet_it(OutputIterator& _output, Filter& _filter): output(_output), filter(_filter) {}
      Facet_it& operator*() {return *this;};
      Facet_it& operator++() {return *this;};
      Facet_it operator++(int) {return *this;};
      template<class T>
      Facet_it& operator=(const T& e) {
        if(filter(e))
          return *this;
        *output++ = e;
        return *this;
      }
      Facet_it& operator=(const Facet_it& f) {
        output = f.output;
        filter = f.filter;
        return *this;
      }  
      Facet_it(const Facet_it&)=default;
    };
    Facet_it facet_it() {
      return Facet_it(output, filter);
    }
  };

        // Visitor for visit_incident_cells:
        // outputs the cells
  template <class OutputIterator, class Filter>
  class Cell_extractor {
    OutputIterator output;
    Filter filter;
    public:
    Cell_extractor(Vertex_handle, OutputIterator _output, const Tds*, Filter _filter):
    output(_output), filter(_filter) {}

    void operator()(Cell_handle c) {
      if(filter(c))
        return;
      *output++ = c;
    }
    CGAL::Emptyset_iterator facet_it() {return CGAL::Emptyset_iterator();}
    OutputIterator result() {
      return output;
    }
  };

  // Visitor for visit_incident_cells:
  // WARNING: 2D ONLY
  // outputs the faces obtained as degenerated cells
  template <class OutputIterator, class Filter>
  class DegCell_as_Facet_extractor {
    OutputIterator output;
    Filter filter;
    public:
    DegCell_as_Facet_extractor(Vertex_handle, OutputIterator _output, const Tds*, Filter _filter):
    output(_output), filter(_filter) {}

    void operator()(Cell_handle c) {
      Facet f = Facet(c,3);
      if(filter(f))
        return;
      *output++ = f;
    }
    CGAL::Emptyset_iterator facet_it() {return CGAL::Emptyset_iterator();}
    OutputIterator result() {
      return output;
    }
  };

  template<class Treatment, class OutputIterator, class Filter, bool hasVisited>
  class Vertex_extractor;

        // Visitor for visit_incident_cells:
        // outputs the result of Treatment applied to the vertices
  template<class Treatment, class OutputIterator, class Filter>
  class Vertex_extractor<Treatment,OutputIterator,Filter,false> {
    Vertex_handle v;

    boost::unordered_set<Vertex_handle, Handle_hash_function> tmp_vertices;

    Treatment treat;
    const Tds* t;
    Filter filter;
  public:
    Vertex_extractor(Vertex_handle _v, OutputIterator _output, const Tds* _t, Filter _filter):
      v(_v), treat(_output, _t), t(_t), filter(_filter) 
    {
#if ( BOOST_VERSION >= 105000 )
      tmp_vertices.reserve(64);
#endif
    }

    void operator()(Cell_handle c) {
      for (int j=0; j<= t->dimension(); ++j) {
        Vertex_handle w = vertex(c, j);
        if(filter(w))
          continue;
        if (w != v) {
          if(tmp_vertices.insert(w).second) {
            treat(c, v, j);
          }

        }
      }
    }


    CGAL::Emptyset_iterator facet_it() {return CGAL::Emptyset_iterator();}
    OutputIterator result() {
      return treat.result();
    }
  };

  template<class Treatment, class OutputIterator, class Filter>
  class Vertex_extractor<Treatment,OutputIterator,Filter,true> {
    Vertex_handle v;
    std::vector<Vertex_handle> tmp_vertices;

    Treatment treat;
    Tds* t;
    Filter filter;
  public:
    Vertex_extractor(Vertex_handle _v, OutputIterator _output, const Tds* _t, Filter _filter):
      v(_v), treat(_output, _t), t(const_cast<Tds*>(_t)), filter(_filter) {
      tmp_vertices.reserve(64);
    }

    void operator()(Cell_handle c) {
      for (int j=0; j<= t->dimension(); ++j) {
	Vertex_handle w = t->vertex(c, j);
	if(filter(w))
	  continue;
	if (w != v){

          if(! w->visited_for_vertex_extractor){
            t->visited_for_vertex_extractor(w) = true;
            tmp_vertices.push_back(w);
	    treat(c, v, j);
          }
        }
      }
    }

    ~Vertex_extractor()
    {
      for(std::size_t i=0; i < tmp_vertices.size(); ++i){
        t->visited_for_vertex_extractor(tmp_vertices[i]) = false;
      }
    }


    CGAL::Emptyset_iterator facet_it() {return CGAL::Emptyset_iterator();}
    OutputIterator result() {
      return treat.result();
    }
  };



  // Treatment for Vertex_extractor:
  // outputs the vertices
  template<class OutputIterator>
  class Vertex_feeder_treatment {
    OutputIterator output;
    const Tds* t;
  public:
    Vertex_feeder_treatment(OutputIterator _output, const Tds* t): output(_output), t(t) {};
    void operator()(Cell_handle c, Vertex_handle, int index) {
      *output++ = t->vertex(c,index);
    }
    OutputIterator result() {
      return output;
    }
  };

  // Treatment for Vertex_extractor:
  // outputs the edges corresponding to the vertices
  template<class OutputIterator>
  class Edge_feeder_treatment {
    OutputIterator output;
    const Tds* t;
  public:
    Edge_feeder_treatment(OutputIterator _output, const Tds* t): output(_output), t(t) {};
    void operator()(Cell_handle c, Vertex_handle v, int index) {
      *output++ = Edge(c, t->index(c,v), index);
    }
    OutputIterator result() {
      return output;
    }
  };

  template <class Filter, class OutputIterator>
  OutputIterator
  incident_cells(Vertex_handle v, OutputIterator cells, Filter f = Filter()) const
  {
    return visit_incident_cells<Cell_extractor<OutputIterator, Filter>,
      OutputIterator>(v, cells, f);
  }

  template <class OutputIterator>
  OutputIterator
  incident_cells(Vertex_handle v, OutputIterator cells) const
  {
    return incident_cells<False_filter>(v, cells);
  }

  // This version only works for vectors and only in 3D
  void incident_cells_3(Vertex_handle v,
                        std::vector<Cell_handle>& cells) const
  {
    just_incident_cells_3(v, cells);  
    typename std::vector<Cell_handle>::iterator cit,end;
    for(cit = cells.begin(), end = cells.end();
              cit != end;
              ++cit)
    {
      const_cast<Tds*>(this)->tds_data(*cit).clear(); // AF:  why did that work before
    }
  }
  
  template <class Filter, class OutputIterator>
  OutputIterator
  incident_cells_threadsafe(Vertex_handle v, OutputIterator cells, Filter f = Filter()) const
  {
    return visit_incident_cells_threadsafe<Cell_extractor<OutputIterator, Filter>,
      OutputIterator>(v, cells, f);
  }

  template <class OutputIterator>
  OutputIterator
  incident_cells_threadsafe(Vertex_handle v, OutputIterator cells) const
  {
    return incident_cells_threadsafe<False_filter>(v, cells);
  }

  template <class Filter, class OutputIterator>
  OutputIterator
  incident_facets(Vertex_handle v, OutputIterator facets, Filter f = Filter()) const
  {
    CGAL_triangulation_precondition( dimension() > 1 );
    if(dimension() == 3)
        return visit_incident_cells<Facet_extractor<OutputIterator, Filter>, OutputIterator>(v, facets, f);
    else
        return visit_incident_cells<DegCell_as_Facet_extractor<OutputIterator, Filter>, OutputIterator>(v, facets, f);
  }

  template <class OutputIterator>
  OutputIterator
  incident_facets(Vertex_handle v, OutputIterator facets) const
  {
    return incident_facets<False_filter>(v, facets);
  }
  
  template <class Filter, class OutputIterator>
  OutputIterator
  incident_facets_threadsafe(Vertex_handle v, OutputIterator facets, Filter f = Filter()) const
  {
    CGAL_triangulation_precondition( dimension() > 1 );
    if(dimension() == 3)
        return visit_incident_cells_threadsafe<Facet_extractor<OutputIterator, Filter>, OutputIterator>(v, facets, f);
    else
        return visit_incident_cells_threadsafe<DegCell_as_Facet_extractor<OutputIterator, Filter>, OutputIterator>(v, facets, f);
  }

  template <class OutputIterator>
  OutputIterator
  incident_facets_threadsafe(Vertex_handle v, OutputIterator facets) const
  {
    return incident_facets_threadsafe<False_filter>(v, facets);
  }

  template <class Filter, class OutputIterator>
  OutputIterator
  incident_edges_1d(Vertex_handle v, OutputIterator edges, Filter f = Filter()) const
  {
    CGAL_assertion (dimension() == 1);
    CGAL_triangulation_assertion( number_of_vertices() >= 3);
    Cell_handle n0 = cell(v);
    const int index_v_in_n0 = index(n0, v);
    CGAL_assume(index_v_in_n0 <= 1);
    Cell_handle n1 = neighbor(n0, 1-index_v_in_n0);
    const int index_v_in_n1 = index(n1, v);
    CGAL_assume(index_v_in_n1 <= 1);
    if(!f(vertex(n0, 1-index_v_in_n0)))
      *edges++ = Edge(n0, index(n0, v), 1-index_v_in_n0);
    if(!f(vertex(n1, 1-index_v_in_n1)))
      *edges++ = Edge(n1, index(n1, v), 1-index_v_in_n1);
    return edges;
  }

  template <class Filter, class OutputIterator>
  OutputIterator
  incident_edges(Vertex_handle v, OutputIterator edges, Filter f = Filter()) const
  {
    CGAL_triangulation_precondition( v != Vertex_handle() );
    CGAL_triangulation_precondition( dimension() >= 1 );
    CGAL_triangulation_expensive_precondition( is_vertex(v) );
    CGAL_triangulation_expensive_precondition( is_valid() );

    if (dimension() == 1) {
      return incident_edges_1d(v, edges, f);
    }
    return visit_incident_cells<Vertex_extractor<Edge_feeder_treatment<OutputIterator>,
                                                 OutputIterator, Filter,
                                                 internal::Has_member_visited<Vertex>::value>,
    OutputIterator>(v, edges, f);
  }

  template <class Filter, class OutputIterator>
  OutputIterator
  incident_edges_threadsafe(Vertex_handle v, OutputIterator edges,
			    Filter f = Filter()) const
  {
    CGAL_triangulation_precondition( v != Vertex_handle() );
    CGAL_triangulation_precondition( dimension() >= 1 );
    CGAL_triangulation_expensive_precondition( is_vertex(v) );
    CGAL_triangulation_expensive_precondition( is_valid() );

    if (dimension() == 1) {
      return incident_edges_1d(v, edges, f);
    }
    return visit_incident_cells_threadsafe<
      Vertex_extractor<Edge_feeder_treatment<OutputIterator>,
                       OutputIterator, Filter,
                       internal::Has_member_visited<Vertex>::value>,
      OutputIterator>(v, edges, f);
  }

  template <class OutputIterator>
  OutputIterator
  incident_edges(Vertex_handle v, OutputIterator edges) const
  {
    return incident_edges<False_filter>(v, edges);
  }

  template <class OutputIterator>
  OutputIterator
  incident_edges_threadsafe(Vertex_handle v, OutputIterator edges) const
  {
    return incident_edges_threadsafe<False_filter>(v, edges);
  }

  template <class Filter, class OutputIterator>
  OutputIterator
  adjacent_vertices(Vertex_handle v, OutputIterator vertices, Filter f = Filter()) const
  {
    CGAL_triangulation_precondition( v != Vertex_handle() );
    CGAL_triangulation_precondition( dimension() >= -1 );
    CGAL_triangulation_expensive_precondition( is_vertex(v) );
    CGAL_triangulation_expensive_precondition( is_valid() );

    if (dimension() == -1)
    return vertices;

    if (dimension() == 0) {
      Vertex_handle v1 = vertex(neighbor(cell(v),0), 0);
      if(!f(v1)) *vertices++ = v1;
      return vertices;
    }

    if (dimension() == 1) {
      CGAL_triangulation_assertion( number_of_vertices() >= 3);
      Cell_handle n0 = cell(v);
      const int index_v_in_n0 = index(n0, v);
      CGAL_assume(index_v_in_n0 <= 1);
      Cell_handle n1 = neighbor(n0, 1-index_v_in_n0);
      const int index_v_in_n1 = index(n1, v);
      CGAL_assume(index_v_in_n1 <= 1);
      Vertex_handle v1 = vertex(n0, 1-index_v_in_n0);
      Vertex_handle v2 = vertex(n1, 1-index_v_in_n1);
      if(!f(v1)) *vertices++ = v1;
      if(!f(v2)) *vertices++ = v2;
      return vertices;
    }
    return visit_incident_cells<Vertex_extractor<Vertex_feeder_treatment<OutputIterator>,
                                OutputIterator, Filter,
                                internal::Has_member_visited<Vertex>::value>,
    OutputIterator>(v, vertices, f);
  }

  // old name - kept for backwards compatibility but not documented
  template <class OutputIterator>
  OutputIterator
  incident_vertices(Vertex_handle v, OutputIterator vertices) const
  {
    return adjacent_vertices<False_filter>(v, vertices);
  }

  // correct name
  template <class OutputIterator>
  OutputIterator
  adjacent_vertices(Vertex_handle v, OutputIterator vertices) const
  {
    return adjacent_vertices<False_filter>(v, vertices);
  }

  template <class Visitor, class OutputIterator, class Filter>
  OutputIterator
  visit_incident_cells(Vertex_handle v, OutputIterator output, Filter f) const
  {
    CGAL_triangulation_precondition( v != Vertex_handle() );
    CGAL_triangulation_expensive_precondition( is_vertex(v) );

    if ( dimension() < 2 )
    return output;

    Visitor visit(v, output, this, f);

    std::vector<Cell_handle> tmp_cells;
    tmp_cells.reserve(64);
    if ( dimension() == 3 )
    incident_cells_3(v, cell(v), std::make_pair(std::back_inserter(tmp_cells), visit.facet_it()));
    else
    incident_cells_2(v, cell(v), std::back_inserter(tmp_cells));

    typename std::vector<Cell_handle>::iterator cit;
    for(cit = tmp_cells.begin();
        cit != tmp_cells.end();
        ++cit)
    {
      const_cast<Tds*>(this)->tds_data(*cit).clear();
      visit(*cit);
    } 

    return visit.result();
  }
  
  template <class Visitor, class OutputIterator, class Filter>
  OutputIterator
  visit_incident_cells_threadsafe(
    Vertex_handle v, OutputIterator output, Filter f) const
  {
    CGAL_triangulation_precondition( v != Vertex_handle() );
    CGAL_triangulation_expensive_precondition( is_vertex(v) );

    if ( dimension() < 2 )
    return output;

    Visitor visit(v, output, this, f);

    std::vector<Cell_handle> tmp_cells;
    tmp_cells.reserve(64);
    if ( dimension() == 3 )
      incident_cells_3_threadsafe(
        v, cell(v), tmp_cells, visit.facet_it());
    else
      incident_cells_2(v, cell(v), std::back_inserter(tmp_cells));

    typename std::vector<Cell_handle>::iterator cit;
    for(cit = tmp_cells.begin();
        cit != tmp_cells.end();
        ++cit)
    {
      visit(*cit);
    } 

    return visit.result();
  }
  
  template <class Visitor, class OutputIterator, class Filter>
  OutputIterator
  visit_incident_cells(Vertex_handle v, OutputIterator output, 
                       std::vector<Cell_handle> &cells, Filter f) const
  {
    CGAL_triangulation_precondition( v != Vertex_handle() );
    CGAL_triangulation_expensive_precondition( is_vertex(v) );

    if ( dimension() < 2 )
    return output;

    Visitor visit(v, output, this, f);

    if ( dimension() == 3 )
    incident_cells_3(v, cell(v), std::make_pair(std::back_inserter(cells), visit.facet_it()));
    else
    incident_cells_2(v, cell(v), std::back_inserter(cells));

    typename std::vector<Cell_handle>::iterator cit;
    for(cit = cells.begin();
        cit != cells.end();
        ++cit)
    {
      const_cast<Tds*>(this)->tds_data(*cit).clear();
      visit(*cit);
    }
    return visit.result();
  }

  template <class Visitor, class OutputIterator, class Filter>
  OutputIterator
  visit_just_incident_cells(Vertex_handle v, OutputIterator output, Filter f) const
  {
    CGAL_triangulation_precondition( v != Vertex_handle() );
    CGAL_triangulation_expensive_precondition( is_vertex(v) );

    if ( dimension() < 2 )
    return output;

    Visitor visit(v, output, this, f);

    std::vector<Cell_handle> tmp_cells;
    tmp_cells.reserve(64);

    if ( dimension() == 3 )
      just_incident_cells_3(v, tmp_cells);
    else
    incident_cells_2(v, cell(v), std::back_inserter(tmp_cells));

    typename std::vector<Cell_handle>::iterator cit;
    for(cit = tmp_cells.begin();
        cit != tmp_cells.end();
        ++cit)
    {
      const_cast<Tds*>(this)->tds_data(*cit).clear();
      visit(*cit);
    }
    return visit.result();
  }
  
  // For dimension 3 only
  template <class VertexFilter, class OutputVertexIterator>
  OutputVertexIterator
  adjacent_vertices_and_cells_3(Vertex_handle v, OutputVertexIterator vertices,
                                std::vector<Cell_handle> &cells,
                                VertexFilter f = VertexFilter()) const
  {
    CGAL_triangulation_precondition( v != Vertex_handle() );
    CGAL_triangulation_precondition( dimension() == 3 );
    CGAL_triangulation_expensive_precondition( is_vertex(v) );
    CGAL_triangulation_expensive_precondition( is_valid() );

    return 
      visit_incident_cells
      <
        Vertex_extractor<Vertex_feeder_treatment<OutputVertexIterator>,
                         OutputVertexIterator, 
                         VertexFilter, 
                         internal::Has_member_visited<Vertex>::value>,
        OutputVertexIterator
      >(v, vertices, cells, f);
  }

  // For dimension 3 only
  template <class OutputVertexIterator>
  OutputVertexIterator
  adjacent_vertices_and_cells_3(Vertex_handle v, OutputVertexIterator vertices,
                                std::vector<Cell_handle> &cells) const
  {
    return adjacent_vertices_and_cells_3<False_filter>(v, vertices, cells);
  }

  size_type degree(Vertex_handle v) const;

  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const;
  bool is_valid(Vertex_handle v, bool verbose = false, int level = 0) const;
  bool is_valid(Cell_handle c, bool verbose = false, int level = 0) const;

  // Helping functions
  template <class TDS_src>
  Vertex_handle copy_tds(const TDS_src & tds,
                         typename TDS_src::Vertex_handle vert);

  template <class TDS_src>
  Vertex_handle copy_tds(const TDS_src & tds)
  {
    return copy_tds(tds, typename TDS_src::Vertex_handle());
  }
    // returns the new vertex corresponding to vert in the new tds

  template <class TDS_src,class ConvertVertex,class ConvertCell>
  Vertex_handle copy_tds(const TDS_src&, typename TDS_src::Vertex_handle,const ConvertVertex&,const ConvertCell&);

  
  void swap(Tds & tds);

  void clear();

  void set_adjacency(Cell_handle c0, int i0,
                     Cell_handle c1, int i1)
  {
      CGAL_triangulation_assertion(i0 >= 0 && i0 <= dimension());
      CGAL_triangulation_assertion(i1 >= 0 && i1 <= dimension());
      CGAL_triangulation_assertion(c0 != c1);
      set_neighbor(c0,i0,c1);
      set_neighbor(c1,i1,c0);
  }

  int mirror_index(Cell_handle c, int i) const
  {
      CGAL_triangulation_precondition ( i>=0 && i<4 );
      return index(neighbor(c,i), c);
  }

  Vertex_handle mirror_vertex(Cell_handle c, int i) const
  {
    return vertex(neighbor(c,i), mirror_index(c, i));
  }

  Facet mirror_facet(Facet f) const
  {
    Cell_handle neighbor_cell = neighbor(f.first, f.second);
    const int opposite_index = index(neighbor_cell, f.first);
    return Facet(neighbor_cell, opposite_index);
  }

  // We need the const_cast<>s because TDS is not const-correct.
  Cell_range & cells() { return _cells; }
  Cell_range & cells() const
  { return const_cast<Tds*>(this)->_cells; }

  Vertex_range & vertices() {return _vertices;}
  Vertex_range & vertices() const
  { return const_cast<Tds*>(this)->_vertices; }

private:

  // Change the orientation of the cell by swapping indices 0 and 1.
  void change_orientation(Cell_handle c)
  {
    Vertex_handle tmp_v = vertex(c, 0);
    set_vertex(c, 0, vertex(c, 1));
    set_vertex(c, 1, tmp_v);
    Cell_handle tmp_c = neighbor(c,0);
    set_neighbor(c, 0, neighbor(c,1));
    set_neighbor(c, 1, tmp_c);
  }

  // in dimension i, number of vertices >= i+2
  // ( the boundary of a simplex in dimension i+1 has i+2 vertices )
  int _dimension;

  Cell_range   _cells;
  Vertex_range _vertices;

  // used by is-valid :
  bool count_vertices(size_type &i, bool verbose = false, int level = 0) const;
  // counts AND checks the validity
  bool count_facets(size_type &i, bool verbose = false, int level = 0) const;
  // counts but does not check
  bool count_edges(size_type &i, bool verbose = false, int level = 0) const;
  // counts but does not check
  bool count_cells(size_type &i, bool verbose = false, int level = 0) const;
  // counts AND checks the validity
};

#ifdef CGAL_TDS_USE_RECURSIVE_CREATE_STAR_3
template <class Vb, class Cb, class Ct>
typename Triangulation_data_structure_3<Vb,Cb,Ct>::Cell_handle
Triangulation_data_structure_3<Vb,Cb,Ct>::
create_star_3(Vertex_handle v, Cell_handle c, int li, int prev_ind2)
{
    CGAL_triangulation_precondition( dimension() == 3);
    CGAL_triangulation_precondition( tds_data(c).is_in_conflict() );
    CGAL_triangulation_precondition( ! tds_data(neighbor(c, li)).is_in_conflict() );

    Cell_handle cnew = create_cell(vertex(c, 0),
                                   vertex(c, 1),
                                   vertex(c, 2),
                                   vertex(c, 3));
    set_vertex(cnew, li, v);
    Cell_handle c_li = neighbor(c, li);
    set_adjacency(cnew, li, c_li, index(c_li, c));

    // Look for the other neighbors of cnew.
    for (int ii=0; ii<4; ++ii) {
      if (ii == prev_ind2 || neighbor(cnew, ii) != Cell_handle())
          continue;
      set_cell(vertex(cnew, ii), cnew);

      // Indices of the vertices of cnew such that ii,vj1,vj2,li positive.
      Vertex_handle vj1 = vertex(c, next_around_edge(ii, li));
      Vertex_handle vj2 = vertex(c, next_around_edge(li, ii));
      Cell_handle cur = c;
      int zz = ii;
      Cell_handle n = neighbor(cur, zz);
      // turn around the oriented edge vj1 vj2
      while ( tds_data(n).is_in_conflict() ) {
        CGAL_triangulation_assertion( n != c );
        cur = n;
        zz = next_around_edge(index(n, vj1), index(n, vj2));
        n = neighbor(cur, zz);
      }
      // Now n is outside region, cur is inside.
      tds_data(n).clear(); // Reset the flag for boundary cells.

      int jj1 = index(n, vj1);
      int jj2 = index(n, vj2);
      Vertex_handle vvv = vertex(n, next_around_edge(jj1, jj2));
      Cell_handle nnn = neighbor(n, next_around_edge(jj2, jj1));
      int zzz = index(nnn, vvv);
      if (nnn == cur) {
        // Neighbor relation is reciprocal, ie
        // the cell we are looking for is not yet created.
        nnn = create_star_3(v, nnn, zz, zzz);
      }

      set_adjacency(nnn, zzz, cnew, ii);
    }

    return cnew;
}
#else
template <class Vb, class Cb, class Ct>
typename Triangulation_data_structure_3<Vb,Cb,Ct>::Cell_handle
Triangulation_data_structure_3<Vb,Cb,Ct>::
recursive_create_star_3(Vertex_handle v, Cell_handle c, int li,
                        int prev_ind2, int depth)
{
    if ( depth == 100 ) return non_recursive_create_star_3(v,c,li,prev_ind2);
    CGAL_triangulation_precondition( dimension() == 3);
    CGAL_triangulation_precondition( tds_data(c).is_in_conflict() );
    CGAL_triangulation_precondition( ! tds_data(neighbor(c,li)).is_in_conflict() );

    Cell_handle cnew = create_cell(vertex(c, 0),
                                   vertex(c, 1),
                                   vertex(c, 2),
                                   vertex(c, 3));
    set_vertex(cnew, li, v);
    Cell_handle c_li = neighbor(c,li);
    set_adjacency(cnew, li, c_li, index(c_li, c));

    // Look for the other neighbors of cnew.
    for (int ii=0; ii<4; ++ii) {
      if (ii == prev_ind2 || neighbor(cnew, ii) != Cell_handle())
          continue;
      set_cell(vertex(cnew, ii), cnew);

      // Indices of the vertices of cnew such that ii,vj1,vj2,li positive.
      Vertex_handle vj1 = vertex(c, next_around_edge(ii, li));
      Vertex_handle vj2 = vertex(c, next_around_edge(li, ii));
      Cell_handle cur = c;
      int zz = ii;
      Cell_handle n = neighbor(cur, zz);
      // turn around the oriented edge vj1 vj2
      while ( tds_data(n).is_in_conflict() ) {
        CGAL_triangulation_assertion( n != c );
        cur = n;
        zz = next_around_edge(index(n, vj1), index(n, vj2));
        n = neighbor(cur, zz);
      }
      // Now n is outside region, cur is inside.
      tds_data(n).clear(); // Reset the flag for boundary cells.

      int jj1 = index(n, vj1);
      int jj2 = index(n, vj2);
      Vertex_handle vvv = vertex(n, next_around_edge(jj1, jj2));
      Cell_handle nnn = neighbor(n, next_around_edge(jj2, jj1));
      int zzz = index(nnn, vvv);
      if (nnn == cur) {
        // Neighbor relation is reciprocal, ie
        // the cell we are looking for is not yet created.
        nnn = recursive_create_star_3(v, nnn, zz, zzz,depth+1);
      }

      set_adjacency(nnn, zzz, cnew, ii);
    }

    return cnew;
}

//We use the class iAdjacency_info instead of a tuple because
//at the moment we made the change it was faster like this.
template <class Vb, class Cb, class Ct>
typename Triangulation_data_structure_3<Vb,Cb,Ct>::Cell_handle
Triangulation_data_structure_3<Vb,Cb,Ct>::
non_recursive_create_star_3(Vertex_handle v, Cell_handle c, int li, int prev_ind2)
{
    CGAL_triangulation_precondition( dimension() == 3);
    CGAL_triangulation_precondition( tds_data(c).is_in_conflict() );
    CGAL_triangulation_precondition( ! tds_data(neighbor(c, li)).is_in_conflict() );

    Cell_handle cnew = create_cell(vertex(c, 0),
                                   vertex(c, 1),
                                   vertex(c, 2),
                                   vertex(c, 3));
    set_vertex(cnew, li, v);
    Cell_handle c_li = neighbor(c, li);
    set_adjacency(cnew, li, c_li, index(c_li, c));
    
    std::stack<iAdjacency_info> adjacency_info_stack;
  
    int ii=0;
    do
    {
      // Look for the other neighbors of cnew.
      if ( ! (ii == prev_ind2 || neighbor(cnew, ii) != Cell_handle()) ){
        set_cell(vertex(cnew, ii), cnew);

        // Indices of the vertices of cnew such that ii,vj1,vj2,li positive.
        Vertex_handle vj1 = vertex(c, next_around_edge(ii, li));
        Vertex_handle vj2 = vertex(c, next_around_edge(li, ii));
        Cell_handle cur = c;
        int zz = ii;
        Cell_handle n = neighbor(cur, zz);
        // turn around the oriented edge vj1 vj2
        while ( tds_data(n).is_in_conflict() ) {
          CGAL_triangulation_assertion( n != c );
          cur = n;
          zz = next_around_edge(index(n, vj1), index(n, vj2));
          n = neighbor(cur, zz);
        }
        // Now n is outside region, cur is inside.
        tds_data(n).clear(); // Reset the flag for boundary cells.

        int jj1 = index(n, vj1);
        int jj2 = index(n, vj2);
        Vertex_handle vvv = vertex(n, next_around_edge(jj1, jj2));
        Cell_handle nnn = neighbor(n, next_around_edge(jj2, jj1));
        int zzz = index(nnn, vvv);
        if (nnn == cur) {
          // Neighbor relation is reciprocal, ie
          // the cell we are looking for is not yet created.
          //re-run the loop
          adjacency_info_stack.push( iAdjacency_info(zzz,cnew,ii,c,li,prev_ind2) );
          c=nnn;  li=zz; prev_ind2=zzz; ii=0;
          //copy-pasted from beginning to avoid if ii==0
          CGAL_triangulation_precondition( tds_data(c).is_in_conflict() );
          CGAL_triangulation_precondition( ! tds_data(neighbor(c, li)).is_in_conflict() );
          cnew = create_cell(vertex(c, 0),vertex(c, 1),vertex(c, 2),vertex(c, 3));
          set_vertex(cnew, li, v);
          c_li = neighbor(c, li);
          set_adjacency(cnew, li, c_li, index(c_li, c));          
          continue;
        }
        set_adjacency(nnn, zzz, cnew, ii);
      }
      while (++ii==4){
        if ( adjacency_info_stack.empty() ) return cnew;
        Cell_handle nnn=cnew; 
        int zzz;
        adjacency_info_stack.top().update_variables(zzz,cnew,ii,c,li,prev_ind2);
        adjacency_info_stack.pop();
        set_adjacency(nnn, zzz, cnew, ii);
      }
    }
    while (true);
}
#endif

template <class Vb, class Cb, class Ct>
typename Triangulation_data_structure_3<Vb,Cb,Ct>::Cell_handle
Triangulation_data_structure_3<Vb,Cb,Ct>::
create_star_2(Vertex_handle v, Cell_handle c, int li )
{
  CGAL_triangulation_assertion( dimension() == 2 );
  Cell_handle cnew;

  // i1 i2 such that v,i1,i2 positive
  int i1=ccw(li);
  // traversal of the boundary of region in ccw order to create all
  // the new facets
  Cell_handle bound = c;
  Vertex_handle v1 = vertex(c, i1);
  int ind = index(neighbor(c, li), c); // to be able to find the
                                       // first cell that will be created
  Cell_handle cur;
  Cell_handle pnew = Cell_handle();
  do {
    cur = bound;
    // turn around v2 until we reach the boundary of region
    while ( tds_data(neighbor(cur, cw(i1))).is_in_conflict() ) {
      // neighbor in conflict
      cur = neighbor(cur, cw(i1));
      i1 = index(cur,  v1 );
    }
    tds_data(neighbor(cur, cw(i1))).clear();
    // here cur has an edge on the boundary of region
    cnew = create_face( v, v1, vertex(cur,  ccw(i1) ) );
    set_adjacency(cnew, 0, neighbor(cur, cw(i1)),
                  index(neighbor(cur, cw(i1)), cur));
    set_neighbor(cnew, 1, Cell_handle());
    set_neighbor(cnew, 2, pnew);
    // pnew is null at the first iteration
    set_cell(v1, cnew);
    //pnew->set_neighbor( cw(pnew->index(v1)), cnew );
    if (pnew != Cell_handle()) { set_neighbor(pnew, 1, cnew );}

    bound = cur;
    i1 = ccw(i1);
    v1 = vertex(bound, i1);
    pnew = cnew;
    //} while ( ( bound != c ) || ( li != cw(i1) ) );
  } while ( v1 != vertex(c, ccw(li)) );
  // missing neighbors between the first and the last created cells
  cur = neighbor(neighbor(c, li), ind); // first created cell
  set_adjacency(cnew, 1, cur, 2);
  return cnew;
}

template <class Vb, class Cb, class Ct>
std::istream&
operator>>(std::istream& is, Triangulation_data_structure_3<Vb,Cb,Ct>& tds)
  // reads :
  // the dimension
  // the number of vertices
  // the number of cells
  // the cells by the indices of their vertices
  // the neighbors of each cell by their index in the preceding list of cells
  // when dimension < 3 : the same with faces of maximal dimension
{
  typedef Triangulation_data_structure_3<Vb,Cb,Ct> Tds;
  typedef typename Tds::Vertex_handle  Vertex_handle;
  typedef typename Tds::Cell_handle    Cell_handle;

  tds.clear();

  std::size_t n;
  int d;
  if(is_ascii(is))
     is >> d >> n;
  else {
    read(is, n);
    read(is,d);
  }
  tds.set_dimension(d);

  if(n == 0)
    return is;

  std::vector<Vertex_handle > V(n);

  // creation of the vertices
  for (std::size_t i=0; i < n; i++) {
    //    is >> p;
    //    V[i] = tds.create_vertex();
    //    V[i]->set_point(p);
    V[i] = tds.create_vertex();
  }

  std::vector< Cell_handle > C;
  std::size_t m;

  tds.read_cells(is, V, m, C);
  CGAL_triangulation_assertion( tds.is_valid() );
  return is;
}

template <class Vb, class Cb, class Ct>
std::ostream&
operator<<(std::ostream& os, const Triangulation_data_structure_3<Vb,Cb,Ct> &tds)
  // writes :
  // the dimension
  // the number of vertices
  // the number of cells
  // the cells by the indices of their vertices
  // the neighbors of each cell by their index in the preceding list of cells
  // when dimension < 3 : the same with faces of maximal dimension
{
  typedef Triangulation_data_structure_3<Vb,Cb,Ct> Tds;
  typedef typename Tds::size_type               size_type;
  typedef typename Tds::Vertex_handle           Vertex_handle;
  typedef typename Tds::Vertex_iterator         Vertex_iterator;


  Unique_hash_map<Vertex_handle, size_type> V;

  // outputs dimension and number of vertices
  size_type n = tds.number_of_vertices();

  if (is_ascii(os))
      os << tds.dimension() << std::endl << n << std::endl;
  else
  {
    write(os,tds.dimension());
    write(os,n);
  }

  if (n == 0)
    return os;

  // index the vertices
  size_type i = 0;
  for (Vertex_iterator it=tds.vertices_begin(); it != tds.vertices_end(); ++it)
    V[it] = i++;

  CGAL_triangulation_assertion( i == n );

  tds.print_cells(os, V);

  return os;
}

template <class Vb, class Cb, class Ct>
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
is_simplex( Cell_handle c ) const
{
  switch(dimension()) {
    case 3 : return is_cell(c);
    case 2 : return is_facet(c, 3);
    case 1 : return is_edge(c, 0, 1);
  case 0 : return is_vertex(vertex(c, 0));
    case -1 : return c == cells().begin();
  }
  return false;
}

template <class Vb, class Cb, class Ct>
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
is_vertex(Vertex_handle v) const
{
    return vertices().owns_dereferencable(v);
}

template <class Vb, class Cb, class Ct>
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
is_edge(Vertex_handle u, Vertex_handle v,
        Cell_handle &c, int &i, int &j) const
  // returns false when dimension <1 or when indices wrong
{
    CGAL_triangulation_expensive_precondition( is_vertex(u) && is_vertex(v) );

    if (u==v)
        return false;

    std::vector<Cell_handle> cells;
    cells.reserve(64);
    incident_cells(u, std::back_inserter(cells));

    for (typename std::vector<Cell_handle>::iterator cit = cells.begin();
              cit != cells.end(); ++cit)
      if (has_vertex(*cit, v, j)) {
            c = *cit;
            i = index(c, u);
            return true;
        }
    return false;
}

template <class Vb, class Cb, class Ct>
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
is_edge(Vertex_handle u, Vertex_handle v) const
{
    Cell_handle c;
    int i, j;
    return is_edge(u, v, c, i, j);
}

template <class Vb, class Cb, class Ct>
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
is_edge(Cell_handle c, int i, int j) const
{
  if (dimension() < 1)
    return false;

  if ( i==j ) return false;
  if ( (i<0) || (j<0) ) return false;
  if ( (dimension() == 1) && ((i>1) || (j>1)) ) return false;
  if ( (dimension() == 2) && ((i>2) || (j>2)) ) return false;
  if ((i>3) || (j>3)) return false;

  return cells().owns_dereferencable(c);
}

template <class Vb, class Cb, class Ct>
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
is_facet(Vertex_handle u, Vertex_handle v,
         Vertex_handle w,
         Cell_handle & c, int & i, int & j, int & k) const
  // returns false when dimension <2 or when indices wrong
{
    CGAL_triangulation_expensive_precondition( is_vertex(u) &&
                                               is_vertex(v) &&
                                               is_vertex(w) );

    if ( u==v || u==w || v==w )
        return false;
    if (dimension() < 2)
        return false;

    std::vector<Cell_handle> cells;
    cells.reserve(64);
    incident_cells(u, std::back_inserter(cells));

    for (typename std::vector<Cell_handle>::iterator cit = cells.begin();
              cit != cells.end(); ++cit)
      if (has_vertex(*cit, v, j) && has_vertex(*cit, w, k)) {
            c = *cit;
            i = index(c, u);
            return true;
        }
    return false;
}

template <class Vb, class Cb, class Ct>
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
is_facet(Cell_handle c, int i) const
{
    CGAL_triangulation_precondition(i>=0 && i<4);

    if ( dimension() < 2 )
        return false;

    if ( (dimension() == 2) && (i!=3) )
        return false;

    return cells().owns_dereferencable(c);
}

template <class Vb, class Cb, class Ct>
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
is_cell( Cell_handle c ) const
  // returns false when dimension <3
{
    if (dimension() < 3)
        return false;

    return cells().owns_dereferencable(c);
}

template <class Vb, class Cb, class Ct>
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
is_cell(Vertex_handle u, Vertex_handle v,
        Vertex_handle w, Vertex_handle t,
        Cell_handle & c, int & i, int & j, int & k, int & l) const
  // returns false when dimension <3
{
    CGAL_triangulation_expensive_precondition( is_vertex(u) &&
                                               is_vertex(v) &&
                                               is_vertex(w) &&
                                               is_vertex(t) );

    if ( u==v || u==w || u==t || v==w || v==t || w==t )
        return false;

    std::vector<Cell_handle> cells;
    cells.reserve(64);
    incident_cells(u, std::back_inserter(cells));

    for (typename std::vector<Cell_handle>::iterator cit = cells.begin();
              cit != cells.end(); ++cit)
      if (has_vertex(*cit, v, j) && has_vertex(*cit, w, k) &&
          has_vertex(*cit, t, l)) {
            c = *cit;
            i = index(c, u);
            return true;
        }
    return false;
}

template <class Vb, class Cb, class Ct>
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
is_cell(Vertex_handle u, Vertex_handle v,
        Vertex_handle w, Vertex_handle t)
    const
  // returns false when dimension <3
{
    Cell_handle c;
    int i, j, k, l;
    return is_cell(u, v, w, t, c, i, j, k, l);
}

template <class Vb, class Cb, class Ct>
inline
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
has_vertex(Cell_handle c, int i, Vertex_handle v, int & j) const
  // computes the index j of the vertex in the cell c giving the query
  // facet (c,i)
  // j has no meaning if false is returned
{
  CGAL_triangulation_precondition( dimension() == 3 );
  return ( has_vertex(c,v,j) && (j != i) );
}

template <class Vb, class Cb, class Ct>
inline
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
has_vertex(Cell_handle c, int i, Vertex_handle v) const
  // checks whether the query facet (c,i) has vertex v
{
  CGAL_triangulation_precondition( dimension() == 3 );
  int j;
  return ( has_vertex(c,v,j) && (j != i) );
}

template <class Vb, class Cb, class Ct>
inline
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
has_vertex(const Facet & f, Vertex_handle v, int & j) const
{
  return has_vertex(f.first, f.second, v, j);
}

template <class Vb, class Cb, class Ct>
inline
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
has_vertex(const Facet & f, Vertex_handle v) const
{
  return has_vertex(f.first, f.second, v);
}

template <class Vb, class Cb, class Ct>
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
are_equal(Cell_handle c, int i, Cell_handle n, int j) const
  // tests whether facets c,i and n,j, have the same 3 vertices
  // the triangulation is supposed to be valid, the orientation of the
  // facets is not checked here
  // the neighbor relations between c and  n are not tested either,
  // which allows to use this method before setting these relations
  // (see remove in Delaunay_3)
  //   if ( neighbor(c, i) != n ) return false;
  //   if ( neighbor(n, j) != c ) return false;
{
  CGAL_triangulation_precondition( dimension() == 3 );

  if ( (c==n) && (i==j) ) return true;

  int j1,j2,j3;
  return( has_vertex(n, vertex(c, (i+1)&3), j1 ) &&
          has_vertex(n, vertex(c, (i+2)&3), j2 ) &&
          has_vertex(n, vertex(c, (i+3)&3), j3 ) &&
          ( j1+j2+j3+j == 6 ) );
}

template <class Vb, class Cb, class Ct>
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
are_equal(const Facet & f, const Facet & g) const
{
  return are_equal(f.first, f.second, g.first, g.second);
}

template <class Vb, class Cb, class Ct>
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
are_equal(const Facet & f, Cell_handle n, int j) const
{
  return are_equal(f.first, f.second, n, j);
}

template <class Vb, class Cb, class Ct>
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
flip( Cell_handle c, int i )
  // returns false if the facet is not flippable
  // true other wise and
  // flips facet i of cell c
  // c will be replaced by one of the new cells
{
  CGAL_triangulation_precondition( (dimension() == 3) && (0<=i) && (i<4)
                                   && (number_of_vertices() >= 6) );
  CGAL_triangulation_expensive_precondition( is_cell(c) );

  Cell_handle n = neighbor(c, i);
  int in = index(n, c);

  // checks that the facet is flippable,
  // ie the future edge does not already exist
  if (is_edge(vertex(c, i), vertex(n, in)))
      return false;

  flip_really(c,i,n,in);
  return true;
}

template <class Vb, class Cb, class Ct>
void
Triangulation_data_structure_3<Vb,Cb,Ct>::
flip_flippable(Cell_handle c, int i )
  // flips facet i of cell c
  // c will be replaced by one of the new cells
{
  CGAL_triangulation_precondition( (dimension() == 3) && (0<=i) && (i<4)
                                   && (number_of_vertices() >= 6) );
  CGAL_triangulation_expensive_precondition( is_cell(c) );

  Cell_handle n = neighbor(c, i);
  int in = index(n, c);

  // checks that the facet is flippable,
  // ie the future edge does not already exist
  CGAL_triangulation_expensive_precondition( !is_edge(vertex(c, i),
                                                      vertex(n, in)));
  flip_really(c,i,n,in);
}

template <class Vb, class Cb, class Ct>
inline
void
Triangulation_data_structure_3<Vb,Cb,Ct>::
flip_really( Cell_handle c, int i, Cell_handle n, int in )
  // private - used by flip and flip_flippable
{
  int i1 = (i+1)&3;
  int i2 = (i+2)&3;
  int i3 = (i+3)&3;

  int in1 = index(n, vertex(c, i1));
  int in2 = index(n, vertex(c, i2));
  int in3 = index(n, vertex(c, i3));

  set_adjacency(c, i, neighbor(n, in3), index(neighbor(n, in3), n));
  set_vertex(c, i3, vertex(n, in) );

  set_adjacency(n, in, neighbor(c, i1), index(neighbor(c, i1), c));
  set_vertex(n, in1, vertex(c, i) );

  Cell_handle cnew = create_cell(vertex(c, i), vertex(c, i1),
                                 vertex(n, in), vertex(n, in3));

  set_adjacency(cnew, 0, neighbor(n, in2), index(neighbor(n, in2), n));
  set_adjacency(cnew, 1, n, in2);
  set_adjacency(cnew, 2, neighbor(c, i2), index(neighbor(c, i2), c));
  set_adjacency(cnew, 3, c, i2);
  set_adjacency(c, i1, n, in3);

  if ((i&1) != 0)
      change_orientation(cnew);

  set_cell(vertex(c, i1), cnew);
  set_cell(vertex(c, i2), c);
  set_cell(vertex(n, in3), n);
  // to be implemented : 2d case
  // CGAL_triangulation_precondition( (0<=i) && (i<3) );
}

template <class Vb, class Cb, class Ct>
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
flip( Cell_handle c, int i, int j )
  // returns false if the edge is not flippable
  // true otherwise and
  // flips edge i,j of cell c
  // c will be deleted
{
  CGAL_triangulation_precondition( (dimension() == 3)
                                   && (0<=i) && (i<4)
                                   && (0<=j) && (j<4)
                                   && ( i != j )
                                   && (number_of_vertices() >= 6) );
  CGAL_triangulation_expensive_precondition( is_cell(c) );

  // checks that the edge is flippable ie degree 3
  int degree = 0;
  Cell_circulator ccir = incident_cells(c,i,j);
  Cell_circulator cdone = ccir;
  do {
    ++degree;
    ++ccir;
  } while ( ccir != cdone );

  if ( degree != 3 ) return false;

  int next = next_around_edge(i,j);
  Cell_handle c1 = neighbor(c,  next );
  Vertex_handle v1 = vertex(c,  next ); // will become vertex of c1
  int i1 = index(c1,  vertex(c, i) );
  int j1 = index(c1,  vertex(c, j) );

  int next1 = next_around_edge(i1,j1);
  Cell_handle c2  = neighbor(c1,  next1 );
  Vertex_handle v2 = vertex(c1,  next1 ); // will become vertex of c2
  int i2 = index(c2,  vertex(c, i) );
  int j2 = index(c2,  vertex(c, j) );

  int next2 = next_around_edge(i2,j2);
  Vertex_handle v3 = vertex(c2,  next2 );

  // checks that the edge is flippable,
  // is the future cells do not already exist
  if ( is_cell(v1,v2,v3,vertex(c, i)) ) return false;
  if ( is_cell(v1,v2,v3,vertex(c, j)) ) return false;

  flip_really(c,i,j,c1,v1,i1,j1,next1,c2,v2,i2,j2,next2,v3);

  return true;
}

template <class Vb, class Cb, class Ct>
void
Triangulation_data_structure_3<Vb,Cb,Ct>::
flip_flippable( Cell_handle c, int i, int j )
  // flips edge i,j of cell c
  // c will be deleted
{
  CGAL_triangulation_precondition( (dimension() == 3)
                                   && (0<=i) && (i<4)
                                   && (0<=j) && (j<4)
                                   && ( i != j )
                                   && (number_of_vertices() >= 6) );
  CGAL_triangulation_expensive_precondition( is_cell(c) );

  // checks that the edge is flippable ie degree 3
  CGAL_triangulation_precondition_code( int degree = 0; );
  CGAL_triangulation_precondition_code
    ( Cell_circulator ccir = incident_cells(c,i,j); );
  CGAL_triangulation_precondition_code( Cell_circulator cdone = ccir; );
  CGAL_triangulation_precondition_code( do {
                                          ++degree;
                                          ++ccir;
                                        } while ( ccir != cdone ); );

  CGAL_triangulation_precondition( degree == 3 );

  int next = next_around_edge(i,j);
  Cell_handle c1 = neighbor(c,  next );
  Vertex_handle v1 = vertex(c,  next ); // will become vertex of c1
  int i1 = index(c1,  vertex(c, i) );
  int j1 = index(c1,  vertex(c, j) );

  int next1 = next_around_edge(i1,j1);
  Cell_handle c2  = neighbor(c1,  next1 );
  Vertex_handle v2 = vertex(c1,  next1 ); // will become vertex of c2
  int i2 = index(c2,  vertex(c, i) );
  int j2 = index(c2,  vertex(c, j) );

  int next2 = next_around_edge(i2,j2);
  Vertex_handle v3 = vertex(c2,  next2 );

  // checks that the edge is flippable,
  // is the future cells do not already exist
  CGAL_triangulation_expensive_precondition( !is_cell(v1,v2,v3,vertex(c, i)) );
  CGAL_triangulation_expensive_precondition( !is_cell(v1,v2,v3,vertex(c, j)) );

  flip_really(c,i,j,c1,v1,i1,j1,next1,c2,v2,i2,j2,next2,v3);
}

template <class Vb, class Cb, class Ct>
inline
void
Triangulation_data_structure_3<Vb,Cb,Ct>::
flip_really( Cell_handle c, int i, int j,
             Cell_handle c1, Vertex_handle v1,
             int i1, int j1, int next1,
             Cell_handle c2, Vertex_handle v2,
             int i2, int j2, int next2,
             Vertex_handle v3 )
{
  set_cell(vertex(c, i), c1);
  set_cell(vertex(c, j), c2);

  set_vertex(c1, j1, v1 );
  set_cell(v1, c1);
  set_vertex(c2, i2, v2 );
  set_cell(v2, c2);

  set_adjacency(c1, next1, neighbor(c2, j2), index(neighbor(c2, j2), c2));
  set_adjacency(c2,index(c2, v1), neighbor(c1, i1),index(neighbor(c1, i1), c1));

  set_adjacency(c1, i1, c2, j2);

  set_adjacency(c1, 6-i1-j1-next1, neighbor(c, j), index(neighbor(c, j), c));
  set_adjacency(c2, next2, neighbor(c, i), index(neighbor(c, i), c));

  set_cell(v3, c2 );

  delete_cell( c );
}

template <class Vb, class Cb, class Ct>
void
Triangulation_data_structure_3<Vb,Cb,Ct>::
read_cells(std::istream& is, const std::vector< Vertex_handle > &V,
           std::size_t & m, std::vector< Cell_handle > &C)
{
  // creation of the cells and neighbors
  switch (dimension()) {
  case 3:
  case 2:
  case 1:
    {
      if(is_ascii(is))
        is >> m;
      else
        read(is, m);

      C.resize(m);

      for(std::size_t i = 0; i < m; i++) {
        Cell_handle c = create_cell();
        for (int k=0; k<=dimension(); ++k) {
          std::size_t ik;
            if(is_ascii(is))
               is >> ik;
            else
              read(is, ik);
            set_vertex(c, k, V[ik]);
            set_cell(V[ik], c);
        }
        C[i] = c;
      }
      for(std::size_t j = 0; j < m; j++) {
        Cell_handle c = C[j];
        for (int k=0; k<=dimension(); ++k) {
          std::size_t ik;
            if(is_ascii(is))
              is >> ik;
            else
              read(is, ik);
            set_neighbor(c, k, C[ik]);
        }
      }
      break;
    }
  case 0:
    {
      m = 2;
      C.resize(m);
      //      CGAL_triangulation_assertion( n == 2 );
      for (int i=0; i < 2; i++) {
        Cell_handle c = create_face(V[i], Vertex_handle(), Vertex_handle());
        C[i] = c;
        set_cell(V[i], c);
      }
      for (int j=0; j < 2; j++) {
        Cell_handle c = C[j];
        set_neighbor(c, 0, C[1-j]);
      }
      break;
    }
  case -1:
    {
      m = 1;
      C.resize(m);
      //      CGAL_triangulation_assertion( n == 1 );
      Cell_handle c = create_face(V[0], Vertex_handle(), Vertex_handle());
      C[0] = c;
      set_cell(V[0], c);
      break;
    }
  }
}

template <class Vb, class Cb, class Ct>
void
Triangulation_data_structure_3<Vb,Cb,Ct>::
print_cells(std::ostream& os, const Unique_hash_map<Vertex_handle, std::size_t> &V ) const
{
  Unique_hash_map<Cell_handle, std::size_t > C;
  std::size_t i = 0;

  switch ( dimension() ) {
  case 3:
    {
      std::size_t m = number_of_cells();
      if(is_ascii(os))
        os << m << std::endl;
      else
        write(os, m);

      // write the cells
      Cell_iterator it;
      for(it = cells_begin(); it != cells_end(); ++it) {
        C[it] = i++;
        for(int j = 0; j < 4; j++){
          if(is_ascii(os)) {
            os << V[vertex(it, j)];
            if ( j==3 )
              os << '\n';
            else
              os << ' ';
          }
          else
            write(os, V[vertex(it, j)]);
        }
      }
      CGAL_triangulation_assertion( i == m );

      // write the neighbors
      for(it = cells_begin(); it != cells_end(); ++it) {
        for (int j = 0; j < 4; j++) {
          if(is_ascii(os)){
            os << C[neighbor(it, j)];
            if(j==3)
              os << '\n';
            else
              os <<  ' ';
          }
          else
            write(os, C[neighbor(it, j)]);
        }
      }
      break;
    }
  case 2:
    {
      size_type m = number_of_facets();
      if(is_ascii(os))
        os << m << '\n';
      else
        write(os, m);

      // write the facets
      Facet_iterator it;
      for(it = facets_begin(); it != facets_end(); ++it) {
        C[(*it).first] = i++;
        for(int j = 0; j < 3; j++){
          if(is_ascii(os)) {
            os << V[vertex((*it).first, j)];
            if ( j==2 )
              os << '\n';
            else
              os <<  ' ';
          }
          else {
            write(os,  V[vertex((*it).first, j)]);
          }
        }
      }
      CGAL_triangulation_assertion( i == m );

      // write the neighbors
      for(it = facets_begin(); it != facets_end(); ++it) {
        for (int j = 0; j < 3; j++) {
          if(is_ascii(os)){
            os << C[neighbor((*it).first, j)];
            if(j==2)
              os << '\n';
            else
              os <<  ' ';
          }
          else {
            write(os, C[neighbor((*it).first,j)]);
          }
        }
      }
      break;
    }
  case 1:
    {
      size_type m = number_of_edges();
      if(is_ascii(os))
        os << m << '\n';
      else
        write(os, m);
      // write the edges
      Edge_iterator it;
      for(it = edges_begin(); it != edges_end(); ++it) {
        C[(*it).first] = i++;
        for(int j = 0; j < 2; j++){
          if(is_ascii(os)) {
            os << V[vertex((*it).first, j)];
            if ( j==1 )
              os << '\n';
            else
              os <<  ' ';
          }
          else {
            write(os, V[vertex((*it).first, j)]);
          }
        }
      }
      CGAL_triangulation_assertion( i == m );

      // write the neighbors
      for(it = edges_begin(); it != edges_end(); ++it) {
        for (int j = 0; j < 2; j++) {
          if(is_ascii(os)){
            os << C[neighbor((*it).first, j)];
            if(j==1)
              os << '\n';
            else
              os <<  ' ';
          }
          else {
            write(os, C[neighbor((*it).first, j)]);
          }
        }
      }
      break;
    }
  }
}


template <class Vb, class Cb, class Ct>
typename Triangulation_data_structure_3<Vb,Cb,Ct>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb,Ct>::insert_first_finite_cell(
  Vertex_handle &v0, Vertex_handle &v1, Vertex_handle &v2, Vertex_handle &v3,
  Vertex_handle v_infinite)
{
  CGAL_triangulation_precondition( 
    (v_infinite == Vertex_handle() && dimension() == -2)
    || (v_infinite != Vertex_handle() && dimension() == -1));

  if (v_infinite == Vertex_handle())
    v_infinite = create_vertex();

  set_dimension(3);

  v0 = create_vertex();
  v1 = create_vertex();
  v2 = create_vertex();
  v3 = create_vertex();

  Cell_handle c0123 = create_cell(v0,         v1,   v2,   v3);
  Cell_handle ci012 = create_cell(v_infinite, v0,   v1,   v2);
  Cell_handle ci103 = create_cell(v_infinite, v1,   v0,   v3);
  Cell_handle ci023 = create_cell(v_infinite, v0,   v2,   v3);
  Cell_handle ci132 = create_cell(v_infinite, v1,   v3,   v2);

  set_cell(v_infinite, ci012);
  set_cell(v0, c0123);
  set_cell(v1, c0123);
  set_cell(v2, c0123);
  set_cell(v3, c0123);

  set_adjacency(c0123, 0, ci132, 0);
  set_adjacency(c0123, 1, ci023, 0);
  set_adjacency(c0123, 2, ci103, 0);
  set_adjacency(c0123, 3, ci012, 0);

  set_adjacency(ci012, 3, ci103, 3);
  set_adjacency(ci012, 2, ci023, 3);
  set_adjacency(ci012, 1, ci132, 2);
  set_adjacency(ci103, 1, ci023, 2);
  set_adjacency(ci023, 1, ci132, 1);
  set_adjacency(ci103, 2, ci132, 3);

  return v_infinite;
}

template <class Vb, class Cb, class Ct>
typename Triangulation_data_structure_3<Vb,Cb,Ct>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb,Ct>::
insert_in_cell(Cell_handle c)
{
  CGAL_triangulation_precondition( dimension() == 3 );
  CGAL_triangulation_precondition( c != Cell_handle() );
  CGAL_triangulation_expensive_precondition( is_cell(c) );

  Vertex_handle v = create_vertex();

  Vertex_handle v0 = vertex(c, 0);
  Vertex_handle v1 = vertex(c, 1);
  Vertex_handle v2 = vertex(c, 2);
  Vertex_handle v3 = vertex(c, 3);

  Cell_handle n1 = neighbor(c, 1);
  Cell_handle n2 = neighbor(c, 2);
  Cell_handle n3 = neighbor(c, 3);

  // c will be modified to have v,v1,v2,v3 as vertices
  Cell_handle c3 = create_cell(v0,v1,v2,v);
  Cell_handle c2 = create_cell(v0,v1,v,v3);
  Cell_handle c1 = create_cell(v0,v,v2,v3);

  set_adjacency(c3, 0, c, 3);
  set_adjacency(c2, 0, c, 2);
  set_adjacency(c1, 0, c, 1);

  set_adjacency(c2, 3, c3, 2);
  set_adjacency(c1, 3, c3, 1);
  set_adjacency(c1, 2, c2, 1);

  set_adjacency(n1, index(n1, c), c1, 1);
  set_adjacency(n2, index(n2, c), c2, 2);
  set_adjacency(n3, index(n3, c), c3, 3);

  set_vertex(c, 0,v);

  set_cell(v0, c1);
  set_cell(v, c);

  return v;
}

template <class Vb, class Cb, class Ct>
typename Triangulation_data_structure_3<Vb,Cb,Ct>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb,Ct>::
insert_in_facet(Cell_handle c, int i)
{ // inserts v in the facet opposite to vertex i of cell c

  CGAL_triangulation_precondition( c != Cell_handle() );
  CGAL_triangulation_precondition( dimension() >= 2 );

  Vertex_handle v = create_vertex();

  switch ( dimension() ) {

  case 3:
    {
      CGAL_triangulation_expensive_precondition( is_cell(c) );
      CGAL_triangulation_precondition( i == 0 || i == 1 ||
                                       i == 2 || i == 3 );
      // c will be modified to have v replacing vertex(i+3)
      int i1,i2,i3;

      if ( (i&1) == 0 ) {
        i1=(i+1)&3; i2=(i+2)&3; i3=6-i-i1-i2;
      }
      else {
        i1=(i+1)&3; i2=(i+3)&3; i3=6-i-i1-i2;
      }
      // i,i1,i2,i3 is well oriented
      // so v will "replace" the vertices in this order
      // when creating the new cells one after another from c

      Vertex_handle vi=vertex(c, i);
      Vertex_handle v1=vertex(c, i1);
      Vertex_handle v2=vertex(c, i2);
      Vertex_handle v3=vertex(c, i3);

      // new cell with v in place of i1
      Cell_handle nc = neighbor(c, i1);
      Cell_handle cnew1 = create_cell(vi,v,v2,v3);
      set_adjacency(cnew1, 1, nc, index(nc, c));
      set_adjacency(cnew1, 3, c, i1);

      set_cell(v3, cnew1);

      // new cell with v in place of i2
      nc = neighbor(c, i2);
      Cell_handle cnew2 = create_cell(vi,v1,v,v3);
      set_adjacency(cnew2, 2, nc, index(nc, c));
      set_adjacency(cnew2, 3, c, i2);
      set_adjacency(cnew1, 2, cnew2, 1);

      // v replaces i3 in c
      set_vertex(c,i3,v);

      // other side of facet containing v
      Cell_handle d = neighbor(c, i);
      int j = index(d, c);
      int j1=index(d, v1);// triangulation supposed to be valid
      int j2=index(d, v2);
      int j3=6-j-j1-j2;
      // then the orientation of j,j1,j2,j3 depends on the parity
      // of i-j

      // new cell with v in place of j1
      Cell_handle nd = neighbor(d, j1);
      Cell_handle dnew1 = create_cell(vertex(d, j),v,v3,v2);
      set_adjacency(dnew1, 1, nd, index(nd, d));
      set_adjacency(dnew1, 2, d, j1);
      set_adjacency(dnew1, 0, cnew1, 0);

      // new cell with v in place of j2
      nd = neighbor(d, j2);
      Cell_handle dnew2 = create_cell(vertex(d, j),v1,v3,v);

      set_adjacency(dnew2, 3, nd, index(nd, d));
      set_adjacency(dnew2, 2, d, j2);
      set_adjacency(dnew2, 0, cnew2, 0);
      set_adjacency(dnew1, 3, dnew2, 1);

      // v replaces i3 in d
      set_vertex(d, j3,v);
      set_cell(v, d);

      break;
    }
  case 2:
    {
      CGAL_triangulation_expensive_precondition( is_facet(c,i) );
      Cell_handle n = neighbor(c, 2);
      Cell_handle cnew = create_face(vertex(c, 0),vertex(c, 1),v);
      set_adjacency(cnew, 2, n, index(n, c));
      set_adjacency(cnew, 0, c, 2);
      set_cell(vertex(c, 0), cnew);

      n = neighbor(c, 1);
      Cell_handle dnew = create_face(vertex(c, 0),v,vertex(c, 2));
      set_adjacency(dnew, 1, n, index(n, c));
      set_adjacency(dnew, 0, c, 1);
      set_adjacency(dnew, 2, cnew, 1);

      set_vertex(c,0,v);
      set_cell(v,c);
      break;
    }
  }

  return v;
}

template <class Vb, class Cb, class Ct>
typename Triangulation_data_structure_3<Vb,Cb,Ct>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb,Ct>::
insert_in_edge(Cell_handle c, int i, int j)
  // inserts a vertex in the edge of cell c with vertices i and j
{
  CGAL_triangulation_precondition( c != Cell_handle() );
  CGAL_triangulation_precondition( i != j );
  CGAL_triangulation_precondition( dimension() >= 1 );

  switch ( dimension() ) {
  case 3:
    {
      CGAL_triangulation_expensive_precondition( is_cell(c) );
      CGAL_triangulation_precondition( i>=0 && i<=3 && j>=0 && j<=3 );

      std::vector<Cell_handle > cells;
      cells.reserve(32);
      Cell_circulator ccir = incident_cells(c, i, j);
      do {
          Cell_handle cc = ccir;
          cells.push_back(cc);
          tds_data(cc).mark_in_conflict();
          ++ccir;
      } while (c != ccir);

      return _insert_in_hole(cells.begin(), cells.end(), c, i);
    }
  case 2:
    {
      CGAL_triangulation_expensive_precondition( is_edge(c,i,j) );

      Vertex_handle v = create_vertex();
      int k=3-i-j; // index of the third vertex of the facet
      Cell_handle d = neighbor(c, k);
      int kd = index(d, c);
      int id = index(d, vertex(c, i));
      int jd = index(d, vertex(c, j));

      Cell_handle cnew = create_cell();
      set_vertex(cnew, i, vertex(c, i));
      set_cell(vertex(c, i), cnew);
      set_vertex(cnew, j, v);
      set_vertex(cnew, k, vertex(c, k));
      set_vertex(c, i, v);

      Cell_handle dnew = create_cell();
      set_vertex(dnew, id, vertex(d, id));
      // vertex(d, id)->cell() is cnew OK
      set_vertex(dnew, jd,v);
      set_vertex(dnew, kd,vertex(d, kd));
      set_vertex(d, id,v);

      Cell_handle nj = neighbor(c, j);
      set_adjacency(cnew, i, c, j);
      set_adjacency(cnew, j, nj, index(nj, c));

      nj = neighbor(d, jd);
      set_adjacency(dnew, id, d, jd);
      set_adjacency(dnew, jd, nj, index(nj, d));

      set_adjacency(cnew, k, dnew, kd);

      set_cell(v, cnew);
      return v;
    }
  default: // case 1:
    {
      Vertex_handle v = create_vertex();
      CGAL_triangulation_expensive_precondition( is_edge(c,i,j) );
      Cell_handle cnew = create_face(v, vertex(c, 1), Vertex_handle());
      set_cell(vertex(c, 1), cnew);
      set_vertex(c, 1, v);
      set_adjacency(cnew, 0, neighbor(c, 0), 1);
      set_adjacency(cnew, 1, c, 0);

      set_cell(v, cnew);
      return v;
    }
  }
}

template <class Vb, class Cb, class Ct>
typename Triangulation_data_structure_3<Vb,Cb,Ct>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb,Ct>::
insert_increase_dimension(Vertex_handle star)
  // star = vertex from which we triangulate the facet of the
  // incremented dimension
  // ( geometrically : star = infinite vertex )
  // = nullptr only used to insert the 1st vertex (dimension -2 to dimension -1)
  // changes the dimension
{
  CGAL_triangulation_precondition( dimension() < 3);

  Vertex_handle v = create_vertex();

  int dim = dimension();
  if (dim != -2) {
      CGAL_triangulation_precondition( star != Vertex_handle() );
      // In this case, this precondition is not relatively expensive.
      CGAL_triangulation_precondition( is_vertex(star) );
  }

  // this is set now, so that it becomes allowed to reorient
  // new facets or cells by iterating on them (otherwise the
  // dimension is too small)
  set_dimension( dimension()+1 );

  switch ( dim ) {
  case -2:
      // insertion of the first vertex
      // ( geometrically : infinite vertex )
    {
      Cell_handle c = create_face(v, Vertex_handle(), Vertex_handle());
      set_cell(v,c);
      break;
    }

  case -1:
    // insertion of the second vertex
    // ( geometrically : first finite vertex )
    {
      Cell_handle d = create_face(v, Vertex_handle(), Vertex_handle());
      set_cell(v, d);
      set_adjacency(d, 0, cell(star), 0);
      break;
    }

  case 0:
    // insertion of the third vertex
    // ( geometrically : second finite vertex )
    {
      Cell_handle c = cell(star);
      Cell_handle d = neighbor(c, 0);

      set_vertex(c, 1,vertex(d, 0));
      set_vertex(d,1,v);
      set_neighbor(d,1,c);
      Cell_handle e = create_face( v, star, Vertex_handle() );
      set_adjacency(e, 0, c, 1);
      set_adjacency(e, 1, d, 0);

      set_cell(v,d);
      break;
    }

  case 1:
    // general case : 4th vertex ( geometrically : 3rd finite vertex )
    // degenerate cases geometrically : 1st non collinear vertex
    {
      Cell_handle c = cell(star);
      int i = index(c, star); // i== 0 or 1
      CGAL_assertion(i==0 || i==1);
      int j = (i == 0) ? 1 : 0;
      Cell_handle d = neighbor(c, j);
        
      set_vertex(c, 2,v);

      Cell_handle e = neighbor(c, i);
      Cell_handle cnew = c;
      Cell_handle enew = Cell_handle();
        
      while( e != d ){
        enew = create_cell();
        set_vertex(enew, i,vertex(e, j));
        set_vertex(enew, j,vertex(e, i));
        set_vertex(enew, 2,star);
        
        set_adjacency(enew, i, cnew, j);
        // false at the first iteration of the loop where it should
        // be neighbor 2
        // it is corrected after the loop
        set_adjacency(enew, 2, e, 2);
        // neighbor j will be set during next iteration of the loop
        
        set_vertex(e, 2, v);

        e = neighbor(e, i);
        cnew = enew;
      }
        
      set_vertex(d, 2, v);
      set_adjacency(enew, j, d, 2);

      // corrections for cell(star) :
      c = cell(star);
      set_neighbor(c, 2, neighbor(neighbor(c, i), 2));
      set_neighbor(c, j,d);

      set_cell(v, d);

      break;
    }

  case 2:
    // general case : 5th vertex ( geometrically : 4th finite vertex )
    // degenerate cases : geometrically 1st non coplanar vertex
    {
      // used to store the new cells, in order to be able to traverse only
      // them to find the missing neighbors.
      std::vector<Cell_handle > new_cells;
      new_cells.reserve(16);

      Cell_iterator it = cells_begin();
      // allowed since the dimension has already been set to 3

      set_cell(v, it); // ok since there is at least one ``cell''
      for(; it != cells_end(); ++it) {
        // Here we must be careful since we create_cells in a loop controlled
        // by an iterator.  So we first take care of the cells newly created
        // by the following test :
        if (neighbor(it, 0) == Cell_handle())
          continue;
        set_neighbor(it, 3, Cell_handle());
        set_vertex(it, 3, v);
        if ( ! has_vertex(it, star) ) {
          Cell_handle cnew = create_cell( vertex(it, 0), vertex(it, 2),
                                          vertex(it, 1), star);
          // The Intel compiler has a problem with passing "it" directly to
          // function "set_adjacency": the adjacency is not changed.
          Cell_handle ch_it = it;
          set_adjacency(cnew, 3, ch_it, 3);
          set_neighbor(cnew, 0, Cell_handle());
          new_cells.push_back(cnew);
        }
      }

      // traversal of the new cells only, to add missing neighbors
      for(typename std::vector<Cell_handle>::iterator ncit = new_cells.begin();
           ncit != new_cells.end(); ++ncit) {
        Cell_handle n = neighbor(*ncit,3); // opposite to star
        for ( int i=0; i<3; i++ ) {
          int j;
          if ( i==0 ) j=0;
          else j=3-i; // vertex 1 and vertex 2 are always switched when
          // creating a new cell (see above)
          Cell_handle c = neighbor(neighbor(n, i), 3);
          if ( c != Cell_handle() ) {
            // i.e. star is not a vertex of neighbor(n, i)
            set_neighbor(*ncit, j, c);
            // opposite relation will be set when ncit arrives on c
            // this avoids to look for the correct index
            // and to test whether *ncit already has neighbor i
          }
          else {
            // star is a vertex of neighbor(n, i)
            set_adjacency(*ncit, j, neighbor(n, i), 3);//neighbor opposite to v
          }
        }
      }
    }
  }// end switch

  return v;
}

template <class Vb, class Cb, class Ct>
void
Triangulation_data_structure_3<Vb,Cb,Ct>::
remove_decrease_dimension(Vertex_handle v, Vertex_handle w)
{
    CGAL_triangulation_expensive_precondition( is_valid() );
    CGAL_triangulation_precondition( dimension() >= -1 );
    CGAL_triangulation_precondition( dimension() != 1 ||
                                     number_of_vertices() == 3);
    CGAL_triangulation_precondition( number_of_vertices() >
                                     (size_type) dimension() + 1 );
    CGAL_triangulation_precondition( degree(v) == number_of_vertices()-1 );

    if (dimension() <= 0) {
        delete_cell(cell(v));
    }
    else {
        // the cells incident to w are down graded one dimension
        // the other cells are deleted
        std::vector<Cell_handle> to_delete, to_downgrade;

        for (Cell_iterator ib = cells().begin();
            ib != cells().end(); ++ib) {
          if ( has_vertex(ib, w) )
                to_downgrade.push_back(ib);
            else
                to_delete.push_back(ib);
        }

        typename std::vector<Cell_handle>::iterator lfit=to_downgrade.begin();
        for( ; lfit != to_downgrade.end(); ++lfit) {
            Cell_handle f = *lfit;
            int j = index(f, w);
            int k; if (has_vertex(f, v, k)) set_vertex(f, k, w);
            if (j != dimension()) {
              set_vertex(f, j, vertex(f, dimension()));
              set_neighbor(f, j, neighbor(f, dimension()));
                if (dimension() >= 1)
                    change_orientation(f);
            }
            set_vertex(f, dimension(), Vertex_handle());
            set_neighbor(f, dimension(), Cell_handle());

            // Update cell(vertex) pointers.
            for (int i = 0; i < dimension(); ++i)
              set_cell(vertex(f, i), f);
        }

        delete_cells(to_delete.begin(), to_delete.end());
    }
    delete_vertex(v);
    set_dimension(dimension()-1);
    CGAL_triangulation_postcondition(is_valid());
}

template <class Vb, class Cb, class Ct>
typename Triangulation_data_structure_3<Vb,Cb,Ct>::Cell_handle
Triangulation_data_structure_3<Vb,Cb,Ct>::
remove_from_maximal_dimension_simplex(Vertex_handle v)
{
    CGAL_triangulation_precondition(dimension() >= 1);
    CGAL_triangulation_precondition(degree(v) == (size_type) dimension() + 1);
    CGAL_triangulation_precondition(number_of_vertices() >
                                    (size_type) dimension() + 1);

    if (number_of_vertices() == (size_type) dimension() + 2) {
        remove_decrease_dimension(v);
        return Cell_handle();
    }

    if (dimension() == 3)
        return remove_degree_4(v);
    if (dimension() == 2)
        return remove_degree_3(v);

    // dimension() == 1
    return remove_degree_2(v);
}

template <class Vb, class Cb, class Ct>
typename Triangulation_data_structure_3<Vb,Cb,Ct>::Cell_handle
Triangulation_data_structure_3<Vb,Cb,Ct>::
remove_degree_2(Vertex_handle v)
{
    CGAL_triangulation_precondition(dimension() == 1);
    CGAL_triangulation_precondition(degree(v) == 2);
    CGAL_triangulation_precondition(number_of_vertices() >= 4);

    // Cells to be killed.
    Cell_handle c0, c1;
    // Indices of v in these cells.
    int i0, i1;

    c0 = cell(v);
    i0 = index(c0, v);
    c1 = neighbor(c0, (i0 == 0) ? 1 : 0);
    i1 = index(c1, v);

    // New cell : we copy the content of c0, so we keep the orientation.
    Cell_handle newc = create_face(vertex(c0, 0),
                                   vertex(c0, 1),
                                   Vertex_handle());

    set_vertex(newc, i0, vertex(c1, index(c1, c0)));

    set_adjacency(newc,    i0, neighbor(c0, i0), mirror_index(c0, i0));
    set_adjacency(newc,  1-i0, neighbor(c1, i1), mirror_index(c1, i1));

    set_cell(vertex(newc, 0), newc);
    set_cell(vertex(newc, 1), newc);

    delete_cell(c0);
    delete_cell(c1);
    delete_vertex(v);

    return newc;
}

template <class Vb, class Cb, class Ct>
typename Triangulation_data_structure_3<Vb,Cb,Ct>::Cell_handle
Triangulation_data_structure_3<Vb,Cb,Ct>::
remove_degree_3(Vertex_handle v)
{
    CGAL_triangulation_precondition(dimension() == 2);
    CGAL_triangulation_precondition(degree(v) == 3);
    CGAL_triangulation_precondition(number_of_vertices() >= 5);

    // Cells to be killed.
    Cell_handle c0, c1, c2;
    // Indices of v in these cells.
    int i0, i1, i2;

    c0 = cell(v);
    i0 = index(c0, v);
    c1 = neighbor(c0, cw(i0));
    i1 = index(c1, v);
    c2 = neighbor(c0, ccw(i0));
    i2 = index(c2, v);

    // New cell : we copy the content of c0, so we keep the orientation.
    Cell_handle newc = create_face(vertex(c0, 0),
                                   vertex(c0, 1),
                                   vertex(c0, 2));

    set_vertex(newc, i0, vertex(c1, index(c1, c0)));

    set_adjacency(newc,      i0, neighbor(c0, i0), mirror_index(c0, i0));
    set_adjacency(newc,  cw(i0), neighbor(c1, i1), mirror_index(c1, i1));
    set_adjacency(newc, ccw(i0), neighbor(c2, i2), mirror_index(c2, i2));

    set_cell(vertex(newc, 0), newc);
    set_cell(vertex(newc, 1), newc);
    set_cell(vertex(newc, 2), newc);

    delete_cell(c0);
    delete_cell(c1);
    delete_cell(c2);
    delete_vertex(v);

    return newc;
}

template <class Vb, class Cb, class Ct>
typename Triangulation_data_structure_3<Vb,Cb,Ct>::Cell_handle
Triangulation_data_structure_3<Vb,Cb,Ct>::
remove_degree_4(Vertex_handle v)
{
    CGAL_triangulation_precondition(dimension() == 3);
    CGAL_triangulation_precondition(degree(v) == 4);
    CGAL_triangulation_precondition(number_of_vertices() >= 6);

    // Cells to be killed.
    Cell_handle c0, c1, c2, c3;
    // Indices of v in these cells.
    int i0, i1, i2, i3;

    c0 = cell(v);
    i0 = index(c0, v);
    c1 = neighbor(c0, i0^1);
    i1 = index(c1, v);
    c2 = neighbor(c0, i0^2);
    i2 = index(c2, v);
    c3 = neighbor(c0, i0^3);
    i3 = index(c3, v);

    // New cell : we copy the content of c0, so we keep the orientation.
    Cell_handle newc = create_cell(vertex(c0, 0),
                                   vertex(c0, 1),
                                   vertex(c0, 2),
                                   vertex(c0, 3));

    set_vertex(newc, i0, vertex(c1, index(c1, c0)));

    set_adjacency(newc,   i0, neighbor(c0, i0), mirror_index(c0, i0));
    set_adjacency(newc, i0^1, neighbor(c1, i1), mirror_index(c1, i1));
    set_adjacency(newc, i0^2, neighbor(c2, i2), mirror_index(c2, i2));
    set_adjacency(newc, i0^3, neighbor(c3, i3), mirror_index(c3, i3));

    set_cell(vertex(newc, 0), newc);
    set_cell(vertex(newc, 1), newc);
    set_cell(vertex(newc, 2), newc);
    set_cell(vertex(newc, 3), newc);

    delete_cell(c0);
    delete_cell(c1);
    delete_cell(c2);
    delete_cell(c3);
    delete_vertex(v);

    return newc;
}

template <class Vb, class Cb, class Ct>
void
Triangulation_data_structure_3<Vb,Cb,Ct>::
decrease_dimension(Cell_handle c, int i)
{
  CGAL_triangulation_expensive_precondition( is_valid() );;
  CGAL_triangulation_precondition( dimension() >= 2);
  CGAL_triangulation_precondition( number_of_vertices() >
                                   (size_type) dimension() + 1 );
  CGAL_triangulation_precondition( degree(vertex(c, i)) == number_of_vertices()-1 );

  Vertex_handle v = vertex(c, i);
  Vertex_handle w = vertex(c, i);

  // the cells incident to w are down graded one dimension
  // the other cells are deleted
  std::vector<Cell_handle> to_delete, to_downgrade;

  for (Cell_iterator ib = cells().begin();
       ib != cells().end(); ++ib) {
    if ( has_vertex(ib, w) )
      to_downgrade.push_back(ib);
    else
      to_delete.push_back(ib);
  }

  typename std::vector<Cell_handle>::iterator lfit=to_downgrade.begin();
  for( ; lfit != to_downgrade.end(); ++lfit) {
    Cell_handle f = *lfit;
    int j = index(f, w);
    int k; 
    if (has_vertex(f,v, k)) set_vertex(f, k, w);
    if (j != dimension()) {
      set_vertex(f, j, vertex(f, dimension()));
      set_neighbor(f, j, neighbor(f, dimension()));
      if (dimension() >= 1)
        change_orientation(f);
    }
    set_vertex(f, dimension(), Vertex_handle());
    set_neighbor(f, dimension(), Cell_handle());

    // Update cell(vertex) pointers.
    for (int i = 0; i < dimension(); ++i)
      set_cell(vertex(f, i), f);
  }

  delete_cells(to_delete.begin(), to_delete.end());

  //delete_vertex(v);
  set_dimension(dimension()-1);

  if(dimension() == 2)
  {
    Cell_handle n0 = neighbor(c, 0);
    Cell_handle n1 = neighbor(c, 1);
    Cell_handle n2 = neighbor(c, 2);
    Vertex_handle v0 = vertex(c, 0);
    Vertex_handle v1 = vertex(c, 1);
    Vertex_handle v2 = vertex(c, 2);
                
    int i0 = 0, i1 = 0, i2 = 0;
                
    for(int i=0; i<3; i++) if(neighbor(n0, i) == c) { i0 = i; break; }
    for(int i=0; i<3; i++) if(neighbor(n1, i) == c) { i1 = i; break; }
    for(int i=0; i<3; i++) if(neighbor(n2, i) == c) { i2 = i; break; }
                
    Cell_handle c1 = create_cell(v, v0, v1, Vertex_handle());
    Cell_handle c2 = create_cell(v, v1, v2, Vertex_handle());

    set_vertex(c, 0, v);
    set_vertex(c, 1, v2);
    set_vertex(c, 2, v0);
    set_vertex(c, 3, Vertex_handle());

    //Cell_handle c3 = create_cell(v, v2, v0, Vertex_handle());
    Cell_handle c3 = c;
                
    set_neighbor(c1, 0, n2); set_neighbor(n2, i2, c1);
    set_neighbor(c1, 1, c2); 
    set_neighbor(c1, 2, c3);
    set_neighbor(c1, 3, Cell_handle());
                
    set_neighbor(c2, 0, n0); set_neighbor(n0, i0, c2);
    set_neighbor(c2, 1, c3); 
    set_neighbor(c2, 2, c1);
    set_neighbor(c2, 3, Cell_handle());
                
    set_neighbor(c3, 0, n1); set_neighbor(n1, i1, c3);
    set_neighbor(c3, 1, c1); 
    set_neighbor(c3, 2, c2);
    set_neighbor(c3, 3, Cell_handle());
                
    set_cell(v, c1);
    set_cell(v0, c1);
    set_cell(v1, c1);
    set_cell(v2, c2);
  }
        
  if(dimension() == 1)
  {
    Cell_handle n0 = neighbor(c, 0);
    Cell_handle n1 = neighbor(c, 1);
    Vertex_handle v0 = vertex(c, 0);
    Vertex_handle v1 = vertex(c, 1);
                
    int i0 = 0 , i1 = 0;
                
    for(int i=0; i<2; i++) if(neighbor(n0, i) == c) { i0 = i; break; }
    for(int i=0; i<2; i++) if(neighbor(n1, i) == c) { i1 = i; break; }
                
    Cell_handle c1 = create_cell(v0, v, Vertex_handle(), Vertex_handle());
                
    set_vertex(c, 0, v);
    set_vertex(c, 1, v1);
    set_vertex(c, 2,Vertex_handle());
    set_vertex(c, 3, Vertex_handle());

    //Cell_handle c2 = create_cell(v, v1, Vertex_handle(), Vertex_handle());
    Cell_handle c2 = c;
                
    set_neighbor(c1, 0, c2); 
    set_neighbor(c1, 1, n1); set_neighbor(n1, i1, c1);
    set_neighbor(c1, 2, Cell_handle());
    set_neighbor(c1, 3, Cell_handle());
                
    set_neighbor(c2, 0, n0); set_neighbor(n0, i0, c2);
    set_neighbor(c2, 1, c1); 
    set_neighbor(c2, 2, Cell_handle());
    set_neighbor(c2, 3, Cell_handle());
                
    set_cell(v, c1);
    set_cell(v0, c1);
    set_cell(v1, c2);
  }
        
  CGAL_triangulation_postcondition(is_valid());
}


template <class Vb, class Cb, class Ct>
typename Triangulation_data_structure_3<Vb,Cb,Ct>::size_type
Triangulation_data_structure_3<Vb,Cb,Ct>::
degree(Vertex_handle v) const
{
    std::size_t res;
    adjacent_vertices(v, Counting_output_iterator(&res));
    return res;
}

template <class Vb, class Cb, class Ct>
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
is_valid(bool verbose, int level ) const
{
  switch ( dimension() ) {
  case 3:
    {
        
      if(number_of_vertices() <= 4) {
        if (verbose)
          std::cerr << "wrong number of vertices" << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }
        
      size_type vertex_count;
      if ( ! count_vertices(vertex_count,verbose,level) )
        return false;
      if ( number_of_vertices() != vertex_count ) {
        if (verbose)
          std::cerr << "wrong number of vertices" << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }

      size_type cell_count;
      if ( ! count_cells(cell_count,verbose,level) )
        return false;
      size_type edge_count;
      if ( ! count_edges(edge_count,verbose,level) )
          return false;
      size_type facet_count;
      if ( ! count_facets(facet_count,verbose,level) )
          return false;

      // Euler relation
      if ( cell_count - facet_count + edge_count - vertex_count != 0 ) {
        if (verbose)
            std::cerr << "Euler relation unsatisfied" << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }

      break;
    }
  case 2:
    {
        
      if(number_of_vertices() <= 3) {
        if (verbose)
          std::cerr << "wrong number of vertices" << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }
        
      size_type vertex_count;
      
      if ( ! count_vertices(vertex_count,verbose,level) )
        return false;
      if ( number_of_vertices() != vertex_count ) {
        if (verbose)
            std::cerr << "false number of vertices" << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }

      size_type edge_count;
      if ( ! count_edges(edge_count,verbose,level) )
          return false;
      // Euler for edges
      if ( edge_count != 3 * vertex_count - 6 ) {
        if (verbose)
            std::cerr << "Euler relation unsatisfied - edges/vertices"
                      << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }

      size_type facet_count;
      if ( ! count_facets(facet_count,verbose,level) )
          return false;
      // Euler for facets
      if ( facet_count != 2 * vertex_count - 4 ) {
        if (verbose)
            std::cerr << "Euler relation unsatisfied - facets/vertices"
                      << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }
      break;
    }
  case 1:
    {
        
      if(number_of_vertices() <= 1) {
        if (verbose)
          std::cerr << "wrong number of vertices" << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }
        
      size_type vertex_count;
      if ( ! count_vertices(vertex_count,verbose,level) )
          return false;
      if ( number_of_vertices() != vertex_count ) {
        if (verbose)
            std::cerr << "false number of vertices" << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }
      size_type edge_count;
      if ( ! count_edges(edge_count,verbose,level) )
          return false;
      // Euler for edges
      if ( edge_count != vertex_count ) {
        if (verbose)
            std::cerr << "false number of edges" << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }
      break;
    }
  case 0:
    {
      if ( number_of_vertices() < 2 ) {
        if (verbose)
            std::cerr << "less than 2 vertices but dimension 0" << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }
      CGAL_FALLTHROUGH;
    }
  case -1:
    {
      if ( number_of_vertices() < 1 ) {
        if (verbose)
          std::cerr << "no vertex but dimension -1" << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }
      // vertex count
      size_type vertex_count;
      if ( ! count_vertices(vertex_count,verbose,level) )
        return false;
      if ( number_of_vertices() != vertex_count ) {
        if (verbose)
          std::cerr << "false number of vertices" << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }
    }
  } // end switch
  if (verbose)
      std::cerr << "valid data structure" << std::endl;
  return true;
}

template <class Vb, class Cb, class Ct>
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
is_valid(Vertex_handle v, bool verbose, int level) const
{
  bool result = v->is_valid(verbose,level);
  result = result && has_vertex(cell(v), v);
  if ( ! result ) {
    if ( verbose )
      std::cerr << "invalid vertex" << std::endl;
    CGAL_triangulation_assertion(false);
  }
  return result;
}

template <class Vb, class Cb, class Ct>
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
is_valid(Cell_handle c, bool verbose, int level) const
{
  if ( ! c->iis_valid(verbose, level) )
        return false;

    switch (dimension()) {
    case -2:
    case -1:
    {
      if ( vertex(c, 0) == Vertex_handle() ) {
        if (verbose)
            std::cerr << "vertex 0 nullptr" << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }
      is_valid(vertex(c, 0),verbose,level);
      if ( vertex(c, 1) != Vertex_handle() || vertex(c, 2) != Vertex_handle()) {
        if (verbose)
            std::cerr << "vertex 1 or 2 != nullptr" << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }
      if ( neighbor(c, 0) != Cell_handle() ||
           neighbor(c, 1) != Cell_handle() ||
           neighbor(c, 2) != Cell_handle()) {
        if (verbose)
            std::cerr << "one neighbor != nullptr" << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }
      break;
    }

    case 0:
      {
        if ( vertex(c, 0) == Vertex_handle() ) {
        if (verbose)
            std::cerr << "vertex 0 nullptr" << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }
        is_valid(vertex(c, 0),verbose,level);
        if ( neighbor(c, 0) == Cell_handle() ) {
        if (verbose)
            std::cerr << "neighbor 0 nullptr" << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }
      if ( vertex(c, 1) != Vertex_handle() ||
           vertex(c, 2) != Vertex_handle() ) {
        if (verbose)
            std::cerr << "vertex 1 or 2 != nullptr" << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }
      if ( neighbor(c, 1) != Cell_handle() ||
           neighbor(c, 2) != Cell_handle() ) {
        if (verbose)
            std::cerr << "neighbor 1 or 2 != nullptr" << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }

      if ( ! has_vertex(neighbor(c, 0), vertex(c, 0)) ) {
        if (verbose)
            std::cerr << "neighbor 0 does not have vertex 0" << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }
      break;
      }

    case 1:
      {
        Vertex_handle v0 = vertex(c, 0);
        Vertex_handle v1 = vertex(c, 1);
        Cell_handle n0 = neighbor(c, 0);
        Cell_handle n1 = neighbor(c, 1);

      if ( v0 == Vertex_handle() || v1 == Vertex_handle() ) {
        if (verbose)
            std::cerr << "vertex 0 or 1 nullptr" << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }
      is_valid(vertex(c, 0),verbose,level);
      is_valid(vertex(c, 1),verbose,level);
      if ( n0 == Cell_handle() || n1 == Cell_handle() ) {
        if (verbose)
            std::cerr << "neighbor 0 or 1 nullptr" << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }

      if ( v0 !=  vertex(n1, 1) ) {
        if (verbose)
            std::cerr << "neighbor 1 does not have vertex 0 as vertex 1"
                      << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }
      if ( v1 != vertex(n0, 0) ) {
        if (verbose)
            std::cerr << "neighbor 0 does not have vertex 1 as vertex 0"
                      << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }

      if ( neighbor(n0, 1) != c ) {
        if (verbose)
            std::cerr << "neighbor 0 does not have this as neighbor 1"
                      << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }
      if ( neighbor(n1, 0) != c ) {
        if (verbose)
            std::cerr << "neighbor 1 does not have this as neighbor 0"
                      << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }

      break;
      }

    case 2:
      {
        if ( vertex(c, 0) == Vertex_handle() ||
             vertex(c, 1) == Vertex_handle() ||
             vertex(c, 2) == Vertex_handle() ) {
        if (verbose)
            std::cerr << "vertex 0, 1, or 2 nullptr" << std::endl;
        CGAL_triangulation_assertion(false);
        return false;
      }
        is_valid(vertex(c, 0),verbose,level);
        is_valid(vertex(c, 1),verbose,level);
        is_valid(vertex(c, 2),verbose,level);
      int in;
      Cell_handle n;
      for(int i = 0; i < 3; i++) {
        n = neighbor(c, i);
        if ( n == Cell_handle() ) {
          if (verbose)
              std::cerr << "neighbor " << i << " nullptr" << std::endl;
          CGAL_triangulation_assertion(false);
          return false;
        }
        if ( ! has_vertex(n, vertex(c, cw(i)),in ) ) {
          if (verbose)
              std::cerr << "vertex " << cw(i)
                        << " not vertex of neighbor " << i << std::endl;
          CGAL_triangulation_assertion(false);
          return false;
        }
        in = cw(in);
        if ( neighbor(n, in) != c ) {
          if (verbose)
              std::cerr << "neighbor " << i
                        << " does not have this as neighbor "
                        << in << std::endl;
          CGAL_triangulation_assertion(false);
          return false;
        }
        if ( vertex(c, ccw(i)) != vertex(n, cw(in)) ) {
          if (verbose)
              std::cerr << "vertex " << ccw(i)
                        << " is not vertex " << cw(in)
                        << " of neighbor " << i << std::endl;
          CGAL_triangulation_assertion(false);
          return false;
        }
      }
      break;
      }

    case 3:
      {
        int i;
        for(i = 0; i < 4; i++) {
          if ( vertex(c, i) == Vertex_handle() ) {
            if (verbose)
                std::cerr << "vertex " << i << " nullptr" << std::endl;
            CGAL_triangulation_assertion(false);
            return false;
          }
          is_valid(vertex(c, i),verbose,level);
        }

        for(i = 0; i < 4; i++) {
          Cell_handle n = neighbor(c, i);
          if ( n == Cell_handle() ) {
            if (verbose)
              std::cerr << "neighbor " << i << " nullptr" << std::endl;
            CGAL_triangulation_assertion(false);
            return false;
          }

          int in = 5;
          // if ( ! n->has_neighbor(handle(), in) ) {
          if ( neighbor(n, 0) == c) in = 0;
          if ( neighbor(n, 1) == c) in = 1;
          if ( neighbor(n, 2) == c) in = 2;
          if ( neighbor(n, 3) == c) in = 3;
          if (in == 5) {
            if (verbose)
              std::cerr << "neighbor of c has not c as neighbor" << std::endl;
            CGAL_triangulation_assertion(false);
            return false;
          }
        
          int j1n=4,j2n=4,j3n=4;
          if ( ! has_vertex(n, vertex(c, (i+1)&3),j1n) ) {
            if (verbose) { std::cerr << "vertex " << ((i+1)&3)
                                     << " not vertex of neighbor "
                                     << i << std::endl; }
            CGAL_triangulation_assertion(false);
            return false;
          }
          if ( ! has_vertex(n, vertex(c, (i+2)&3),j2n) ) {
            if (verbose) { std::cerr << "vertex " << ((i+2)&3)
                                     << " not vertex of neighbor "
                                     << i << std::endl; }
            CGAL_triangulation_assertion(false);
            return false;
          }
          if ( ! has_vertex(n, vertex(c, (i+3)&3),j3n) ) {
            if (verbose) { std::cerr << "vertex " << ((i+3)&3)
                                     << " not vertex of neighbor "
                                     << i << std::endl; }
            CGAL_triangulation_assertion(false);
            return false;
          }
        
          if ( in+j1n+j2n+j3n != 6) {
            if (verbose) { std::cerr << "sum of the indices != 6 "
                                     << std::endl; }
            CGAL_triangulation_assertion(false);
            return false;
          }
        
          // tests whether the orientations of this and n are consistent
          if ( ((i+in)&1) == 0 ) { // i and in have the same parity
            if ( j1n == ((in+1)&3) ) {
              if ( ( j2n != ((in+3)&3) ) || ( j3n != ((in+2)&3) ) ) {
                if (verbose)
                  std::cerr << " pb orientation with neighbor "
                            << i << std::endl;
                CGAL_triangulation_assertion(false);
                return false;
              }
            }
            if ( j1n == ((in+2)&3) ) {
              if ( ( j2n != ((in+1)&3) ) || ( j3n != ((in+3)&3) ) ) {
                if (verbose)
                  std::cerr << " pb orientation with neighbor "
                            << i << std::endl;
                CGAL_triangulation_assertion(false);
                return false;
              }
            }
            if ( j1n == ((in+3)&3) ) {
              if ( ( j2n != ((in+2)&3) ) || ( j3n != ((in+1)&3) ) ) {
                if (verbose)
                  std::cerr << " pb orientation with neighbor "
                            << i << std::endl;
                CGAL_triangulation_assertion(false);
                return false;
              }
            }
          }
          else { // i and in do not have the same parity
            if ( j1n == ((in+1)&3) ) {
              if ( ( j2n != ((in+2)&3) ) || ( j3n != ((in+3)&3) ) ) {
                if (verbose)
                  std::cerr << " pb orientation with neighbor "
                            << i << std::endl;
                CGAL_triangulation_assertion(false);
                return false;
              }
            }
            if ( j1n == ((in+2)&3) ) {
              if ( ( j2n != ((in+3)&3) ) || ( j3n != ((in+1)&3) ) ) {
                if (verbose)
                  std::cerr << " pb orientation with neighbor "
                            << i << std::endl;
                CGAL_triangulation_assertion(false);
                return false;
              }
            }
            if ( j1n == ((in+3)&3) ) {
              if ( ( j2n != ((in+1)&3) ) || ( j3n != ((in+2)&3) ) ) {
                if (verbose)
                  std::cerr << " pb orientation with neighbor "
                            << i << std::endl;
                CGAL_triangulation_assertion(false);
                return false;
              }
            }
          }
        } // end looking at neighbors
      }// end case dim 3
    } // end switch
    return true;
}

template <class Vb, class Cb, class Ct>
template <class TDS_src,class ConvertVertex,class ConvertCell>
typename Triangulation_data_structure_3<Vb,Cb,Ct>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb,Ct>::
copy_tds(const TDS_src& tds,
        typename TDS_src::Vertex_handle vert,
        const ConvertVertex& convert_vertex,
        const ConvertCell& convert_cell)
{
  CGAL_triangulation_expensive_precondition( vert == Vertex_handle()
                                          || tds.is_vertex(vert) );

  clear();

  size_type n = tds.number_of_vertices();
  set_dimension(tds.dimension());

  if (n == 0)  return Vertex_handle(); 

  // Number of pointers to cell/vertex to copy per cell.
  int dim = (std::max)(1, dimension() + 1);

  // Create the vertices.
  std::vector<typename TDS_src::Vertex_handle> TV(n);
  size_type i = 0;

  for (typename TDS_src::Vertex_iterator vit = tds.vertices_begin();
       vit != tds.vertices_end(); ++vit)
    TV[i++] = vit;

  CGAL_triangulation_assertion( i == n );

  Unique_hash_map< typename TDS_src::Vertex_handle,Vertex_handle > V;
  Unique_hash_map< typename TDS_src::Cell_handle,Cell_handle > F;
  
  for (i=0; i <= n-1; ++i){
    Vertex_handle vh=create_vertex( convert_vertex(*TV[i]) );
    V[ TV[i] ] = vh;
    convert_vertex(*TV[i],*vh);
  }

  // Create the cells.
  for (typename TDS_src::Cell_iterator cit = tds.cells().begin();
          cit != tds.cells_end(); ++cit) {
      Cell_handle ch=create_cell(convert_cell(*cit));
      F[cit]=ch;
      for (int j = 0; j < dim; j++)
        set_vertex(ch, j, V[vertex(cit,j)] );
      convert_cell(*cit,*ch);
  }

  // Link the vertices to a cell.
  for (typename TDS_src::Vertex_iterator vit2 = tds.vertices_begin();
       vit2 != tds.vertices_end(); ++vit2)
    set_cell(V[vit2], F[cell(vit2)] );

  // Hook neighbor pointers of the cells.
  for (typename TDS_src::Cell_iterator cit2 = tds.cells().begin();
          cit2 != tds.cells_end(); ++cit2) {
    for (int j = 0; j < dim; j++)
      set_neighbor(F[cit2], j, F[neighbor(cit2, j)] );
  }

  CGAL_triangulation_postcondition( is_valid() );

  return (vert != typename TDS_src::Vertex_handle()) ? V[vert] : Vertex_handle();
}

//utilities for copy_tds
namespace internal { namespace TDS_3{
  template <class Vertex_src,class Vertex_tgt>
  struct Default_vertex_converter
  {
    Vertex_tgt operator()(const Vertex_src& src) const {
      return Vertex_tgt(src.point());
    }
    
    void operator()(const Vertex_src&,Vertex_tgt&) const {}
  };

  template <class Cell_src,class Cell_tgt>
  struct Default_cell_converter
  {
    Cell_tgt operator()(const Cell_src&) const {
      return Cell_tgt();
    }
    
    void operator()(const Cell_src&,Cell_tgt&) const {}
  };
  
  template <class Vertex>
  struct Default_vertex_converter<Vertex,Vertex>
  {
    const Vertex& operator()(const Vertex& src) const {
      return src;
    }
    
    void operator()(const Vertex&,Vertex&) const {}
  };
  
  template <class Cell>
  struct Default_cell_converter<Cell,Cell>{
    const Cell& operator()(const Cell& src) const {
      return src;
    } 
    
    void operator()(const Cell&,Cell&) const {}
  };
} } //namespace internal::TDS_3

template <class Vb, class Cb, class Ct>
template<class TDS_src>
typename Triangulation_data_structure_3<Vb,Cb,Ct>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb,Ct>::
copy_tds(const TDS_src& src,typename TDS_src::Vertex_handle vert)
{
  internal::TDS_3::Default_vertex_converter<typename TDS_src::Vertex,Vertex> setv;
  internal::TDS_3::Default_cell_converter<typename TDS_src::Cell,Cell>  setc;
  return copy_tds(src,vert,setv,setc);
}

template <class Vb, class Cb, class Ct>
void
Triangulation_data_structure_3<Vb,Cb,Ct>::
swap(Tds & tds)
{
  CGAL_triangulation_expensive_precondition(tds.is_valid() && is_valid());

  std::swap(_dimension, tds._dimension);
  cells().swap(tds.cells());
  vertices().swap(tds.vertices());
}

template <class Vb, class Cb, class Ct>
void
Triangulation_data_structure_3<Vb,Cb,Ct>::
clear()
{
  cells().clear();
  vertices().clear();
  set_dimension(-2);
}

template <class Vb, class Cb, class Ct>
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
count_vertices(size_type & i, bool verbose, int level) const
  // counts AND checks the validity
{
  i = 0;

  for (Vertex_iterator it = vertices_begin(); it != vertices_end(); ++it) {
    if ( ! is_valid(it,verbose,level) ) {
      if (verbose)
          std::cerr << "invalid vertex" << std::endl;
      CGAL_triangulation_assertion(false);
      return false;
    }
    ++i;
  }
  return true;
}

template <class Vb, class Cb, class Ct>
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
count_facets(size_type & i, bool verbose, int level) const
  // counts but does not check
{
  i = 0;

  for (Facet_iterator it = facets_begin(); it != facets_end(); ++it) {
    if ( ! is_valid((*it).first,verbose, level) ) {
      if (verbose)
          std::cerr << "invalid facet" << std::endl;
      CGAL_triangulation_assertion(false);
      return false;
    }
    ++i;
  }
  return true;
}

template <class Vb, class Cb, class Ct>
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
count_edges(size_type & i, bool verbose, int level) const
  // counts but does not check
{
  i = 0;

  for (Edge_iterator it = edges_begin(); it != edges_end(); ++it) {
    if ( ! is_valid((*it).first,verbose, level) ) {
      if (verbose)
          std::cerr << "invalid edge" << std::endl;
      CGAL_triangulation_assertion(false);
      return false;
    }
    ++i;
  }
  return true;
}

template <class Vb, class Cb, class Ct>
bool
Triangulation_data_structure_3<Vb,Cb,Ct>::
count_cells(size_type & i, bool verbose, int level) const
  // counts AND checks the validity
{
  i = 0;

  for (Cell_iterator it = cells_begin(); it != cells_end(); ++it) {
    if ( ! is_valid(it,verbose, level) ) {
      if (verbose)
          std::cerr << "invalid cell" << std::endl;
      CGAL_triangulation_assertion(false);
      return false;
    }
    ++i;
  }
  return true;
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_TRIANGULATION_DATA_STRUCTURE_3_H
