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
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Sylvain Pion

// combinatorial triangulation of the boundary of a polytope
// of dimension d in dimension d+1
// for -1 <= d <= 3

#ifndef CGAL_TRIANGULATION_DATA_STRUCTURE_3_H
#define CGAL_TRIANGULATION_DATA_STRUCTURE_3_H

#include <CGAL/basic.h>

#include <utility>
#include <map>
#include <set>
#include <vector>
#include <stack>

#include <CGAL/utility.h>
#include <CGAL/iterator.h>

#include <CGAL/Unique_hash_map.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_utils_3.h>

#include <CGAL/Compact_container.h>

#include <CGAL/Triangulation_ds_cell_base_3.h>
#include <CGAL/Triangulation_ds_vertex_base_3.h>
#include <CGAL/Triangulation_simplex_3.h>

#include <CGAL/internal/Triangulation_ds_iterators_3.h>
#include <CGAL/internal/Triangulation_ds_circulators_3.h>

#ifdef CGAL_HAS_THREADS
#  include <boost/thread/tss.hpp>
#endif


namespace CGAL {

// TODO : noms : Vb != Vertex_base : clarifier.

template < class Vb = Triangulation_ds_vertex_base_3<>,
           class Cb = Triangulation_ds_cell_base_3<> >
class Triangulation_data_structure_3
  : public Triangulation_utils_3
{
  typedef Triangulation_data_structure_3<Vb,Cb>         Tds;

public:

  // Tools to change the Vertex and Cell types of the TDS.
  template < typename Vb2 >
  struct Rebind_vertex {
    typedef Triangulation_data_structure_3<Vb2, Cb>  Other;
  };

  template < typename Cb2 >
  struct Rebind_cell {
    typedef Triangulation_data_structure_3<Vb, Cb2>  Other;
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

  typedef Compact_container<Cell>                  Cell_range;
  typedef Compact_container<Vertex>                Vertex_range;

  typedef typename Cell_range::size_type       size_type;
  typedef typename Cell_range::difference_type difference_type;

  typedef typename Cell_range::iterator        Cell_iterator;
  typedef typename Vertex_range::iterator      Vertex_iterator;

  typedef internal::Triangulation_ds_facet_iterator_3<Tds>   Facet_iterator;
  typedef internal::Triangulation_ds_edge_iterator_3<Tds>    Edge_iterator;

  typedef internal::Triangulation_ds_cell_circulator_3<Tds>  Cell_circulator;
  typedef internal::Triangulation_ds_facet_circulator_3<Tds> Facet_circulator;

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
      Cell_handle newf = create_face(f0->vertex(cw(i0)),
			             f1->vertex(cw(i1)),
			             f2->vertex(cw(i2)));
      set_adjacency(newf, 2, f0, i0);
      set_adjacency(newf, 0, f1, i1);
      set_adjacency(newf, 1, f2, i2);
      return newf;
    }

  Cell_handle create_face(Cell_handle f0, int i0,
	                  Cell_handle f1, int i1)
    {
      CGAL_triangulation_precondition(dimension() <= 2);
      Cell_handle newf = create_face(f0->vertex(cw(i0)),
			             f1->vertex(cw(i1)),
			             f1->vertex(ccw(i1)));
      set_adjacency(newf, 2, f0, i0);
      set_adjacency(newf, 0, f1, i1);
      return newf;
    }

  Cell_handle create_face(Cell_handle f, int i, Vertex_handle v)
    {
      CGAL_triangulation_precondition(dimension() <= 2);
      Cell_handle newf = create_face(f->vertex(cw(i)),
			             f->vertex(ccw(i)),
				     v);
      set_adjacency(newf, 2, f, i);
      return newf;
    }

  // not documented
  void read_cells(std::istream& is, std::map< std::size_t, Vertex_handle > &V,
                  std::size_t & m, std::map< std::size_t, Cell_handle > &C );
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
      // if begin == NULL (default arg), we could compute one by walking in
      // CellIt.  At the moment, the functionality is not available, you have
      // to specify a starting facet.

      Cell_handle cnew;
      if (dimension() == 3)
	  cnew = create_star_3(newv, begin, i);
      else
	  cnew = create_star_2(newv, begin, i);

      newv->set_cell(cnew);
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
	  (*cit)->tds_data().mark_in_conflict();

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

  Vertex_iterator vertices_begin() const
  {
    return vertices().begin();
  }

  Vertex_iterator vertices_end() const
  {
    return vertices().end();
  }

  // CIRCULATOR METHODS

  // cells around an edge
  Cell_circulator incident_cells(const Edge & e) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    return Cell_circulator(e);
  }
  Cell_circulator incident_cells(Cell_handle ce, int i, int j) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    return Cell_circulator(ce, i, j);
  }

  Cell_circulator incident_cells(const Edge &e, Cell_handle start) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    return Cell_circulator(e, start);
  }
  Cell_circulator incident_cells(Cell_handle ce, int i, int j,
	                         Cell_handle start) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    return Cell_circulator(ce, i, j, start);
  }

  //facets around an edge
  Facet_circulator incident_facets(const Edge & e) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    return Facet_circulator(e);
  }
  Facet_circulator incident_facets(Cell_handle ce, int i, int j) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    return Facet_circulator(ce, i, j);
  }
  Facet_circulator incident_facets(const Edge & e, const Facet & start) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    return Facet_circulator(e, start);
  }
  Facet_circulator incident_facets(Cell_handle ce, int i, int j,
				   const Facet & start) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    return Facet_circulator(ce, i, j, start);
  }
  Facet_circulator incident_facets(const Edge & e,
	                           Cell_handle start, int f) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    return Facet_circulator(e, start, f);
  }
  Facet_circulator incident_facets(Cell_handle ce, int i, int j,
				   Cell_handle start, int f) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    return Facet_circulator(ce, i, j, start, f);
  }

  // 2D : circulates on the faces adjacent to a vertex.
  Face_circulator incident_faces(Vertex_handle v) const
  {
    CGAL_triangulation_precondition( dimension() == 2 );
    return Face_circulator(v, v->cell());
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
	d->tds_data().mark_in_conflict();
	*it.first++ = d;
	
	do {
		Cell_handle c = cell_stack.top();
		cell_stack.pop();
		
		for (int i=0; i<4; ++i) {
			if (c->vertex(i) == v)
				continue;
			Cell_handle next = c->neighbor(i);
			if (c < next)
				*it.second++ = Facet(c, i); // Incident facet.
			if (! next->tds_data().is_clear())
				continue;			
			cell_stack.push(next);
			next->tds_data().mark_in_conflict();
			*it.first++ = next;
		}
	} while(!cell_stack.empty());

	return it;
  }


  template <class OutputIterator>
  void
  incident_cells_2(Vertex_handle v, Cell_handle c,
	           OutputIterator cells) const
  {
      CGAL_triangulation_precondition(dimension() == 2);

      // TODO : in 2D, there's no real need for tds_data, we could use
      // a smarter algorithm.  We could use the 2D Face_circulator.
      // Should we just have this Face_circulator ?

      // Flag values :
      // 1 : incident cell already visited
      // 0 : unknown
      c->tds_data().mark_in_conflict();
      *cells++ = c;

      for (int i=0; i<3; ++i) {
	  if (c->vertex(i) == v)
	      continue;
	  Cell_handle next = c->neighbor(i);
	  if (! next->tds_data().is_clear())
	      continue;
	  incident_cells_2(v, next, cells);
      }
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


	// Visitor for visit_incident_cells:
	// outputs the result of Treatment applied to the vertices
  template<class Treatment, class OutputIterator, class Filter>
  class Vertex_extractor {
    Vertex_handle v;
    std::set<Vertex_handle> tmp_vertices;
    Treatment treat;
    const Tds* t;
    Filter filter;
  public:
    Vertex_extractor(Vertex_handle _v, OutputIterator _output, const Tds* _t, Filter _filter):
    v(_v), treat(_output), t(_t), filter(_filter) {}
    void operator()(Cell_handle c) {
      for (int j=0; j<= t->dimension(); ++j) {
	Vertex_handle w = c->vertex(j);
	if(filter(w))
	  continue;
	if (w != v)
	  if(tmp_vertices.insert(w).second) {
	    treat(c, v, j);
	  }
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
  public:
    Vertex_feeder_treatment(OutputIterator _output): output(_output) {};
    void operator()(Cell_handle c, Vertex_handle, int index) {
      *output++ = c->vertex(index);
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
  public:
    Edge_feeder_treatment(OutputIterator _output): output(_output) {};
    void operator()(Cell_handle c, Vertex_handle v, int index) {
      *output++ = Edge(c, c->index(v), index);
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
  incident_edges(Vertex_handle v, OutputIterator edges, Filter f = Filter()) const
  {
    CGAL_triangulation_precondition( v != Vertex_handle() );
    CGAL_triangulation_precondition( dimension() >= 1 );
    CGAL_triangulation_expensive_precondition( is_vertex(v) );
    CGAL_triangulation_expensive_precondition( is_valid() );

    if (dimension() == 1) {
      CGAL_triangulation_assertion( number_of_vertices() >= 3);
      Cell_handle n0 = v->cell();
      const int index_v_in_n0 = n0->index(v);
      CGAL_assume(index_v_in_n0 <= 1);
      Cell_handle n1 = n0->neighbor(1-index_v_in_n0);
      const int index_v_in_n1 = n1->index(v);
      CGAL_assume(index_v_in_n1 <= 1);
      if(!f(n0->vertex(1-index_v_in_n0))) *edges++ = Edge(n0, n0->index(v), 1-index_v_in_n0);
      if(!f(n1->vertex(1-index_v_in_n1))) *edges++ = Edge(n1, n1->index(v), 1-index_v_in_n1);
      return edges;
    }
    return visit_incident_cells<Vertex_extractor<Edge_feeder_treatment<OutputIterator>,
    OutputIterator, Filter>,
    OutputIterator>(v, edges, f);
  }

  template <class OutputIterator>
  OutputIterator
  incident_edges(Vertex_handle v, OutputIterator edges) const
  {
    return incident_edges<False_filter>(v, edges);
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
      *vertices++ = v->cell()->neighbor(0)->vertex(0);
      return vertices;
    }

    if (dimension() == 1) {
      CGAL_triangulation_assertion( number_of_vertices() >= 3);
      Cell_handle n0 = v->cell();
      const int index_v_in_n0 = n0->index(v);
      CGAL_assume(index_v_in_n0 <= 1);
      Cell_handle n1 = n0->neighbor(1-index_v_in_n0);
      const int index_v_in_n1 = n1->index(v);
      CGAL_assume(index_v_in_n1 <= 1);
      Vertex_handle v1 = n0->vertex(1-index_v_in_n0);
      Vertex_handle v2 = n1->vertex(1-index_v_in_n1);
      if(!f(v1)) *vertices++ = v1;
      if(!f(v2)) *vertices++ = v2;
      return vertices;
    }
    return visit_incident_cells<Vertex_extractor<Vertex_feeder_treatment<OutputIterator>,
    OutputIterator, Filter>,
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
    incident_cells_3(v, v->cell(), std::make_pair(std::back_inserter(tmp_cells), visit.facet_it()));
    else
    incident_cells_2(v, v->cell(), std::back_inserter(tmp_cells));

    typename std::vector<Cell_handle>::iterator cit;
    for(cit = tmp_cells.begin();
	cit != tmp_cells.end();
	++cit)
    {
      (*cit)->tds_data().clear();
      visit(*cit);
    }
    return visit.result();
  }

  size_type degree(Vertex_handle v) const;

  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const;
  bool is_valid(Vertex_handle v, bool verbose = false, int level = 0) const;
  bool is_valid(Cell_handle c, bool verbose = false, int level = 0) const;

  // Helping functions
  Vertex_handle copy_tds(const Tds & tds,
	                 Vertex_handle vert = Vertex_handle() );
    // returns the new vertex corresponding to vert in the new tds

  void swap(Tds & tds);

  void clear();

  void set_adjacency(Cell_handle c0, int i0,
	             Cell_handle c1, int i1) const
  {
      CGAL_triangulation_assertion(i0 >= 0 && i0 <= dimension());
      CGAL_triangulation_assertion(i1 >= 0 && i1 <= dimension());
      CGAL_triangulation_assertion(c0 != c1);
      c0->set_neighbor(i0,c1);
      c1->set_neighbor(i1,c0);
  }

  int mirror_index(Cell_handle c, int i) const
  {
      CGAL_triangulation_precondition ( i>=0 && i<4 );
      return c->neighbor(i)->index(c);
  }

  Vertex_handle mirror_vertex(Cell_handle c, int i) const
  {
      return c->neighbor(i)->vertex(mirror_index(c, i));
  }

  Facet mirror_facet(Facet f) const
  {
    Cell_handle neighbor_cell = f.first->neighbor(f.second);
    const int opposite_index = neighbor_cell->index(f.first);
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
  void change_orientation(Cell_handle c) const
  {
      Vertex_handle tmp_v = c->vertex(0);
      c->set_vertex(0, c->vertex(1));
      c->set_vertex(1, tmp_v);
      Cell_handle tmp_c = c->neighbor(0);
      c->set_neighbor(0, c->neighbor(1));
      c->set_neighbor(1, tmp_c);
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
template <class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Cell_handle
Triangulation_data_structure_3<Vb,Cb>::
create_star_3(Vertex_handle v, Cell_handle c, int li, int prev_ind2)
{
    CGAL_triangulation_precondition( dimension() == 3);
    CGAL_triangulation_precondition( c->tds_data().is_in_conflict() );
    CGAL_triangulation_precondition( ! c->neighbor(li)->tds_data().is_in_conflict() );

    Cell_handle cnew = create_cell(c->vertex(0),
				   c->vertex(1),
				   c->vertex(2),
				   c->vertex(3));
    cnew->set_vertex(li, v);
    Cell_handle c_li = c->neighbor(li);
    set_adjacency(cnew, li, c_li, c_li->index(c));

    // Look for the other neighbors of cnew.
    for (int ii=0; ii<4; ++ii) {
      if (ii == prev_ind2 || cnew->neighbor(ii) != Cell_handle())
	  continue;
      cnew->vertex(ii)->set_cell(cnew);

      // Indices of the vertices of cnew such that ii,vj1,vj2,li positive.
      Vertex_handle vj1 = c->vertex(next_around_edge(ii, li));
      Vertex_handle vj2 = c->vertex(next_around_edge(li, ii));
      Cell_handle cur = c;
      int zz = ii;
      Cell_handle n = cur->neighbor(zz);
      // turn around the oriented edge vj1 vj2
      while ( n->tds_data().is_in_conflict() ) {
	CGAL_triangulation_assertion( n != c );
	cur = n;
	zz = next_around_edge(n->index(vj1), n->index(vj2));
	n = cur->neighbor(zz);
      }
      // Now n is outside region, cur is inside.

      n->tds_data().clear(); // Reset the flag for boundary cells.

      int jj1 = n->index(vj1);
      int jj2 = n->index(vj2);
      Vertex_handle vvv = n->vertex(next_around_edge(jj1, jj2));
      Cell_handle nnn = n->neighbor(next_around_edge(jj2, jj1));
      int zzz = nnn->index(vvv);
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
template <class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Cell_handle
Triangulation_data_structure_3<Vb,Cb>::
recursive_create_star_3(Vertex_handle v, Cell_handle c, int li,
                        int prev_ind2, int depth)
{
    if ( depth == 100 ) return non_recursive_create_star_3(v,c,li,prev_ind2);
    CGAL_triangulation_precondition( dimension() == 3);
    CGAL_triangulation_precondition( c->tds_data().is_in_conflict() );
    CGAL_triangulation_precondition( ! c->neighbor(li)->tds_data().is_in_conflict() );

    Cell_handle cnew = create_cell(c->vertex(0),
				   c->vertex(1),
				   c->vertex(2),
				   c->vertex(3));
    cnew->set_vertex(li, v);
    Cell_handle c_li = c->neighbor(li);
    set_adjacency(cnew, li, c_li, c_li->index(c));

    // Look for the other neighbors of cnew.
    for (int ii=0; ii<4; ++ii) {
      if (ii == prev_ind2 || cnew->neighbor(ii) != Cell_handle())
	  continue;
      cnew->vertex(ii)->set_cell(cnew);

      // Indices of the vertices of cnew such that ii,vj1,vj2,li positive.
      Vertex_handle vj1 = c->vertex(next_around_edge(ii, li));
      Vertex_handle vj2 = c->vertex(next_around_edge(li, ii));
      Cell_handle cur = c;
      int zz = ii;
      Cell_handle n = cur->neighbor(zz);
      // turn around the oriented edge vj1 vj2
      while ( n->tds_data().is_in_conflict() ) {
	CGAL_triangulation_assertion( n != c );
	cur = n;
	zz = next_around_edge(n->index(vj1), n->index(vj2));
	n = cur->neighbor(zz);
      }
      // Now n is outside region, cur is inside.

      n->tds_data().clear(); // Reset the flag for boundary cells.

      int jj1 = n->index(vj1);
      int jj2 = n->index(vj2);
      Vertex_handle vvv = n->vertex(next_around_edge(jj1, jj2));
      Cell_handle nnn = n->neighbor(next_around_edge(jj2, jj1));
      int zzz = nnn->index(vvv);
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
template <class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Cell_handle
Triangulation_data_structure_3<Vb,Cb>::
non_recursive_create_star_3(Vertex_handle v, Cell_handle c, int li, int prev_ind2)
{
    CGAL_triangulation_precondition( dimension() == 3);
    CGAL_triangulation_precondition( c->tds_data().is_in_conflict() );
    CGAL_triangulation_precondition( ! c->neighbor(li)->tds_data().is_in_conflict() );

    Cell_handle cnew = create_cell(c->vertex(0),
				   c->vertex(1),
				   c->vertex(2),
				   c->vertex(3));
    cnew->set_vertex(li, v);
    Cell_handle c_li = c->neighbor(li);
    set_adjacency(cnew, li, c_li, c_li->index(c));
    
    std::stack<iAdjacency_info> adjacency_info_stack;
  
    int ii=0;
    do
    {
      // Look for the other neighbors of cnew.
      if ( ! (ii == prev_ind2 || cnew->neighbor(ii) != Cell_handle()) ){
        cnew->vertex(ii)->set_cell(cnew);

        // Indices of the vertices of cnew such that ii,vj1,vj2,li positive.
        Vertex_handle vj1 = c->vertex(next_around_edge(ii, li));
        Vertex_handle vj2 = c->vertex(next_around_edge(li, ii));
        Cell_handle cur = c;
        int zz = ii;
        Cell_handle n = cur->neighbor(zz);
        // turn around the oriented edge vj1 vj2
        while ( n->tds_data().is_in_conflict() ) {
          CGAL_triangulation_assertion( n != c );
          cur = n;
          zz = next_around_edge(n->index(vj1), n->index(vj2));
          n = cur->neighbor(zz);
        }
        // Now n is outside region, cur is inside.

        n->tds_data().clear(); // Reset the flag for boundary cells.

        int jj1 = n->index(vj1);
        int jj2 = n->index(vj2);
        Vertex_handle vvv = n->vertex(next_around_edge(jj1, jj2));
        Cell_handle nnn = n->neighbor(next_around_edge(jj2, jj1));
        int zzz = nnn->index(vvv);
        if (nnn == cur) {
          // Neighbor relation is reciprocal, ie
          // the cell we are looking for is not yet created.
          //re-run the loop
          adjacency_info_stack.push( iAdjacency_info(zzz,cnew,ii,c,li,prev_ind2) );
          c=nnn;  li=zz; prev_ind2=zzz; ii=0;
          //copy-pasted from beginning to avoid if ii==0
          CGAL_triangulation_precondition( c->tds_data().is_in_conflict() );
          CGAL_triangulation_precondition( ! c->neighbor(li)->tds_data().is_in_conflict() );
          cnew = create_cell(c->vertex(0),c->vertex(1),c->vertex(2),c->vertex(3));
          cnew->set_vertex(li, v);
          c_li = c->neighbor(li);
          set_adjacency(cnew, li, c_li, c_li->index(c));          
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

template <class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Cell_handle
Triangulation_data_structure_3<Vb,Cb>::
create_star_2(Vertex_handle v, Cell_handle c, int li )
{
  CGAL_triangulation_assertion( dimension() == 2 );
  Cell_handle cnew;

  // i1 i2 such that v,i1,i2 positive
  int i1=ccw(li);
  // traversal of the boundary of region in ccw order to create all
  // the new facets
  Cell_handle bound = c;
  Vertex_handle v1 = c->vertex(i1);
  int ind = c->neighbor(li)->index(c); // to be able to find the
                                       // first cell that will be created
  Cell_handle cur;
  Cell_handle pnew = Cell_handle();
  do {
    cur = bound;
    // turn around v2 until we reach the boundary of region
    while ( cur->neighbor(cw(i1))->tds_data().is_in_conflict() ) {
      // neighbor in conflict
      cur = cur->neighbor(cw(i1));
      i1 = cur->index( v1 );
    }
    cur->neighbor(cw(i1))->tds_data().clear();
    // here cur has an edge on the boundary of region
    cnew = create_face( v, v1, cur->vertex( ccw(i1) ) );
    set_adjacency(cnew, 0, cur->neighbor(cw(i1)),
	                   cur->neighbor(cw(i1))->index(cur));
    cnew->set_neighbor(1, Cell_handle());
    cnew->set_neighbor(2, pnew);
    // pnew is null at the first iteration
    v1->set_cell(cnew);
    //pnew->set_neighbor( cw(pnew->index(v1)), cnew );
    if (pnew != Cell_handle()) { pnew->set_neighbor( 1, cnew );}

    bound = cur;
    i1 = ccw(i1);
    v1 = bound->vertex(i1);
    pnew = cnew;
    //} while ( ( bound != c ) || ( li != cw(i1) ) );
  } while ( v1 != c->vertex(ccw(li)) );
  // missing neighbors between the first and the last created cells
  cur = c->neighbor(li)->neighbor(ind); // first created cell
  set_adjacency(cnew, 1, cur, 2);
  return cnew;
}

template < class Vb, class Cb>
std::istream&
operator>>(std::istream& is, Triangulation_data_structure_3<Vb,Cb>& tds)
  // reads :
  // the dimension
  // the number of vertices
  // the number of cells
  // the cells by the indices of their vertices
  // the neighbors of each cell by their index in the preceding list of cells
  // when dimension < 3 : the same with faces of maximal dimension
{
  typedef Triangulation_data_structure_3<Vb,Cb> Tds;
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

  std::map<std::size_t , Vertex_handle > V;

  // creation of the vertices
  for (std::size_t i=0; i < n; i++) {
    //    is >> p;
    //    V[i] = tds.create_vertex();
    //    V[i]->set_point(p);
    V[i] = tds.create_vertex();
  }

  std::map< std::size_t, Cell_handle > C;
  std::size_t m;

  tds.read_cells(is, V, m, C);
  CGAL_triangulation_assertion( tds.is_valid() );
  return is;
}

template < class Vb, class Cb>
std::ostream&
operator<<(std::ostream& os, const Triangulation_data_structure_3<Vb,Cb> &tds)
  // writes :
  // the dimension
  // the number of vertices
  // the number of cells
  // the cells by the indices of their vertices
  // the neighbors of each cell by their index in the preceding list of cells
  // when dimension < 3 : the same with faces of maximal dimension
{
  typedef Triangulation_data_structure_3<Vb,Cb> Tds;
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

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
is_simplex( Cell_handle c ) const
{
  switch(dimension()) {
    case 3 : return is_cell(c);
    case 2 : return is_facet(c, 3);
    case 1 : return is_edge(c, 0, 1);
    case 0 : return is_vertex(c->vertex(0));
    case -1 : return c == cells().begin();
  }
  return false;
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
is_vertex(Vertex_handle v) const
{
    return vertices().owns_dereferencable(v);
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
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
        if ((*cit)->has_vertex(v, j)) {
	    c = *cit;
	    i = c->index(u);
	    return true;
	}
    return false;
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
is_edge(Vertex_handle u, Vertex_handle v) const
{
    Cell_handle c;
    int i, j;
    return is_edge(u, v, c, i, j);
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
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

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
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
        if ((*cit)->has_vertex(v, j) && (*cit)->has_vertex(w, k)) {
	    c = *cit;
	    i = c->index(u);
	    return true;
	}
    return false;
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
is_facet(Cell_handle c, int i) const
{
    CGAL_triangulation_precondition(i>=0 && i<4);

    if ( dimension() < 2 )
        return false;

    if ( (dimension() == 2) && (i!=3) )
        return false;

    return cells().owns_dereferencable(c);
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
is_cell( Cell_handle c ) const
  // returns false when dimension <3
{
    if (dimension() < 3)
        return false;

    return cells().owns_dereferencable(c);
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
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
        if ((*cit)->has_vertex(v, j) && (*cit)->has_vertex(w, k) &&
	    (*cit)->has_vertex(t, l)) {
	    c = *cit;
	    i = c->index(u);
	    return true;
	}
    return false;
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
is_cell(Vertex_handle u, Vertex_handle v,
	Vertex_handle w, Vertex_handle t)
    const
  // returns false when dimension <3
{
    Cell_handle c;
    int i, j, k, l;
    return is_cell(u, v, w, t, c, i, j, k, l);
}

template < class Vb, class Cb>
inline
bool
Triangulation_data_structure_3<Vb,Cb>::
has_vertex(Cell_handle c, int i, Vertex_handle v, int & j) const
  // computes the index j of the vertex in the cell c giving the query
  // facet (c,i)
  // j has no meaning if false is returned
{
  CGAL_triangulation_precondition( dimension() == 3 );
  return ( c->has_vertex(v,j) && (j != i) );
}

template < class Vb, class Cb>
inline
bool
Triangulation_data_structure_3<Vb,Cb>::
has_vertex(Cell_handle c, int i, Vertex_handle v) const
  // checks whether the query facet (c,i) has vertex v
{
  CGAL_triangulation_precondition( dimension() == 3 );
  int j;
  return ( c->has_vertex(v,j) && (j != i) );
}

template < class Vb, class Cb>
inline
bool
Triangulation_data_structure_3<Vb,Cb>::
has_vertex(const Facet & f, Vertex_handle v, int & j) const
{
  return has_vertex(f.first, f.second, v, j);
}

template < class Vb, class Cb>
inline
bool
Triangulation_data_structure_3<Vb,Cb>::
has_vertex(const Facet & f, Vertex_handle v) const
{
  return has_vertex(f.first, f.second, v);
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
are_equal(Cell_handle c, int i, Cell_handle n, int j) const
  // tests whether facets c,i and n,j, have the same 3 vertices
  // the triangulation is supposed to be valid, the orientation of the
  // facets is not checked here
  // the neighbor relations between c and  n are not tested either,
  // which allows to use this method before setting these relations
  // (see remove in Delaunay_3)
  //   if ( c->neighbor(i) != n ) return false;
  //   if ( n->neighbor(j) != c ) return false;
{
  CGAL_triangulation_precondition( dimension() == 3 );

  if ( (c==n) && (i==j) ) return true;

  int j1,j2,j3;
  return( n->has_vertex( c->vertex((i+1)&3), j1 ) &&
	  n->has_vertex( c->vertex((i+2)&3), j2 ) &&
	  n->has_vertex( c->vertex((i+3)&3), j3 ) &&
	  ( j1+j2+j3+j == 6 ) );
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
are_equal(const Facet & f, const Facet & g) const
{
  return are_equal(f.first, f.second, g.first, g.second);
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
are_equal(const Facet & f, Cell_handle n, int j) const
{
  return are_equal(f.first, f.second, n, j);
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
flip( Cell_handle c, int i )
  // returns false if the facet is not flippable
  // true other wise and
  // flips facet i of cell c
  // c will be replaced by one of the new cells
{
  CGAL_triangulation_precondition( (dimension() == 3) && (0<=i) && (i<4)
				   && (number_of_vertices() >= 6) );
  CGAL_triangulation_expensive_precondition( is_cell(c) );

  Cell_handle n = c->neighbor(i);
  int in = n->index(c);

  // checks that the facet is flippable,
  // ie the future edge does not already exist
  if (is_edge(c->vertex(i), n->vertex(in)))
      return false;

  flip_really(c,i,n,in);
  return true;
}

template < class Vb, class Cb>
void
Triangulation_data_structure_3<Vb,Cb>::
flip_flippable(Cell_handle c, int i )
  // flips facet i of cell c
  // c will be replaced by one of the new cells
{
  CGAL_triangulation_precondition( (dimension() == 3) && (0<=i) && (i<4)
				   && (number_of_vertices() >= 6) );
  CGAL_triangulation_expensive_precondition( is_cell(c) );

  Cell_handle n = c->neighbor(i);
  int in = n->index(c);

  // checks that the facet is flippable,
  // ie the future edge does not already exist
  CGAL_triangulation_expensive_precondition( !is_edge(c->vertex(i),
	                                              n->vertex(in)));
  flip_really(c,i,n,in);
}

template < class Vb, class Cb>
inline
void
Triangulation_data_structure_3<Vb,Cb>::
flip_really( Cell_handle c, int i, Cell_handle n, int in )
  // private - used by flip and flip_flippable
{
  int i1 = (i+1)&3;
  int i2 = (i+2)&3;
  int i3 = (i+3)&3;

  int in1 = n->index(c->vertex(i1));
  int in2 = n->index(c->vertex(i2));
  int in3 = n->index(c->vertex(i3));

  set_adjacency(c, i, n->neighbor(in3), n->neighbor(in3)->index(n));
  c->set_vertex( i3, n->vertex(in) );

  set_adjacency(n, in, c->neighbor(i1), c->neighbor(i1)->index(c));
  n->set_vertex( in1, c->vertex(i) );

  Cell_handle cnew = create_cell(c->vertex(i), c->vertex(i1),
			         n->vertex(in), n->vertex(in3));

  set_adjacency(cnew, 0, n->neighbor(in2), n->neighbor(in2)->index(n));
  set_adjacency(cnew, 1, n, in2);
  set_adjacency(cnew, 2, c->neighbor(i2), c->neighbor(i2)->index(c));
  set_adjacency(cnew, 3, c, i2);
  set_adjacency(c, i1, n, in3);

  if ((i&1) != 0)
      change_orientation(cnew);

  c->vertex(i1)->set_cell(cnew);
  c->vertex(i2)->set_cell(c);
  n->vertex(in3)->set_cell(n);
  // to be implemented : 2d case
  // CGAL_triangulation_precondition( (0<=i) && (i<3) );
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
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
  Cell_handle c1 = c->neighbor( next );
  Vertex_handle v1 = c->vertex( next ); // will become vertex of c1
  int i1 = c1->index( c->vertex(i) );
  int j1 = c1->index( c->vertex(j) );

  int next1 = next_around_edge(i1,j1);
  Cell_handle c2  = c1->neighbor( next1 );
  Vertex_handle v2 = c1->vertex( next1 ); // will become vertex of c2
  int i2 = c2->index( c->vertex(i) );
  int j2 = c2->index( c->vertex(j) );

  int next2 = next_around_edge(i2,j2);
  Vertex_handle v3 = c2->vertex( next2 );

  // checks that the edge is flippable,
  // is the future cells do not already exist
  if ( is_cell(v1,v2,v3,c->vertex(i)) ) return false;
  if ( is_cell(v1,v2,v3,c->vertex(j)) ) return false;

  flip_really(c,i,j,c1,v1,i1,j1,next1,c2,v2,i2,j2,next2,v3);

  return true;
}

template < class Vb, class Cb>
void
Triangulation_data_structure_3<Vb,Cb>::
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
  Cell_handle c1 = c->neighbor( next );
  Vertex_handle v1 = c->vertex( next ); // will become vertex of c1
  int i1 = c1->index( c->vertex(i) );
  int j1 = c1->index( c->vertex(j) );

  int next1 = next_around_edge(i1,j1);
  Cell_handle c2  = c1->neighbor( next1 );
  Vertex_handle v2 = c1->vertex( next1 ); // will become vertex of c2
  int i2 = c2->index( c->vertex(i) );
  int j2 = c2->index( c->vertex(j) );

  int next2 = next_around_edge(i2,j2);
  Vertex_handle v3 = c2->vertex( next2 );

  // checks that the edge is flippable,
  // is the future cells do not already exist
  CGAL_triangulation_expensive_precondition( !is_cell(v1,v2,v3,c->vertex(i)) );
  CGAL_triangulation_expensive_precondition( !is_cell(v1,v2,v3,c->vertex(j)) );

  flip_really(c,i,j,c1,v1,i1,j1,next1,c2,v2,i2,j2,next2,v3);
}

template < class Vb, class Cb>
inline
void
Triangulation_data_structure_3<Vb,Cb>::
flip_really( Cell_handle c, int i, int j,
	     Cell_handle c1, Vertex_handle v1,
	     int i1, int j1, int next1,
	     Cell_handle c2, Vertex_handle v2,
	     int i2, int j2, int next2,
	     Vertex_handle v3 )
{
  c->vertex(i)->set_cell(c1);
  c->vertex(j)->set_cell(c2);

  c1->set_vertex( j1, v1 );
  v1->set_cell(c1);
  c2->set_vertex( i2, v2 );
  v2->set_cell(c2);

  set_adjacency(c1, next1,c2->neighbor(j2), c2->neighbor(j2)->index(c2));
  set_adjacency(c2,c2->index(v1),c1->neighbor(i1),c1->neighbor(i1)->index(c1));

  set_adjacency(c1, i1, c2, j2);

  set_adjacency(c1, 6-i1-j1-next1, c->neighbor(j), c->neighbor(j)->index(c));
  set_adjacency(c2, next2, c->neighbor(i), c->neighbor(i)->index(c));

  v3->set_cell( c2 );

  delete_cell( c );
}

template < class Vb, class Cb >
void
Triangulation_data_structure_3<Vb,Cb>::
read_cells(std::istream& is, std::map< std::size_t, Vertex_handle > &V,
	   std::size_t & m, std::map< std::size_t, Cell_handle > &C)
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

      for(std::size_t i = 0; i < m; i++) {
	Cell_handle c = create_cell();
	for (int k=0; k<=dimension(); ++k) {
          std::size_t ik;
            if(is_ascii(is))
               is >> ik;
            else
              read(is, ik);
	    c->set_vertex(k, V[ik]);
	    V[ik]->set_cell(c);
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
	    c->set_neighbor(k, C[ik]);
	}
      }
      break;
    }
  case 0:
    {
      m = 2;

      //      CGAL_triangulation_assertion( n == 2 );
      for (int i=0; i < 2; i++) {
	Cell_handle c = create_face(V[i], Vertex_handle(), Vertex_handle());
	C[i] = c;
	V[i]->set_cell(c);
      }
      for (int j=0; j < 2; j++) {
	Cell_handle c = C[j];
        c->set_neighbor(0, C[1-j]);
      }
      break;
    }
  case -1:
    {
      m = 1;
      //      CGAL_triangulation_assertion( n == 1 );
      Cell_handle c = create_face(V[0], Vertex_handle(), Vertex_handle());
      C[0] = c;
      V[0]->set_cell(c);
      break;
    }
  }
}

template < class Vb, class Cb>
void
Triangulation_data_structure_3<Vb,Cb>::
print_cells(std::ostream& os, const Unique_hash_map<Vertex_handle, std::size_t> &V ) const
{
  std::map<Cell_handle, std::size_t > C;
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
            os << V[it->vertex(j)];
	    if ( j==3 )
	      os << std::endl;
	    else
	      os << ' ';
	  }
          else
            write(os, V[it->vertex(j)]);
	}
      }
      CGAL_triangulation_assertion( i == m );

      // write the neighbors
      for(it = cells_begin(); it != cells_end(); ++it) {
	for (int j = 0; j < 4; j++) {
	  if(is_ascii(os)){
            os << C[it->neighbor(j)];
	    if(j==3)
	      os << std::endl;
	    else
	      os <<  ' ';
          }
          else
            write(os, C[it->neighbor(j)]);
	}
      }
      break;
    }
  case 2:
    {
      size_type m = number_of_facets();
      if(is_ascii(os))
        os << m << std::endl;
      else
        write(os, m);

      // write the facets
      Facet_iterator it;
      for(it = facets_begin(); it != facets_end(); ++it) {
	C[(*it).first] = i++;
	for(int j = 0; j < 3; j++){
	  if(is_ascii(os)) {
	    os << V[(*it).first->vertex(j)];
	    if ( j==2 )
	      os << std::endl;
	    else
	      os <<  ' ';
	  }
	  else {
	    write(os,  V[(*it).first->vertex(j)]);
	  }
	}
      }
      CGAL_triangulation_assertion( i == m );

      // write the neighbors
      for(it = facets_begin(); it != facets_end(); ++it) {
	for (int j = 0; j < 3; j++) {
	  if(is_ascii(os)){
	    os << C[(*it).first->neighbor(j)];
	    if(j==2)
	      os << std::endl;
	    else
	      os <<  ' ';
	  }
	  else {
	    write(os, C[(*it).first->neighbor(j)]);
	  }
	}
      }
      break;
    }
  case 1:
    {
      size_type m = number_of_edges();
      if(is_ascii(os))
        os << m << std::endl;
      else
        write(os, m);
      // write the edges
      Edge_iterator it;
      for(it = edges_begin(); it != edges_end(); ++it) {
	C[(*it).first] = i++;
	for(int j = 0; j < 2; j++){
	  if(is_ascii(os)) {
	    os << V[(*it).first->vertex(j)];
	    if ( j==1 )
	      os << std::endl;
	    else
	      os <<  ' ';
	  }
	  else {
	    write(os, V[(*it).first->vertex(j)]);
	  }
	}
      }
      CGAL_triangulation_assertion( i == m );

      // write the neighbors
      for(it = edges_begin(); it != edges_end(); ++it) {
	for (int j = 0; j < 2; j++) {
	  if(is_ascii(os)){
	    os << C[(*it).first->neighbor(j)];
	    if(j==1)
	      os << std::endl;
	    else
	      os <<  ' ';
	  }
	  else {
	    write(os, C[(*it).first->neighbor(j)]);
	  }
	}
      }
      break;
    }
  }
}

template <class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
insert_in_cell(Cell_handle c)
{
  CGAL_triangulation_precondition( dimension() == 3 );
  CGAL_triangulation_precondition( c != Cell_handle() );
  CGAL_triangulation_expensive_precondition( is_cell(c) );

  Vertex_handle v = create_vertex();

  Vertex_handle v0 = c->vertex(0);
  Vertex_handle v1 = c->vertex(1);
  Vertex_handle v2 = c->vertex(2);
  Vertex_handle v3 = c->vertex(3);

  Cell_handle n1 = c->neighbor(1);
  Cell_handle n2 = c->neighbor(2);
  Cell_handle n3 = c->neighbor(3);

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

  set_adjacency(n1, n1->index(c), c1, 1);
  set_adjacency(n2, n2->index(c), c2, 2);
  set_adjacency(n3, n3->index(c), c3, 3);

  c->set_vertex(0,v);

  v0->set_cell(c1);
  v->set_cell(c);

  return v;
}

template <class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
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

      Vertex_handle vi=c->vertex(i);
      Vertex_handle v1=c->vertex(i1);
      Vertex_handle v2=c->vertex(i2);
      Vertex_handle v3=c->vertex(i3);

      // new cell with v in place of i1
      Cell_handle nc = c->neighbor(i1);
      Cell_handle cnew1 = create_cell(vi,v,v2,v3);
      set_adjacency(cnew1, 1, nc, nc->index(c));
      set_adjacency(cnew1, 3, c, i1);

      v3->set_cell(cnew1);

      // new cell with v in place of i2
      nc = c->neighbor(i2);
      Cell_handle cnew2 = create_cell(vi,v1,v,v3);
      set_adjacency(cnew2, 2, nc, nc->index(c));
      set_adjacency(cnew2, 3, c, i2);
      set_adjacency(cnew1, 2, cnew2, 1);

      // v replaces i3 in c
      c->set_vertex(i3,v);

      // other side of facet containing v
      Cell_handle d = c->neighbor(i);
      int j = d->index(c);
      int j1=d->index(v1);// triangulation supposed to be valid
      int j2=d->index(v2);
      int j3=6-j-j1-j2;
      // then the orientation of j,j1,j2,j3 depends on the parity
      // of i-j

      // new cell with v in place of j1
      Cell_handle nd = d->neighbor(j1);
      Cell_handle dnew1 = create_cell(d->vertex(j),v,v3,v2);
      set_adjacency(dnew1, 1, nd, nd->index(d));
      set_adjacency(dnew1, 2, d, j1);
      set_adjacency(dnew1, 0, cnew1, 0);

      // new cell with v in place of j2
      nd = d->neighbor(j2);
      Cell_handle dnew2 = create_cell(d->vertex(j),v1,v3,v);

      set_adjacency(dnew2, 3, nd, nd->index(d));
      set_adjacency(dnew2, 2, d, j2);
      set_adjacency(dnew2, 0, cnew2, 0);
      set_adjacency(dnew1, 3, dnew2, 1);

      // v replaces i3 in d
      d->set_vertex(j3,v);
      v->set_cell(d);

      break;
    }
  case 2:
    {
      CGAL_triangulation_expensive_precondition( is_facet(c,i) );
      Cell_handle n = c->neighbor(2);
      Cell_handle cnew = create_face(c->vertex(0),c->vertex(1),v);
      set_adjacency(cnew, 2, n, n->index(c));
      set_adjacency(cnew, 0, c, 2);
      c->vertex(0)->set_cell(cnew);

      n = c->neighbor(1);
      Cell_handle dnew = create_face(c->vertex(0),v,c->vertex(2));
      set_adjacency(dnew, 1, n, n->index(c));
      set_adjacency(dnew, 0, c, 1);
      set_adjacency(dnew, 2, cnew, 1);

      c->set_vertex(0,v);
      v->set_cell(c);
      break;
    }
  }

  return v;
}

template <class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
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
	  cc->tds_data().mark_in_conflict();
	  ++ccir;
      } while (c != ccir);

      return _insert_in_hole(cells.begin(), cells.end(), c, i);
    }
  case 2:
    {
      CGAL_triangulation_expensive_precondition( is_edge(c,i,j) );

      Vertex_handle v = create_vertex();
      int k=3-i-j; // index of the third vertex of the facet
      Cell_handle d = c->neighbor(k);
      int kd = d->index(c);
      int id = d->index(c->vertex(i));
      int jd = d->index(c->vertex(j));

      Cell_handle cnew = create_cell();
      cnew->set_vertex(i,c->vertex(i));
      c->vertex(i)->set_cell(cnew);
      cnew->set_vertex(j,v);
      cnew->set_vertex(k,c->vertex(k));
      c->set_vertex(i,v);

      Cell_handle dnew = create_cell();
      dnew->set_vertex(id,d->vertex(id));
      // d->vertex(id)->cell() is cnew OK
      dnew->set_vertex(jd,v);
      dnew->set_vertex(kd,d->vertex(kd));
      d->set_vertex(id,v);

      Cell_handle nj = c->neighbor(j);
      set_adjacency(cnew, i, c, j);
      set_adjacency(cnew, j, nj, nj->index(c));

      nj = d->neighbor(jd);
      set_adjacency(dnew, id, d, jd);
      set_adjacency(dnew, jd, nj, nj->index(d));

      set_adjacency(cnew, k, dnew, kd);

      v->set_cell(cnew);
      return v;
    }
  default: // case 1:
    {
      Vertex_handle v = create_vertex();
      CGAL_triangulation_expensive_precondition( is_edge(c,i,j) );
      Cell_handle cnew = create_face(v, c->vertex(1), Vertex_handle());
      c->vertex(1)->set_cell(cnew);
      c->set_vertex(1,v);
      set_adjacency(cnew, 0, c->neighbor(0), 1);
      set_adjacency(cnew, 1, c, 0);

      v->set_cell(cnew);
      return v;
    }
  }
}

template <class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
insert_increase_dimension(Vertex_handle star)
  // star = vertex from which we triangulate the facet of the
  // incremented dimension
  // ( geometrically : star = infinite vertex )
  // = NULL only used to insert the 1st vertex (dimension -2 to dimension -1)
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
      v->set_cell(c);
      break;
    }

  case -1:
    // insertion of the second vertex
    // ( geometrically : first finite vertex )
    {
      Cell_handle d = create_face(v, Vertex_handle(), Vertex_handle());
      v->set_cell(d);
      set_adjacency(d, 0, star->cell(), 0);
      break;
    }

  case 0:
    // insertion of the third vertex
    // ( geometrically : second finite vertex )
    {
      Cell_handle c = star->cell();
      Cell_handle d = c->neighbor(0);

      c->set_vertex(1,d->vertex(0));
      d->set_vertex(1,v);
      d->set_neighbor(1,c);
      Cell_handle e = create_face( v, star, Vertex_handle() );
      set_adjacency(e, 0, c, 1);
      set_adjacency(e, 1, d, 0);

      v->set_cell(d);
      break;
    }

  case 1:
    // general case : 4th vertex ( geometrically : 3rd finite vertex )
    // degenerate cases geometrically : 1st non collinear vertex
    {
      Cell_handle c = star->cell();
      int i = c->index(star); // i== 0 or 1
      CGAL_assertion(i==0 || i==1);
      int j = (i == 0) ? 1 : 0;
      Cell_handle d = c->neighbor(j);
	
      c->set_vertex(2,v);

      Cell_handle e = c->neighbor(i);
      Cell_handle cnew = c;
      Cell_handle enew = Cell_handle();
	
      while( e != d ){
	enew = create_cell();
	enew->set_vertex(i,e->vertex(j));
	enew->set_vertex(j,e->vertex(i));
	enew->set_vertex(2,star);
	
	set_adjacency(enew, i, cnew, j);
	// false at the first iteration of the loop where it should
	// be neighbor 2
	// it is corrected after the loop
	set_adjacency(enew, 2, e, 2);
	// neighbor j will be set during next iteration of the loop
	
	e->set_vertex(2,v);

	e = e->neighbor(i);
	cnew = enew;
      }
	
      d->set_vertex(2,v);
      set_adjacency(enew, j, d, 2);

      // corrections for star->cell() :
      c = star->cell();
      c->set_neighbor(2,c->neighbor(i)->neighbor(2));
      c->set_neighbor(j,d);

      v->set_cell(d);

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

      v->set_cell(it); // ok since there is at least one ``cell''
      for(; it != cells_end(); ++it) {
	// Here we must be careful since we create_cells in a loop controlled
	// by an iterator.  So we first take care of the cells newly created
	// by the following test :
	if (it->neighbor(0) == Cell_handle())
	  continue;
	it->set_neighbor(3, Cell_handle());
	it->set_vertex(3, v);
	if ( ! it->has_vertex(star) ) {
	  Cell_handle cnew = create_cell( it->vertex(0), it->vertex(2),
			                  it->vertex(1), star);
	  // The Intel compiler has a problem with passing "it" directly to
	  // function "set_adjacency": the adjacency is not changed.
	  Cell_handle ch_it = it;
	  set_adjacency(cnew, 3, ch_it, 3);
	  cnew->set_neighbor(0, Cell_handle());
	  new_cells.push_back(cnew);
	}
      }

      // traversal of the new cells only, to add missing neighbors
      for(typename std::vector<Cell_handle>::iterator ncit = new_cells.begin();
           ncit != new_cells.end(); ++ncit) {
	Cell_handle n = (*ncit)->neighbor(3); // opposite to star
	for ( int i=0; i<3; i++ ) {
	  int j;
	  if ( i==0 ) j=0;
	  else j=3-i; // vertex 1 and vertex 2 are always switched when
	  // creating a new cell (see above)
          Cell_handle c = n->neighbor(i)->neighbor(3);
	  if ( c != Cell_handle() ) {
	    // i.e. star is not a vertex of n->neighbor(i)
	    (*ncit)->set_neighbor(j, c);
	    // opposite relation will be set when ncit arrives on c
	    // this avoids to look for the correct index
	    // and to test whether *ncit already has neighbor i
	  }
	  else {
	    // star is a vertex of n->neighbor(i)
	    set_adjacency(*ncit, j, n->neighbor(i), 3);//neighbor opposite to v
	  }
	}
      }
    }
  }// end switch

  return v;
}

template <class Vb, class Cb >
void
Triangulation_data_structure_3<Vb,Cb>::
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
	delete_cell(v->cell());
    }
    else {
        // the cells incident to w are down graded one dimension
        // the other cells are deleted
        std::vector<Cell_handle> to_delete, to_downgrade;

        for (Cell_iterator ib = cells().begin();
            ib != cells().end(); ++ib) {
            if ( ib->has_vertex(w) )
	        to_downgrade.push_back(ib);
            else
	        to_delete.push_back(ib);
        }

        typename std::vector<Cell_handle>::iterator lfit=to_downgrade.begin();
        for( ; lfit != to_downgrade.end(); ++lfit) {
	    Cell_handle f = *lfit;
	    int j = f->index(w);
	    int k; if (f->has_vertex(v, k)) f->set_vertex(k, w);
            if (j != dimension()) {
	        f->set_vertex(j, f->vertex(dimension()));
	        f->set_neighbor(j, f->neighbor(dimension()));
	        if (dimension() >= 1)
		    change_orientation(f);
	    }
	    f->set_vertex(dimension(), Vertex_handle());
	    f->set_neighbor(dimension(), Cell_handle());

	    // Update vertex->cell() pointers.
	    for (int i = 0; i < dimension(); ++i)
	        f->vertex(i)->set_cell(f);
        }

        delete_cells(to_delete.begin(), to_delete.end());
    }
    delete_vertex(v);
    set_dimension(dimension()-1);
    CGAL_triangulation_postcondition(is_valid());
}

template <class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Cell_handle
Triangulation_data_structure_3<Vb,Cb>::
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

template <class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Cell_handle
Triangulation_data_structure_3<Vb,Cb>::
remove_degree_2(Vertex_handle v)
{
    CGAL_triangulation_precondition(dimension() == 1);
    CGAL_triangulation_precondition(degree(v) == 2);
    CGAL_triangulation_precondition(number_of_vertices() >= 4);

    // Cells to be killed.
    Cell_handle c0, c1;
    // Indices of v in these cells.
    int i0, i1;

    c0 = v->cell();
    i0 = c0->index(v);
    c1 = c0->neighbor((i0 == 0) ? 1 : 0);
    i1 = c1->index(v);

    // New cell : we copy the content of c0, so we keep the orientation.
    Cell_handle newc = create_face(c0->vertex(0),
                                   c0->vertex(1),
				   Vertex_handle());

    newc->set_vertex(i0, c1->vertex(c1->index(c0)));

    set_adjacency(newc,    i0, c0->neighbor(i0), mirror_index(c0, i0));
    set_adjacency(newc,  1-i0, c1->neighbor(i1), mirror_index(c1, i1));

    newc->vertex(0)->set_cell(newc);
    newc->vertex(1)->set_cell(newc);

    delete_cell(c0);
    delete_cell(c1);
    delete_vertex(v);

    return newc;
}

template <class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Cell_handle
Triangulation_data_structure_3<Vb,Cb>::
remove_degree_3(Vertex_handle v)
{
    CGAL_triangulation_precondition(dimension() == 2);
    CGAL_triangulation_precondition(degree(v) == 3);
    CGAL_triangulation_precondition(number_of_vertices() >= 5);

    // Cells to be killed.
    Cell_handle c0, c1, c2;
    // Indices of v in these cells.
    int i0, i1, i2;

    c0 = v->cell();
    i0 = c0->index(v);
    c1 = c0->neighbor(cw(i0));
    i1 = c1->index(v);
    c2 = c0->neighbor(ccw(i0));
    i2 = c2->index(v);

    // New cell : we copy the content of c0, so we keep the orientation.
    Cell_handle newc = create_face(c0->vertex(0),
                                   c0->vertex(1),
                                   c0->vertex(2));

    newc->set_vertex(i0, c1->vertex(c1->index(c0)));

    set_adjacency(newc,      i0, c0->neighbor(i0), mirror_index(c0, i0));
    set_adjacency(newc,  cw(i0), c1->neighbor(i1), mirror_index(c1, i1));
    set_adjacency(newc, ccw(i0), c2->neighbor(i2), mirror_index(c2, i2));

    newc->vertex(0)->set_cell(newc);
    newc->vertex(1)->set_cell(newc);
    newc->vertex(2)->set_cell(newc);

    delete_cell(c0);
    delete_cell(c1);
    delete_cell(c2);
    delete_vertex(v);

    return newc;
}

template <class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Cell_handle
Triangulation_data_structure_3<Vb,Cb>::
remove_degree_4(Vertex_handle v)
{
    CGAL_triangulation_precondition(dimension() == 3);
    CGAL_triangulation_precondition(degree(v) == 4);
    CGAL_triangulation_precondition(number_of_vertices() >= 6);

    // Cells to be killed.
    Cell_handle c0, c1, c2, c3;
    // Indices of v in these cells.
    int i0, i1, i2, i3;

    c0 = v->cell();
    i0 = c0->index(v);
    c1 = c0->neighbor(i0^1);
    i1 = c1->index(v);
    c2 = c0->neighbor(i0^2);
    i2 = c2->index(v);
    c3 = c0->neighbor(i0^3);
    i3 = c3->index(v);

    // New cell : we copy the content of c0, so we keep the orientation.
    Cell_handle newc = create_cell(c0->vertex(0),
                                   c0->vertex(1),
                                   c0->vertex(2),
                                   c0->vertex(3));

    newc->set_vertex(i0, c1->vertex(c1->index(c0)));

    set_adjacency(newc,   i0, c0->neighbor(i0), mirror_index(c0, i0));
    set_adjacency(newc, i0^1, c1->neighbor(i1), mirror_index(c1, i1));
    set_adjacency(newc, i0^2, c2->neighbor(i2), mirror_index(c2, i2));
    set_adjacency(newc, i0^3, c3->neighbor(i3), mirror_index(c3, i3));

    newc->vertex(0)->set_cell(newc);
    newc->vertex(1)->set_cell(newc);
    newc->vertex(2)->set_cell(newc);
    newc->vertex(3)->set_cell(newc);

    delete_cell(c0);
    delete_cell(c1);
    delete_cell(c2);
    delete_cell(c3);
    delete_vertex(v);

    return newc;
}

template <class Vb, class Cb >
void
Triangulation_data_structure_3<Vb,Cb>::
decrease_dimension(Cell_handle c, int i)
{
  CGAL_triangulation_expensive_precondition( is_valid() );;
  CGAL_triangulation_precondition( dimension() >= 2);
  CGAL_triangulation_precondition( number_of_vertices() >
                                   (size_type) dimension() + 1 );
  CGAL_triangulation_precondition( degree(c->vertex(i)) == number_of_vertices()-1 );

  Vertex_handle v = c->vertex(i);
  Vertex_handle w = c->vertex(i);

  // the cells incident to w are down graded one dimension
  // the other cells are deleted
  std::vector<Cell_handle> to_delete, to_downgrade;

  for (Cell_iterator ib = cells().begin();
       ib != cells().end(); ++ib) {
    if ( ib->has_vertex(w) )
      to_downgrade.push_back(ib);
    else
      to_delete.push_back(ib);
  }

  typename std::vector<Cell_handle>::iterator lfit=to_downgrade.begin();
  for( ; lfit != to_downgrade.end(); ++lfit) {
    Cell_handle f = *lfit;
    int j = f->index(w);
    int k; 
    if (f->has_vertex(v, k)) f->set_vertex(k, w);
    if (j != dimension()) {
      f->set_vertex(j, f->vertex(dimension()));
      f->set_neighbor(j, f->neighbor(dimension()));
      if (dimension() >= 1)
        change_orientation(f);
    }
    f->set_vertex(dimension(), Vertex_handle());
    f->set_neighbor(dimension(), Cell_handle());

    // Update vertex->cell() pointers.
    for (int i = 0; i < dimension(); ++i)
      f->vertex(i)->set_cell(f);
  }

  delete_cells(to_delete.begin(), to_delete.end());

  //delete_vertex(v);
  set_dimension(dimension()-1);

  if(dimension() == 2)
  {
    Cell_handle n0 = c->neighbor(0);
    Cell_handle n1 = c->neighbor(1);
    Cell_handle n2 = c->neighbor(2);
    Vertex_handle v0 = c->vertex(0);
    Vertex_handle v1 = c->vertex(1);
    Vertex_handle v2 = c->vertex(2);
		
    int i0 = 0, i1 = 0, i2 = 0;
		
    for(int i=0; i<3; i++) if(n0->neighbor(i) == c) { i0 = i; break; }
    for(int i=0; i<3; i++) if(n1->neighbor(i) == c) { i1 = i; break; }
    for(int i=0; i<3; i++) if(n2->neighbor(i) == c) { i2 = i; break; }				
		
    Cell_handle c1 = create_cell(v, v0, v1, Vertex_handle());
    Cell_handle c2 = create_cell(v, v1, v2, Vertex_handle());

    c->set_vertex(0, v);
    c->set_vertex(1, v2);
    c->set_vertex(2, v0);
    c->set_vertex(3, Vertex_handle());

    //Cell_handle c3 = create_cell(v, v2, v0, Vertex_handle());		
    Cell_handle c3 = c;
		
    c1->set_neighbor(0, n2); n2->set_neighbor(i2, c1);
    c1->set_neighbor(1, c2); 
    c1->set_neighbor(2, c3);
    c1->set_neighbor(3, Cell_handle());
		
    c2->set_neighbor(0, n0); n0->set_neighbor(i0, c2);
    c2->set_neighbor(1, c3); 
    c2->set_neighbor(2, c1);
    c2->set_neighbor(3, Cell_handle());	
		
    c3->set_neighbor(0, n1); n1->set_neighbor(i1, c3);
    c3->set_neighbor(1, c1); 
    c3->set_neighbor(2, c2);
    c3->set_neighbor(3, Cell_handle());
		
    v->set_cell(c1);
    v0->set_cell(c1);
    v1->set_cell(c1);
    v2->set_cell(c2);	
  }
	
  if(dimension() == 1)
  {
    Cell_handle n0 = c->neighbor(0);
    Cell_handle n1 = c->neighbor(1);
    Vertex_handle v0 = c->vertex(0);
    Vertex_handle v1 = c->vertex(1);
		
    int i0 = 0 , i1 = 0;
		
    for(int i=0; i<2; i++) if(n0->neighbor(i) == c) { i0 = i; break; }
    for(int i=0; i<2; i++) if(n1->neighbor(i) == c) { i1 = i; break; }
		
    Cell_handle c1 = create_cell(v0, v, Vertex_handle(), Vertex_handle());
		
    c->set_vertex(0, v);
    c->set_vertex(1, v1);
    c->set_vertex(2, Vertex_handle());
    c->set_vertex(3, Vertex_handle());

    //Cell_handle c2 = create_cell(v, v1, Vertex_handle(), Vertex_handle());	
    Cell_handle c2 = c;
		
    c1->set_neighbor(0, c2); 
    c1->set_neighbor(1, n1); n1->set_neighbor(i1, c1);
    c1->set_neighbor(2, Cell_handle());
    c1->set_neighbor(3, Cell_handle());
		
    c2->set_neighbor(0, n0); n0->set_neighbor(i0, c2);
    c2->set_neighbor(1, c1); 
    c2->set_neighbor(2, Cell_handle());
    c2->set_neighbor(3, Cell_handle());
		
    v->set_cell(c1);
    v0->set_cell(c1);
    v1->set_cell(c2);
  }
	
  CGAL_triangulation_postcondition(is_valid());
}


template <class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::size_type
Triangulation_data_structure_3<Vb,Cb>::
degree(Vertex_handle v) const
{
    std::size_t res;
    adjacent_vertices(v, Counting_output_iterator(&res));
    return res;
}

template <class Vb, class Cb >
bool
Triangulation_data_structure_3<Vb,Cb>::
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
      // no break; continue
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

template <class Vb, class Cb >
bool
Triangulation_data_structure_3<Vb,Cb>::
is_valid(Vertex_handle v, bool verbose, int level) const
{
  bool result = v->is_valid(verbose,level);
  result = result && v->cell()->has_vertex(v);
  if ( ! result ) {
    if ( verbose )
      std::cerr << "invalid vertex" << std::endl;
    CGAL_triangulation_assertion(false);
  }
  return result;
}

template <class Vb, class Cb >
bool
Triangulation_data_structure_3<Vb,Cb>::
is_valid(Cell_handle c, bool verbose, int level) const
{
    if ( ! c->is_valid(verbose, level) )
	return false;

    switch (dimension()) {
    case -2:
    case -1:
    {
      if ( c->vertex(0) == Vertex_handle() ) {
	if (verbose)
	    std::cerr << "vertex 0 NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      is_valid(c->vertex(0),verbose,level);
      if ( c->vertex(1) != Vertex_handle() || c->vertex(2) != Vertex_handle()) {
	if (verbose)
	    std::cerr << "vertex 1 or 2 != NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      if ( c->neighbor(0) != Cell_handle() ||
	   c->neighbor(1) != Cell_handle() ||
	   c->neighbor(2) != Cell_handle()) {
	if (verbose)
	    std::cerr << "one neighbor != NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      break;
    }

    case 0:
      {
      if ( c->vertex(0) == Vertex_handle() ) {
	if (verbose)
	    std::cerr << "vertex 0 NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      is_valid(c->vertex(0),verbose,level);
      if ( c->neighbor (0) == Cell_handle() ) {
	if (verbose)
	    std::cerr << "neighbor 0 NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      if ( c->vertex(1) != Vertex_handle() ||
           c->vertex(2) != Vertex_handle() ) {
	if (verbose)
	    std::cerr << "vertex 1 or 2 != NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      if ( c->neighbor(1) != Cell_handle() ||
           c->neighbor(2) != Cell_handle() ) {
	if (verbose)
	    std::cerr << "neighbor 1 or 2 != NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }

      if ( ! c->neighbor(0)->has_vertex(c->vertex(0)) ) {
	if (verbose)
	    std::cerr << "neighbor 0 does not have vertex 0" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      break;
      }

    case 1:
      {
      Vertex_handle v0 = c->vertex(0);
      Vertex_handle v1 = c->vertex(1);
      Cell_handle n0 = c->neighbor(0);
      Cell_handle n1 = c->neighbor(1);

      if ( v0 == Vertex_handle() || v1 == Vertex_handle() ) {
	if (verbose)
	    std::cerr << "vertex 0 or 1 NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      is_valid(c->vertex(0),verbose,level);
      is_valid(c->vertex(1),verbose,level);
      if ( n0 == Cell_handle() || n1 == Cell_handle() ) {
	if (verbose)
	    std::cerr << "neighbor 0 or 1 NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }

      if ( v0 !=  n1->vertex(1) ) {
	if (verbose)
	    std::cerr << "neighbor 1 does not have vertex 0 as vertex 1"
		      << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      if ( v1 != n0->vertex(0) ) {
	if (verbose)
	    std::cerr << "neighbor 0 does not have vertex 1 as vertex 0"
		      << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }

      if ( n0->neighbor(1) != c ) {
	if (verbose)
	    std::cerr << "neighbor 0 does not have this as neighbor 1"
		      << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      if ( n1->neighbor(0) != c ) {
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
      if ( c->vertex(0) == Vertex_handle() ||
	   c->vertex(1) == Vertex_handle() ||
	   c->vertex(2) == Vertex_handle() ) {
	if (verbose)
	    std::cerr << "vertex 0, 1, or 2 NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      is_valid(c->vertex(0),verbose,level);
      is_valid(c->vertex(1),verbose,level);
      is_valid(c->vertex(2),verbose,level);
      int in;
      Cell_handle n;
      for(int i = 0; i < 3; i++) {
	n = c->neighbor(i);
	if ( n == Cell_handle() ) {
	  if (verbose)
	      std::cerr << "neighbor " << i << " NULL" << std::endl;
	  CGAL_triangulation_assertion(false);
	  return false;
	}
	if ( ! n->has_vertex(c->vertex(cw(i)),in ) ) {
	  if (verbose)
	      std::cerr << "vertex " << cw(i)
		        << " not vertex of neighbor " << i << std::endl;
	  CGAL_triangulation_assertion(false);
	  return false;
	}
	in = cw(in);
	if ( n->neighbor(in) != c ) {
	  if (verbose)
	      std::cerr << "neighbor " << i
		        << " does not have this as neighbor "
		        << in << std::endl;
	  CGAL_triangulation_assertion(false);
	  return false;
	}
	if ( c->vertex(ccw(i)) != n->vertex(cw(in)) ) {
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
	  if ( c->vertex(i) == Vertex_handle() ) {
	    if (verbose)
		std::cerr << "vertex " << i << " NULL" << std::endl;
	    CGAL_triangulation_assertion(false);
	    return false;
	  }
	  is_valid(c->vertex(i),verbose,level);
	}

	for(i = 0; i < 4; i++) {
	  Cell_handle n = c->neighbor(i);
	  if ( n == Cell_handle() ) {
	    if (verbose)
	      std::cerr << "neighbor " << i << " NULL" << std::endl;
	    CGAL_triangulation_assertion(false);
	    return false;
	  }

	  int in = 5;
	  // if ( ! n->has_neighbor(handle(), in) ) {
          if ( n->neighbor(0) == c) in = 0;
          if ( n->neighbor(1) == c) in = 1;
          if ( n->neighbor(2) == c) in = 2;
          if ( n->neighbor(3) == c) in = 3;
          if (in == 5) {
	    if (verbose)
              std::cerr << "neighbor of c has not c as neighbor" << std::endl;
	    CGAL_triangulation_assertion(false);
	    return false;
	  }
	
	  int j1n,j2n,j3n;
	  if ( ! n->has_vertex(c->vertex((i+1)&3),j1n) ) {
	    if (verbose) { std::cerr << "vertex " << ((i+1)&3)
				     << " not vertex of neighbor "
				     << i << std::endl; }
	    CGAL_triangulation_assertion(false);
	    return false;
	  }
	  if ( ! n->has_vertex(c->vertex((i+2)&3),j2n) ) {
	    if (verbose) { std::cerr << "vertex " << ((i+2)&3)
				     << " not vertex of neighbor "
				     << i << std::endl; }
	    CGAL_triangulation_assertion(false);
	    return false;
	  }
	  if ( ! n->has_vertex(c->vertex((i+3)&3),j3n) ) {
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

template <class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
copy_tds(const Tds & tds, Vertex_handle vert )
  // returns the new vertex corresponding to vert in the new tds
{
  CGAL_triangulation_expensive_precondition( vert == Vertex_handle()
	                                  || tds.is_vertex(vert) );

  clear();

  size_type n = tds.number_of_vertices();
  set_dimension(tds.dimension());

  // Number of pointers to cell/vertex to copy per cell.
  int dim = (std::max)(1, dimension() + 1);

  if (n == 0)
      return vert;

  // Create the vertices.
  std::vector<Vertex_handle> TV(n);
  size_type i = 0;

  for (Vertex_iterator vit = tds.vertices_begin();
       vit != tds.vertices_end(); ++vit)
    TV[i++] = vit;

  CGAL_triangulation_assertion( i == n );

  Unique_hash_map<Vertex_handle, Vertex_handle> V;
  Unique_hash_map<Cell_handle, Cell_handle > F;

  for (size_type i=0; i <= n-1; ++i)
    V[ TV[i] ] = create_vertex(TV[i]);

  // Create the cells.
  for (Cell_iterator cit = tds.cells().begin();
	  cit != tds.cells_end(); ++cit) {
      F[cit] = create_cell(cit);
      for (int j = 0; j < dim; j++)
        F[cit]->set_vertex(j, V[cit->vertex(j)] );
  }

  // Link the vertices to a cell.
  for (Vertex_iterator vit2 = tds.vertices_begin();
       vit2 != tds.vertices_end(); ++vit2)
    V[vit2]->set_cell( F[vit2->cell()] );

  // Hook neighbor pointers of the cells.
  for (Cell_iterator cit2 = tds.cells().begin();
	  cit2 != tds.cells_end(); ++cit2) {
    for (int j = 0; j < dim; j++)
      F[cit2]->set_neighbor(j, F[cit2->neighbor(j)] );
  }

  CGAL_triangulation_postcondition( is_valid() );

  return (vert != Vertex_handle()) ? V[vert] : vert;
}

template <class Vb, class Cb >
void
Triangulation_data_structure_3<Vb,Cb>::
swap(Tds & tds)
{
  CGAL_triangulation_expensive_precondition(tds.is_valid() && is_valid());

  std::swap(_dimension, tds._dimension);
  cells().swap(tds.cells());
  vertices().swap(tds.vertices());
}

template <class Vb, class Cb >
void
Triangulation_data_structure_3<Vb,Cb>::
clear()
{
  cells().clear();
  vertices().clear();
  set_dimension(-2);
}

template <class Vb, class Cb >
bool
Triangulation_data_structure_3<Vb,Cb>::
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

template <class Vb, class Cb >
bool
Triangulation_data_structure_3<Vb,Cb>::
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

template <class Vb, class Cb >
bool
Triangulation_data_structure_3<Vb,Cb>::
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

template <class Vb, class Cb >
bool
Triangulation_data_structure_3<Vb,Cb>::
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

#endif // CGAL_TRIANGULATION_DATA_STRUCTURE_3_H
