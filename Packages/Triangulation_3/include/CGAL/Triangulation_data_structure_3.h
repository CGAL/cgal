// ============================================================================
//
// Copyright (c) 1999,2000,2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Triangulation_data_structure_3.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================
//
// combinatorial triangulation of the boundary of a polytope
// of dimension d in dimension d+1
// for -1 <= d <= 3
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_DATA_STRUCTURE_3_H
#define CGAL_TRIANGULATION_DATA_STRUCTURE_3_H

#include <CGAL/basic.h>

#include <utility>
#include <map>
#include <set>
#include <vector>

#include <CGAL/triple.h>

#include <CGAL/Triangulation_short_names_3.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_utils_3.h>

#include <CGAL/Triangulation_cell_base_3.h>

#include <CGAL/Triangulation_ds_cell_3.h>
#include <CGAL/Triangulation_ds_vertex_3.h>

#include <CGAL/Triangulation_ds_iterators_3.h>
#include <CGAL/Triangulation_ds_circulators_3.h>
#include <CGAL/Triangulation_handles_3.h>

#include <CGAL/DS_Container.h>

CGAL_BEGIN_NAMESPACE

template <class Vb, class Cb>
class Triangulation_data_structure_3
  : public Triangulation_utils_3
{
  typedef Triangulation_data_structure_3<Vb, Cb>   Self;
public:

  typedef Triangulation_data_structure_3<Vb,Cb>    Tds;
  typedef Vb                                       Vertex_base;
  typedef Cb                                       Cell_base;

  typedef Triangulation_ds_vertex_3<Tds>           Vertex;
  typedef Triangulation_ds_cell_3<Tds>             Cell;

  typedef Triangulation_cell_handle_3<Self>        Cell_handle;
  typedef Triangulation_vertex_handle_3<Self>      Vertex_handle;
  // TODO : The following (simpler), should work.  One day...
  //typedef Pointer<Cell>        Cell_handle;
  //typedef Pointer<Vertex>      Vertex_handle;
  //typedef Cell *        Cell_handle;
  //typedef Vertex *      Vertex_handle;

  typedef std::pair<Cell_handle, int>              Facet;
  typedef triple<Cell_handle, int, int>            Edge;

  friend class Triangulation_ds_facet_iterator_3<Tds>;
  friend class Triangulation_ds_edge_iterator_3<Tds>;

  friend class Triangulation_ds_cell_circulator_3<Tds>;
  friend class Triangulation_ds_facet_circulator_3<Tds>;

  typedef DS_Container<Cell>                       Cell_container;
  typedef DS_Container<Vertex>                     Vertex_container;
  typedef typename Cell_container::iterator        Cell_iterator;
  typedef typename Vertex_container::iterator      Vertex_iterator;
  typedef Triangulation_ds_facet_iterator_3<Tds>   Facet_iterator;
  typedef Triangulation_ds_edge_iterator_3<Tds>    Edge_iterator;

  typedef Triangulation_ds_cell_circulator_3<Tds>  Cell_circulator;
  typedef Triangulation_ds_facet_circulator_3<Tds> Facet_circulator;

  Triangulation_data_structure_3() 
    : _dimension(-2), _number_of_vertices(0)
  {}

  Triangulation_data_structure_3(const Tds & tds)
    : _number_of_vertices(0)
    // _number_of_vertices is set to 0 so that clear() in copy_tds() works
  {
    copy_tds(tds);
  }

  ~Triangulation_data_structure_3()
  {
    clear();
  }

  Tds & operator= (const Tds & tds)
  {
    copy_tds(tds);
    return *this;
  }  

  int number_of_vertices() const {return _number_of_vertices;}
  
  int dimension() const {return _dimension;}

  int number_of_cells() const 
    { 
      if ( dimension() < 3 ) return 0;
      return cell_container().size();
    }
  
  int number_of_facets() const
    {
      if ( dimension() < 2 ) return 0;
      return std::distance(facets_begin(), facets_end());
    }

  int number_of_edges() const
    {
      if ( dimension() < 1 ) return 0;
      return std::distance(edges_begin(), edges_end());
    }

  // USEFUL CONSTANT TIME FUNCTIONS

  // SETTING
  // to be protected ?

  void set_number_of_vertices(int n) { _number_of_vertices = n; }

  void set_dimension(int n) { _dimension = n; }

  Vertex_handle create_vertex()
  {
      return vertex_container().get_new_element();
  }

  Cell_handle create_cell() 
    { 
      return get_new_cell();
    }

  Cell_handle create_cell(Cell_handle c)
    {
      Cell_handle cnew = get_new_cell();
      *cnew = *c;
      cnew->init();
      return cnew; 
    }

  Cell_handle create_cell(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3)
    {
      Cell_handle c = get_new_cell();
      c->set_vertices(v0,v1,v2,v3);
      return c; 
    }

  Cell_handle create_cell(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3,
		    Cell_handle n0, Cell_handle n1, Cell_handle n2, Cell_handle n3)
    {
      Cell_handle c = get_new_cell();
      c->set_vertices(v0,v1,v2,v3);
      c->set_neighbors(n0,n1,n2,n3);
      return c; 
    }

private:

  Cell_handle get_new_cell()
  {
      Cell_handle r = cell_container().get_new_element();
      r->init();
      return r;
  }

public:
  // not documented
  void read_cells(std::istream& is, std::map< int, Vertex_handle > &V,
			   int & m, std::map< int, Cell_handle > &C );
  // not documented
  void print_cells(std::ostream& os, std::map<Vertex_handle, int> &V ) const;

  // ACCESS FUNCTIONS

  void delete_vertex( Vertex_handle v )
  {
      CGAL_triangulation_expensive_precondition( is_vertex(v) );
      vertex_container().release_element(&*v);
  }

  void delete_cell( Cell_handle c )
  {
      CGAL_triangulation_expensive_precondition( dimension() != 3 ||
                                                 is_cell(c) );
      CGAL_triangulation_expensive_precondition( dimension() != 2 ||
                                                 is_facet(c,3) );
      CGAL_triangulation_expensive_precondition( dimension() != 1 ||
                                                 is_edge(c,0,1) );
      CGAL_triangulation_expensive_precondition( dimension() != 0 ||
                                                 is_vertex(c->vertex(0)) );

      cell_container().release_element(&*c);
  }

  template <class It>
  void delete_cells(It begin, It end)
  {
      for(It i = begin; i != end; ++i)
	  delete_cell(&**i);
  }

  // QUERIES

  bool is_vertex(Vertex_handle v) const;
  bool is_edge(Cell_handle c, int i, int j) const;
  bool is_edge(Vertex_handle u, Vertex_handle v, Cell_handle & c, int & i, int & j) const;
  bool is_facet(Cell_handle c, int i) const;
  bool is_facet(Vertex_handle u, Vertex_handle v, Vertex_handle w, 
		Cell_handle & c, int & i, int & j, int & k) const;
  bool is_cell(Cell_handle c) const;
  bool is_cell(Vertex_handle u, Vertex_handle v, Vertex_handle w, Vertex_handle t, 
	       Cell_handle & c, int & i, int & j, int & k, int & l) const;
  bool is_cell(Vertex_handle u, Vertex_handle v, Vertex_handle w, Vertex_handle t) const; 

  bool has_vertex(const Facet & f, Vertex_handle v, int & j) const;
  bool has_vertex(Cell_handle c, int i, Vertex_handle v, int & j) const;
  bool has_vertex(const Facet & f, Vertex_handle v) const;
  bool has_vertex(Cell_handle c, int i, Vertex_handle v) const;

  bool are_equal(Cell_handle c, int i, Cell_handle n, int j) const;
  bool are_equal(const Facet & f, const Facet & g) const;
  bool are_equal(const Facet & f, Cell_handle n, int j) const;

  // MODIFY

  bool flip(Cell_handle c, int i);
  bool flip(const Facet &f)
  { return flip( f.first, f.second); }

  void flip_flippable(Cell_handle c, int i);
  void flip_flippable(const Facet &f)
  { return flip_flippable( f.first, f.second ); }

  bool flip(Cell_handle c, int i, int j);
  bool flip(const Edge &e)
  { return flip( e.first, e.second, e.third ); }

  void flip_flippable(Cell_handle c, int i, int j);
  void flip_flippable(const Edge &e)
  { return flip_flippable( e.first, e.second, e.third ); }

private:
  // common to flip and flip_flippable
  void flip_really(Cell_handle c, int i, Cell_handle n, int in);
  void flip_really(Cell_handle c, int i, int j,
		   Cell_handle c1, Vertex_handle v1, int i1, int j1, int next1,
		   Cell_handle c2, Vertex_handle v2, int i2, int j2, int next2,
		   Vertex_handle v3);
public:

  //INSERTION

  Vertex_handle insert_in_cell(Vertex_handle v, Cell_handle c);

  Vertex_handle insert_in_facet(Vertex_handle v, const Facet & f)
    { return insert_in_facet(w,f.first,f.second); }
  
  Vertex_handle insert_in_facet(Vertex_handle v, Cell_handle c, int i);
  
  Vertex_handle insert_in_edge(Vertex_handle v, const Edge & e)   
    { return insert_in_edge(w, e.first, e.second, e.third); }
  
  Vertex_handle insert_in_edge(Vertex_handle v, Cell_handle c, int i, int j);   

  Vertex_handle insert_increase_dimension(Vertex_handle v, // new vertex
				     Vertex_handle star = NULL,
				     bool reorient = false);

  // I think the find_conflicts_[23] should be in the TDS, as they do not
  // depend on geometry, although they may, depending on the actual
  // Conflict_test.  But e.g. for incident_cells(), the tester is not
  // geometric, so it makes sense to have it in the TDS.

  template <class FacetIt, class CellIt>
  void star_hole_3(Vertex_handle newv, FacetIt facet_begin, FacetIt facet_end,
	           CellIt cell_begin, CellIt cell_end)
  {
      star_hole_3(newv, facet_begin, facet_end);
      delete_cells(cell_begin, cell_end);
  }

  template <class FacetIt, class CellIt>
  void star_hole_2(Vertex_handle newv, FacetIt facet_begin, FacetIt facet_end,
	           CellIt cell_begin, CellIt cell_end)
  {
      star_hole_2(newv, facet_begin, facet_end);
      delete_cells(cell_begin, cell_end);
  }

  // Facets->first is in conflict, and we walk inside the hole.
  //
  // Note #1 : If we can merge the FacetIt loops, then maybe we can get rid
  // of the corresponding container by just walking over the boundary ?
  // Thinking a bit more about that : I think it's either we have a container
  // of facets, or we have a recursive function over the boundary... (?)
  template <class FacetIt>
  void star_hole_3(Vertex_handle newv, FacetIt facet_begin, FacetIt facet_end)
  {
      CGAL_triangulation_precondition(dimension()==3);

      // Would be nice if there were already room reserved in the facet
      // vector.
      std::vector<Cell_handle > V;
      V.reserve(std::distance(facet_begin, facet_end));

      // For each facet on the boundary :
      // - create a new cell, link its vertices and one cell pointer.
      for (FacetIt fit = facet_begin; fit != facet_end; ++fit) {
	  Cell_handle old = fit->first;
	  Cell_handle bound = old->neighbor(fit->second);
	  // Note that the initial orientation of the new cells is positive,
	  // as we copy it from an existing one.
	  Cell_handle newc = create_cell(old->vertex(0),
		                   old->vertex(1),
		                   old->vertex(2),
		                   old->vertex(3));
	  newc->set_vertex(fit->second, newv);
	  set_adjacency(newc, bound, fit->second, bound->index(old));
	  newc->vertex(0)->set_cell(newc);
	  newc->vertex(1)->set_cell(newc);
	  newc->vertex(2)->set_cell(newc);
	  newc->vertex(3)->set_cell(newc);
	  V.push_back(newc);
      }

      // For each facet on the boundary, for each of the 3 edges :
      // - we must find the neighbor facet
      // - link the 2 corresponding new cells.
      int zz = -1;
      for (FacetIt fit2 = facet_begin; fit2 != facet_end; ++fit2) {
	  ++zz;
	  Cell_handle old = fit2->first;
	  Cell_handle newc = V[zz];

	  for (int i=0; i<=3; ++i) {
	      // We must avoid i == fit->second, but the following
	      // test will avoid it too.
	      if (newc->neighbor(i) != NULL)
		  continue;

	      // Now we turn around the edge inside the hole.
	      // To recognize when we hit the boundary, we look at the
	      // neighbor, and see if it doesn't point back to us, in which
	      // case it's the boundary cell we are looking for.
	      Cell_handle t = old;
	      Vertex_handle k = t->vertex(fit2->second);
	      int j = i;
	      Cell_handle newt = t->neighbor(j);
	      int z;
	      while (newt->has_neighbor(t, z)) {
		  j = newt->index(k);
		  k = newt->vertex(z);
		  t = newt;
		  newt = t->neighbor(j);
	      };

	      // Compute the address of the corresponding new cell.
	      Cell_handle back;
	      for (int l=0;; ++l) {
	          back = newt->neighbor(l);
		  if (l==3)
		      break;
	          // The vertices are those of t except vertex(j) = newv.
		  if (back->vertex(j)   == newv &&
		      back->vertex(j^1) == t->vertex(j^1) &&
		      back->vertex(j^2) == t->vertex(j^2) &&
		      back->vertex(j^3) == t->vertex(j^3))
		      break;
	      }

	      set_adjacency(newc, back, i, back->index(k));
	  }
      }
  }

  // Note : the code is almost entirely duplicated, just to have dimension() a
  // constant for performance.  That's not perfect...
  template <class FacetIt>
  void star_hole_2(Vertex_handle newv, FacetIt facet_begin, FacetIt facet_end)
  {
      CGAL_triangulation_precondition(dimension()==2);

      // Would be nice if there were already room reserved in the facet
      // vector.
      std::vector<Cell_handle > V;
      V.reserve(std::distance(facet_begin, facet_end));

      // For each facet on the boundary :
      // - create a new cell, link its vertices and one cell pointer.
      for (FacetIt fit = facet_begin; fit != facet_end; ++fit) {
	  Cell_handle old = fit->first;
	  Cell_handle bound = old->neighbor(fit->second);
	  // Note that the initial orientation of the new cells is positive,
	  // as we copy it from an existing one.
	  Cell_handle newc = create_cell(old->vertex(0),
		                   old->vertex(1),
		                   old->vertex(2),
		                   NULL);
	  newc->set_vertex(fit->second, newv);
	  set_adjacency(newc, bound, fit->second, bound->index(old));
	  newc->vertex(0)->set_cell(newc);
	  newc->vertex(1)->set_cell(newc);
	  newc->vertex(2)->set_cell(newc);
	  V.push_back(newc);
      }

      // For each facet on the boundary, for each of the 3 edges :
      // - we must find the neighbor facet
      // - link the 2 corresponding new cells.
      int zz = -1;
      for (FacetIt fit2 = facet_begin; fit2 != facet_end; ++fit2) {
	  ++zz;
	  Cell_handle old = fit2->first;
	  Cell_handle newc = V[zz];

	  for (int i=0; i<=2; ++i) {
	      // We must avoid i == fit->second, but the following
	      // test will avoid it too.
	      if (newc->neighbor(i) != NULL)
		  continue;

	      // Now we turn around the edge inside the hole.
	      // To recognize when we hit the boundary, we look at the
	      // neighbor, and see if it doesn't point back to us, in which
	      // case it's the boundary cell we are looking for.
	      Cell_handle t = old;
	      Vertex_handle k = t->vertex(fit2->second);
	      int j = i;
	      Cell_handle newt = t->neighbor(j);
	      int z;
	      while (newt->has_neighbor(t, z)) {
		  j = newt->index(k);
		  k = newt->vertex(z);
		  t = newt;
		  newt = t->neighbor(j);
	      };

	      // Compute the address of the corresponding new cell.
	      Cell_handle back;
	      for (int l=0;; ++l) {
	          back = newt->neighbor(l);
		  if (l==2)
		      break;
	          // The vertices are those of t except vertex(j) = newv.
		  if (back->vertex(j)      == newv &&
		      back->vertex(cw(j))  == t->vertex(cw(j)) &&
		      back->vertex(ccw(j)) == t->vertex(ccw(j)))
		      break;
	      }

	      set_adjacency(newc, back, i, back->index(k));
	  }
      }
  }

  // ITERATOR METHODS

  Cell_iterator cells_begin() const
  {
    if ( dimension() < 3 )
	return cells_end();
    return cell_container().begin();
  }

  Cell_iterator cells_end() const
  {
    return cell_container().end();
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
    return vertex_container().begin();
  }

  Vertex_iterator vertices_end() const
  {
    return vertex_container().end();
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

  Cell_circulator incident_cells(const Edge & e, Cell_handle start) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    return Cell_circulator(e, start);
  }
  Cell_circulator incident_cells(Cell_handle ce, int i, int j, Cell_handle start) const
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
  Facet_circulator incident_facets(const Edge & e, Cell_handle start, int f) const
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

  // around a vertex
private:
  class Incident_tester {
      Vertex_handle _v;
  public :
      Incident_tester(Vertex_handle v) : _v(v) {}
      bool operator()(Cell_handle c) const
      { return c->has_vertex(v); }
  };

private:
  // TODO : This should be reimplemented.  Using find_conflict is overkill...
  template <class OutIt>
  void
  incident_cells(Vertex_handle v, OutIt outcells) const
  {
      if ( dimension() < 3 )
          return;

      std::vector<Facet> facets;
      std::vector<Cell_handle > cells;
      facets.reserve(64);
      cells.reserve(64);
      find_conflicts_3(v->cell(), Incident_tester(v), facets, cells);
      for(typename std::vector<Facet>::iterator fit = facets.begin();
	      fit != facets.end(); ++fit)
	  fit->first->set_in_conflict_flag(0);
      for(typename std::vector<Cell_handle >::iterator cit = cells.begin();
	      cit != cells.end(); ++cit) {
	  *cit->set_in_conflict_flag(0);
	  *outcells++ = *cit;
      }
  }
  
public:
  void
  incident_cells(Vertex_handle v, std::set<Cell_handle> & cells,
	         Cell_handle c = NULL ) const;

  void
  incident_vertices(Vertex_handle v, std::set<Vertex_handle> & vertices,
		    Cell_handle c = NULL ) const;

  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const;


  // Helping functions
  Vertex_handle copy_tds(const Tds & tds, Vertex_handle vert = NULL);
    // returns the new vertex corresponding to vert in the new tds 

  void swap(Tds & tds);

  void clear();

  void clear_cells_only();

  void set_adjacency(Cell_handle c0, Cell_handle c1, int i0, int i1) const
  {
      CGAL_triangulation_assertion(i0 >= 0 && i0 <= dimension());
      CGAL_triangulation_assertion(i1 >= 0 && i1 <= dimension());
      CGAL_triangulation_assertion(c0 != c1);
      c0->set_neighbor(i0,c1);
      c1->set_neighbor(i1,c0);
  }

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

private:

  Cell_container & cell_container()             { return _cell_container; }
  const Cell_container & cell_container() const { return _cell_container; }

  Vertex_container & vertex_container()             {return _vertex_container;}
  const Vertex_container & vertex_container() const {return _vertex_container;}

  // in dimension i, number of vertices >= i+2 
  // ( the boundary of a simplex in dimension i+1 has i+2 vertices )
  int _dimension;
  int _number_of_vertices;

  Cell_container   _cell_container;
  Vertex_container _vertex_container;

  // used by is-valid :
  bool count_vertices(int & i, bool verbose = false, int level = 0) const;
  // counts AND checks the validity
  bool count_facets(int & i, bool verbose = false, int level = 0) const;
  // counts but does not check
  bool count_edges(int & i, bool verbose = false, int level = 0) const;
  // counts but does not check
  bool count_cells(int & i, bool verbose = false, int level = 0) const;
  // counts AND checks the validity
};

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

  int n, d;
  is >> d >> n;
  tds.set_dimension(d);
  tds.set_number_of_vertices(n);

  if(n == 0)
    return is;

  std::map< int, Vertex_handle > V;
  
  // creation of the vertices    
  for (int i=0; i < n; i++) {
    //    is >> p;
    //    V[i] = tds.create_vertex();
    //    V[i]->set_point(p);
    V[i] = tds.create_vertex();
  }

  std::map< int, Cell_handle > C;
  int m;
 
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
  typedef typename Tds::Vertex_handle           Vertex_handle;
  typedef typename Tds::Vertex_iterator         Vertex_iterator;

  std::map<Vertex_handle, int> V;

  // outputs dimension and number of vertices
  int n = tds.number_of_vertices();

  if (is_ascii(os))
      os << tds.dimension() << std::endl << n << std::endl;
  else
      os << tds.dimension() << n;

  if (n == 0)
    return os;
  
  // index the vertices
  int i = 0;
  for (Vertex_iterator it=tds.vertices_begin(); it != tds.vertices_end(); ++it)
    V[&(*it)] = i++;

  CGAL_triangulation_assertion( i == n );

  tds.print_cells(os, V);

  return os;
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
is_vertex(Vertex_handle v) const
{
    return vertex_container().is_element(&*v);
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
is_edge(Vertex_handle u, Vertex_handle v, Cell_handle &c, int &i, int &j) const
  // returns false when dimension <1 or when indices wrong
{
  if (u==v)
      return false;

  for(Cell_iterator cit = cell_container().begin(); cit != cells_end(); ++cit)
    if (cit->has_vertex(u,i) && cit->has_vertex(v,j)) {
      c = &*cit;
      return true; 
    }
  return false;
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
is_edge(Cell_handle c, int i, int j) const
  // returns false when dimension <1
{
  if ( i==j ) return false;
  if ( (i<0) || (j<0) ) return false;
  if ( (dimension() == 1) && ((i>1) || (j>1)) ) return false;
  if ( (dimension() == 2) && ((i>2) || (j>2)) ) return false;
  if ((i>3) || (j>3)) return false;

  for(Cell_iterator cit = cell_container().begin(); cit != cells_end(); ++cit)
    if (&*cit == &*c)
	return true;
  return false;
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
is_facet(Vertex_handle u, Vertex_handle v, Vertex_handle w, 
	 Cell_handle & c, int & i, int & j, int & k) const
  // returns false when dimension <2 or when indices wrong
{
      if ( (u==v) || (u==w) || (v==w) ) return false;
      Facet_iterator it = facets_begin();
      while ( it != facets_end() ) {
	if ( ( ((*it).first)->has_vertex(u,i) )
	     && ( ((*it).first)->has_vertex(v,j) )
	     && ( ((*it).first)->has_vertex(w,k) ) ) {
	  c = (*it).first;
	  return true;
	}
	++it;
      }
      return false;
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
is_facet(Cell_handle c, int i) const
  // returns false when dimension <2
{
  CGAL_triangulation_precondition(i>=0 && i<4);
  if ( (dimension() == 2) && (i!=3) )
      return false;
  return cell_container().is_element(&*c);
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
is_cell( Cell_handle c ) const
  // returns false when dimension <3
{
  if (dimension() < 3)
      return false;
  return cell_container().is_element(&*c);
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
is_cell(Vertex_handle u, Vertex_handle v, Vertex_handle w, Vertex_handle t, 
	Cell_handle & c, int & i, int & j, int & k, int & l) const
  // returns false when dimension <3
{
  if ( (u==v) || (u==w) || (u==t) || (v==w) || (v==t) || (w==t) )
    return false;
  for(Cell_iterator it = cells_begin(); it != cells_end(); ++it) {
    if ( ( it->has_vertex(u,i) )
	 && ( it->has_vertex(v,j) )
	 && ( it->has_vertex(w,k) ) 
	 && ( it->has_vertex(t,l) ) ) {
      c = &*it;
      return true;
    }
  }
  return false;
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
is_cell(Vertex_handle u, Vertex_handle v, Vertex_handle w, Vertex_handle t)
    const
  // returns false when dimension <3
{
  if ( (u==v) || (u==w) || (u==t) || (v==w) || (v==t) || (w==t) )
    return false;
  for(Cell_iterator it = cells_begin(); it != cells_end(); ++it) {
    if ( it->has_vertex(u) &&
	 it->has_vertex(v) &&
	 it->has_vertex(w) &&
	 it->has_vertex(t) ) {
      return true;
    }
  }
  return false;
}

template < class Vb, class Cb>
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
bool
Triangulation_data_structure_3<Vb,Cb>::
has_vertex(const Facet & f, Vertex_handle v, int & j) const
{
  return has_vertex(f.first, f.second, v, j);
}

template < class Vb, class Cb>
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
				   && (number_of_vertices() > 6) );
  CGAL_triangulation_expensive_precondition( is_cell(c) );

  Cell_handle n = c->neighbor(i);
  int in = n->index(c);

  // checks that the facet is flippable,
  // ie the future edge does not already exist
  std::set<Vertex_handle> setc;
  incident_vertices( c->vertex(i), setc );
  if ( setc.find( n->vertex(in) ) != setc.end() ) return false;

  flip_really(c,i,n,in);
  return true;
}

template < class Vb, class Cb>
void
Triangulation_data_structure_3<Vb,Cb>::
flip_flippable( Cell_handle c, int i )
  // flips facet i of cell c
  // c will be replaced by one of the new cells
{
  CGAL_triangulation_precondition( (dimension() == 3) && (0<=i) && (i<4) 
				   && (number_of_vertices() > 6) );
  CGAL_triangulation_expensive_precondition( is_cell(c) );

  Cell_handle n = c->neighbor(i);
  int in = n->index(c);

  // checks that the facet is flippable,
  // ie the future edge does not already exist
  typedef std::set<Vertex_handle> set_of_vertices;
  CGAL_triangulation_expensive_precondition_code( set_of_vertices setc; );
  CGAL_triangulation_expensive_precondition_code
    ( incident_vertices( c->vertex(i), setc ); );
  CGAL_triangulation_expensive_precondition
    ( ( setc.find( n->vertex(in) ) == setc.end() ) );

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

  set_adjacency(c, n->neighbor(in3), i, n->neighbor(in3)->index(n));
  c->set_vertex( i3, n->vertex(in) );

  set_adjacency(n, c->neighbor(i1), in, c->neighbor(i1)->index(c));
  n->set_vertex( in1, c->vertex(i) );

  Cell_handle cnew = create_cell(c->vertex(i), c->vertex(i1),
			         n->vertex(in), n->vertex(in3));

  set_adjacency(cnew, n->neighbor(in2), 0, n->neighbor(in2)->index(n));
  set_adjacency(cnew, n, 1, in2);
  set_adjacency(cnew, c->neighbor(i2), 2, c->neighbor(i2)->index(c));
  set_adjacency(cnew, c, 3, i2);
  set_adjacency(c, n, i1, in3);

  if (i&1 != 0)
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
				   && (number_of_vertices() > 6) );
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
				   && (number_of_vertices() > 6) );
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
	     Cell_handle c1, Vertex_handle v1, int i1, int j1, int next1,
	     Cell_handle c2, Vertex_handle v2, int i2, int j2, int next2,
	     Vertex_handle v3 )
{
  c->vertex(i)->set_cell(c1);
  c->vertex(j)->set_cell(c2);

  c1->set_vertex( j1, v1 );
  v1->set_cell(c1);
  c2->set_vertex( i2, v2 );
  v2->set_cell(c2);

  set_adjacency(c1,c2->neighbor(j2), next1, c2->neighbor(j2)->index(c2));
  set_adjacency(c2,c1->neighbor(i1),c2->index(v1),c1->neighbor(i1)->index(c1));

  set_adjacency(c1, c2, i1, j2);

  set_adjacency(c1, c->neighbor(j), 6-i1-j1-next1, c->neighbor(j)->index(c));
  set_adjacency(c2, c->neighbor(i), next2, c->neighbor(i)->index(c));

  v3->set_cell( c2 );

  delete_cell( c );
}

template < class Vb, class Cb >
void
Triangulation_data_structure_3<Vb,Cb>::
read_cells(std::istream& is, std::map< int, Vertex_handle > &V,
	   int & m, std::map< int, Cell_handle > &C)
{
  // creation of the cells and neighbors
  switch (dimension()) {
  case 3:
  case 2:
  case 1:
    {
      is >> m;

      for(int i = 0; i < m; i++) {
	Cell_handle c = create_cell();
	for (int k=0; k<=dimension(); ++k) {
	    int ik;
	    is >> ik;
	    c->set_vertex(k, V[ik]);
	    V[ik]->set_cell(c);
	}
	C[i] = c;
      }
      for(int j = 0; j < m; j++) {
        Cell_handle c = C[j];
	for (int k=0; k<=dimension(); ++k) {
	    int ik;
	    is >> ik;
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
	Cell_handle c = create_cell(V[i], NULL, NULL, NULL);
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
      Cell_handle c = create_cell(V[0], NULL, NULL, NULL);
      C[0] = c;
      V[0]->set_cell(c);
      break;
    }
  }
}

template < class Vb, class Cb>
void
Triangulation_data_structure_3<Vb,Cb>::
print_cells(std::ostream& os, std::map<Vertex_handle, int> &V ) const
{
  std::map<Cell_handle, int > C;

  int i = 0;
  int j;
  int m;
  
  switch ( dimension() ) {
  case 3:
    {
      m = number_of_cells();
      os << m;
      if(is_ascii(os))
	  os << std::endl;

      // write the cells
      Cell_iterator it;
      for(it = cells_begin(); it != cells_end(); ++it) {
	C[&(*it)] = i++;
	for(j = 0; j < 4; j++){
	  os << V[it->vertex(j)];
	  if(is_ascii(os)) {
	    if ( j==3 )
	      os << std::endl;
	    else
	      os << ' ';
	  }
	}
      }
      CGAL_triangulation_assertion( i == m );
      
      // write the neighbors
      for(it = cells_begin(); it != cells_end(); ++it) {
	for (j = 0; j < 4; j++) {
	  os << C[it->neighbor(j)];
	  if(is_ascii(os)){
	    if(j==3)
	      os << std::endl;
	    else
	      os <<  ' ';
	  }
	}
      }
      break;
    }
  case 2:
    {
      m = number_of_facets();
      os << m;
      if(is_ascii(os))
	  os << std::endl;

      // write the facets
      Facet_iterator it;
      for(it = facets_begin(); it != facets_end(); ++it) {
	C[(*it).first] = i++;
	for(j = 0; j < 3; j++){
	  os << V[(*it).first->vertex(j)];
	  if(is_ascii(os)) {
	    if ( j==2 )
	      os << std::endl;
	    else
	      os <<  ' ';
	  }
	}
      }
      CGAL_triangulation_assertion( i == m );
      
      // write the neighbors
      for(it = facets_begin(); it != facets_end(); ++it) {
	for (j = 0; j < 3; j++) {
	  os << C[(*it).first->neighbor(j)];
	  if(is_ascii(os)){
	    if(j==2)
	      os << std::endl;
	    else
	      os <<  ' ';
	  }
	}
      }
      break;
    }
  case 1:
    {
      m = number_of_edges();
      os << m;
      if(is_ascii(os))
	  os << std::endl;

      // write the edges
      Edge_iterator it;
      for(it = edges_begin(); it != edges_end(); ++it) {
	C[(*it).first] = i++;
	for(j = 0; j < 2; j++){
	  os << V[(*it).first->vertex(j)];
	  if(is_ascii(os)) {
	    if ( j==1 )
	      os << std::endl;
	    else
	      os <<  ' ';
	  }
	}
      }
      CGAL_triangulation_assertion( i == m );
      
      // write the neighbors
      for(it = edges_begin(); it != edges_end(); ++it) {
	for (j = 0; j < 2; j++) {
	  os << C[(*it).first->neighbor(j)];
	  if(is_ascii(os)){
	    if(j==1)
	      os << std::endl;
	    else
	      os <<  ' ';
	  }
	}
      }
      break;
    }
  }
}

template <class Vb, class Cb >
Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
insert_in_cell( Vertex_handle v, Cell_handle c )
{
  CGAL_triangulation_precondition( dimension() == 3 );
  CGAL_triangulation_precondition( (c != NULL) );
  CGAL_triangulation_expensive_precondition( is_cell(c) );

  if ( v == NULL )
    v = create_vertex();

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

  set_adjacency(c3, c, 0, 3);
  set_adjacency(c2, c, 0, 2);
  set_adjacency(c1, c, 0, 1);

  set_adjacency(c2, c3, 3, 2);
  set_adjacency(c1, c3, 3, 1);
  set_adjacency(c1, c2, 2, 1);

  set_adjacency(n1, c1, n1->index(c), 1);
  set_adjacency(n2, c2, n2->index(c), 2);
  set_adjacency(n3, c3, n3->index(c), 3);

  c->set_vertex(0,v);

  v0->set_cell(c1);
  v->set_cell(c);
  set_number_of_vertices(number_of_vertices() +1);

  return v;
}

template <class Vb, class Cb >
Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
insert_in_facet(Vertex_handle v, Cell_handle c, int i)
{ // inserts v in the facet opposite to vertex i of cell c

  CGAL_triangulation_precondition( (c != NULL)); 
  CGAL_triangulation_precondition( dimension() >= 2 );

  if ( v == NULL )
    v = create_vertex();

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
      set_adjacency(cnew1, nc, 1, nc->index(c));
      set_adjacency(cnew1, c, 3, i1);

      v3->set_cell(cnew1);

      // new cell with v in place of i2
      nc = c->neighbor(i2);
      Cell_handle cnew2 = create_cell(vi,v1,v,v3);
      set_adjacency(cnew2, nc, 2, nc->index(c));
      set_adjacency(cnew2, c, 3, i2);
      set_adjacency(cnew1, cnew2, 2, 1);

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
      set_adjacency(dnew1, nd, 1, nd->index(d));
      set_adjacency(dnew1, d, 2, j1);
      set_adjacency(dnew1, cnew1, 0, 0);

      // new cell with v in place of j2
      nd = d->neighbor(j2);
      Cell_handle dnew2 = create_cell(d->vertex(j),v1,v3,v);

      set_adjacency(dnew2, nd, 3, nd->index(d));
      set_adjacency(dnew2, d, 2, j2);
      set_adjacency(dnew2, cnew2, 0, 0);
      set_adjacency(dnew1, dnew2, 3, 1);

      // v replaces i3 in d
      d->set_vertex(j3,v);
      v->set_cell(d);

      break;
    }
  case 2:
    {
      CGAL_triangulation_expensive_precondition( is_facet(c,i) );
      Cell_handle n = c->neighbor(2);
      Cell_handle cnew = create_cell(c->vertex(0),c->vertex(1),v,NULL);
      set_adjacency(cnew, n, 2, n->index(c));
      set_adjacency(cnew, c, 0, 2);
      c->vertex(0)->set_cell(cnew);

      n = c->neighbor(1);
      Cell_handle dnew = create_cell(c->vertex(0),v,c->vertex(2),NULL);
      set_adjacency(dnew, n, 1, n->index(c));
      set_adjacency(dnew, c, 0, 1);
      set_adjacency(dnew, cnew, 2, 1);

      c->set_vertex(0,v);
      v->set_cell(c);
      break;
    }
  }
  set_number_of_vertices(number_of_vertices() +1);

  return v;
}
// end insert_in_facet

template <class Vb, class Cb >
Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
insert_in_edge(Vertex_handle v, Cell_handle c, int i, int j)   
  // inserts v in the edge of cell c with vertices i and j
{ 
  CGAL_triangulation_precondition( c != NULL ); 
  CGAL_triangulation_precondition( i != j );
  CGAL_triangulation_precondition( dimension() >= 1 );

  if ( v == NULL )
    v = create_vertex();

  switch ( dimension() ) {
  case 3:
    {
      CGAL_triangulation_expensive_precondition( is_cell(c) );
      CGAL_triangulation_precondition( i>=0 && i<=3 && j>=0 && j<=3 );

      std::vector<Cell_handle > cells;
      std::vector<Facet> facets;
      cells.reserve(32);
      facets.reserve(64);
      const Vertex_handle vi=c->vertex(i);
      const Vertex_handle vj=c->vertex(j);
      Cell_circulator ccir = incident_cells(c, i, j);
      do {
	  Cell_handle cc = &*ccir;
	  cells.push_back(cc);
	  facets.push_back(Facet(cc, cc->index(vi)));
	  facets.push_back(Facet(cc, cc->index(vj)));
	  ++ccir;
      } while (&*ccir != &*c);

      star_hole_3(v, facets.begin(), facets.end(), cells.begin(), cells.end());
      break;
    }
  case 2:
    {
      CGAL_triangulation_expensive_precondition( is_edge(c,i,j) );
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
      set_adjacency(cnew, c, i, j);
      set_adjacency(cnew, nj, j, nj->index(c));

      nj = d->neighbor(jd);
      set_adjacency(dnew, d, id, jd);
      set_adjacency(dnew, nj, jd, nj->index(d));

      set_adjacency(cnew, dnew, k, kd);

      v->set_cell(cnew);
      break;
    }

  case 1:
    {
      CGAL_triangulation_expensive_precondition( is_edge(c,i,j) );
      Cell_handle cnew = create_cell(v,c->vertex(1),NULL,NULL);
      c->vertex(1)->set_cell(cnew);
      c->set_vertex(1,v);
      set_adjacency(cnew, c->neighbor(0), 0, 1);
      set_adjacency(cnew, c, 1, 0);

      v->set_cell(cnew); 
      break;
    }
  }
  set_number_of_vertices(number_of_vertices() +1);

  return v;
}// end insert_in_edge

template <class Vb, class Cb >
Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
insert_increase_dimension(Vertex_handle v, // new vertex
			  Vertex_handle star,
			  bool reorient) 
  // star = vertex from which we triangulate the facet of the
  // incremented dimension  
  // ( geometrically : star = infinite vertex )
  // = Null only used to insert the 1st vertex (dimension -2 to dimension -1)
  // changes the dimension
  // if (reorient) the orientation of the cells is modified
{
  CGAL_triangulation_precondition( dimension() < 3);

  if ( v == NULL ) 
    v = create_vertex();

  int dim = dimension();
  if (dim != -2) {
      CGAL_triangulation_precondition( star != NULL );
      // In this case, this precondition is not relatively expensive.
      CGAL_triangulation_precondition( is_vertex(star) ); 
  }

  // this is set now, so that it becomes allowed to reorient
  // new facets or cells by iterating on them (otherwise the
  // dimension is too small)
  set_number_of_vertices( number_of_vertices()+1 );
  set_dimension( dimension()+1 );

  switch ( dim ) {
  case -2:
      // insertion of the first vertex
      // ( geometrically : infinite vertex )
    {
      Cell_handle c = create_cell( v, NULL, NULL, NULL);
      v->set_cell(c);
      break;
    }

  case -1:
    // insertion of the second vertex
    // ( geometrically : first finite vertex )
    {
      Cell_handle d = create_cell( v, NULL, NULL, NULL);
      v->set_cell(d);
      set_adjacency(d, star->cell(), 0, 0);
      break;
    }

  case 0:
    // insertion of the third vertex
    // ( geometrically : second finite vertex )
    {
      Cell_handle c = star->cell();
      Cell_handle d = c->neighbor(0);

      if (reorient) {
	c->set_vertex(0,d->vertex(0));
	c->set_vertex(1,star);
	d->set_vertex(1,d->vertex(0));
	d->set_vertex(0,v);
	set_adjacency(c, d, 1, 0);
	Cell_handle e = create_cell( star, v, NULL, NULL);
	set_adjacency(e, d, 0, 1);
	set_adjacency(e, c, 1, 0);
      }
      else {
	c->set_vertex(1,d->vertex(0));
	d->set_vertex(1,v);
	d->set_neighbor(1,c);
	Cell_handle e = create_cell( v, star, NULL, NULL);
	set_adjacency(e, c, 0, 1);
	set_adjacency(e, d, 1, 0);
      }
	
      v->set_cell(d);
      break;
    }

  case 1:
    // general case : 4th vertex ( geometrically : 3rd finite vertex )
    // degenerate cases geometrically : 1st non collinear vertex
    {
      Cell_handle c = star->cell();
      int i = c->index(star); // i== 0 or 1
      int j = (1-i);
      Cell_handle d = c->neighbor(j);
	
      c->set_vertex(2,v);

      Cell_handle e = c->neighbor(i);
      Cell_handle cnew = c;
      Cell_handle enew=NULL;
	
      while( e != d ){
	enew = create_cell( );
	enew->set_vertex(i,e->vertex(j));
	enew->set_vertex(j,e->vertex(i));
	enew->set_vertex(2,star);
	  
	set_adjacency(enew, cnew, i, j);
	// false at the first iteration of the loop where it should
	// be neighbor 2 
	// it is corrected after the loop
	set_adjacency(enew, e, 2, 2);
	// neighbor j will be set during next iteration of the loop
	  
	e->set_vertex(2,v);

	e = e->neighbor(i);
	cnew = enew;
      }
	
      d->set_vertex(2,v);
      set_adjacency(enew, d, j, 2);

      // corrections for star->cell() :
      c = star->cell();
      c->set_neighbor(2,c->neighbor(i)->neighbor(2));
      c->set_neighbor(j,d);

      v->set_cell(d);

      if (reorient) {
	// reorientation of all the cells
	for (Facet_iterator fit = facets_begin(); fit != facets_end(); ++fit)
	  change_orientation((*fit).first);
      }
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

      v->set_cell(&(*it)); // ok since there is at least one ``cell''
      for(; it != cells_end(); ++it) {
	// Here we must be careful since we create_cells in a loop controlled
	// by an iterator.  So we first take care of the cells newly created
	// by the following test :
	if (it->neighbor(0) == NULL)
	  continue;
	it->set_vertex(3,v);
	if ( ! it->has_vertex(star) ) {
	  Cell_handle cnew = create_cell( it->vertex(0), it->vertex(2),
			            it->vertex(1), star);
	  set_adjacency(cnew, &*it, 3, 3);
	  new_cells.push_back(cnew);
	}
      }

      // traversal of the new cells only, to add missing neighbors
      for (typename std::vector<Cell_handle >::iterator ncit = new_cells.begin(); 
           ncit != new_cells.end(); ++ncit) {
	Cell_handle n = (*ncit)->neighbor(3); // opposite to star
	for ( int i=0; i<3; i++ ) {
	  int j;
	  if ( i==0 ) j=0;
	  else j=3-i; // vertex 1 and vertex 2 are always switched when
	  // creating a new cell (see above)
          Cell_handle c = n->neighbor(i)->neighbor(3);
	  if ( c != NULL ) {
	    // i.e. star is not a vertex of n->neighbor(i)
	    (*ncit)->set_neighbor(j, c);
	    // opposite relation will be set when ncit arrives on c
	    // this avoids to look for the correct index 
	    // and to test whether *ncit already has neighbor i
	  }
	  else {
	    // star is a vertex of n->neighbor(i)
	    set_adjacency(*ncit, n->neighbor(i), j, 3);//neighbor opposite to v
	  }
	}
      }
	
      // reorientation of all the cells
      if (reorient)
	  for(it = cells_begin(); it != cells_end(); ++it)
	      change_orientation(&*it);
    }
  }// end switch
    
  return v;
}

template <class Vb, class Cb >
void
Triangulation_data_structure_3<Vb,Cb>::
incident_cells(Vertex_handle v, std::set<Cell_handle> & cells, Cell_handle c) const
{
  CGAL_triangulation_precondition( v != NULL );
  CGAL_triangulation_expensive_precondition( is_vertex(v) );

  if ( dimension() < 3 )
      return;

  if ( c == NULL )
    c = v->cell();
  else
    CGAL_triangulation_precondition( c->has_vertex(v) );

  if ( cells.find( c ) != cells.end() )
    return; // c was already found

  cells.insert( c );
      
  for ( int j=0; j<4; j++ )
    if ( j != c->index(v) )
      incident_cells( v, cells, c->neighbor(j) );
}
  
template <class Vb, class Cb >
void
Triangulation_data_structure_3<Vb,Cb>::
incident_vertices(Vertex_handle v, std::set<Vertex_handle> & vertices, Cell_handle c) const
{
  CGAL_triangulation_precondition( v != NULL );
  CGAL_triangulation_expensive_precondition( is_vertex(v) );
      
  if ( number_of_vertices() < 2 )
      return;

  if ( c == NULL )
    c = v->cell();
  else
    CGAL_triangulation_precondition( c->has_vertex(v) );

  int d = dimension();
  int j;
  int found = 0;
  for ( j=0; j <= d; j++ ) {
    if ( j != c->index(v) ) {
      if ( vertices.find( c->vertex(j) ) == vertices.end() )
	vertices.insert( c->vertex(j) );
      else
	found++; // c->vertex(j) was already found 
    }
  }
  if ( found == 3 )
      return; // c was already visited
      
  for ( j=0; j <= d; j++ )
    if ( j != c->index(v) )
      incident_vertices( v, vertices, c->neighbor(j) );
}

template <class Vb, class Cb >
bool
Triangulation_data_structure_3<Vb,Cb>::
is_valid(bool verbose, int level ) const
{
  switch ( dimension() ) {
  case 3:
    {
      int vertex_count;
      if ( ! count_vertices(vertex_count,verbose,level) )
	  return false;
      if ( number_of_vertices() != vertex_count ) {
	if (verbose)
	    std::cerr << "wrong number of vertices" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }

      int cell_count;
      if ( ! count_cells(cell_count,verbose,level) )
	  return false;
      int edge_count;
      if ( ! count_edges(edge_count,verbose,level) )
	  return false;
      int facet_count;
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
      int vertex_count;
      if ( ! count_vertices(vertex_count,verbose,level) )
	  return false;
      if ( number_of_vertices() != vertex_count ) {
	if (verbose)
	    std::cerr << "false number of vertices" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }

      int edge_count;
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

      int facet_count;
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
      int vertex_count;
      if ( ! count_vertices(vertex_count,verbose,level) )
	  return false;
      if ( number_of_vertices() != vertex_count ) {
	if (verbose)
	    std::cerr << "false number of vertices" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      int edge_count;
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
      int vertex_count;
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
Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
copy_tds(const Tds & tds, Vertex_handle vert )
  // returns the new vertex corresponding to vert in the new tds 
{
  CGAL_triangulation_expensive_precondition( vert == NULL
	                                  || tds.is_vertex(vert) );

  clear();

  int n = tds.number_of_vertices();
  set_number_of_vertices(n);
  set_dimension(tds.dimension());

  if (n == 0)
      return vert;

  // Create the vertices.
  // the vertices must be indexed by their order of creation so
  // that when reread from file, the orders of vertices are the
  // same - important for remove 
  std::vector<Vertex_handle> TV(n);
  int i = 0;

  for (Vertex_iterator vit = tds.vertices_begin();
       vit != tds.vertices_end(); ++vit)
    TV[i++] = &*vit; 
  
  CGAL_triangulation_assertion( i == n ); 
  std::sort(TV.begin(), TV.end(), 
	    Vertex_tds_compare_order_of_creation<Vertex_handle>()); 

  std::map< Vertex_handle, Vertex_handle > V;
  std::map< Cell_handle, Cell_handle > F;

  for (i=0; i <= n-1; i++) {
    V[ TV[i] ] = create_vertex();
    *V[ TV[i] ] = *TV[i];
  }

  // Create the cells.
  for (Cell_iterator cit = tds.cell_container().begin();
	  cit != tds.cells_end(); ++cit) {
      F[&(*cit)] = create_cell(&*cit);
      F[&(*cit)]->set_vertices(V[cit->vertex(0)],
			       V[cit->vertex(1)],
			       V[cit->vertex(2)],
			       V[cit->vertex(3)]);
  }

  // Link the vertices to a cell.
  for (Vertex_iterator vit2 = tds.vertices_begin();
       vit2 != tds.vertices_end(); ++vit2)
    V[&(*vit2)]->set_cell( F[vit2->cell()] );

  // Hook neighbor pointers of the cells.
  for (Cell_iterator cit2 = tds.cell_container().begin();
	  cit2 != tds.cells_end(); ++cit2) {
    for (int j = 0; j < 4; j++)
      F[&(*cit2)]->set_neighbor(j, F[cit2->neighbor(j)] );
  }

  CGAL_triangulation_postcondition( is_valid() );

  return (vert != NULL) ? V[vert] : (Vertex_handle) NULL;
}

template <class Vb, class Cb >
void
Triangulation_data_structure_3<Vb,Cb>::
swap(Tds & tds)
{
  // tds and *this are supposed to be valid
  std::swap(_dimension, tds._dimension);
  std::swap(_number_of_vertices, tds._number_of_vertices);
  cell_container().swap(tds.cell_container());
  vertex_container().swap(tds.vertex_container());
}

template <class Vb, class Cb >
void
Triangulation_data_structure_3<Vb,Cb>::
clear()
{
  cell_container().clear();
  vertex_container().clear();

  set_number_of_vertices(0);
  set_dimension(-2);
}

template <class Vb, class Cb >
void
Triangulation_data_structure_3<Vb,Cb>::
clear_cells_only()
{
  cell_container().clear();
}

template <class Vb, class Cb >
bool
Triangulation_data_structure_3<Vb,Cb>::
count_vertices(int & i, bool verbose, int level) const
  // counts AND checks the validity
{
  i = 0;

  for (Vertex_iterator it = vertices_begin(); it != vertices_end(); ++it) {
    if ( ! it->is_valid(verbose,level) ) {
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
count_facets(int & i, bool verbose, int level) const
  // counts but does not check
{
  i = 0;

  for (Facet_iterator it = facets_begin(); it != facets_end(); ++it) {
    if ( ! (*it).first->is_valid(dimension(),verbose, level) ) {
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
count_edges(int & i, bool verbose, int level) const
  // counts but does not check
{
  i = 0;

  for (Edge_iterator it = edges_begin(); it != edges_end(); ++it) {
    if ( ! (*it).first->is_valid(dimension(),verbose, level) ) {
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
count_cells(int & i, bool verbose, int level) const
  // counts AND checks the validity
{
  i = 0;

  for (Cell_iterator it = cells_begin(); it != cells_end(); ++it) {
    if ( ! it->is_valid(dimension(),verbose, level) ) {
      if (verbose)
	  std::cerr << "invalid cell" << std::endl;
      CGAL_triangulation_assertion(false);
      return false;
    }
    ++i;
  }
  return true;
}

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_DATA_STRUCTURE_3_H
