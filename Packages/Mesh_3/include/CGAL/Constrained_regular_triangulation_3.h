#ifndef CGAL_CONSTRAINED_REGULAR_TRIANGULATION_3_H
#define CGAL_CONSTRAINED_REGULAR_TRIANGULATION_3_H

#include <iostream>
#include <CGAL/IO/File_header_OFF.h>
#include <CGAL/IO/File_scanner_OFF.h>

CGAL_BEGIN_NAMESPACE

template <class Tr>
class Constrained_regular_triangulation_3 : public Tr
{
public:
  typedef Constrained_regular_triangulation_3 Self;
  typedef Tr Triangulation;
  typedef typename Triangulation::Geom_traits Geom_traits;

  typedef typename Triangulation::Vertex_handle Vertex_handle;
  typedef typename Triangulation::Cell_handle Cell_handle;
  typedef typename Geom_traits::Bare_point Bare_point;
  typedef typename Triangulation::Weighted_point Weighted_point;
  typedef Weighted_point Point;
  typedef typename Triangulation::Locate_type Locate_type;

  typedef typename Triangulation::Facet Facet;


  Vertex_handle off_file_input( std::istream& is, bool verbose = false);

  // --- CONSTUCTORS ---
  Constrained_regular_triangulation_3(const Geom_traits& gt = Geom_traits())
    : Triangulation(gt)
  {}


  // --- INSERTIONS ---

  // -- points insertions --
  Vertex_handle insert(const Weighted_point & p, Cell_handle start = NULL);

  Vertex_handle insert(const Weighted_point & p, Locate_type lt,
	               Cell_handle c, int li, int);

  bool insert_constrained_edge(const Vertex_handle& va,
			       const Vertex_handle& vb)
  {
    Cell_handle ch;
    int i, j;

    if( !is_edge(va, vb, ch, i, j) )
      return false; /** \todo Conform the edge. */
    else
      {
	va->set_is_adjacent_by_constraint(vb, true);
	vb->set_is_adjacent_by_constraint(va, true);
      }
    return true;
  }

  bool insert_constrained_facet(const Vertex_handle& va,
				const Vertex_handle& vb,
				const Vertex_handle& vc)
  {
    Cell_handle c;
    int i, j, k;
    if( !is_facet(va, vb, vc,
		  c, i, j, k) )
      return false; /** \todo Force the facet into the triangulation. */
    else
      {
	const int l = 6-i-j-k;
	const Cell_handle& n = c->neighbor(l);
	c->set_constrained(l, true);
	n->set_constrained(n->index(c), true);
	return true;
      }
  }
private:

  Oriented_side
  power_test(const Weighted_point &p, const Weighted_point &q) const
  {
      CGAL_precondition(equal(p, q));
      return geom_traits().power_test_3_object()(p, q);
  }

  Oriented_side
  power_test(const Weighted_point &p, const Weighted_point &q,
	     const Weighted_point &r) const
  {
      CGAL_precondition(collinear(p, q, r));
      return geom_traits().power_test_3_object()(p, q, r);
  }

  Oriented_side
  power_test(const Weighted_point &p, const Weighted_point &q,
	     const Weighted_point &r, const Weighted_point &s) const
  {
      CGAL_precondition(coplanar(p, q, r, s));
      return geom_traits().power_test_3_object()(p, q, r, s);
  }

  Oriented_side
  power_test(const Weighted_point &p, const Weighted_point &q,
	     const Weighted_point &r, const Weighted_point &s,
	     const Weighted_point &t) const
  {
      return geom_traits().power_test_3_object()(p, q, r, s, t);
  }

  bool in_conflict_3(const Weighted_point &p, const Cell_handle c) const
  {
      return side_of_power_sphere(c, p) == ON_BOUNDED_SIDE;
  }

  bool in_conflict_2(const Weighted_point &p, const Cell_handle c, int i) const
  {
      return side_of_power_circle(c, i, p) == ON_BOUNDED_SIDE;
  }

  bool in_conflict_1(const Weighted_point &p, const Cell_handle c) const
  {
      return side_of_power_segment(c, p) == ON_BOUNDED_SIDE;
  }

  bool in_conflict_0(const Weighted_point &p, const Cell_handle c) const
  {
      return power_test(c->vertex(0)->point(), p) == ON_POSITIVE_SIDE;
  }

  class Conflict_tester_3
  {
      const Weighted_point &p;
      const Self *t;
      mutable std::vector<Vertex_handle> cv;

  public:

      Conflict_tester_3(const Weighted_point &pt, const Self *tr)
	  : p(pt), t(tr) {}

      bool operator()(const Cell_handle c) const
      {
	  // We mark the vertices so that we can find the deleted ones easily.
	  if (t->in_conflict_3(p, c))
	  {
	      for (int i=0; i<4; i++)
	      {
		  Vertex_handle v = c->vertex(i);
		  if (v->cell() != NULL)
		  {
		      cv.push_back(v);
		      v->set_cell(NULL);
		  }
	      }
	      return true;
	  }
	  return false;
      }

      std::vector<Vertex_handle> & conflict_vector()
      {
	  return cv;
      }
  };

  class Conflict_tester_2
  {
      const Weighted_point &p;
      const Self *t;
      mutable std::vector<Vertex_handle> cv;

  public:

      Conflict_tester_2(const Weighted_point &pt, const Self *tr)
	  : p(pt), t(tr) {}

      bool operator()(const Cell_handle c) const
      {
	  if (t->in_conflict_2(p, c, 3))
	  {
	      for (int i=0; i<3; i++)
	      {
		  Vertex_handle v = c->vertex(i);
		  if (v->cell() != NULL)
		  {
		      cv.push_back(v);
		      v->set_cell(NULL);
		  }
	      }
	      return true;
	  }
	  return false;
      }

      std::vector<Vertex_handle> & conflict_vector()
      {
	  return cv;
      }
  };

  friend class Conflict_tester_3;
  friend class Conflict_tester_2;


  template <class Conflict_test,
            class OutputIteratorBoundaryFacets,
            class OutputIteratorCells,
            class OutputIteratorInternalFacets>
  Triple<OutputIteratorBoundaryFacets,
         OutputIteratorCells,
         OutputIteratorInternalFacets>
  find_conflicts_2(Cell_handle c, const Conflict_test &tester,
	           Triple<OutputIteratorBoundaryFacets,
                          OutputIteratorCells,
		          OutputIteratorInternalFacets> it) const
  {
    CGAL_triangulation_precondition( dimension()==2 );
    CGAL_triangulation_precondition( tester(c) );

    c->set_in_conflict_flag(1);
    *it.second++ = c;

    for (int i=0; i<3; ++i) {
      Cell_handle test = c->neighbor(i);
      if (test->get_in_conflict_flag() == 1) {
	  if (c < test)
	      *it.third++ = Facet(c, i); // Internal facet.
          continue; // test was already in conflict.
      }
      if (test->get_in_conflict_flag() == 0) {
	  if (tester(test)) {
	      if (c < test)
		  *it.third++ = Facet(c, i); // Internal facet.
              it = find_conflicts_2(test, tester, it);
	      continue;
	  }
	  test->set_in_conflict_flag(2); // test is on the boundary.
      }
      *it.first++ = Facet(c, i);
    }
    return it;
  }

  // Note: the code duplication between _2 and _3 should be avoided one day.
  template <class Conflict_test,
            class OutputIteratorBoundaryFacets,
            class OutputIteratorCells,
            class OutputIteratorInternalFacets>
  Triple<OutputIteratorBoundaryFacets,
         OutputIteratorCells,
         OutputIteratorInternalFacets>
  find_conflicts_3(Cell_handle c, const Conflict_test &tester,
	           Triple<OutputIteratorBoundaryFacets,
                          OutputIteratorCells,
		          OutputIteratorInternalFacets> it) const
  {
    CGAL_triangulation_precondition( dimension()==3 );
    CGAL_triangulation_precondition( tester(c) );

    c->set_in_conflict_flag(1);
    *it.second++ = c;

    for (int i=0; i<4; ++i) {
      Cell_handle test = c->neighbor(i);
      if (test->get_in_conflict_flag() == 1) { // test was already in conflict.
	  if (c < test)
	      *it.third++ = Facet(c, i); // Internal facet.
          continue;
      }
      if ( /*! c->is_constrained(i) && */test->get_in_conflict_flag() == 0) {
	  if (tester(test)) {
	      if (c < test)
		  *it.third++ = Facet(c, i); // Internal facet.
              it = find_conflicts_3(test, tester, it);
	      continue;
	  }
	  test->set_in_conflict_flag(2); // test is on the boundary.
      }
      *it.first++ = Facet(c, i);
    }
    return it;
  }

  // This one takes a function object to recursively determine the cells in
  // conflict, then calls _tds._insert_in_hole().
  template < class Conflict_test >
  Vertex_handle
  insert_conflict_2(Cell_handle c, const Conflict_test &tester)
  {
    CGAL_triangulation_precondition( dimension() == 2 );
    CGAL_triangulation_precondition( c != NULL );
    CGAL_triangulation_precondition( tester(c) );

    std::vector<Cell_handle> cells;
    cells.reserve(32);

    Facet facet;

    // Find the cells in conflict
    find_conflicts_2(c, tester, make_triple(Oneset_iterator<Facet>(facet),
		                            std::back_inserter(cells),
				            Emptyset_iterator()));

    // Create the new cells and delete the old.
    return _tds._insert_in_hole(cells.begin(), cells.end(),
	                        facet.first, facet.second);
  }

  // This one takes a function object to recursively determine the cells in
  // conflict, then calls _tds._insert_in_hole().
  template < class Conflict_test >
  Vertex_handle
  insert_conflict_3(Cell_handle c, const Conflict_test &tester)
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    CGAL_triangulation_precondition( c != NULL );
    CGAL_triangulation_precondition( tester(c) );

    std::vector<Cell_handle> cells;
    cells.reserve(32);

    Facet facet;

    // Find the cells in conflict
    find_conflicts_3(c, tester, make_triple(Oneset_iterator<Facet>(facet),
		                            std::back_inserter(cells),
				            Emptyset_iterator()));

    // Create the new cells and delete the old.
    return _tds._insert_in_hole(cells.begin(), cells.end(),
	                        facet.first, facet.second);
  }

};

template <class Tr>
typename Constrained_regular_triangulation_3<Tr>::Vertex_handle
Constrained_regular_triangulation_3<Tr>::
insert(const Weighted_point & p, Cell_handle start) 
{
    Locate_type lt;
    int li, lj;
    Cell_handle c = locate(p, lt, li, lj, start);
    return insert(p, lt, c, li, lj);
}

template <class Tr>
typename Constrained_regular_triangulation_3<Tr>::Vertex_handle
Constrained_regular_triangulation_3<Tr>::
insert(const Weighted_point & p, Locate_type lt, Cell_handle c, int li, int)
{
  switch (dimension()) {
  case 3:
    {
      // TODO :
      // In case the point is completely equal (including weight), then we need
      // to discard it (don't update the triangulation, nor hide it), right ?
      if (! in_conflict_3(p, c)) {  // new point is hidden
          if (lt == Tr::VERTEX)
              return c->vertex(li); // by coinciding point
          else
              return NULL;          // by cell
      }

      // Should I mark c's vertices too ?
      Conflict_tester_3 tester(p, this);
      Vertex_handle v = insert_conflict_3(c, tester);
      v->set_point(p);
      for( typename std::vector<Vertex_handle>::iterator
		it = tester.conflict_vector().begin();
		it != tester.conflict_vector().end(); ++it)
      {
        if ((*it)->cell() == NULL)
	{
          // vertex has to be deleted
          tds().delete_vertex(*it);
	}
      }
      // TODO : manage the hidden points.
      return v;
    }
  case 2:
    {
      switch (lt) {
      case Tr::OUTSIDE_CONVEX_HULL:
      case Tr::CELL:
      case Tr::FACET:
      case Tr::EDGE:
      case Tr::VERTEX:
	{
          if (! in_conflict_2(p, c, 3)) {  // new point is hidden
              if (lt == Tr::VERTEX)
                  return c->vertex(li); // by coinciding point
              else
                  return NULL;          // by face
          }

	  Conflict_tester_2 tester(p, this);
	  Vertex_handle v = insert_conflict_2(c, tester);
	  v->set_point(p);
          for( typename std::vector<Vertex_handle>::iterator
		it = tester.conflict_vector().begin();
		it != tester.conflict_vector().end(); ++it)
	  {
            if ((*it)->cell() == NULL)
	    {
              // vertex has to be deleted
              tds().delete_vertex(*it);
	    }
	  }
	  return v;
	}
      case Tr::OUTSIDE_AFFINE_HULL:
	{
	  // if the 2d triangulation is Regular, the 3d
	  // triangulation will be Regular
	  return Tr::insert_outside_affine_hull(p);
	}
      }
    }//dim 2
  case 1:
    {
      switch (lt) {
      case Tr::OUTSIDE_CONVEX_HULL:
      case Tr::EDGE:
      case Tr::VERTEX:
	{
          if (! in_conflict_1(p, c)) {  // new point is hidden
              if (lt == Tr::VERTEX)
                  return c->vertex(li); // by coinciding point
              else
                  return NULL;          // by edge
          }

	  Cell_handle bound[2];
          // corresponding index: bound[j]->neighbor(1-j) is in conflict.
	  std::vector<Vertex_handle>  hidden_vertices;
	  std::vector<Cell_handle>    conflicts;
          conflicts.push_back(c);

          // We get all cells in conflict,
          // and remember the 2 external boundaries.

	  for (int j = 0; j<2; ++j) {
	    Cell_handle n = c->neighbor(j);
	    while ( in_conflict_1( p, n) ) {
	      conflicts.push_back(n);
              hidden_vertices.push_back(n->vertex(j));
	      n = n->neighbor(j);
	    }
	    bound[j] = n;
	  }

          // We preserve the order (like the orientation in 2D-3D).

	  Vertex_handle v = tds().create_vertex();
	  v->set_point(p);
          Cell_handle c0 = tds().create_face(v, bound[0]->vertex(0), NULL);
          Cell_handle c1 = tds().create_face(bound[1]->vertex(1), v, NULL);
          tds().set_adjacency(c0, 1, c1, 0);
          tds().set_adjacency(bound[0], 1, c0, 0);
          tds().set_adjacency(c1, 1, bound[1], 0);
          bound[0]->vertex(0)->set_cell(bound[0]);
          bound[1]->vertex(1)->set_cell(bound[1]);
          v->set_cell(c0);

	  tds().delete_cells(conflicts.begin(), conflicts.end());
	  tds().delete_vertices(hidden_vertices.begin(), hidden_vertices.end());
	  return v;
	}
      case Tr::OUTSIDE_AFFINE_HULL:
	return Tr::insert_outside_affine_hull(p);
      case Tr::FACET:
      case Tr::CELL:
	// impossible in dimension 1
        CGAL_assertion(false);
	return NULL;
      }
    }
  case 0:
    {
        // We need to compare the weights when the points are equal.
        if (lt == Tr::VERTEX && in_conflict_0(p, c)) {
            CGAL_assertion(li == 0);
            c->vertex(li)->set_point(p); // replace by heavier point
        }
        else
            return Tr::insert(p, c);
    }
  default :
    {
      return Tr::insert(p, c);
    }
  }
}

template <class Tr>
typename Tr::Vertex_handle
Constrained_regular_triangulation_3<Tr>::
off_file_input(std::istream& is, bool verbose)
{
  Vertex_handle vinf(0);
  File_scanner_OFF scanner(is, verbose);
  if (! is) {
    if (scanner.verbose()) {
         std::cerr 
	   << " " << std::endl
	   << "Constrained_regular_triangulation_3::off_file_input"
	   << std::endl
	   << " input error: file format is not OFF." << std::endl;
    }
    return vinf;
  }
  
  clear();
  
  std::vector<Vertex_handle> vvh(scanner.size_of_vertices());
  //  std::map<Vh_pair, Edge> edge_map;

  // insert points
  int i;
  for ( i = 0; i < scanner.size_of_vertices(); i++) {
    Bare_point p;
    file_scan_vertex( scanner, p);
    vvh[i] = insert(p); // insert the point in the triangulation
    scanner.skip_to_next_vertex( i);
  }

  if ( ! is ) {
    is.clear( std::ios::badbit);
    return vinf;
  }

  // inserts constrained edges and facets
  for ( i = 0; i < scanner.size_of_facets(); i++) {
    Integer32 no;
    scanner.scan_facet( no, i);
    if( ! is || (no != 3 && no!= 2) ) {
      if ( scanner.verbose()) {
	std::cerr 
	  << " " << std::endl
	  << "Constrained_regular_triangulation_3::off_file_input"
	  << "edge or facet " << i << "does not have 2 or 3 vertices." 
	  << std::endl;
      }
      is.clear( std::ios::badbit);
      return vinf;
    }

    Integer32 index0;
    Integer32 index1;
    scanner.scan_facet_vertex_index( index0, i);
    scanner.scan_facet_vertex_index( index1, i);
    if( no == 3 ) // facet
      {
	Integer32 index2;
	scanner.scan_facet_vertex_index( index2, i);
	bool r = insert_constrained_facet(vvh[index0], vvh[index1],
					  vvh[index2]);
	CGAL_assertion(r);
      }
    else // edge
      {
	bool r = insert_constrained_edge(vvh[index0], vvh[index1]);
	CGAL_assertion(r);
      }
  }
  return vinf;
}

CGAL_END_NAMESPACE

#endif // CGAL_CONSTRAINED_REGULAR_TRIANGULATION_3_H
