#ifndef CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_3_H
#define CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_3_H

#include <iostream>
#include <CGAL/IO/File_header_OFF.h>
#include <CGAL/IO/File_scanner_OFF.h>

CGAL_BEGIN_NAMESPACE

template <class Tr>
class Constrained_Delaunay_triangulation_3 : public Tr
{
public:
  typedef Constrained_Delaunay_triangulation_3 Self;
  typedef Tr Triangulation;
  typedef typename Triangulation::Geom_traits Geom_traits;

  typedef typename Triangulation::Vertex_handle Vertex_handle;
  typedef typename Triangulation::Cell_handle Cell_handle;

  Vertex_handle off_file_input( std::istream& is, bool verbose = false);
};

template <class Tr>
typename Tr::Vertex_handle
Constrained_Delaunay_triangulation_3<Tr>::
off_file_input(std::istream& is, bool verbose)
{
  Vertex_handle vinf;
  File_header_OFF scanner(is, verbose);
  if (! is) {
    if (scanner.verbose()) {
         std::cerr 
	   << " " << std::endl
	   << "Constrained_Delaunay_triangulation_3::off_file_input"
	   << std::endl
	   << " input error: file format is not OFF." << std::endl;
    }
    return vinf;
  }
  
  clear();
  
  std::vector<Vertex_handle > vvh(scanner.size_of_vertices());
  std::map<Vh_pair, Edge> edge_map;

  // insert points
  for ( i = 0; i < scanner.size_of_vertices(); i++) {
    Point p;
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
	  << "Constrained_Delaunay_triangulation_3::off_file_input"
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
	insert_constrained_facet(vvh[index0], vvh[index1],
				 vvh[index2]);
      }
    else // edge
      insert_constrained_edge(vvh[index0], vvh[index1]);
  }
}

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
      const l = 6-i-j-k;
      const Cell_handle& n = c->neighbor(l);
      c->set_constraint(l, true);
      n->set_constraint(n->index(c), true);
    }
}

CGAL_END_NAMESPACE

#endif // CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_3_H
