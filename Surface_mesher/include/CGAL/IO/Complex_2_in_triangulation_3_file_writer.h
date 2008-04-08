// Copyright (c) 2003-2007  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
//
//
// Author(s)     : Steve Oudot, Laurent Rineau

#define CGAL_C2T3_USE_FILE_WRITER_OFF

#ifndef CGAL_IO_COMPLEX_2_IN_TRIANGULATION_3_FILE_WRITER_H
#define CGAL_IO_COMPLEX_2_IN_TRIANGULATION_3_FILE_WRITER_H

#ifdef CGAL_C2T3_USE_FILE_WRITER_OFF
#  include <CGAL/IO/File_writer_OFF.h>
#endif

#ifdef CGAL_C2T3_USE_POLYHEDRON
#  include <CGAL/Polyhedron_3.h>
#  include <CGAL/Polyhedron_incremental_builder_3.h>
#endif

#include <iomanip>
#include <stack>

namespace Surface_mesher_io {

// using namespace ::CGAL::Surface_mesher;

template <class C2t3>
void
output_surface_facets_to_off (std::ostream& os, const C2t3& c2t3)
{
  typedef typename C2t3::Triangulation Tr;
  typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Facet Facet;
  typedef typename Tr::Edge Edge;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Point Point;
  typedef typename Tr::Geom_traits Gt;
#ifdef CGAL_C2T3_USE_FILE_WRITER_OFF
  CGAL::File_writer_OFF off;
  const Tr& tr = c2t3.triangulation();
  const typename Tr::size_type number_of_facets = c2t3.number_of_facets();

  off.header().set_no_comments(true);
  off.write_header(os,
		   tr.number_of_vertices(),
		   0, // fake number of halfedges, not used.
		   number_of_facets);

  std::map<Vertex_handle, int> V;
  int inum = 0;
  for(Finite_vertices_iterator vit = tr.finite_vertices_begin();
      vit != tr.finite_vertices_end();
      ++vit)
  {
    V[vit] = inum++;
    Point p = static_cast<Point>(vit->point());
    off.write_vertex(p.x(), p.y(), p.z());
  }

  off.write_facet_header();

  Finite_facets_iterator fit = tr.finite_facets_begin();
  std::set<Facet> oriented_set;
  std::stack<Facet> stack;

  CGAL_assertion_code(typename Tr::size_type nb_facets = 0; )

    while (oriented_set.size() != number_of_facets) {
      while ( fit->first->is_facet_on_surface(fit->second) == false ||
	      oriented_set.find(*fit) != oriented_set.end() ||

	      oriented_set.find(c2t3.opposite_facet(*fit)) !=
	      oriented_set.end() ) {
	++fit;
      }
      oriented_set.insert(*fit);
      stack.push(*fit);
      while(! stack.empty() ) {
	Facet f = stack.top();
	stack.pop();
	for(int ih = 0 ; ih < 3 ; ++ih) {
	  const int i1  = tr.vertex_triple_index(f.second, tr. cw(ih));
	  const int i2  = tr.vertex_triple_index(f.second, tr.ccw(ih));
	  if( c2t3.face_status(Edge(f.first, i1, i2)) == C2t3::REGULAR ) {
	    Facet fn = c2t3.neighbor(f, ih);
	    if (oriented_set.find(fn) == oriented_set.end() &&
		oriented_set.find(c2t3.opposite_facet(fn)) == oriented_set.end())
	    {
	      oriented_set.insert(fn);
	      stack.push(fn);
	    }
	  } // end "if the edge is regular"
	} // end "for each neighbor of f"
      } // end "stack non empty"
    } // end "oriented_set not full"

  for(typename std::set<Facet>::const_iterator fit = 
	oriented_set.begin();
      fit != oriented_set.end();
      ++fit)
  {
    off.write_facet_begin(3);
    for (int i=0; i<3; i++)
      off.write_facet_vertex_index(V[fit->first->vertex(tr.vertex_triple_index(fit->second, i))]);
    CGAL_assertion_code(++nb_facets);
    off.write_facet_end();
  }

  CGAL_assertion(nb_facets == number_of_facets);
  off.write_footer();

#elif defined(CGAL_C2T3_USE_POLYHEDRON)
  typedef CGAL::Polyhedron_3<Gt> Polyhedron_3;
  typedef typename Polyhedron_3::HalfedgeDS HalfedgeDS;

  Polyhedron_3 p;

  struct Off_builder : public CGAL::Modifier_base<HalfedgeDS>{
    const C2t3& c2t3;
    const Tr& tr;
    Off_builder(const C2t3& c2t3) : c2t3(c2t3), tr(c2t3.triangulation()) {};
    void operator()( HalfedgeDS& hds) {
      CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> builder(hds, true);
      const typename Tr::size_type number_of_facets = c2t3.number_of_facets();
      builder.begin_surface(tr.number_of_vertices(), 
			    number_of_facets);
      {
	// Finite vertices coordinates.
	std::map<Vertex_handle, int> V;
	int inum = 0;
	for(Finite_vertices_iterator vit = tr.finite_vertices_begin();
	    vit != tr.finite_vertices_end();
	    ++vit)
	{
	  V[vit] = inum++;
	  Point p = static_cast<Point>(vit->point());
	  builder.add_vertex(p);
	}
	Finite_facets_iterator fit = tr.finite_facets_begin();
	std::set<Facet> oriented_set;
	std::stack<Facet> stack;

	CGAL_assertion_code(typename Tr::size_type nb_facets = 0; )

	while (oriented_set.size() != number_of_facets) {
	  while ( fit->first->is_facet_on_surface(fit->second) == false ||
		  oriented_set.find(*fit) != oriented_set.end() ||

		  oriented_set.find(c2t3.opposite_facet(*fit)) !=
		  oriented_set.end() ) {
	    ++fit;
	  }
	  oriented_set.insert(*fit);
	  stack.push(*fit);
	  while(! stack.empty() ) {
	    Facet f = stack.top();
	    stack.pop();
	    for(int ih = 0 ; ih < 3 ; ++ih) {
	      const int i1  = tr.vertex_triple_index(f.second, tr. cw(ih));
	      const int i2  = tr.vertex_triple_index(f.second, tr.ccw(ih));
	      if( c2t3.face_status(Edge(f.first, i1, i2)) == C2t3::REGULAR ) {
		Facet fn = c2t3.neighbor(f, ih);
		if (oriented_set.find(fn) == oriented_set.end() &&
		    oriented_set.find(c2t3.opposite_facet(fn)) == oriented_set.end())
		{
		  oriented_set.insert(fn);
		  stack.push(fn);
		}
	      } // end "if the edge is regular"
	    } // end "for each neighbor of f"
	  } // end "stack non empty"
	} // end "oriented_set not full"

	for(typename std::set<Facet>::const_iterator fit = 
	      oriented_set.begin();
	    fit != oriented_set.end();
	    ++fit)
	{
	  int indices[3];
	  int index = 0;
	  for (int i=0; i<3; i++)
	    indices[index++] = 
	      V[fit->first->vertex(tr.vertex_triple_index(fit->second, i))];
	  builder.add_facet(indices+0, indices+3);
	  CGAL_assertion_code(++nb_facets);
	}
	CGAL_assertion(nb_facets == number_of_facets);
// 	for( Finite_facets_iterator fit = tr.finite_facets_begin();
// 	     fit != tr.finite_facets_end(); ++fit)
// 	  if ((*fit).first->is_facet_on_surface((*fit).second)==true)
// 	  {
// 	    int indices[3];
// 	    int index = 0;
// 	    for (int i=0; i<3; i++)
// 	      std::cerr << ( indices[index++] = V[(*fit).first->vertex(tr.vertex_triple_index(fit->second, i))] ) << ", ";
// 	    std::cerr << "\n";
// 	    if( builder.test_facet(indices+0, indices+3) )
// 	      builder.add_facet(indices+0, indices+3);
// 	    else
// 	    {
// 	      builder.begin_facet();
// 	      builder.add_vertex_to_facet(indices[2]);
// 	      builder.add_vertex_to_facet(indices[1]);
// 	      builder.add_vertex_to_facet(indices[0]);
// 	      builder.end_facet();
// 	    }
// 	    CGAL_assertion_code(++nb_facets);
// 	  }
      }
      builder.end_surface();
    }
  } off_builder(c2t3);

  p.delegate( off_builder );
  os << p;
#else
  // Header.
  const Tr& tr = c2t3.triangulation();

  os << "OFF \n"
     << tr.number_of_vertices() << " "
     << c2t3.number_of_facets()
     << " " << 0 << "\n";

  CGAL_assertion(c2t3.number_of_facets() == number_of_facets_on_surface(tr));

  os << std::setprecision(20);

  // Finite vertices coordinates.
  std::map<Vertex_handle, int> V;
  int inum = 0;
  for(Finite_vertices_iterator vit = tr.finite_vertices_begin();
      vit != tr.finite_vertices_end();
      ++vit)
  {
    V[vit] = inum++;
    Point p = static_cast<Point>(vit->point());
    os << p.x() << " " << p.y() << " " << p.z() << "\n";
  }

  // Finite facets indices.
  for( Finite_facets_iterator fit = tr.finite_facets_begin();
       fit != tr.finite_facets_end(); ++fit)
    if ((*fit).first->is_facet_on_surface((*fit).second)==true)
    {
      os << "3 ";
      for (int i=0; i<4; i++)
        if (i != (*fit).second)
          os << V[(*fit).first->vertex(i)] << " ";
      
      os << "\n"; // without color.
    }
#endif
}

// only if cells have is_in_domain() method.
template < class Tr>
void
output_oriented_surface_facets_to_off (std::ostream& os, const Tr & T) {
  typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Point Point;

  // Header.


  os << "OFF \n"
     << T.number_of_vertices() << " " <<
    number_of_facets_on_surface (T) <<
    " " << 0 << "\n";

  os << std::setprecision(20);

  // Finite vertices coordinates.
  std::map<Vertex_handle, int> V;
  int inum = 0;
  for( Finite_vertices_iterator
      vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit) {
    V[vit] = inum++;
    Point p = static_cast<Point>(vit->point());
    os << p.x() << " " << p.y() << " " << p.z() << "\n";
  }

  // Finite facets indices.
  for( Finite_facets_iterator fit = T.finite_facets_begin();
       fit != T.finite_facets_end(); ++fit)
    if ((*fit).first->is_facet_on_surface((*fit).second)==true)
      {
	typename Tr::Facet f = *fit;
	typename Tr::Facet opposite = T.mirror_facet(f);
	CGAL_assertion (f.first->is_in_domain() !=
			opposite.first->is_in_domain());
	if(!f.first->is_in_domain())
	  f = T.mirror_facet(f);
	os << "3 "
	   << V[f.first->vertex(T.vertex_triple_index(f.second,0))] << " "
	   << V[f.first->vertex(T.vertex_triple_index(f.second,1))] << " "
	   << V[f.first->vertex(T.vertex_triple_index(f.second,2))] << " "
	   << "\n"; // without color.
      }
}

// only if cells have is_in_domain() method.
template < class Tr>
void
output_surface_facets_to_ghs   (std::ostream& os_points,
				std::ostream& os_faces,
				const Tr & T) {
  typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Point Point;
  typedef typename Tr::Facet Facet;

  // Header.


  os_points << T.number_of_vertices() << "\n";

  os_points << std::setprecision(20);

  // Finite vertices coordinates.
  std::map<Vertex_handle, int> V;
  int inum = 1;
  for( Finite_vertices_iterator
      vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit) {
    V[vit] = inum++;
    Point p = static_cast<Point>(vit->point());
    os_points << p.x() << " " << p.y() << " " << p.z() << " 0\n";
  }

  os_faces << number_of_facets_on_surface(T) << "\n";
  // Finite facets indices.
  for( Finite_facets_iterator fit = T.finite_facets_begin();
       fit != T.finite_facets_end(); ++fit)
    if ((*fit).first->is_facet_on_surface((*fit).second)==true)
      {
	Facet f = *fit;
	if(!f.first->is_in_domain())
	  f = T.mirror_facet(f);
	os_faces
	  << "3 "
	  << V[f.first->vertex(T.vertex_triple_index(f.second,0))] << " "
	  << V[f.first->vertex(T.vertex_triple_index(f.second,1))] << " "
	  << V[f.first->vertex(T.vertex_triple_index(f.second,2))] << " "
	  << "0 0 0 0\n";
      }
}

template < class Tr>
int number_of_facets_in_domain(const Tr& T) {
  int result=0;
  for (typename Tr::Finite_facets_iterator fit = T.finite_facets_begin();
       fit != T.finite_facets_end(); ++fit) {
    typename Tr::Cell_handle neighb = fit->first->neighbor (fit->second);
    if ((fit->first->is_in_domain () || neighb->is_in_domain()) &&
	!fit->first->is_facet_on_surface (fit->second))
      ++result;
  }
  return result;
}

template < class Tr>
void
output_interior_facets_to_off (std::ostream& os, const Tr & T) {
  typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;

  // Header.


  os << "OFF \n"
     << T.number_of_vertices() << " " <<
    number_of_facets_in_domain (T) <<
    " " << 0 << "\n";

  os << std::setprecision(20);

  // Finite vertices coordinates.
  std::map<Vertex_handle, int> V;
  int inum = 0;
  for( Finite_vertices_iterator
      vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit) {
    V[vit] = inum++;
    os << vit->point() << "\n";
  }

  // Finite facets indices.
  for( Finite_facets_iterator fit = T.finite_facets_begin();
       fit != T.finite_facets_end(); ++fit){
    typename Tr::Cell_handle neighb = fit->first->neighbor (fit->second);
    if ((fit->first->is_in_domain () || neighb->is_in_domain()) &&
	!fit->first->is_facet_on_surface (fit->second))
      {
	os << "3 ";
	for (int i=0; i<4; i++)
          if (i != (*fit).second)
	    os << V[(*fit).first->vertex(i)] << " ";

	os << "\n"; // without color.
      }
  }
}

} // end namespace Surface_mesher_io

// backward compatibility: the CGAL namespace was first forgotten
// TODO: fix this
namespace CGAL {
  using namespace Surface_mesher_io;
}
using namespace Surface_mesher_io;

#endif  // CGAL_IO_COMPLEX_2_IN_TRIANGULATION_3_FILE_WRITER_H
