// Copyright (c) 2003-2007  INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Steve Oudot, Laurent Rineau

#ifndef CGAL_IO_COMPLEX_2_IN_TRIANGULATION_3_FILE_WRITER_H
#define CGAL_IO_COMPLEX_2_IN_TRIANGULATION_3_FILE_WRITER_H

#include <CGAL/license/Surface_mesher.h>


#define CGAL_C2T3_USE_FILE_WRITER_OFF

#include <CGAL/IO/File_writer_OFF.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <iomanip>
#include <stack>
#include <set>

namespace CGAL { namespace Surface_mesher {
  template < class Tr>
  typename Tr::size_type number_of_facets_on_surface(const Tr& T);

  template <class Triangulation>
  class Write_to_OFF_file 
  {
    CGAL::File_writer_OFF off;
    std::ostream& os;

    typedef Triangulation Tr;
  public:
    Write_to_OFF_file(std::ostream& os,  bool verbose)
      : off(verbose), os(os)
    {
    }

    bool write_header(const typename Tr::size_type number_of_vertices,
		      const typename Tr::size_type number_of_facets)
    {
      off.header().set_no_comments(true);
      off.write_header(os,
		       number_of_vertices,
		       0, // fake number of halfedges, not used.
		       number_of_facets);
      return os.good();
    }

    bool write_vertex(const typename Tr::Vertex_handle& v)
    {
      const typename Tr::Point& p = v->point();
      off.write_vertex(p.x(), p.y(), p.z());
      return os.good();
    }

    bool begin_facets()
    {
      off.write_facet_header();
      return os.good();
    }

    bool write_facet(const int index1,
		     const int index2,
		     const int index3)
    {
      off.write_facet_begin(3);
      off.write_facet_vertex_index(index1);
      off.write_facet_vertex_index(index2);
      off.write_facet_vertex_index(index3);
      off.write_facet_end();
      return os.good();
    }

    bool write_footer()
    {
      off.write_footer();
      return os.good();
    }
  }; // end class Write_to_OFF_file

  template <class Triangulation, class HDS>
  class Write_to_HDS 
  {
    CGAL::Polyhedron_incremental_builder_3<HDS> builder;

    typedef Triangulation Tr;
  public:
    Write_to_HDS(HDS& hds)
      : builder(hds)
    {
    }

    bool write_header(const typename Tr::size_type number_of_vertices,
		      const typename Tr::size_type number_of_facets)
    {
      builder.begin_surface(number_of_vertices,
			    number_of_facets);
      return !builder.error();
    }

    bool write_vertex(const typename Tr::Vertex_handle& v)
    {
      const typename Tr::Point& p = v->point();
      builder.add_vertex(p);
      return !builder.error();
    }

    bool begin_facets()
    {
      return !builder.error();
    }

    bool write_facet(const int index1,
		     const int index2,
		     const int index3)
    {
      int indices[3];
      indices[0]=index1;
      indices[1]=index2;
      indices[2]=index3;
      builder.add_facet(indices+0, indices+3);
      return !builder.error();
    }

    bool write_footer()
    {
      builder.end_surface();
      return builder.error();
    }
  }; // end class Write_to_HDS

  enum IO_option { NO_OPTION = 0,
		   IO_ORIENT_SURFACE = 1,
		   IO_VERBOSE = 2 };

  } // end namespace Surface_mesher

template <class C2t3>
bool output_surface_facets_to_off (std::ostream& os,
				   const C2t3& c2t3,
				   int options = 
				   Surface_mesher::IO_ORIENT_SURFACE)
{
  using CGAL::Surface_mesher::number_of_facets_on_surface;

  typedef typename C2t3::Triangulation Tr;
  typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Facet Facet;
  typedef typename Tr::Edge Edge;
  typedef typename Tr::Vertex_handle Vertex_handle;

  // Header.
  const Tr& tr = c2t3.triangulation();

  bool success = true;

  Surface_mesher::Write_to_OFF_file<Tr> 
    off(os, (options & Surface_mesher::IO_VERBOSE) != 0);

  success &= off.write_header(tr.number_of_vertices(),
			      c2t3.number_of_facets());
  
  CGAL_assertion(c2t3.number_of_facets() == number_of_facets_on_surface(tr));
  
  // Finite vertices coordinates.
  std::map<Vertex_handle, int> V;
  int inum = 0;
  for(Finite_vertices_iterator vit = tr.finite_vertices_begin();
      vit != tr.finite_vertices_end();
      ++vit)
  {
    V[vit] = inum++;
    success &= off.write_vertex(vit);
  }

  success &= off.begin_facets();

  if((options & Surface_mesher::IO_ORIENT_SURFACE) == 0) 
  {
    for( Finite_facets_iterator fit = tr.finite_facets_begin();
	 fit != tr.finite_facets_end(); ++fit)
    {
      const typename Tr::Cell_handle cell = fit->first;
      const int& index = fit->second;
      if (cell->is_facet_on_surface(index)==true)
      {
	const int index1 = V[cell->vertex(tr.vertex_triple_index(index, 0))];
	const int index2 = V[cell->vertex(tr.vertex_triple_index(index, 1))];
	const int index3 = V[cell->vertex(tr.vertex_triple_index(index, 2))];
	success &= off.write_facet(index1, index2, index3);
      }
    }
  }
  else // if facets must be oriented
  {

    Finite_facets_iterator fit = tr.finite_facets_begin();
    std::set<Facet> oriented_set;
    std::stack<Facet> stack;

    typename Tr::size_type number_of_facets = c2t3.number_of_facets();

    CGAL_assertion_code(typename Tr::size_type nb_facets = 0; )

      while (oriented_set.size() != number_of_facets) 
      {
	while ( fit->first->is_facet_on_surface(fit->second) == false ||
		oriented_set.find(*fit) != oriented_set.end() ||
		
		oriented_set.find(c2t3.opposite_facet(*fit)) !=
		oriented_set.end() ) 
	{
	  ++fit;
	}
	oriented_set.insert(*fit);
	stack.push(*fit);
	while(! stack.empty() )
	{
	  Facet f = stack.top();
	  stack.pop();
	  for(int ih = 0 ; ih < 3 ; ++ih) {
	    const int i1  = tr.vertex_triple_index(f.second, tr. cw(ih));
	    const int i2  = tr.vertex_triple_index(f.second, tr.ccw(ih));

	    const typename C2t3::Face_status face_status
	      = c2t3.face_status(Edge(f.first, i1, i2));
	    if(face_status == C2t3::REGULAR) {
	      Facet fn = c2t3.neighbor(f, ih);
	      if (oriented_set.find(fn) == oriented_set.end()) {
		if(oriented_set.find(c2t3.opposite_facet(fn)) == oriented_set.end())
		{
		  oriented_set.insert(fn);
		  stack.push(fn);
		}
		else {
		  success = false; // non-orientable
		}
	      }
	    }
	    else if(face_status != C2t3::BOUNDARY) {
	      success = false; // non manifold, thus non-orientable
	    }
	  } // end "for each neighbor of f"
	} // end "stack non empty"
      } // end "oriented_set not full"
    
    for(typename std::set<Facet>::const_iterator fit = 
	  oriented_set.begin();
	fit != oriented_set.end();
	++fit)
    {
      const typename Tr::Cell_handle cell = fit->first;
      const int& index = fit->second;
      const int index1 = V[cell->vertex(tr.vertex_triple_index(index, 0))];
      const int index2 = V[cell->vertex(tr.vertex_triple_index(index, 1))];
      const int index3 = V[cell->vertex(tr.vertex_triple_index(index, 2))];
      success &= off.write_facet(index1, index2, index3);
      CGAL_assertion_code(++nb_facets);
    }

    CGAL_assertion(nb_facets == number_of_facets);
  } // end if(facets must be oriented)

  success &= off.write_footer();
  return success;
}

//   typedef CGAL::Polyhedron_3<Gt> Polyhedron_3;
//   typedef typename Polyhedron_3::HalfedgeDS HalfedgeDS;

//   Polyhedron_3 p;

//   struct Build_polyhedron : public CGAL::Modifier_base<HalfedgeDS>{
//     const C2t3& c2t3;
//     const Tr& tr;
//     Build_polyhedron(const C2t3& c2t3) : c2t3(c2t3), tr(c2t3.triangulation()) {};
//     void operator()( HalfedgeDS& hds) {
//       CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> builder(hds, true);
//       const typename Tr::size_type number_of_facets = c2t3.number_of_facets();
//       builder.begin_surface(tr.number_of_vertices(), 
// 			    number_of_facets);
//       {
// 	// Finite vertices coordinates.
// 	std::map<Vertex_handle, int> V;
// 	int inum = 0;
// 	for(Finite_vertices_iterator vit = tr.finite_vertices_begin();
// 	    vit != tr.finite_vertices_end();
// 	    ++vit)
// 	{
// 	  V[vit] = inum++;
// 	  Point p = static_cast<Point>(vit->point());
// 	  builder.add_vertex(p);
// 	}
// 	Finite_facets_iterator fit = tr.finite_facets_begin();
// 	std::set<Facet> oriented_set;
// 	std::stack<Facet> stack;

// 	CGAL_assertion_code(typename Tr::size_type nb_facets = 0; )

// 	while (oriented_set.size() != number_of_facets) {
// 	  while ( fit->first->is_facet_on_surface(fit->second) == false ||
// 		  oriented_set.find(*fit) != oriented_set.end() ||

// 		  oriented_set.find(c2t3.opposite_facet(*fit)) !=
// 		  oriented_set.end() ) {
// 	    ++fit;
// 	  }
// 	  oriented_set.insert(*fit);
// 	  stack.push(*fit);
// 	  while(! stack.empty() ) {
// 	    Facet f = stack.top();
// 	    stack.pop();
// 	    for(int ih = 0 ; ih < 3 ; ++ih) {
// 	      const int i1  = tr.vertex_triple_index(f.second, tr. cw(ih));
// 	      const int i2  = tr.vertex_triple_index(f.second, tr.ccw(ih));
// 	      if( c2t3.face_status(Edge(f.first, i1, i2)) == C2t3::REGULAR ) {
// 		Facet fn = c2t3.neighbor(f, ih);
// 		if (oriented_set.find(fn) == oriented_set.end() &&
// 		    oriented_set.find(c2t3.opposite_facet(fn)) == oriented_set.end())
// 		{
// 		  oriented_set.insert(fn);
// 		  stack.push(fn);
// 		}
// 	      } // end "if the edge is regular"
// 	    } // end "for each neighbor of f"
// 	  } // end "stack non empty"
// 	} // end "oriented_set not full"

// 	for(typename std::set<Facet>::const_iterator fit = 
// 	      oriented_set.begin();
// 	    fit != oriented_set.end();
// 	    ++fit)
// 	{
// 	  int indices[3];
// 	  int index = 0;
// 	  for (int i=0; i<3; i++)
// 	    indices[index++] = 
// 	      V[fit->first->vertex(tr.vertex_triple_index(fit->second, i))];
// 	  builder.add_facet(indices+0, indices+3);
// 	  CGAL_assertion_code(++nb_facets);
// 	}
// 	CGAL_assertion(nb_facets == number_of_facets);
// // 	for( Finite_facets_iterator fit = tr.finite_facets_begin();
// // 	     fit != tr.finite_facets_end(); ++fit)
// // 	  if ((*fit).first->is_facet_on_surface((*fit).second)==true)
// // 	  {
// // 	    int indices[3];
// // 	    int index = 0;
// // 	    for (int i=0; i<3; i++)
// // 	      std::cerr << ( indices[index++] = V[(*fit).first->vertex(tr.vertex_triple_index(fit->second, i))] ) << ", ";
// // 	    std::cerr << "\n";
// // 	    if( builder.test_facet(indices+0, indices+3) )
// // 	      builder.add_facet(indices+0, indices+3);
// // 	    else
// // 	    {
// // 	      builder.begin_facet();
// // 	      builder.add_vertex_to_facet(indices[2]);
// // 	      builder.add_vertex_to_facet(indices[1]);
// // 	      builder.add_vertex_to_facet(indices[0]);
// // 	      builder.end_facet();
// // 	    }
// // 	    CGAL_assertion_code(++nb_facets);
// // 	  }
//       }
//       builder.end_surface();
//     }
//   } build_polyhedron(c2t3);

//   p.delegate( build_polyhedron );
//   os << p;
// }

// only if cells have is_in_domain() method.
template < class Tr>
void
output_oriented_surface_facets_to_off (std::ostream& os, const Tr & T) {
  using CGAL::Surface_mesher::number_of_facets_on_surface;

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

} // end namespace CGAL

#endif  // CGAL_IO_COMPLEX_2_IN_TRIANGULATION_3_FILE_WRITER_H
