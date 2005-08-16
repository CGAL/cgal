// Copyright (c) 2003-2005  INRIA Sophia-Antipolis (France).
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
// $Source: 
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Steve OUDOT

#ifndef _POLYHEDRAL_H
#define _POLYHEDRAL_H


#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Data_structure_using_octree_3.h>
#include <CGAL/iterator.h>
#include <CGAL/Surface_mesher/Oracles/Null_oracle_visitor.h>

namespace CGAL {
  namespace Surface_mesher {

template <class DT_3,
          class Visitor = Null_oracle_visitor>
class Polyhedral
{
public:

  typedef typename DT_3::Geom_traits GT_3;
  typedef GT_3 Kernel_GT_3;

  typedef typename GT_3::Point_3 Point;
  typedef typename GT_3::Segment_3 Segment;
  typedef typename GT_3::Ray_3 Ray;
  typedef typename GT_3::Line_3 Line;
  typedef typename GT_3::Triangle_3 Triangle;

  typedef typename DT_3::Vertex_iterator Vertex_iterator;

  typedef std::list<Point> Points;

  // Private members

private:
  CGAL::Data_structure_using_octree_3<GT_3>  data_struct;
  DT_3 dt;
  Visitor visitor;

private: // private types
  typedef typename DT_3::Cell_handle Cell_handle;
  

public:

  // Public members

  // Default constructor
  Polyhedral()
    {}

  // Surface constructor
  Polyhedral(std::ifstream& is,
             Visitor visitor_ = Visitor() )
    : visitor(visitor_)
    {
      is.seekg(0,std::ios::beg);
      dt.clear();
      // The data structure for testing intersections is set
      std::cerr << "Creating data structure for intersections detection... ";
      data_struct.input(is, CGAL::Insert_iterator<DT_3>(dt));
      std::cerr << "done\n\n";
    }

  Vertex_iterator finite_vertices_begin()
  {
    return dt.finite_vertices_begin();
  }

  
  Vertex_iterator finite_vertices_end()
  {
    return dt.finite_vertices_end();
  }


  // @todo random_points is ad hoc for the implicit oracle only
  Points random_points (int n) {
    Points res;

    for (Vertex_iterator vit = finite_vertices_begin(); 
	 vit != finite_vertices_end() && n > 0;
	 ++vit, --n)
      res.push_back (vit->point());

    return res;
  }

  bool is_in_volume(const Point& p)
  {
    Cell_handle c = dt.locate(p);
    return c->is_in_domain();
  }


  // Basic intersection function for segments/rays/lines with the polyhedron 
  template <class Elt>
  CGAL::Object intersect_with_surface (Elt e) {
    typedef CGAL::Data_structure_using_octree_3<GT_3> Octree;
    for ( typename Octree::Constraint_map_iterator cit = data_struct.c_m.begin();
	  cit != data_struct.c_m.end(); ++cit ) {
      if (cit->second->does_intersect (e))
	return cit->second->intersection (e);
    }
    
    return CGAL::Object();
  }

  CGAL::Object intersect_segment_surface(Segment s)
    {      
      // debug: test if segment is degenerate 
      // (can happen, because of rounding in circumcenter computations)
      if (s.vertex(0)==s.vertex(1)) {
	std::cerr << "Warning: degenerate segment (" << s << ")\n";
	return CGAL::Object();
      }

      // debug: for detecting whether Marie's code works
      // (we compare with our basic intersection function)
//       CGAL::Object oun, odeux;
//       Point p;
//       oun = data_struct.intersect(s.vertex(0), s.vertex(1));
//       odeux = intersect_with_surface (s);
//       odeux = oun;


//       if ((assign(p, oun) && !assign(p,odeux)) ||
// 	  !assign(p, oun) && assign(p,odeux))
// 	std::cout << "s " << s 
// 		  << " " << (assign(p, odeux))
// 		  << std::endl;

      Object o = data_struct.intersect(s.vertex(0), s.vertex(1));
      Point p;
      if( assign(p, o) )
      {
        visitor.new_point(p);
      }
      return o;

/*       return  intersect_with_surface (s);  // basic intersection function */
/*       return data_struct.intersect (s.vertex(0), s.vertex(1));  // Marie */
    }

  CGAL::Object intersect_ray_surface(Ray &r)
    {
      // debug: for detecting whether Marie's code works
      // (we compare with our basic intersection function)
//       CGAL::Object oun, odeux;
//       Point p;
//       oun = data_struct.intersect (r);
//       //      odeux = intersect_with_surface (r);
//       odeux = oun;

//       if ((assign(p, oun) && !assign(p,odeux)) ||
// 	  !assign(p, oun) && assign(p,odeux))
// 	std::cout << "r " << r 
// 		  << " " << (assign(p, odeux)) 
// 		  << std::endl;

//       return odeux;

      Object o = data_struct.intersect (r);
      Point p;
      if( assign(p, o) )
      {
        visitor.new_point(p);
      }
      return o;
     

//       return intersect_with_surface (r);  // basic intersection function
//       return data_struct.intersect (r);   // Marie's code
    }


  CGAL::Object intersect_line_surface(Line &)
    {      
      CGAL_assertion(false);
      return CGAL::Object();
    }
private:


}; // end class Polyhedral

  } // end namespace Surface_mesher
} // end namespace CGAL

#endif // _POLYHEDRAL_H
