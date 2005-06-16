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

namespace CGAL {
  namespace Surface_mesher {

template <class DT_3>
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

private: // private types
  typedef typename DT_3::Cell_handle Cell_handle;
  

public:

  // Public members

  // Default constructor
  Polyhedral()
    {}

  // Surface constructor
  Polyhedral(std::ifstream& is)
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
      CGAL::Object oun, odeux;
      Point p;
      oun = data_struct.intersect(s.vertex(0), s.vertex(1));
      odeux = intersect_with_surface (s);

      if ((assign(p, oun) && !assign(p,odeux)) ||
	  !assign(p, oun) && assign(p,odeux))
	std::cout << "s " << s 
		  << " " << (assign(p, odeux))
		  << std::endl;

      return odeux;

/*       return  intersect_with_surface (s); */
/*       return data_struct.intersect (s.vertex(0), s.vertex(1)); */
    }

    


  CGAL::Object intersect_ray_surface(Ray &r)
    {      
      CGAL::Object oun, odeux;
      Point p;
      oun = data_struct.intersect (r);
      odeux = intersect_with_surface (r);

      if ((assign(p, oun) && !assign(p,odeux)) ||
	  !assign(p, oun) && assign(p,odeux))
	std::cout << "r " << r 
		  << " " << (assign(p, odeux)) 
		  << std::endl;

      return odeux;

/*       return intersect_with_surface (r); */
/*       return data_struct.intersect (r); */
    }


  CGAL::Object intersect_line_surface(Line &l)
    {      
      CGAL_assertion(false);
      return CGAL::Object();
    }
private:


}; // end class Polyhedral

  } // end namespace Surface_mesher
} // end namespace CGAL

#endif // _POLYHEDRAL_H
