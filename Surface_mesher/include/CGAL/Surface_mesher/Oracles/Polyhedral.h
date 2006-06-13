// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Steve OUDOT, Laurent Rineau

#ifndef _POLYHEDRAL_H
#define _POLYHEDRAL_H


#include <CGAL/Polyhedral_surface_3.h>
#include <CGAL/iterator.h>
#include <CGAL/Surface_mesher/Oracles/Null_oracle_visitor.h>

namespace CGAL {
  namespace Surface_mesher {

template <class Tr,
          class Visitor = Null_oracle_visitor
         >
class Polyhedral
{
public:
  typedef typename Tr::Geom_traits Geom_traits;

  typedef typename Geom_traits::Point_3 Point;
  typedef typename Kernel_traits<Point>::Kernel::Point_3 Kernel_point;
  typedef typename Geom_traits::Segment_3 Segment_3;
  typedef typename Geom_traits::Ray_3 Ray_3;
  typedef typename Geom_traits::Line_3 Line_3;
  typedef typename Geom_traits::Triangle_3 Triangle_3;

  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;

  typedef Polyhedral<Tr, Visitor> Self;

  typedef Data_structure_using_octree_3<Geom_traits> Subfacets_octree;
  typedef Subfacets_octree Surface_3;

  // Private members

private:
  Tr tr;
  Visitor visitor;

private: // private types
  typedef typename Tr::Cell_handle Cell_handle;

public:

  // Public members

  // Default constructor
  Polyhedral()
  {}

  // Surface constructor
  Polyhedral(Tr& tr,
             Visitor visitor_ = Visitor() )
    : tr(tr), visitor(visitor_)
  {
//       is.seekg(0,std::ios::beg);
//       tr.clear();
//       // The data structure for testing intersections is set
//       std::cerr << "Creating data structure for intersections detection... ";
//       data_struct.input(is, CGAL::Insert_iterator<Tr>(tr));
//       std::cerr << "done\n\n";
  }

//   Finite_vertices_iterator finite_vertices_begin()
//   {
//     return tr.finite_vertices_begin();
//   }


//   Finite_vertices_iterator finite_vertices_end()
//   {
//     return tr.finite_vertices_end();
//   }

  class Intersect_3;

  friend class Intersect_3;

  class Intersect_3 {
    Self& self;
  public:
    Intersect_3(Self& self) : self(self)
    {
    }

    Object operator()(Surface_3& subfacets_octree, Segment_3 s) const
    {
      return self.intersect_segment_surface(subfacets_octree, s);
    }
    
    Object operator()(Surface_3& subfacets_octree, const Ray_3& r) const {
      return self.intersect_ray_surface(subfacets_octree, r);
    }
      
    Object operator()(Surface_3& subfacets_octree, const Line_3& l) const {
      return self.intersect_line_surface(subfacets_octree, l);
    }
  };

  Intersect_3 intersect_3_object()
  {
    return Intersect_3(*this);
  }

  class Construct_initial_points
  {
  public:
    template <typename OutputIteratorPoints>
    OutputIteratorPoints operator() (const Surface_3& subfacets_octree, 
                                     OutputIteratorPoints out, 
                                     int n = 20) // WARNING: why 20?    
    {
      for (typename Surface_3::Finite_vertices_iterator vit =
             subfacets_octree.finite_vertices_begin();
           vit != subfacets_octree.finite_vertices_end() && n > 0;
           ++vit, --n)
        *out++= *vit;
      
      return out;
    }
  };

  Construct_initial_points construct_initial_points_object() 
  {
    return Construct_initial_points();
  }

  bool is_in_volume(Surface_3& subfacets_octree, const Point& p)
  {
    return subfacets_octree.number_of_intersection(p) & 1;
  }

//   // Basic intersection function for segments/rays/lines with the polyhedron
//   template <class Elt>
//   CGAL::Object intersect_with_surface (Octree data_struct, Elt e) {
//     typedef CGAL::Data_structure_using_octree_3<Geom_traits> Octree;
//     for ( typename Octree::Constraint_map_iterator cit = data_struct.c_m.begin();
// 	  cit != data_struct.c_m.end(); ++cit ) {
//       if (cit->second->does_intersect (e))
// 	return cit->second->intersection (e);
//     }

//     return CGAL::Object();
//   }

  CGAL::Object intersect_segment_surface(Subfacets_octree& data_struct, Segment_3 s)
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
      Kernel_point kp;
      if( assign(kp, o) )
      {
        Point p = kp;
        visitor.new_point(p);
        return make_object(p);
      }
      else
        return Object();

/*       return  intersect_with_surface (s);  // basic intersection function */
/*       return data_struct.intersect (s.vertex(0), s.vertex(1));  // Marie */
    }

  CGAL::Object intersect_ray_surface(Subfacets_octree& data_struct, Ray_3 &r)
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
      Kernel_point kp;
      if( assign(kp, o) )
      {
        Point p = kp;
        visitor.new_point(p);
        return make_object(p);
      }
      else
        return Object();


//       return intersect_with_surface (r);  // basic intersection function
//       return data_struct.intersect (r);   // Marie's code
    }


  CGAL::Object intersect_line_surface(Subfacets_octree& data_struct, Line_3 &)
    {
      CGAL_assertion(false);
      return CGAL::Object();
    }
private:


}; // end class Polyhedral

  } // end namespace Surface_mesher
} // end namespace CGAL

#endif // _POLYHEDRAL_H
