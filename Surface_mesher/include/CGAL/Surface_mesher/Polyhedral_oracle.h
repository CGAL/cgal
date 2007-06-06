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

#ifndef CGAL_SURFACE_MESHER_POLYHEDRAL_ORACLE_H
#define CGAL_SURFACE_MESHER_POLYHEDRAL_ORACLE_H

#include <boost/static_warning.hpp>
#include <utility>

#include <CGAL/iterator.h>
#include <CGAL/Surface_mesher/Null_oracle_visitor.h>

namespace CGAL {
  namespace Surface_mesher {

template <
  class Surface,
  class Point_creator = Creator_uniform_3<typename Surface::Geom_traits::FT,
    typename Surface::Geom_traits::Point_3>,
  class Visitor = Null_oracle_visitor
  >
class Polyhedral_oracle
{
public:
  typedef typename Surface::Geom_traits GT;
  typedef GT Geom_traits;
  typedef typename GT::FT FT;

  typedef typename Geom_traits::Point_3 Point;
  typedef typename Kernel_traits<Point>::Kernel::Point_3 Kernel_point;
  typedef typename Geom_traits::Segment_3 Segment_3;
  typedef typename Geom_traits::Ray_3 Ray_3;
  typedef typename Geom_traits::Line_3 Line_3;
  typedef typename Geom_traits::Triangle_3 Triangle_3;

  typedef Polyhedral_oracle<Surface, Point_creator, Visitor> Self;

  typedef Surface Surface_3;

  typedef typename Surface::Subfacets_octree Subfacets_octree;

  // Private members

private:
  Visitor visitor;

public:

  // Public members

  // Surface constructor
  Polyhedral_oracle(Visitor visitor_ = Visitor() )
    : visitor(visitor_)
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
    const Self& self;
  public:
    Intersect_3(const Self& self) : self(self)
    {
    }

    Object operator()(const Surface_3& surface, const Segment_3& s) const
    {
      return self.intersect_segment_surface(surface.subfacets_octree, s);
    }
    
    Object operator()(const Surface_3& surface, const Ray_3& r) const {
      return self.intersect_ray_surface(surface.subfacets_octree, r);
    }
      
    Object operator()(const Surface_3& surface, const Line_3& l) const {
      return self.intersect_line_surface(surface.subfacets_octree, l);
    }
  };

  Intersect_3 intersect_3_object() const
  {
    return Intersect_3(*this);
  }

  class Construct_initial_points;

  friend class Construct_initial_points;

  class Construct_initial_points
  {
    const Self& self;
  public:
    Construct_initial_points(const Self& self) : self(self)
    {
    }

    template <typename OutputIteratorPoints>
    OutputIteratorPoints operator() (const Surface_3& surface, 
                                     OutputIteratorPoints out, 
                                     int n = 20) const // WARNING: why 20?    
    {
      for (typename std::vector<Point>::const_iterator vit =
             surface.corner_points.begin();
           vit != surface.corner_points.end();
           ++vit)
      {
        Point p = *vit;
        self.visitor.new_point(p);
        *out++= p;
      }
      for (typename std::vector<Point>::const_iterator vit =
             surface.edges_points.begin();
           vit != surface.edges_points.end() && n > 0;
           ++vit, --n)
      {
        Point p = *vit;
        self.visitor.new_point(p);
        *out++= p;
      }
      for (typename std::set<Point>::const_iterator vit =
             surface.input_points.begin();
           vit != surface.input_points.end() && n > 0;
           ++vit, --n)
      {
        Point p = *vit;
        self.visitor.new_point(p);
        *out++= p;
      }
      return out;
    }
  };

  Construct_initial_points construct_initial_points_object() const
  {
    return Construct_initial_points(*this);
  }

  bool is_in_volume(const Surface_3& surface, const Point& p)
  {
    typename CGAL::Random_points_on_sphere_3<Point,
      Point_creator> random_point(FT(1));
    typename Geom_traits::Construct_vector_3 vector =
      Geom_traits().construct_vector_3_object();
    typename Geom_traits::Construct_ray_3 ray =
      Geom_traits().construct_ray_3_object();

    std::pair<bool, int> result = std::make_pair(false, 0);
    while(! result.first)
    {
      result = surface.subfacets_octree.
        number_of_intersections(ray(p, vector(CGAL::ORIGIN,
                                              *random_point++)));
    }
    return (result.second % 2) == 1;
  }

  Object intersect_curves_with_triangle(const Surface_3& surface,
					const Triangle_3& t) const
  {
    Object o = surface.subsegments_octree.intersection(t);
    Kernel_point kp;
    if( assign(kp, o) )
    {
      Point p = kp;
      visitor.new_point(p);
      return make_object(p);
    }
    else
      return Object();
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

  CGAL::Object intersect_segment_surface(const Subfacets_octree& data_struct,
                                         const Segment_3& s) const
    {
      typename Geom_traits::Is_degenerate_3  is_degenerate;
      // debug: test if segment is degenerate
      // (can happen, because of rounding in circumcenter computations)
      if (is_degenerate(s)) {
	std::cerr << "Warning: degenerate segment (" << s << ")\n";
	return CGAL::Object();
      }

      // debug: for detecting whether Marie's code works
      // (we compare with our basic intersection function)
//       CGAL::Object oun, odeux;
//       Point p;
//       oun = data_struct.intersection(s.vertex(0), s.vertex(1));
//       odeux = intersect_with_surface (s);
//       odeux = oun;


//       if ((assign(p, oun) && !assign(p,odeux)) ||
// 	  !assign(p, oun) && assign(p,odeux))
// 	std::cout << "s " << s
// 		  << " " << (assign(p, odeux))
// 		  << std::endl;

      Object o = data_struct.intersection(s.vertex(0), s.vertex(1));
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
/*       return data_struct.intersection (s.vertex(0), s.vertex(1));  // Marie */
    }

  CGAL::Object intersect_ray_surface(const Subfacets_octree& data_struct,
                                     const Ray_3 &r) const
    {
      typename Geom_traits::Is_degenerate_3  is_degenerate;
      // debug: test if segment is degenerate
      // (can happen, because of rounding in circumcenter computations)
      if (is_degenerate(r)) {
	std::cerr << "Warning: degenerate ray (" << r << ")\n";
	return CGAL::Object();
      }
      // debug: for detecting whether Marie's code works
      // (we compare with our basic intersection function)
//       CGAL::Object oun, odeux;
//       Point p;
//       oun = data_struct.intersection (r);
//       //      odeux = intersect_with_surface (r);
//       odeux = oun;

//       if ((assign(p, oun) && !assign(p,odeux)) ||
// 	  !assign(p, oun) && assign(p,odeux))
// 	std::cout << "r " << r
// 		  << " " << (assign(p, odeux))
// 		  << std::endl;

//       return odeux;

      Object o = data_struct.intersection (r);
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
//       return data_struct.intersection (r);   // Marie's code
    }


  CGAL::Object intersect_line_surface(const Subfacets_octree&,
                                      const Line_3 &) const
    {
      CGAL_assertion(false);
      return CGAL::Object();
    }
private:


}; // end class Polyhedral_oracle

template <class GT,
          class Visitor = Null_oracle_visitor
         >
class Polyhedral : public Polyhedral_oracle<GT, Visitor>
{
  typedef int Deprecated__class__use__Polyhedral_oracle__instead;

  Polyhedral(Visitor visitor = Visitor())
    : Polyhedral_oracle<GT, Visitor>(visitor)
  {
    BOOST_STATIC_WARNING(Deprecated__class__use__Polyhedral_oracle__instead() == 1);
  }
};

  } // end namespace Surface_mesher
} // end namespace CGAL

#endif // CGAL_SURFACE_MESHER_POLYHEDRAL_ORACLE_H
