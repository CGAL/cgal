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
// Author(s)     : Steve OUDOT, Laurent Rineau

#ifndef CGAL_SURFACE_MESHER_POLYHEDRAL_ORACLE_H
#define CGAL_SURFACE_MESHER_POLYHEDRAL_ORACLE_H

#include <utility>

#include <CGAL/config.h> // CGAL_DEPRECATED
#include <CGAL/iterator.h>
#include <CGAL/Surface_mesher/Null_oracle_visitor.h>
#include <CGAL/Surface_mesher/Has_edges.h>

namespace CGAL {
  namespace Surface_mesher {

template <
  class Surface,
  class Point_creator = Creator_uniform_3<typename Surface::Geom_traits::FT,
    typename Surface::Geom_traits::Point_3>,
  class Visitor = Null_oracle_visitor,
  typename Has_edges_tag_ = typename Surface::Has_edges_tag,
  bool mesh_the_whole_bounding_box = false
  >
class Polyhedral_oracle : Has_edges_tag_
{
public:
  typedef typename Surface::Geom_traits GT;
  typedef GT Geom_traits;
  typedef typename GT::FT FT;

  typedef typename Geom_traits::Point_3 Point;
  typedef typename Surface::Subfacets_tree Subfacets_tree;
  typedef typename Surface::Subsegments_tree Subsegments_tree;
  typedef typename Subsegments_tree::Point_with_index Point_with_index;
  typedef typename Geom_traits::Segment_3 Segment_3;
  typedef typename Geom_traits::Ray_3 Ray_3;
  typedef typename Geom_traits::Line_3 Line_3;
  typedef typename Geom_traits::Triangle_3 Triangle_3;

  typedef Polyhedral_oracle<Surface, 
                            Point_creator, 
                            Visitor,
                            Has_edges_tag_,
                            mesh_the_whole_bounding_box> Self;

  typedef Surface Surface_3;

  typedef Point Intersection_point;

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
#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_PINPOLYHEDRON
      double pa[3], pb[3];
      pa[0] = s[0].x();
      pa[1] = s[0].y();
      pa[2] = s[0].z();
      pb[0] = s[1].x();
      pb[1] = s[1].y();
      pb[2] = s[1].z();
      if(surface.pinpolyhedron_ptr->isPinPolyhedron(pa) == 
	 surface.pinpolyhedron_ptr->isPinPolyhedron(pb))
	return Object();
#endif // CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_PINPOLYHEDRON
#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_OCTREE
      return self.intersect_segment_surface(*surface.subfacets_tree_ptr, s);
#endif
#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_INTERSECTION_DATA_STRUCTURE
      return surface.subfacets_tree_ptr->intersection(s);
#endif
    }
    
    Object operator()(const Surface_3& surface, const Ray_3& r) const {
#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_OCTREE
      return self.intersect_ray_surface(*surface.subfacets_tree_ptr, r);
#endif
      return Object();
    }
      
    Object operator()(const Surface_3& surface, const Line_3& l) const {
#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_OCTREE
      return self.intersect_line_surface(*surface.subfacets_tree_ptr, l);
#endif
      return Object();
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
      for (typename Surface_3::Corner_vertices::const_iterator vit =
             surface.corner_vertices_ptr->begin();
           vit != surface.corner_vertices_ptr->end();
           ++vit)
      {
        Point p = (*vit)->point();
        CGAL_assertion((*vit)->tag() >= 0);
        p.set_on_vertex((*vit)->tag());
        self.visitor.new_point(p);
        *out++= p;
      }
      for (typename Surface_3::Edges_vertices::const_iterator vit =
             surface.edges_vertices_ptr->begin();
           vit != surface.edges_vertices_ptr->end() && n > 0;
           ++vit, --n)
      {
        Point p = (*vit)->point();
        CGAL_assertion((*vit)->tag() >= 0);
        p.set_on_curve((*vit)->tag());
        self.visitor.new_point(p);
        *out++= p;
      }
      for (typename Surface_3::Input_vertices::const_iterator vit =
             surface.input_vertices_ptr->begin();
           vit != surface.input_vertices_ptr->end() && n > 0;
           ++vit, --n)
      {
        Point p = (*vit)->point();
        CGAL_assertion((*vit)->tag() >= 0);
        p.set_on_surface((*vit)->tag());
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

  template <class P>
  bool is_in_volume(const Surface_3& surface, const P& p)
  {
    if(mesh_the_whole_bounding_box)
      return CGAL::do_overlap(surface.bbox(),p.bbox());
#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_PINPOLYHEDRON
    double point[3];
    point[0] = p.x();
    point[1] = p.y();
    point[2] = p.z();
    return surface.pinpolyhedron_ptr->isPinPolyhedron(point);
#endif // CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_PINPOLYHEDRON
    typename CGAL::Random_points_on_sphere_3<Point,
      Point_creator> random_point(FT(1));
    typename Geom_traits::Construct_vector_3 vector =
      Geom_traits().construct_vector_3_object();
    typename Geom_traits::Construct_segment_3 segment =
      Geom_traits().construct_segment_3_object();
    typename Geom_traits::Construct_translated_point_3 translate =
      Geom_traits().construct_translated_point_3_object();
    typename Geom_traits::Bounded_side_3 bounded_side =
      Geom_traits().bounded_side_3_object();
    typename Geom_traits::Construct_scaled_vector_3 scale = 
      Geom_traits().construct_scaled_vector_3_object();

    const typename GT::Iso_cuboid_3& cuboid = surface.subfacets_tree_ptr->iso_cuboid();

    if( bounded_side(cuboid,
		     p) == ON_UNBOUNDED_SIDE )
      return false;

#ifdef CGAL_SURFACE_MESHER_DEBUG_INTERSECTION_DATA_STRUCTURE
    std::cerr << "(in volume) ";
#endif
    std::pair<bool, int> result = std::make_pair(false, 0);

    // upper bound of the diameter of the bounding box
    const FT& diameter = 2*FT(surface.subfacets_tree_ptr->max_length());
    while(! result.first)
    {
      result = surface.subfacets_tree_ptr->
        number_of_intersections(segment(p, 
					translate(p, 
						  scale(vector(ORIGIN,
							       *random_point++),
							diameter))));
    }
    return (result.second % 2) == 1;
  }

  Object intersect_curves_with_triangle(const Surface_3& surface,
					const Triangle_3& t) const
  {
    if(! surface.has_edges())
      return Object();

    const Object o = surface.subsegments_tree_ptr->intersection(t);
    if(const Point_with_index* pi = object_cast<Point_with_index>(&o))
    {
      Point p = *pi;
      CGAL_assertion(pi->surface_index()>=0);
      p.set_on_curve(pi->surface_index());
      visitor.new_point(p);
      return make_object(p);
    }
    else
      return Object();
  }
//   // Basic intersection function for segments/rays/lines with the polyhedron
//   template <class Elt>
//   CGAL::Object intersect_with_surface (Octree data_struct, Elt e) {
//     typedef CGAL::Data_structure_using_tree_3<Geom_traits> Octree;
//     for ( typename Octree::Constraint_map_iterator cit = data_struct.c_m.begin();
// 	  cit != data_struct.c_m.end(); ++cit ) {
//       if (cit->second->does_intersect (e))
// 	return cit->second->intersection (e);
//     }

//     return CGAL::Object();
//   }

#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_OCTREE
  CGAL::Object intersect_segment_surface(const Subfacets_tree& data_struct,
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

      const Object o = data_struct.intersection(s.vertex(0), s.vertex(1));
      if(const Point_with_index* pi = object_cast<Point_with_index>(&o))
      {
        Point p = *pi;
        CGAL_assertion(pi->surface_index()>=0);
        p.set_on_surface(pi->surface_index());
        visitor.new_point(p);
        return make_object(p);
      }
      else
        return Object();

/*       return  intersect_with_surface (s);  // basic intersection function */
/*       return data_struct.intersection (s.vertex(0), s.vertex(1));  // Marie */
    }

  CGAL::Object intersect_ray_surface(const Subfacets_tree& data_struct,
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

      const Object o = data_struct.intersection (r);
      if(const Point_with_index* pi = object_cast<Point_with_index>(&o))
      {
        Point p = *pi;
        CGAL_assertion(pi->surface_index()>=0);
        p.set_on_surface(pi->surface_index());
        visitor.new_point(p);
        return make_object(p);
      }
      else
        return Object();


//       return intersect_with_surface (r);  // basic intersection function
//       return data_struct.intersection (r);   // Marie's code
    }


  CGAL::Object intersect_line_surface(const Subfacets_tree&,
                                      const Line_3 &) const
    {
      CGAL_error();
      return CGAL::Object();
    }
#endif // CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_OCTREE
private:


}; // end class Polyhedral_oracle

template <class GT,
          class Visitor = Null_oracle_visitor
         >
class CGAL_DEPRECATED Polyhedral : public Polyhedral_oracle<GT, Visitor>
{
  Polyhedral(Visitor visitor = Visitor())
    : Polyhedral_oracle<GT, Visitor>(visitor)
  {
  }
};

  } // end namespace Surface_mesher
} // end namespace CGAL

#endif // CGAL_SURFACE_MESHER_POLYHEDRAL_ORACLE_H
