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
// Author(s)     : Steve OUDOT, Laurent RINEAU


#ifndef CGAL_SURFACE_MESHER_IMPLICIT_ORACLE_H
#define CGAL_SURFACE_MESHER_IMPLICIT_ORACLE_H

#include <CGAL/Surface_mesher/Oracles/Null_oracle_visitor.h>
#include <CGAL/point_generators_3.h>

// NB: this oracle requires that the user provide a function that can
// compute the value of the potential in any point of space
namespace CGAL {

  namespace Surface_mesher {

  template <
    class GT,
    class Surface,
    class Point_creator = Creator_uniform_3<typename GT::RT,
                                            typename GT::Point_3>,
    class Visitor = Null_oracle_visitor
    >
  class Implicit_surface_oracle
  {
    // private types
    typedef Implicit_surface_oracle<GT, Surface, Point_creator, Visitor> Self;
    
    typedef typename GT::Point_3 Point;
    typedef typename Kernel_traits<Point>::Kernel::Point_3 Kernel_point;

    typedef typename GT::FT FT;
    typedef typename Surface::Sphere_3 Sphere_3;

  public:

    // Public types
    typedef GT Geom_traits;
    typedef typename GT::Point_3 Point_3;
    typedef typename GT::Segment_3 Segment_3;
    typedef typename GT::Ray_3 Ray_3;
    typedef typename GT::Line_3 Line_3;

    typedef Surface Surface_3;

  private:
    // Private members
    bool parity_oracle;  // flag that tells whether the surface has no boundary
    Visitor visitor; // a visitor that can modify a point, before returning it.

  public:

    // Constructors
    Implicit_surface_oracle (bool parity = false,
			     Visitor visitor_ = Visitor() ) :
      parity_oracle (parity),
      visitor(visitor_)
      {}

    // Predicates and Constructions

    bool is_in_volume(const Point& p)
    {
      return surf_equation(p)<0.;
    }

    class Intersect_3 
    {
      Self& oracle;
    public:
      Intersect_3(Self& oracle) : oracle(oracle)
      {
      }

      Object operator()(const Surface_3& surface, Segment_3 s) const
      // s is passed by value, because it is used in a CGAL::assign below.
      {
        GT ker;

        // First rescale segment if necessary
        Object obj=oracle.rescale_seg_bounding_sphere(surface, s);
        if (!assign(s,obj))
          return Object();

        // ATTENTION: optimization for closed surfaces: if both
        // extremities are on the same side of the surface, return no
        // intersection
        if(oracle.parity_oracle &&
           (surf_equation(surface, s.source()) * 
            surf_equation(surface, s.target())>0))
          return Object();


        // Code for surfaces with boundaries

        std::list<Segment_3> f;
        f.push_back(s);


        while(!f.empty()) {
          const Segment_3 sf = f.front();
          f.pop_front();
          const Point& p1 = sf.source();
          const Point& p2 = sf.target();

          if (surf_equation(surface, p1) * surf_equation(surface, p2) < 0) {
            return oracle.intersect_segment_surface_rec(surface, p1, p2);
          }

          if (ker.compute_squared_distance_3_object()
              (p1,p2) >= surface.squared_error_bound())
	  {
	    Point mid=ker.construct_midpoint_3_object()(p1,p2);
	    f.push_back(Segment_3(p1,mid));
	    f.push_back(Segment_3(mid,p2));
	  }
        }
        return Object();
      } // end operator()(Surface_3, Segment_3)

      Object operator()(const Surface_3& surface, const Ray_3& r) const {
        GT ker;
        Point p1,p2;

        const Sphere_3& sphere = surface.bounding_sphere();
        const Point center = ker.construct_center_3_object()(sphere);
        const FT squared_radius = 
          ker.compute_squared_radius_3_object()(sphere);

        p1=r.point(0);

        // The second point is calculated with the radius of the bounding ball
        p2=p1+
          ker.construct_scaled_vector_3_object()
          (ker.construct_vector_3_object()
           (p1,r.point(1)),
           2*approximate_sqrt(std::max(squared_radius,approximate_sqrt
                                       (ker.compute_squared_distance_3_object()
                                        (center,p1))) /
                              ker.compute_squared_distance_3_object()(p1,r.point(1))));


        return( operator()(surface, Segment_3(p1,p2)) );
      } // end operator()(Surface_3, Ray_3)

      Object operator()(const Surface_3& surface, const Line_3& l) const {
        GT ker;

        const Sphere_3& sphere = surface.bounding_sphere();
        const Point center = ker.construct_center_3_object()(sphere);
        const FT squared_radius = 
          ker.compute_squared_radius_3_object()(sphere);

         Point p1=l.point(0);
        Point p2=l.point(1);

        // The other points are calculated with the radius of the bounding ball

        Point p3=p1+
          ker.construct_scaled_vector_3_object()
          (ker.construct_vector_3_object()
           (p1,p2),
           2*approximate_sqrt(std::max(squared_radius,approximate_sqrt
                                       (ker.compute_squared_distance_3_object()
                                        (center,p1))) /
                              ker.compute_squared_distance_3_object()(p1,p2)));

        Point p4=p1+ // BEURK
          ker.construct_scaled_vector_3_object()
          (ker.construct_vector_3_object()
           (p2,p1),
           2*approximate_sqrt(std::max(squared_radius,approximate_sqrt
                                       (ker.compute_squared_distance_3_object()
                                        (center,p1))) /
                              ker.compute_squared_distance_3_object()(p1,p2)));

        Object result_temp = operator()(surface, Segment_3(p1,p3));

        Point result;
        if (assign(result,result_temp))
          return result_temp;
        else
          return( operator()(surface, Segment_3(p1,p4)) );
      } // end operator()(Surface_3, Line_3)

    }; // end nested class Intersect_3

    class Construct_initial_points
    {
      Self& oracle;
    public:
      Construct_initial_points(Self& oracle) : oracle(oracle)
      {
      }
      
      // Random points
      template <typename OutputIteratorPoints>
      OutputIteratorPoints operator() (const Surface_3& surface, 
                                       OutputIteratorPoints out, 
                                       int n = 20) // WARNING: why 20?
      {
        GT geom_traits;

        const Sphere_3& sphere = surface.bounding_sphere();
        const Point center = 
          geom_traits.construct_center_3_object()(sphere);
        const FT squared_radius = 
          geom_traits.compute_squared_radius_3_object()(sphere);
        const FT radius = CGAL::sqrt(squared_radius);

        typename CGAL::Random_points_in_sphere_3<Point,
          Point_creator> random_point_in_sphere(CGAL::to_double(radius));
        typename GT::Construct_line_3 line_3 = 
          geom_traits.construct_line_3_object();
        typename GT::Construct_vector_3 vector_3 =
          geom_traits.construct_vector_3_object();
        typename GT::Construct_translated_point_3 translate =
          geom_traits.construct_translated_point_3_object();

        // the exhaustive oracle is used
        bool save_parity = oracle.parity_oracle;
        //      parity_oracle = false; // !! WHY??

        while (n>0) {
          Point p1 = translate(*random_point_in_sphere++,
                               vector_3(CGAL::ORIGIN, center));
          Point p2 = translate(*random_point_in_sphere++,
                               vector_3(CGAL::ORIGIN, center));

          Object o = oracle.intersect_3_object()(surface, line_3(p1,p2));
          Point p;
          if (assign(p,o)) {
            *out++= p;
            --n;
          }
        }

        // We restore the user-defined oracle
        oracle.parity_oracle = save_parity;

        return out;
      }
    }; // end nested class Construct_initial_points

    Construct_initial_points construct_initial_points_object()
    {
      return Construct_initial_points(*this);
    }

    Intersect_3 intersect_3_object()
    {
      return Intersect_3(*this);
    }

  private:

    // Private functions


    Object intersect_segment_surface_rec(const Surface_3& surface,
                                         Point p1,
                                         Point p2) {
      GT ker;

      while(true){
	Point mid = ker.construct_midpoint_3_object()(p1,p2);
	// If the two points are close, then we must decide
	if (ker.compute_squared_distance_3_object()
	    (p1,p2) < surface.squared_error_bound())
	  {
	    visitor.new_point(mid);
	    return make_object(mid);
	  }

	// Else we must go on
	if (surf_equation(surface, p1) * surf_equation(surface, mid) < 0)
	  p2 = mid;
	else
	  p1 = mid;
      }
    }

  private:
    static FT surf_equation (Surface_3 surface, Point p) {
      return surface(p.x(), p.y(), p.z());
    } // @WARNING: we use x(), y() and z()

  private:
    // Rescale segment according to bounding sphere
    Object rescale_seg_bounding_sphere(const Surface_3& surface,
                                       const Segment_3& s)
      {
	GT ker;

        const Sphere_3& sphere = surface.bounding_sphere();
        const Point center = ker.construct_center_3_object()(sphere);
        const FT squared_radius = 
          ker.compute_squared_radius_3_object()(sphere);
        const FT radius = CGAL::sqrt(squared_radius);

 	Point p1=s.source();
	Point p2=s.target();

	// If both points are too far away from each other, then replace them
	if (ker.compute_squared_distance_3_object()(p1,p2) >
	    17*squared_radius)
	  {

	    // cerr << "[" << p1 << "," << p2 << "] -> [";


	    // If p1 is in the sphere, then replace p2
	    if (ker.compute_squared_distance_3_object()(center,p1)
		<= squared_radius)
	      p2=p1 +
		ker.construct_scaled_vector_3_object()
		(ker.construct_vector_3_object()(p1,p2),
		 2*radius /
		 approximate_sqrt(ker.compute_squared_distance_3_object()(p1,p2)));

	    // If p2 is in the sphere, then replace p1
	    else if (ker.compute_squared_distance_3_object()(center,p2)
		     <= squared_radius)
	      p1=p2 +
		ker.construct_scaled_vector_3_object()
		(ker.construct_vector_3_object()(p2,p1),
		 2*radius /
		 approximate_sqrt(ker.compute_squared_distance_3_object()(p2,p1)));

	    // If both p1 and p2 are outside the sphere, replace both of them
	    else
	      {
		// First find a point on [p1,p2] in the sphere
		Object obj = intersect_segment_sphere(surface, p1, p2);

		// If no point is found, then return
		Kernel_point mid;
		if (!assign(mid,obj))
		  return Object();

		// Then reset p1 and p2
		else
		  {
		    p1=mid+
		      ker.construct_scaled_vector_3_object()
		      (ker.construct_vector_3_object()(mid,p1),
		       2*radius /
		       approximate_sqrt(ker.compute_squared_distance_3_object()(mid,p1)));

		    p2=mid+
		      ker.construct_scaled_vector_3_object()
		      (ker.construct_vector_3_object()(mid,p2),
		       2*radius /
		       approximate_sqrt(ker.compute_squared_distance_3_object()(mid,p2)));
		  }
	      }

	    // cerr << p1 << "," << p2 << "]" << endl;

	  }

	return(make_object(Segment_3(p1,p2)));
      }

    Object intersect_segment_sphere(const Surface_3& surface,
                                    const Point& a,
                                    const Point& b) {
      FT cosine, deltaprime, root1, root2;
      GT ker;

      const Sphere_3& sphere = surface.bounding_sphere();
        const Point center = ker.construct_center_3_object()(sphere);
      const FT squared_radius = 
        ker.compute_squared_radius_3_object()(sphere);

      // Compute the squared cosine of angle ([a,center),[a,b))
      cosine=(center-a)*(b-a);
      cosine=cosine*cosine /
	(ker.compute_squared_distance_3_object()(a,center) *
	 ker.compute_squared_distance_3_object()(a,b));

      // Compute the reduced discriminant
      deltaprime =
	ker.compute_squared_distance_3_object()(a,b) *
	((cosine-1)*
	 ker.compute_squared_distance_3_object()(a,center)
	 + squared_radius);

      //      cerr << endl << "   -> discriminant=" << deltaprime << endl;

      // If the discriminant is negative, then no intersection
      if (deltaprime < 0)
	return Object();

      // Else compute the two roots
      deltaprime=approximate_sqrt(deltaprime);
      root1 = ((center-a)*(b-a)) + deltaprime /
	ker.compute_squared_distance_3_object()(a,b);
      root2 = ((center-a)*(b-a) - deltaprime) /
	ker.compute_squared_distance_3_object()(a,b);

      //      cerr << "   -> roots: " << root1 << " and " << root2 << endl;

      // There is an intersection iff a root is between 0 and 1
      if (root1>=0 && root1<=1)
	return make_object(a + ker.construct_scaled_vector_3_object()
				 (ker.construct_vector_3_object()(a,b),
				  root1));
      else if (root2>=0 && root2<=1)
	return make_object(a + ker.construct_scaled_vector_3_object()
				 (ker.construct_vector_3_object()(a,b),
				  root2));
      else
	return Object();
    }
  };  // end Implicit_surface_oracle

template <typename FT>
FT approximate_sqrt(const FT x) {
  return FT (CGAL_NTS sqrt(CGAL_NTS to_double(x)));
}

  }  // namespace Surface_mesher

} // namespace CGAL


#endif  // CGAL_SURFACE_MESHER_IMPLICIT_ORACLE_H
