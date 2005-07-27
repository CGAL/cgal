// Copyright (c) 2003-2004  INRIA Sophia-Antipolis (France).
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


#ifndef CGAL_SURFACE_MESHER_MULTI_IMPLICIT_ORACLE_H
#define CGAL_SURFACE_MESHER_MULTI_IMPLICIT_ORACLE_H

#include <vector>
#include <iterator>
#include <algorithm>
#include <CGAL/Point_traits.h>


// NB: this oracle requires that the user provide a function that can 
// compute the value of the potential in any point of space
namespace CGAL {

  namespace Surface_mesher {

  template < class GT, class Function > class Multi_implicit_oracle
  {
  public:
    // Public types
    
    typedef typename GT::Point_3 Point;
    typedef Point_traits<Point> Point_traits;
    typedef typename Point_traits::Bare_point Bare_point_3;
    typedef typename GT::Segment_3 Segment;
    typedef typename GT::Ray_3 Ray;
    typedef typename GT::Line_3 Line;
    typedef typename GT::Triangle_3 Triangle;
    typedef typename GT::FT FT;


    typedef typename std::list<Point> Points;
    
    
  private:
    typedef std::vector<Function*> Functions;
    typedef typename Functions::const_iterator Functions_iterator;

    // Private members
    
    int number_of_functions;
    Functions functions;
    Point center;  // center of bounding ball
    FT radius;  // radius of bounding ball
    FT min_squared_length;  // minimal length of a segment for 
    // detecting intersections with surface
    bool parity_oracle;  // flag that tells whether the surface has no boundary
    bool debug;  // flag for debug mode

    FT apply(Functions_iterator fit, const Point& p)
    {
      return (**fit)(p.x(), p.y(), p.z());
    }
    
    FT apply(const int i, const Point& p)
    { // warning: surface indices start from 1, not 0, hence "i-1"!
      return (*(functions[i-1]))(p.x(), p.y(), p.z());
    }
    
  public:
    
    // Constructors
    Multi_implicit_oracle (Function& f, Point emb_center, FT emb_radius,
			   FT precision, 
			   bool parity = false, bool dbg = false) :
      number_of_functions(1),
      center (emb_center),
      radius (emb_radius),
      min_squared_length (precision * precision),
      parity_oracle (parity),
      debug (dbg)
      {
        functions.push_back(&f);
      }

    template <class It>
    Multi_implicit_oracle (It first, It end, Point emb_center, FT emb_radius,
			   FT precision, 
			   bool parity = false, bool dbg = false) :
      center (emb_center),
      radius (emb_radius),
      min_squared_length (precision * precision),
      parity_oracle (parity),
      debug (dbg)
      {
        std::copy(first, end, std::back_inserter(functions));
        number_of_functions = functions.size();
      }

    //Public Queries
    
    Point get_center() {
      return center; 
    }
    
    FT get_radius() {
      return radius;
    }
    
    FT get_precision() {
      return std::sqrt (min_squared_length);
    }
    
    
    
    // Predicates and Constructions
    
    bool is_in_volume(const Point& p)
    {
      for(Functions_iterator fit = functions.begin();
          fit != functions.end();
          ++fit)
        if(apply(fit, p)<0.)
          return true;
      return false;
    }

    Object intersect_segment_surface(Segment s) {      
      GT ker;
      
      // First rescale segment if necessary
      Object obj=rescale_seg_bounding_sphere(s);
      if (!assign(s,obj))
	return Object();
      
      // ATTENTION: optimization for closed surfaces: if both
      // extremities are on the same side of the surface, return no
      // intersection
      if(parity_oracle)
        {
          bool no_intersection = true;
          for(Functions_iterator fit = functions.begin();
              fit != functions.end();
              ++fit)
            if(apply(fit, s.vertex(0)) * apply(fit, s.vertex(1)) <= 0)
              no_intersection = false;
          if(no_intersection)
            return Object();
        }
      
      // Code for surfaces with boundaries
      
      std::list<Segment> f;
      f.push_back(s);
      
      while(!f.empty()) {
	s=f.front();
	f.pop_front();
	Point p1=s.vertex(0);
	Point p2=s.vertex(1);

	// warning: surface indices are from 1 to number_of_functions
	for(int i = 1; i <= number_of_functions; ++i)
          if (apply(i, p1) * apply(i, p2) < 0)
            return intersect_segment_surface_rec(p1,p2,i);

	if (ker.compute_squared_distance_3_object()
	    (p1,p2) >= min_squared_length)
	  {
	    Point mid=ker.construct_midpoint_3_object()(p1,p2);
	    f.push_back(Segment(p1,mid));
	    f.push_back(Segment(mid,p2));
	  }
      }
      
      return Object();
    }




    Object intersect_ray_surface(Ray r) {      
      GT ker;
      Point p1,p2;
      
      p1=r.point(0);
      
      // The second point is calculated with the radius of the bounding ball
      p2=p1+
	ker.construct_scaled_vector_3_object()
	(ker.construct_vector_3_object()
	 (p1,r.point(1)),
	 2*approximate_sqrt(std::max(radius*radius,approximate_sqrt
			 (ker.compute_squared_distance_3_object()
			  (center,p1))) /
		ker.compute_squared_distance_3_object()(p1,r.point(1))));
      
      
      return(intersect_segment_surface(Segment(p1,p2)));
    }
    
    
    Object intersect_line_surface(Line l) {      
      GT ker;

      Point p1=l.point(0);
      Point p2=l.point(1);
      
      
      // The other points are calculated with the radius of the bounding ball
      
      Point p3=p1+
	ker.construct_scaled_vector_3_object()
	(ker.construct_vector_3_object()
	 (p1,p2),
	 2*approximate_sqrt(std::max(radius*radius,approximate_sqrt
			 (ker.compute_squared_distance_3_object()
			  (center,p1))) /
		ker.compute_squared_distance_3_object()(p1,p2)));
      
      Point p4=p1+ // BEURK
	ker.construct_scaled_vector_3_object()
	(ker.construct_vector_3_object()
	 (p2,p1),
	 2*approximate_sqrt(std::max(radius*radius,approximate_sqrt
			 (ker.compute_squared_distance_3_object()
			  (center,p1))) /
		ker.compute_squared_distance_3_object()(p1,p2)));
      
      
      Point result;
      Object result_temp;
      
      result_temp=intersect_segment_surface(Segment(p1,p3));
      if (assign(result,result_temp))
	return result_temp;
      else
	return(intersect_segment_surface(Segment(p1,p4)));
    }
    
    
    
    // Random points
    Points random_points (int n) {
      CGAL_precondition (n > 0);

      // the exhaustive oracle is used
      bool save_parity = parity_oracle;
//       parity_oracle = false;

      Points result;
      while (n>0) {
	Bare_point_3 p1(2*radius*rand()/INT_MAX-radius+center.x(),
                        2*radius*rand()/INT_MAX-radius+center.y(),
                        2*radius*rand()/INT_MAX-radius+center.z());
	Bare_point_3 p2(2*radius*rand()/INT_MAX-radius+center.x(),
                        2*radius*rand()/INT_MAX-radius+center.y(),
                        2*radius*rand()/INT_MAX-radius+center.z());
	
	Object o = intersect_line_surface (Line (p1,p2));
	Point p;
	if (assign(p,o)) {
	  result.push_back(p);
	  --n;
	}
      }
      
      // We restore the user-defined oracle
      parity_oracle = save_parity;

      return result;
    }
    
    
    
  private:
    
    // Private functions
    
    
    Object intersect_segment_surface_rec(Point& p1, Point& p2,
                                         const int surface_index) {
      GT ker;
      Point mid=ker.construct_midpoint_3_object()(p1,p2);
      
      
      // If the two points are close, then we must decide
      if (ker.compute_squared_distance_3_object()
	  (p1,p2) < min_squared_length)
        {
          mid.set_surface_index(surface_index);
          return make_object(mid);
        }

      // Else we must go on
      else {
        if (apply(surface_index, p1) * apply(surface_index, mid) < 0)
          return intersect_segment_surface_rec (p1, mid, surface_index);
        else
          return intersect_segment_surface_rec (mid, p2, surface_index);
      }
    }

  private:
    // Rescale segment according to bounding sphere
    Object rescale_seg_bounding_sphere(Segment s)
      {
	GT ker;
	
	Point p1=s.vertex(0);
	Point p2=s.vertex(1);
	
	// If both points are too far away from each other, then replace them
	if (ker.compute_squared_distance_3_object()(p1,p2) >
	    17*radius*radius)
	  {
	    
	    // cerr << "[" << p1 << "," << p2 << "] -> [";
	    
	    
	    // If p1 is in the sphere, then replace p2
	    if (ker.compute_squared_distance_3_object()(center,p1)
		<= radius*radius)
	      p2=p1 +
		ker.construct_scaled_vector_3_object()
		(ker.construct_vector_3_object()(p1,p2), 
		 2*radius / 
		 approximate_sqrt(ker.compute_squared_distance_3_object()(p1,p2)));
	    
	    // If p2 is in the sphere, then replace p1
	    else if (ker.compute_squared_distance_3_object()(center,p2)
		     <= radius*radius)
	      p1=p2 +
		ker.construct_scaled_vector_3_object()
		(ker.construct_vector_3_object()(p2,p1), 
		 2*radius / 
		 approximate_sqrt(ker.compute_squared_distance_3_object()(p2,p1)));
	    
	    // If both p1 and p2 are outside the sphere, replace both of them
	    else
	      {
		// First find a point on [p1,p2] in the sphere
		Object obj = intersect_segment_sphere(p1, p2);
		
		// If no point is found, then return
		Point mid;
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
	
	return(make_object(Segment(p1,p2)));
      }
    
    Object intersect_segment_sphere
    (Point a, Point b) {
      FT cosine, deltaprime, root1, root2;
      GT ker;
      
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
	 + radius*radius);
      
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

    // Private methods
  private:
    FT approximate_sqrt(const FT x) {
      return FT (CGAL_NTS sqrt(CGAL_NTS to_double(x)));
    }
  };  // end Multi_Implicit_oracle


  }  // namespace Surface_mesher

} // namespace CGAL


#endif  // CGAL_SURFACE_MESHER_MULTI_IMPLICIT_ORACLE_H
