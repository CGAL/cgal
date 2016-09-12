// Copyright (c) 2007  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MIN_ANGLE_H
#define CGAL_MIN_ANGLE_H

#include <cmath>

template<typename Triangulation>
class Compute_min_angle
{  
  public:
    
    typedef typename Triangulation::Cell_handle Cell_handle;
    typedef typename Triangulation::Tetrahedron Tetrahedron;  
    typedef typename Triangulation::Point Point;
  
    // constructor
    Compute_min_angle(Triangulation _tr) : tr(_tr) {}  
    
    // computes the minimum angle between all 6 faces of a tetrahedra
    double 
    operator()(const Cell_handle ch) const
    {
      double min_quotient = compute_quotient(ch, 0, 1, 2, 3);
      min_quotient = std::min(min_quotient,
                              compute_quotient(ch, 0, 2, 1, 3));
      min_quotient = std::min(min_quotient,
                              compute_quotient(ch, 0, 3, 1, 2));
      min_quotient = std::min(min_quotient,
                              compute_quotient(ch, 1, 2, 0, 3));  
      min_quotient = std::min(min_quotient,
                              compute_quotient(ch, 1, 3, 0, 2));  
      min_quotient = std::min(min_quotient,
                              compute_quotient(ch, 2, 3, 0, 1));  
      
      const double volume = CGAL::to_double(tr.tetrahedron(ch).volume());
      
      return asin( 1.5 * volume * min_quotient) * 180 / CGAL_PI; 
    }
  
  
  private:  
    
    Triangulation tr;    
    
    double compute_quotient(const Cell_handle ch,
                            const int i,
                            const int j,
                            const int k,
                            const int l) const
    {
      const Point& pi = ch->vertex(i)->point();
      const Point& pj = ch->vertex(j)->point();
            
      const double edge_lenght = 
        CGAL::to_double(CGAL::sqrt(CGAL::squared_distance(pi, pj)));
        
      const double area_k =
        CGAL::to_double(CGAL::sqrt(tr.triangle(ch, k).squared_area()));   
      
      const double area_l =
        CGAL::to_double(CGAL::sqrt(tr.triangle(ch, l).squared_area()));    
      
      return edge_lenght / area_k / area_l;
    }
};

namespace CGAL {

  namespace details {

    template <typename K>
    typename K::FT
    min_dihedral_angle_aux_compute_quotient(const typename K::Point_3& p0,
					    const typename K::Point_3& p1,
					    const typename K::Point_3& p2,
					    const typename K::Point_3& p3,
					    K k = K())
    {
      typename K::Construct_triangle_3 make_triangle = 
	k.construct_triangle_3_object();
      typename K::Compute_area_3 area = 
	k.compute_area_3_object();
      typename K::Compute_squared_distance_3 sq_distance = 
	k.compute_squared_distance_3_object();

      return CGAL::sqrt(sq_distance(p0, p1))
	/ area(make_triangle(p0, p1, p3))
	/ area(make_triangle(p0, p1, p2));
    }

  } // end namespace details;

  template <typename K>
  typename K::FT
  minimum_dihedral_angle(const typename K::Point_3& p0,
			 const typename K::Point_3& p1,
			 const typename K::Point_3& p2,
			 const typename K::Point_3& p3,
			 K k = K())
  {
    typedef typename K::FT FT;
    typename K::Compute_volume_3 volume = 
      k.compute_volume_3_object();

    using details::min_dihedral_angle_aux_compute_quotient;

    FT min_quotient = 
      min_dihedral_angle_aux_compute_quotient(p0, p1, p2, p3, k);

    min_quotient = 
      std::min(min_quotient,
	       min_dihedral_angle_aux_compute_quotient(p0, p2, p1, p3, k));
    min_quotient = 
      std::min(min_quotient,
	       min_dihedral_angle_aux_compute_quotient(p0, p3, p1, p2, k));
    min_quotient = 
      std::min(min_quotient,
	       min_dihedral_angle_aux_compute_quotient(p1, p2, p0, p3, k));
    min_quotient = 
      std::min(min_quotient,
	       min_dihedral_angle_aux_compute_quotient(p1, p3, p0, p2, k));
    min_quotient = 
      std::min(min_quotient,
	       min_dihedral_angle_aux_compute_quotient(p2, p3, p0, p1, k));

//     std::cerr << CGAL::sqrt(min_quotient) << " - "
// 	      << volume(p0, p1, p2, p3) << " - "
// 	      << FT(1.5) * volume(p0, p1, p2, p3) *
//       min_quotient << "\n";
    return std::asin( FT(1.5) * volume(p0, p1, p2, p3) * min_quotient )
      * FT(180) / FT(CGAL_PI);
  };

} // end namespace CGAL

#endif // CGAL_MIN_ANGLE_H
