// Copyright (c) 2002  Utrecht University (The Netherlands).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Authors       : Hans Tangelder (<hanst@cs.uu.nl>)



#ifndef CGAL_FUZZY_ISO_RECTANGLE_D_H
#define CGAL_FUZZY_ISO_RECTANGLE_D_H
#include <CGAL/Kd_tree_rectangle.h>

namespace CGAL {

  template <class Point, class Iso_box_d>
  class Fuzzy_iso_rectangle_d{

    public:

    typedef typename Kernel_traits<Point>::Kernel::FT NT;
    
    private:

    Iso_box_d *box;
    NT eps;
    unsigned int dim;

    public:

    	// default constructor
    	Fuzzy_iso_rectangle_d() {}
		

	

        // constructor
	Fuzzy_iso_rectangle_d(const Point& p, const Point& q, NT epsilon=NT(0)) :
        eps(epsilon), dim(p.dimension())
        {box= new Iso_box_d(p,q);}
        	
        bool contains(const Point& p) const {	 
		for (unsigned int i = 0; i < dim; ++i) {
			if ( (p[i] < box->min_coord(i)) || 
			   (p[i] >= box->max_coord(i)) ) return false;
		}
		return true; 
        }

        
	bool inner_range_intersects(const Kd_tree_rectangle<NT>* rectangle) 
        const {   
 		for (unsigned int i = 0; i < dim; ++i) {
        		if ( (box->max_coord(i)-eps < rectangle->min_coord(i)) 
			|| (box->min_coord(i)+eps >= rectangle->max_coord(i)) ) 
			return false;
    		}
    		return true;                                     
	}


	bool outer_range_is_contained_by(const Kd_tree_rectangle<NT>* rectangle) 
        const { 
    		for (unsigned int i = 0; i < dim; ++i) {
        		if (  (box->max_coord(i)+eps < rectangle->max_coord(i) ) 
			|| (box->min_coord(i)-eps >= rectangle->min_coord(i)) ) 
			return false;
    		}
    		return true;
  	} 

	

	~Fuzzy_iso_rectangle_d() {delete box;}

	

  }; // class Fuzzy_iso_rectangle_d

} // namespace CGAL
#endif // FUZZY_ISO_RECTANGLE_D_H
