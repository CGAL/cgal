// ======================================================================
//
// Copyright (c) 2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.5-I-99 $
// release_date  : $CGAL_Date: 2003/05/23 $
//
// file          : include/CGAL/Fuzzy_iso_box_d.h
// package       : ASPAS (3.12)
// maintainer    : Hans Tangelder <hanst@cs.uu.nl>
// revision      : 3.0
// revision_date : 2003/07/10 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================



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
			if ( (p[i] < box->min_coord(i)) || (p[i] >= box->max_coord(i)) ) return false;
		}
		return true; 
        }

        
	bool inner_range_intersects(const Kd_tree_rectangle<NT>* rectangle) const {   
 		for (unsigned int i = 0; i < dim; ++i) {
        		if ( (box->max_coord(i)-eps < rectangle->min_coord(i)) 
			|| (box->min_coord(i)+eps >= rectangle->max_coord(i)) ) return false;
    		}
    		return true;                                     
	}


	bool outer_range_is_contained_by(const Kd_tree_rectangle<NT>* rectangle) const { 
    		for (unsigned int i = 0; i < dim; ++i) {
        		if (  (box->max_coord(i)+eps < rectangle->max_coord(i) ) 
			|| (box->min_coord(i)-eps >= rectangle->min_coord(i)) ) return false;
    		}
    		return true;
  	} 

	

	~Fuzzy_iso_rectangle_d() {delete box;}

	

  }; // class Fuzzy_iso_rectangle_d

} // namespace CGAL
#endif // FUZZY_ISO_RECTANGLE_D_H
