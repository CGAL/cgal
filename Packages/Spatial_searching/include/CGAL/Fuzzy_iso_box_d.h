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
// file          : include/CGAL/Iso_box_d.h
// package       : ASPAS (3.12)
// maintainer    : Hans Tangelder <hanst@cs.uu.nl>
// revision      : 2.4 
// revision_date : 2002/16/08 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================



#ifndef CGAL_FUZZY_ISO_BOX_D_H
#define CGAL_FUZZY_ISO_BOX_D_H
#include <CGAL/Kd_tree_rectangle.h>

namespace CGAL {

  template <class Item, class Iso_box_d>
  class Fuzzy_iso_box_d{

    public:

    typedef typename Item::R::FT NT;
    
    private:

    Iso_box_d box;
    NT eps;
    unsigned int dim;

    public:

    	// default constructor
    	Fuzzy_iso_box_d() {}
		

	// constructor
	Fuzzy_iso_box_d(Iso_box_d b, NT epsilon=NT(0)) : 
	box(b), eps(epsilon), dim(b.min().dimension())
	{
	 std::cout << "dim=" << dim << std::endl;
         for (unsigned int i = 0; i < dim; ++i) {
		std::cout << box.min_coord(i) << " " << box.max_coord(i) << std::endl;
         }
	}
        	
        bool contains(const Item& p) const {	 
		for (unsigned int i = 0; i < dim; ++i) {
			if ( (p[i] < box.min_coord(i)) || (p[i] >= box.max_coord(i)) ) return false;
		}
		return true; 
        }

        
	bool inner_range_intersects(const Kd_tree_rectangle<NT>* rectangle) const {   
 		for (unsigned int i = 0; i < dim; ++i) {
        		if ( (box.max_coord(i)-eps < rectangle->min_coord(i)) 
			|| (box.min_coord(i)+eps >= rectangle->max_coord(i)) ) return false;
    		}
    		return true;                                     
	}


	bool outer_range_is_contained_by(const Kd_tree_rectangle<NT>* rectangle) const { 
    		for (unsigned int i = 0; i < dim; ++i) {
        		if (  (box.max_coord(i)+eps < rectangle->max_coord(i) ) 
			|| (box.min_coord(i)-eps >= rectangle->min_coord(i)) ) return false;
    		}
    		return true;
  	} 

	

	~Fuzzy_iso_box_d() {}

	

  }; // class Fuzzy_iso_box_d

} // namespace CGAL
#endif // FUZZY_ISO_BOX_D_H
