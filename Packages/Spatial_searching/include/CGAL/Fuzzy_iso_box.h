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
// file          : include/CGAL/Fuzzy_iso_box.h
// package       : ASPAS (3.12)
// maintainer    : Hans Tangelder <hanst@cs.uu.nl>
// revision      : 3.0
// revision_date : 2003/07/10 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================



#ifndef CGAL_FUZZY_ISO_BOX_H
#define CGAL_FUZZY_ISO_BOX_H
#include <CGAL/Kd_tree_rectangle.h>

namespace CGAL {

  template <class SearchTraits>
  class Fuzzy_iso_box{

    public:

    typedef typename SearchTraits::Point_d Point_d;
    typedef typename SearchTraits::Iso_box_d Iso_box_d;
    typedef typename SearchTraits::FT FT;
    
    private:

    Iso_box_d box;
    FT eps;
    unsigned int dim;

    public:

    	// default constructor
    	Fuzzy_iso_box() {}

        // constructor
	Fuzzy_iso_box(const Point_d& p, const Point_d& q, FT epsilon=FT(0)) 
	  : eps(epsilon)
        {
	  typename SearchTraits::Construct_cartesian_const_iterator_d construct_it;
	  typename SearchTraits::Cartesian_const_iterator_d begin = construct_it(p),
	    end = construct_it(p,1);
	  dim = end - begin;
	  box = SearchTraits::Construct_iso_box_d()(p,q);
	}

        	
        bool contains(const Point_d& p) const {	
	  typename SearchTraits::Construct_cartesian_const_iterator_d construct_it;
	  typename SearchTraits::Cartesian_const_iterator_d pit = construct_it(p); 
		for (unsigned int i = 0; i < dim; ++i, ++pit) {
			if ( ((*pit) < box.min()[i]) || ((*pit) >= box.max()[i]) ) return false;
		}
		return true; 
        }

	bool inner_range_intersects(const Kd_tree_rectangle<SearchTraits>& rectangle) const {   
 		for (unsigned int i = 0; i < dim; ++i) {
        		if ( (box.max()[i]-eps < rectangle.min_coord(i)) 
			|| (box.min()[i]+eps >= rectangle.max_coord(i)) ) return false;
    		}
    		return true;                                     
	}


	bool outer_range_is_contained_by(const Kd_tree_rectangle<SearchTraits>& rectangle) const { 
    		for (unsigned int i = 0; i < dim; ++i) {
        		if (  (box.max()[i]+eps < rectangle.max_coord(i) ) 
			|| (box.min()[i]-eps >= rectangle.min_coord(i)) ) return false;
    		}
    		return true;
  	} 

	

	

  }; // class Fuzzy_iso_box

} // namespace CGAL
#endif // FUZZY_ISO_BOX_H
