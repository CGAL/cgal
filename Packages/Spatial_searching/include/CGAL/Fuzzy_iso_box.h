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
    typedef typename SearchTraits::Construct_vertex_d Construct_vertex_d;
    typedef typename SearchTraits::Cartesian_const_iterator_d Cartesian_const_iterator_d;
    typedef typename SearchTraits::Construct_cartesian_const_iterator_d Construct_cartesian_const_iterator_d;

    private:

    Point_d min, max;
    Cartesian_const_iterator_d min_begin, max_begin;
    FT eps;
    unsigned int dim;

    public:

    	// default constructor
    	Fuzzy_iso_box() {}

        // constructor
	Fuzzy_iso_box(const Point_d& p, const Point_d& q, FT epsilon=FT(0)) 
	  : eps(epsilon)
        {
	  Construct_cartesian_const_iterator_d construct_it;
	  Cartesian_const_iterator_d begin = construct_it(p),
	    end = construct_it(p,1);
	  dim = end - begin;

	  Iso_box_d box = typename SearchTraits::Construct_iso_box_d()(p,q);
	  Construct_vertex_d construct_vertex_d;
	  min = construct_vertex_d(box, 0);
	  max = construct_vertex_d(box, (1<<dim)-1);
	  min_begin = construct_it(min);
	  max_begin = construct_it(max);
	  
	}

        	
        bool contains(const Point_d& p) const {	
	  Construct_cartesian_const_iterator_d construct_it;
	  Cartesian_const_iterator_d pit = construct_it(p);
	  Cartesian_const_iterator_d minit= min_begin, maxit = max_begin; 
		for (unsigned int i = 0; i < dim; ++i, ++pit, ++minit, ++maxit) {
			if ( ((*pit) < (*minit)) || ((*pit) >= (*maxit)) ) return false;
		}
		return true; 
        }

	bool inner_range_intersects(const Kd_tree_rectangle<SearchTraits>& rectangle) const { 
	  Cartesian_const_iterator_d minit= min_begin, maxit = max_begin;   
 		for (unsigned int i = 0; i < dim; ++i, ++minit, ++maxit) {
        		if ( ((*maxit)-eps < rectangle.min_coord(i)) 
			|| ((*minit)+eps >= rectangle.max_coord(i)) ) return false;
    		}
    		return true;                                     
	}


	bool outer_range_is_contained_by(const Kd_tree_rectangle<SearchTraits>& rectangle) const {  
	  Cartesian_const_iterator_d minit= min_begin, maxit = max_begin;   
    		for (unsigned int i = 0; i < dim; ++i, ++minit, ++maxit) {
        		if (  ((*maxit)+eps < rectangle.max_coord(i) ) 
			|| ((*minit)-eps >= rectangle.min_coord(i)) ) return false;
    		}
    		return true;
  	} 

	

	

  }; // class Fuzzy_iso_box

} // namespace CGAL
#endif // FUZZY_ISO_BOX_H
