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
// release       :
// release_date  :
//
// file          : include/CGAL/Splitters.h
// package       : ASPAS
// revision      : 2.4 
// revision_date : 2002/16/08 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================

// Defines rules used for constructing a split node. That is, it implements, 
// in several ways, the concept
// Boxtree_splitter<NT>.

#ifndef CGAL_SPLITTERS_H
#define CGAL_SPLITTERS_H
#include <CGAL/Point_container.h>
#include <CGAL/Plane_separator.h>
namespace CGAL {

template <class Item>
class Median_of_max_spread {
public:
  typedef typename Item::R::FT NT;
  void split(Plane_separator<NT>& sep, Point_container<Item>& c0,
  			     Point_container<Item>& c1, 
  			     NT Aspect_ratio=NT(3)) {        
        sep=Plane_separator<NT>(c0.max_tight_span_coord(),NT(0));
        sep.set_cutting_val(c0.median(sep.cutting_dimension()));
        c0.split_container(c1,sep,true);
  }
};

template <class Item>
class Fair {
public:
  typedef typename Item::R::FT NT;
  void split(Plane_separator<NT>& sep, Point_container<Item>& c0,
  			     Point_container<Item>& c1, 
  			     NT Aspect_ratio=NT(3)) {
	// find legal cut with max spread
        sep=Plane_separator<NT>(c0.max_tight_span_coord_balanced(Aspect_ratio),
				NT(0));
        sep.set_cutting_val(c0.balanced_fair(sep.cutting_dimension(),
				Aspect_ratio));
        c0.split_container(c1,sep);
  }
};

template <class Item>
class Sliding_fair {
public:
  typedef typename Item::R::FT NT;
  void split(Plane_separator<NT>& sep, Point_container<Item>& c0,
  			     Point_container<Item>& c1, 
  			     NT Aspect_ratio=NT(3))  {
    // find legal cut with max spread
   
    sep=Plane_separator<NT>(c0.max_tight_span_coord_balanced(Aspect_ratio),
			    NT(0));
    
    sep.set_cutting_val(c0.balanced_sliding_fair(sep.cutting_dimension(),
			 Aspect_ratio));
    c0.split_container(c1,sep,true);
  }
};


template <class Item>
class Sliding_midpoint {
public:
  typedef typename Item::R::FT NT;
  void split(Plane_separator<NT>& sep, Point_container<Item>& c0,
  			     Point_container<Item>& c1, 
  			     NT Aspect_ratio=NT(3))
  {
        sep=Plane_separator<NT>(c0.max_span_coord(),
              (c0.max_span_upper() + c0.max_span_lower())/NT(2));
	NT max_span_lower = 
	c0.tight_bounding_box().min_coord(c0.max_span_coord());
	NT max_span_upper = 
	c0.tight_bounding_box().max_coord(c0.max_span_coord());
	if (max_span_upper <= sep.cutting_value()) {
		sep.set_cutting_val(max_span_upper); 
	};
	if (max_span_lower >= sep.cutting_value()) {
		sep.set_cutting_val(max_span_lower); 
	};
	c0.split_container(c1,sep,true);
  }
};

template <class Item>
class Median_of_rectangle {
public:
  typedef typename Item::R::FT NT;
  void split(Plane_separator<NT>& sep, Point_container<Item>& c0,
  			     Point_container<Item>& c1, 
  			     NT Aspect_ratio=NT(3))
  {
    sep=Plane_separator<NT>(c0.max_span_coord(),NT(0));
    sep.set_cutting_val(c0.median(sep.cutting_dimension()));
    c0.split_container(c1,sep,true);
  }
};

template <class Item>
class Midpoint_of_max_spread {
public:
  typedef typename Item::R::FT NT;
  void split(Plane_separator<NT>& sep, Point_container<Item>& c0,
  			     Point_container<Item>& c1, 
  			     NT Aspect_ratio=NT(3))
  {
    sep= Plane_separator<NT>(c0.max_tight_span_coord(),
    (c0.max_tight_span_upper() + c0.max_tight_span_lower())/NT(2));
    c0.split_container(c1,sep);
  }
};

template <class Item>
class Midpoint_of_rectangle {
public:
  typedef typename Item::R::FT NT;
  void split(Plane_separator<NT>& sep, Point_container<Item>& c0,
  			     Point_container<Item>& c1, 
  			     NT Aspect_ratio=NT(3))
  {
    sep = Plane_separator<NT>(c0.max_span_coord(),
              (c0.max_span_upper() + c0.max_span_lower())/NT(2));
    c0.split_container(c1,sep);          
  }
 
};

} // namespace CGAL
#endif // CGAL_SPLITTERS



