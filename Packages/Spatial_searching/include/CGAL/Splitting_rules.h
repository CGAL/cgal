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
// file          : include/CGAL/Splitting_rules.h
// package       : ASPAS
// revision      : 1.4 
// revision_date : 2002/16/08 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================

// Defines rules used for constructing a split node. That is, it implements, 
// in several ways, the concept
// Boxtree_splitter<NT>.

#ifndef CGAL_SPLITTING_RULES_H
#define CGAL_SPLITTING_RULES_H
#include <CGAL/Point_container.h>
#include <CGAL/Plane_separator.h>
namespace CGAL {

class Split_rules {
public:
enum Split_rule {MEDIAN_OF_MAX_SPREAD, MEDIAN_OF_RECTANGLE,
        MIDPOINT_OF_MAX_SPREAD, MIDPOINT_OF_RECTANGLE, FAIR, 
		SLIDING_MIDPOINT, SLIDING_FAIR};
};

template <class Item>
class Splitter {
public:
// empty base class
};

template <class Item>
class Median_of_max_spread : public Splitter<Item> {
public:
  typedef typename Item::R::FT NT;
  Plane_separator<NT>* rule(Point_container<Item>& c) {
        Plane_separator<NT>* sep =
        new Plane_separator<NT>(c.max_tight_span_coord(),NT(0));
        sep->set_cutting_val(c.median(sep->cutting_dimension()));
        return sep;
  }
};

template <class Item>
class Fair : public Splitter<Item> {
public:
  typedef typename Item::R::FT NT;
  Plane_separator<NT>* rule(Point_container<Item>& c, NT Aspect_ratio) {
		// find legal cut with max spread
	    Plane_separator<NT>* sep = 
        new Plane_separator<NT>(c.max_tight_span_coord_balanced(Aspect_ratio),
				NT(0));
        sep->set_cutting_val(c.balanced_fair(sep->cutting_dimension(),
				Aspect_ratio));
        return sep;
  }
};

template <class Item>
class Sliding_fair : public Splitter<Item> {
public:
  typedef typename Item::R::FT NT;
  Plane_separator<NT>* rule(Point_container<Item>& c, NT Aspect_ratio) {
    // find legal cut with max spread
   
    Plane_separator<NT>* sep = 
    new Plane_separator<NT>(c.max_tight_span_coord_balanced(Aspect_ratio),
			    NT(0));
    
    sep->set_cutting_val(c.balanced_sliding_fair(sep->cutting_dimension(),
			 Aspect_ratio));
    
    return sep;
  }
};

template <class Item>
class Sliding_midpoint: public Splitter<Item> {
public:
  typedef typename Item::R::FT NT;
  Plane_separator<NT>* rule(Point_container<Item>& c)
  {
    Plane_separator<NT>* sep = new Plane_separator<NT>(c.max_span_coord(),
              (c.max_span_upper() + c.max_span_lower())/NT(2));
	NT max_span_lower = c.tight_bounding_box().min_coord(c.max_span_coord());
	NT max_span_upper = c.tight_bounding_box().max_coord(c.max_span_coord());
	if (max_span_upper <= sep->cutting_value()) {
		sep->set_cutting_val(max_span_upper); 
	}
	if (max_span_lower >= sep->cutting_value()) {
		sep->set_cutting_val(max_span_lower); 
	}
	return sep;
  }
};

template <class Item>
class Median_of_rectangle : public Splitter<Item> {
public:
  typedef typename Item::R::FT NT;
  Plane_separator<NT>* rule(Point_container<Item>& c)
  {
    Plane_separator<NT>* sep = 
    new Plane_separator<NT>(c.max_span_coord(),NT(0));
    sep->set_cutting_val(c.median(sep->cutting_dimension()));
	return sep;
  }
};

template <class Item>
class Midpoint_of_max_spread : public Splitter<Item> {
public:
  typedef typename Item::R::FT NT;
  Plane_separator<NT>* rule(Point_container<Item>& c)
  {
    Plane_separator<NT>* sep = 
    new Plane_separator<NT>(c.max_tight_span_coord(),
    (c.max_tight_span_upper() + c.max_tight_span_lower())/NT(2));
	return sep;
  }
};

template <class Item>
class Midpoint_of_rectangle: public Splitter<Item> {
public:
  typedef typename Item::R::FT NT;
  Plane_separator<NT>* rule(Point_container<Item>& c)
  {
    Plane_separator<NT>* sep = new Plane_separator<NT>(c.max_span_coord(),
              (c.max_span_upper() + c.max_span_lower())/NT(2));
	return sep;
  }
};

} // namespace CGAL
#endif // CGAL_SPLITTING_RULES_H



