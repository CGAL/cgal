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
#include <CGAL/Points_container.h>
#include <CGAL/Plane_separator.h>

namespace CGAL {

enum Split_rule {MEDIAN_OF_MAX_SPREAD, MEDIAN_OF_BOX,
        MIDPOINT_OF_MAX_SPREAD, MIDPOINT_OF_BOX, FAIR, 
		SLIDING_MIDPOINT, SLIDING_FAIR};

enum Shrink_rule {NONE, SIMPLE, CENTROID};

template <class P>
class Splitter {
public:
  typedef typename Kernel_traits<P>::Kernel K;
  typedef typename K::FT NT;
  // virtual void rule(Points_container<P>& c, Plane_separator<NT>* sep) {}
};

template <class P>
class Median_Of_Max_Spread : public Splitter<P> {
public:
  typedef typename Kernel_traits<P>::Kernel K;
  typedef typename K::FT NT;
  Plane_separator<NT>* rule(Points_container<P>& c) {
        Plane_separator<NT>* sep =
        new Plane_separator<NT>(c.max_tight_span_coord(),0.0);
        sep->set_cutting_val(c.median(sep->cutting_dimension()));
        return sep;
  }
};

template <class P>
class Fair : public Splitter<P> {
public:
  typedef typename Kernel_traits<P>::Kernel K;
  typedef typename K::FT NT;
  Plane_separator<NT>* rule(Points_container<P>& c, NT Aspect_ratio) {
		// find legal cut with max spread
	    Plane_separator<NT>* sep = 
        new Plane_separator<NT>(c.max_tight_span_coord_balanced(Aspect_ratio),0.0);
        sep->set_cutting_val(c.balanced_fair(sep->cutting_dimension(),Aspect_ratio));
        return sep;
  }
};

template <class P>
class Sliding_Fair : public Splitter<P> {
public:
  typedef typename Kernel_traits<P>::Kernel K;
  typedef typename K::FT NT;
  Plane_separator<NT>* rule(Points_container<P>& c, NT Aspect_ratio) {
    // find legal cut with max spread
    Plane_separator<NT>* sep = 
    new Plane_separator<NT>(c.max_tight_span_coord_balanced(Aspect_ratio),0.0);
    sep->set_cutting_val(c.balanced_sliding_fair(sep->cutting_dimension(),Aspect_ratio));
    return sep;
  }
};

template <class P>
class Sliding_MidPoint: public Splitter<P> {
public:
  typedef typename Kernel_traits<P>::Kernel K;
  typedef typename K::FT NT;
  Plane_separator<NT>* rule(Points_container<P>& c)
  {
    Plane_separator<NT>* sep = new Plane_separator<NT>(c.max_span_coord(),
              (c.max_span_upper() + c.max_span_lower())/2.0);
	NT max_span_lower = c.tight_bounding_box().lower(c.max_span_coord());
	NT max_span_upper = c.tight_bounding_box().upper(c.max_span_coord());
	if (max_span_upper <= sep->cutting_value()) {
		sep->set_cutting_val(max_span_upper); 
	}
	if (max_span_lower >= sep->cutting_value()) {
		sep->set_cutting_val(max_span_lower); 
	}
	return sep;
  }
};

template <class P>
class Median_Of_Box : public Splitter<P> {
public:
  typedef typename Kernel_traits<P>::Kernel K;
  typedef typename K::FT NT;
  Plane_separator<NT>* rule(Points_container<P>& c)
  {
    Plane_separator<NT>* sep = new Plane_separator<NT>(c.max_span_coord(),0.0);
    sep->set_cutting_val(c.median(sep->cutting_dimension()));
	return sep;
  }
};

template <class P>
class MidPoint_Of_Max_Spread : public Splitter<P> {
public:
  typedef typename Kernel_traits<P>::Kernel K;
  typedef typename K::FT NT;
  Plane_separator<NT>* rule(Points_container<P>& c)
  {
    Plane_separator<NT>* sep = new Plane_separator<NT>(c.max_tight_span_coord(),
    (c.max_tight_span_upper() + c.max_tight_span_lower())/2.0);
	return sep;
  }
};

template <class P>
class MidPoint_Of_Box: public Splitter<P> {
public:
  typedef typename Kernel_traits<P>::Kernel K;
  typedef typename K::FT NT;
  Plane_separator<NT>* rule(Points_container<P>& c)
  {
    Plane_separator<NT>* sep = new Plane_separator<NT>(c.max_span_coord(),
              (c.max_span_upper() + c.max_span_lower())/2.0);
	return sep;
  }
};

} // namespace CGAL
#endif // CGAL_SPLITTING_RULES_H



