// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
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
// file          : include/CGAL/Extended_internal_node.h
// package       : APSPAS
// revision      : 1.0 
// revision_date : 2001/06/12 
// maintainer : Hans Tangelder (<hanst@cs.uu.nl>)
//
// ======================================================================

#ifndef CGAL_EXTENDED_INTERNAL_NODE_H
#define CGAL_EXTENDED_INTERNAL_NODE_H
#include <CGAL/Internal_node.h>
#include <iomanip>

namespace CGAL {

template <class Traits> // = Kd_tree_traits_2d>
class Extended_internal_node : public Internal_node<Traits> {

public:
  typedef Traits::Item Item;
  typedef typename double NT; // Item::FT NT;
  typedef Traits::Separator Separator;
  typedef Base_node<Traits> Node;
  typedef CGAL::Point_2< CGAL::Cartesian<NT> > Point_2D;
  typedef CGAL::Segment_2< CGAL::Cartesian<NT> > Segment_2D;

private:

  NT low_val;
  NT high_val;

public:

  inline const NT low_value() const { return low_val; }
  inline const NT high_value() const { return high_val; }

  // Extended_internal_node(const  Node* lower, const  Node* upper) :
  //			  lower_ch(lower), upper_ch(upper) {}

  Extended_internal_node(Points_container<Item>& c, Traits& t) {

    Points_container<Item> c_low = Points_container<Item>(c.dimension());

    Box<NT> bbox(c.bounding_box());

    sep = t.split(c, c_low);
	
    int cd  = sep->cutting_dimension();
    low_val = bbox.lower(cd);
    high_val = bbox.upper(cd);

    if (c_low.size() > t.bucket_size())
      lower_ch = new Extended_internal_node<Traits>(c_low,t);
    else
      lower_ch = new Leaf_node<Traits>(c_low);

    // delete *c_low;
    // delete []c_low;

    if (c.size() > t.bucket_size())
      upper_ch = new Extended_internal_node<Traits>(c,t);
    else
      upper_ch = new Leaf_node<Traits>(c);
  }


  ~Extended_internal_node() {
  }
};


} // namespace CGAL
#endif // CGAL_EXTENDED_INTERNAL_NODE_H
