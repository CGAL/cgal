// Copyright (c) 2002 Utrecht University (The Netherlands).
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
// Author(s)     : Hans Tangelder (<hanst@cs.uu.nl>)


// custom point container
#ifndef CGAL_POINT_CONTAINER_H
#define CGAL_POINT_CONTAINER_H

#include <list>
#include <vector>
#include <functional>
#include <algorithm>
#include <CGAL/Kd_tree_rectangle.h>

namespace CGAL {

template <class SearchTraits> 
class Point_container {

private:
  typedef typename SearchTraits::Point_d Point_d;
  typedef std::vector<Point_d*> Point_vector;
  
public:
  typedef typename SearchTraits::FT FT;
  
  typedef typename Point_vector::iterator iterator;
  
private:
  iterator b, e; // the iterator range of the Point_container 
  
  int built_coord;    // a coordinate for which the pointer list is built
  Kd_tree_rectangle<SearchTraits> bbox;       // bounding box, i.e. rectangle of node
  Kd_tree_rectangle<SearchTraits> tbox;       // tight bounding box, 
                                              // i.e. minimal enclosing bounding
	                	              // box of points
	                	    
public:

  inline const Kd_tree_rectangle<SearchTraits>& 
  bounding_box() const 
  { 
    return bbox; 
  }

    inline const Kd_tree_rectangle<SearchTraits>&
    tight_bounding_box() const 
    { 
      return tbox; 
    }

    inline int 
    dimension() const 
  { 
    return bbox.dimension(); 
  } 

    inline int 
    built_coordinate() const 
  { 
    return built_coord; 
  } 

    // coordinate of the maximal span
    inline int 
    max_span_coord() const 
  {
    return bbox.max_span_coord(); 
  }

    // coordinate of the maximal tight span
    inline int 
    max_tight_span_coord() const 
  { 
    return tbox.max_span_coord(); 
  }

  inline FT  
  max_span_lower() const 
  { 
    return bbox.min_coord(max_span_coord());
  }

  inline FT  
  max_tight_span_lower() const 
  {
    return tbox.min_coord(max_tight_span_coord());
  }
  
  inline FT  
  max_span_upper() const 
  { 
    return bbox.max_coord(max_span_coord());
  }
  
  inline FT  
  max_tight_span_upper() const 
  {
    return tbox.max_coord(max_tight_span_coord());
  }

  inline FT 
  max_spread() const 
  { 
    return  max_span_upper() -  max_span_lower(); 
  }

  inline FT 
  max_tight_spread() const 
  {
    return  max_tight_span_upper() -  max_tight_span_lower(); 
  }


  int 
  max_tight_span_coord_balanced(FT Aspect_ratio) const 
  {
    int cut_dim(-1);
    FT max_spread_points(FT(-1));
    FT max_length=max_spread();  // length of longest side of box
    int dim=dimension();
    for (int d=0; d<dim; d++) {
      FT length=bbox.max_coord(d)-bbox.min_coord(d);
      
      if (FT(2)*max_length/length <= Aspect_ratio) {
	FT spread=tbox.max_coord(d)-tbox.min_coord(d);
	
	if (spread > max_spread_points) {
	  max_spread_points = spread;
	  cut_dim = d;
	}
      }
    }
    // assert(cut_dim >= 0);
    return cut_dim;
  }

  FT 
  max_span_upper_without_dim(int d) const 
  {
    FT max_span(FT(0));
    int dim=dimension();
    for (int i=0; i<dim; i++) {
      FT span = bbox.max_coord(i)-bbox.min_coord(i);
      if (d != i && span > max_span) max_span=span;
    }
    return max_span;
  }

  FT 
  balanced_fair(int d, FT Aspect_ratio) 
  {
    FT small_piece = max_span_upper_without_dim(d) / Aspect_ratio;
    FT low_cut = bbox.min_coord(d) + small_piece; // lowest legal cut;
    FT high_cut = bbox.max_coord(d) - small_piece; //highest legal cut;
    // assert (high_cut >= low_cut);
    FT split_value = median(d);
    if (split_value < low_cut) split_value = low_cut;
    if (split_value > high_cut) split_value = high_cut;
    return split_value;
  }

  FT 
  balanced_sliding_fair(int d, FT Aspect_ratio) 
  {
    FT small_piece = max_span_upper_without_dim(d) / Aspect_ratio;
    FT low_cut = bbox.min_coord(d) + small_piece; // lowest legal cut;
    FT high_cut = bbox.max_coord(d) - small_piece; //highest legal cut;
    // assert (high_cut >= low_cut);
    FT split_value = median(d);
    FT max_span_lower = tbox.min_coord(d);
    FT max_span_upper = tbox.max_coord(d);
    if (split_value < low_cut) split_value= max_span_lower; 
    if (split_value > high_cut) split_value = max_span_upper; 
    return split_value;
  }

  //  points
  inline unsigned int 
  size() const 
  {
    return e - b;
  }
  
  inline iterator 
  begin() const {
    return b;
  }
  
  inline iterator 
  end() const 
  {
    return e;
  }
     
  inline bool 
  empty() const
  {
    return b == e;
  }

  // building the container from a sequence of Point_d*
  Point_container(const int d, iterator begin, iterator end) :
    b(begin), e(end), bbox(d), tbox(d)  
  {
    bbox = Kd_tree_rectangle<SearchTraits>(d, begin, end);
    tbox = bbox;
    
    built_coord = max_span_coord();
  }

  void 
  set_range(iterator begin, iterator end)
  {
    b = begin;
    e = end;
  }


  // building an empty container 
  Point_container(const int d) :
    b(NULL), e(NULL), bbox(d), tbox(d)  
  {}
  
  template <class SearchTraits>   
  struct Cmp {
    typedef typename SearchTraits::FT FT;
    typedef typename SearchTraits::Point_d Point_d;
    typedef std::vector<Point_d*> Point_vector;
    
    int split_coord;
    FT value;
    typename SearchTraits::Construct_cartesian_const_iterator_d construct_it;
    
    Cmp(int s, FT c)
      : split_coord(s), value(c)
    {}
    
    bool 
    operator()(Point_d* pt) const
    {
      typename SearchTraits::Cartesian_const_iterator_d ptit;
      ptit = construct_it(*pt);
      return  *(ptit+split_coord) < value; 
    }
  };
  
  void recompute_tight_bounding_box() 
  {
    tbox.update_from_point_pointers(begin(), end(), begin() == end());
  }
  
  
  // note that splitting is restricted to the built coordinate
  template <class Separator>
  void split(Point_container<SearchTraits>& c, Separator& sep,  
	     bool sliding=false) 
  {
    assert(dimension()==c.dimension());
    
    c.bbox=bbox;
    
    const int split_coord = sep.cutting_dimension();
    const FT cutting_value = sep.cutting_value();

    built_coord=split_coord;
    c.built_coord=split_coord;
		
	
    typename SearchTraits::Construct_cartesian_const_iterator_d construct_it;
    typename SearchTraits::Cartesian_const_iterator_d ptit;

    Cmp<SearchTraits> cmp(split_coord, cutting_value);
    iterator it = std::partition(begin(), end(), cmp);
    // now [begin,it) are lower and [it,end) are upper
    if (sliding) { // avoid empty lists 
	  
      if (it == begin()) {
	FT min_value = bbox.max_coord(built_coord);
	bool found_smaller = false;
	for (iterator pt = begin(); (pt != end()); ++pt) {
	  ptit = construct_it((*(*pt)));	
	  if ( *(ptit+split_coord) < min_value) {
	    min_value= *(ptit+split_coord);
	    it = pt;
	    found_smaller = true;
	  }
	}
	if(found_smaller){
	  Cmp<SearchTraits> cmp2(split_coord, min_value);
	  it = std::partition(begin(), end(), cmp2);
	}
      }
      if (it == end()) {
	FT max_value=bbox.min_coord(built_coord);
	bool found_larger = false;
	for (iterator pt = begin(); (pt != end()); ++pt) {
	  ptit = construct_it((*(*pt)));	
	  if ( *(ptit+split_coord) > max_value) {
	    max_value= *(ptit+split_coord) ;
	    it = pt;
	    found_larger = true;
	  }
	}
	if(found_larger){
	  Cmp<SearchTraits> cmp2(split_coord, max_value);
	  it = std::partition(begin(), end(), cmp2);
	}
	    
      }
    }
	
    c.set_range(begin(), it);
    set_range(it, end());

    // adjusting boxes
    bbox.set_lower_bound(split_coord, cutting_value);
    tbox.update_from_point_pointers(begin(),
				    end(),
				    empty());
	
    c.bbox.set_upper_bound(split_coord, cutting_value);
    c.tbox.update_from_point_pointers(c.begin(),
				      c.end(),
				      c.empty());
        
       
  }



  template <class SearchTraits2, class Value>
  struct comp_coord_val {
    
  private:
    Value coord;   
    
    typedef typename SearchTraits2::Point_d Point_d;
  public:
    comp_coord_val (const Value& coordinate) 
      : coord(coordinate) 
    {}
    
    bool 
    operator()(const Point_d *a, const Point_d *b) const
    {
      typename SearchTraits2::Construct_cartesian_const_iterator_d construct_it;
      typename SearchTraits2::Cartesian_const_iterator_d ait = construct_it(*a),
	bit = construct_it(*b);
      return *(ait+coord) < *(bit+coord);
    }
  };
  

  FT 
  median(const int split_coord) 
  {
    Point_vector::iterator mid = begin() + (end() - begin())/2;
    
    std::nth_element(begin(), mid, end(),comp_coord_val<SearchTraits,int>(split_coord));
    
    typename SearchTraits::Construct_cartesian_const_iterator_d construct_it;
    typename SearchTraits::Cartesian_const_iterator_d mpit = construct_it((*(*mid)));
    FT val1= *(mpit+split_coord);
    mid++;
    mpit = construct_it((*(*mid)));
    FT val2= *(mpit+split_coord);
    return (val1+val2)/FT(2); 
  }



private:
  explicit Point_container() 
  {} // disable default constructor
  
};

  template <class Point>
  std::ostream& 
  operator<< (std::ostream& s, Point_container<Point>& c) 
  {
    s << "Points container of size " << c.size() << "\n cell:";
    s << bbox;
    s << "\n minimal box enclosing points:"; s << tbox; 
    return s;
  }

} // namespace CGAL

#endif // CGAL_POINT_CONTAINER_H


