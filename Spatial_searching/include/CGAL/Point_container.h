// Copyright (c) 2002,2011 Utrecht University (The Netherlands).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
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
#include <CGAL/internal/Get_dimension_tag.h>

#include <boost/optional.hpp>

namespace CGAL {

template <class Traits> 
class Point_container {

private:
  typedef typename Traits::Point_d Point_d;
  typedef std::vector<const Point_d*> Point_vector;
  
public:
  typedef typename Traits::FT FT;
  
  typedef typename Point_vector::iterator iterator;
  typedef typename Point_vector::const_iterator const_iterator;
  typedef typename internal::Get_dimension_tag<Traits>::Dimension D;
private:
  Traits traits;
  // the iterator range of the Point_container 
  boost::optional<iterator> m_b ;
  boost::optional<iterator> m_e ;
  
  int built_coord;    // a coordinate for which the pointer list is built
  Kd_tree_rectangle<FT,D> bbox;       // bounding box, i.e. rectangle of node
  Kd_tree_rectangle<FT,D> tbox;       // tight bounding box, 
  // i.e. minimal enclosing bounding
  // box of points
	                	    
public:

  inline const Kd_tree_rectangle<FT,D>& 
  bounding_box() const 
  { 
    return bbox; 
  }
  
  inline const Kd_tree_rectangle<FT,D>&
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
    FT max_length = max_spread();  // length of longest side of box
    int dim = dimension();
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
    // CGAL_assertion(cut_dim >= 0);
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
    // CGAL_assertion (high_cut >= low_cut);
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
    // CGAL_assertion (high_cut >= low_cut);
    FT split_value = median(d);
    FT max_span_lower = tbox.min_coord(d);
    FT max_span_upper = tbox.max_coord(d);
    if (split_value < low_cut) split_value= max_span_lower; 
    if (split_value > high_cut) split_value = max_span_upper; 
    return split_value;
  }

  //  points
  inline std::size_t
  size() const 
  {
    return *m_e - *m_b;
  }
  
  inline const_iterator 
  begin() const {
    return *m_b;
  }
  
  inline const_iterator 
  end() const 
  {
    return *m_e;
  }

  inline iterator 
  begin()
  {
    return *m_b;
  }
  
  inline iterator 
  end()
  {
    return *m_e;
  }
     
  inline bool 
  empty() const
  {
    return !m_b || !m_e || (*m_b == *m_e ) ;
  }

  // building the container from a sequence of Point_d*
  Point_container(const int d, iterator begin, iterator end,const Traits& traits_) :
    traits(traits_),m_b(begin), m_e(end), bbox(d, begin, end,traits.construct_cartesian_const_iterator_d_object()), tbox(bbox)  
  {
    built_coord = max_span_coord();
  }

  void 
  set_range(iterator begin, iterator end)
  {
    m_b = begin;
    m_e = end;
  }


  // building an empty container 
  Point_container(const int d,const Traits& traits_) :
    traits(traits_),bbox(d), tbox(d)  
  {}
  
  template <class Traits2>   
  struct Cmp {
    typedef typename Traits2::FT FT;
    typedef typename Traits2::Point_d Point_d;
    typedef std::vector<const Point_d*> Point_vector;
    
    int split_coord;
    FT value;
    const typename Traits2::Construct_cartesian_const_iterator_d& construct_it;
    
    Cmp(int s, FT c,const typename Traits2::Construct_cartesian_const_iterator_d& cst_it)
      : split_coord(s), value(c), construct_it(cst_it)
    {}
    
    bool 
    operator()(const Point_d* pt) const
    {
      typename Traits2::Cartesian_const_iterator_d ptit;
      ptit = construct_it(*pt);
      return  *(ptit+split_coord) < value; 
    }
  };


  template <class Traits2>   
  struct Between {
    typedef typename Traits2::FT FT;
    typedef typename Traits2::Point_d Point_d;
    typedef std::vector<const Point_d*> Point_vector;
    
    int split_coord;
    FT low, high;
    const typename Traits2::Construct_cartesian_const_iterator_d& construct_it;
    
    Between(int s, FT l, FT h,const typename Traits2::Construct_cartesian_const_iterator_d& cst_it)
      : split_coord(s), low(l), high(h), construct_it(cst_it)
    {}
    
    bool 
    operator()(const Point_d* pt) const
    {
      typename Traits2::Cartesian_const_iterator_d ptit;
      ptit = construct_it(*pt);
      if(! ( *(ptit+split_coord) <= high ) ){
	//	std::cerr << "Point " << *pt << " exceeds " << high << " in dimension " << split_coord << std::endl;
	return false;
      }
      if(! ( *(ptit+split_coord) >= low ) ){
	//std::cerr << "Point " << *pt << " below " << low << " in dimension " << split_coord << std::endl;
	return false;
      }
      return true;
    }
  };


  void recompute_tight_bounding_box() 
  {
    tbox.template update_from_point_pointers<typename Traits::Construct_cartesian_const_iterator_d>(begin(), end(),traits.construct_cartesian_const_iterator_d_object());
  }
  

  bool
  is_valid() const
  {
    if(empty()) return true;
    bool b = true;
    for (int i = 0; i < dimension(); i++){
      CGAL_assertion( b = (b && (bbox.min_coord(i) <= tbox.min_coord(i))));
      CGAL_assertion( b = (b && (bbox.max_coord(i) >= tbox.max_coord(i))));

      typename Traits::Construct_cartesian_const_iterator_d construct_it=traits.construct_cartesian_const_iterator_d_object();
      Between<Traits> between(i,tbox.min_coord(i), tbox.max_coord(i), construct_it);
      for(const_iterator it = begin(); it != end(); it++){
	b = (b && between(*it));
      }
    }
    return b;
  }


  // note that splitting is restricted to the built coordinate
  template <class Separator>
  void split(Point_container<Traits>& c, Separator& sep,  
	     bool sliding=false) 
  {
    CGAL_assertion(dimension()==c.dimension());
    CGAL_assertion(is_valid());
    c.bbox=bbox;
    
    const int split_coord = sep.cutting_dimension();
    FT cutting_value = sep.cutting_value();

    built_coord=split_coord;
    c.built_coord=split_coord;
		
	
    typename Traits::Construct_cartesian_const_iterator_d construct_it=traits.construct_cartesian_const_iterator_d_object();

    Cmp<Traits> cmp(split_coord, cutting_value,construct_it);
    iterator it = std::partition(begin(), end(), cmp);
    // now [begin,it) are lower and [it,end) are upper
    if (sliding) { // avoid empty lists 

      if (it == begin()) {
	iterator minelt = std::min_element(begin(),end(),comp_coord_val<Traits,int>(split_coord,construct_it));
	if(minelt != it){
	  std::iter_swap(minelt,it);
	}
	cutting_value = *(construct_it(**it)+split_coord);
	sep.set_cutting_value(cutting_value);
	it++;
      }
      if (it == end()) {
	iterator maxelt = std::max_element(begin(),end(),comp_coord_val<Traits,int>(split_coord,construct_it));
	it--;
	if(maxelt != it){
	  std::iter_swap(maxelt,it);
	}
	cutting_value = *(construct_it(**it)+split_coord);
	sep.set_cutting_value(cutting_value);
      }
    }

    c.set_range(begin(), it);
    set_range(it, end());
    // adjusting boxes
    bbox.set_lower_bound(split_coord, cutting_value);
    tbox. template update_from_point_pointers<typename Traits::Construct_cartesian_const_iterator_d>(begin(),end(),construct_it);
    c.bbox.set_upper_bound(split_coord, cutting_value);
    c.tbox. template update_from_point_pointers<typename Traits::Construct_cartesian_const_iterator_d>(c.begin(),c.end(),construct_it);
    CGAL_assertion(is_valid());
    CGAL_assertion(c.is_valid());
  }



  template <class Traits2, class Value>
  struct comp_coord_val {
    
  private:
    Value coord;   
    const typename Traits2::Construct_cartesian_const_iterator_d& construct_it;
    
    typedef typename Traits2::Point_d Point_d;
  public:
    comp_coord_val (const Value& coordinate,const typename Traits2::Construct_cartesian_const_iterator_d& cst_it) 
      : coord(coordinate), construct_it(cst_it)
    {}
    
    bool 
    operator()(const Point_d *a, const Point_d *b) const
    {
      typename Traits2::Cartesian_const_iterator_d ait = construct_it(*a),
	bit = construct_it(*b);
      return *(ait+coord) < *(bit+coord);
    }
  };
  

  FT 
  median(const int split_coord)
  {
    typename Traits::Construct_cartesian_const_iterator_d construct_it=traits.construct_cartesian_const_iterator_d_object();
    iterator mid = begin() + (end() - begin())/2;
    std::nth_element(begin(), mid, end(),comp_coord_val<Traits,int>(split_coord,construct_it));
        
    typename Traits::Cartesian_const_iterator_d mpit = construct_it((*(*mid)));
    FT val1 = *(mpit+split_coord);
    mid++;
    mpit = construct_it((*(*mid)));
    FT val2 = *(mpit+split_coord);
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
    s << c.bounding_box();
    s << "\n minimal box enclosing points:"; s << c.tight_bounding_box(); 
    return s;
  }

} // namespace CGAL

#endif // CGAL_POINT_CONTAINER_H
