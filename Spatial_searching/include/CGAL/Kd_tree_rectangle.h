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

#ifndef CGAL_KD_TREE_RECTANGLE_H
#define CGAL_KD_TREE_RECTANGLE_H

#include <functional>
#include <algorithm>
#include <new>
#include <CGAL/assertions.h>
#include <CGAL/array.h>
#include <CGAL/Dimension.h>

namespace CGAL {

  template <class Construct_cartesian_const_iterator_d, class P, class T>
  struct set_bounds_from_pointer : public std::unary_function<P, void> {
    int dim;
    T *lower;
    T *upper;
    Construct_cartesian_const_iterator_d construct_it;
    
    set_bounds_from_pointer(int d, T *l, T *u,Construct_cartesian_const_iterator_d construct_it_) 
      : dim(d), lower(l), upper(u), construct_it(construct_it_)
    {}

    void 
    operator()(P p) 
    {
      T h;
      typename Construct_cartesian_const_iterator_d::result_type pit = construct_it(*p);
      for (int i = 0; i < dim; ++i, ++pit) {
	h=(*pit);
	if (h < lower[i]) lower[i] = h;
	if (h > upper[i]) upper[i] = h;
      }
    }
  };


  template <class FT_, typename D = Dynamic_dimension_tag> 
  class Kd_tree_rectangle {
  public:
    typedef FT_ FT;
    typedef FT T;
    
  private:
    
    //int dim;
    CGAL::cpp11::array<T,D::value> lower_;
    CGAL::cpp11::array<T,D::value> upper_;
    int max_span_coord_;
    
  public:
    
    inline void 
    set_upper_bound(int i, const FT& x) 
    {
      CGAL_assertion(i >= 0 && i < D::value);
      CGAL_assertion(x >= lower_[i]);
      upper_[i] = x;
      set_max_span();
    }

    inline void 
    set_lower_bound(int i, const FT& x) 
    {
      CGAL_assertion(i >= 0 && i < D::value);
      CGAL_assertion(x <= upper_[i]);
      lower_[i] = x;
      set_max_span();
    }

    inline void 
    set_max_span() 
    {
      FT span = upper_[0]-lower_[0];
      max_span_coord_ = 0;
      for (int i = 1; i < D::value; ++i) {
	FT tmp = upper_[i] - lower_[i];
	if (span < tmp) {
	  span = tmp;
	  max_span_coord_ = i;
	}
      }
    }
    
    Kd_tree_rectangle(int)
      : max_span_coord_(0)
    {
      lower_.fill(FT(0));
      upper_.fill(FT(0));
    }

    Kd_tree_rectangle() 
    {}
    
    
    explicit 
    Kd_tree_rectangle(const Kd_tree_rectangle<FT,D>& r) 
    : max_span_coord_(r.max_span_coord_) 
    {
      lower_ = r.lower_;
      upper_ = r.upper_;
    }

    template <class Construct_cartesian_const_iterator_d,class PointPointerIter>
    void update_from_point_pointers(PointPointerIter begin, 
                                    PointPointerIter end,
                                    const Construct_cartesian_const_iterator_d& construct_it
    ) 
    {
      if (begin ==end)
        return;
      // initialize with values of first point
      typename Construct_cartesian_const_iterator_d::result_type bit = construct_it(**begin);
      
      for (int i=0; i < D::value; ++i, ++bit) {
        lower_[i]= *bit; upper_[i]=lower_[i];
      }
      begin++;
      typedef typename std::iterator_traits<PointPointerIter>::value_type P;
      std::for_each(begin, end,set_bounds_from_pointer<Construct_cartesian_const_iterator_d,P,T>(D::value, &(lower_[0]), &(upper_[0]), construct_it));
      set_max_span();
    }
    
    template <class Construct_cartesian_const_iterator_d,class PointPointerIter> // was PointIter
    Kd_tree_rectangle(int,  PointPointerIter begin,  PointPointerIter end,const Construct_cartesian_const_iterator_d& construct_it)
    {
      update_from_point_pointers<Construct_cartesian_const_iterator_d>(begin,end,construct_it);
    }

    inline int 
    max_span_coord() const 
    { 
      return max_span_coord_; 
    }

    inline FT 
    max_span() const 
    {
      return upper_[max_span_coord_] - lower_[max_span_coord_];
    }

    inline FT 
    min_coord(int i) const 
    {
      CGAL_assertion(lower_.size() != 0);
      return lower_[i];
    }
    
    inline FT 
    max_coord(int i) const 
    {
      return upper_[i];
    }
    
    std::ostream& 
    print(std::ostream& s) const 
    {
      s << "Rectangle dimension = " << D::value;
      s << "\n lower: ";
      for (int i=0; i < D::value; ++i)
	s << lower_[i] << " ";
      // std::copy(lower_, lower_ + D,
      // 	      std::ostream_iterator<FT>(s," "));
      s << "\n upper: ";
      for (int j=0; j < D::value; ++j)
	s << upper_[j] << " ";
      // std::copy(upper_, upper_ + D,
      //	      std::ostream_iterator<FT>(s," "));
      s << "\n maximum span " << max_span() <<
	" at coordinate " << max_span_coord() << std::endl;
      return s;
    }
    
    // Splits rectangle by modifying itself to lower half 
    // and returns upper half
    //    Kd_tree_rectangle* 
    void
    split(Kd_tree_rectangle& r, int d, FT value)
    {
      CGAL_assertion(d >= 0 && d < D::value);
      CGAL_assertion(lower_[d] <= value && value <= upper_[d]);
      
      //Kd_tree_rectangle* r = new Kd_tree_rectangle(*this);
      upper_[d]=value;
      r.lower_[d]=value;
      //return r;
    }
    

    int 
    dimension() const 
    {
      return D::value;
    }

    const T* lower() const {return lower_;}
    const T* upper() const {return upper_;}    
 
    Kd_tree_rectangle<FT,D>& 
    operator=(const Kd_tree_rectangle<FT,D>& r) 
    {
      CGAL_assertion(dimension() == r.dimension());
      if (this != &r) {
        lower_ = r.lower_;
        upper_ = r.upper_;
	set_max_span();
      }
      return *this;
    }

  

  }; // of class Kd_tree_rectangle


  // Partial specialization for dynamic dimension, which means dimension at runtime

  template <class FT_> 
  class Kd_tree_rectangle<FT_, Dynamic_dimension_tag> {

  public:
    typedef FT_ FT;
    typedef FT T;
    
  private:
    
    int dim;
    T* lower_;
    T* upper_;
    int max_span_coord_;
    
  public:
    
    inline void 
    set_upper_bound(int i, const FT& x) 
    {
      CGAL_assertion(i >= 0 && i < dim);
      CGAL_assertion(x >= lower_[i]);
      upper_[i] = x;
      set_max_span();
    }

    inline void 
    set_lower_bound(int i, const FT& x) 
    {
      CGAL_assertion(i >= 0 && i < dim);
      CGAL_assertion(x <= upper_[i]);
      lower_[i] = x;
      set_max_span();
    }

    inline void 
    set_max_span() 
    {
      FT span = upper_[0]-lower_[0];
      max_span_coord_ = 0;
      for (int i = 1; i < dim; ++i) {
	FT tmp = upper_[i] - lower_[i];
	if (span < tmp) {
	  span = tmp;
	  max_span_coord_ = i;
	}
      }
    }
    
    Kd_tree_rectangle(int d) 
      : dim(d), lower_(new FT[d]), upper_(new FT[d]), max_span_coord_(0)
    {
      std::fill(lower_, lower_ + dim, FT(0));
      std::fill(upper_, upper_ + dim, FT(0));
    }

    Kd_tree_rectangle() 
      : dim(0), lower_(0), upper_(0) 
    {
}
    
    
    explicit 
    Kd_tree_rectangle(const Kd_tree_rectangle<FT,Dynamic_dimension_tag>& r) 
      : dim(r.dim), lower_(new FT[dim]), upper_(new FT[dim]), 
	max_span_coord_(r.max_span_coord_) 
    {
      std::copy(r.lower_, r.lower_+dim, lower_);
      std::copy(r.upper_, r.upper_+dim, upper_);
    }

    template <class Construct_cartesian_const_iterator_d,class PointPointerIter>
    void update_from_point_pointers(PointPointerIter begin, 
                                    PointPointerIter end,
                                    const Construct_cartesian_const_iterator_d& construct_it
    ) 
    {
      if (begin ==end)
        return;
      // initialize with values of first point
      typename Construct_cartesian_const_iterator_d::result_type bit = construct_it(**begin);
      
      for (int i=0; i < dim; ++i, ++bit) {
        lower_[i]= *bit; upper_[i]=lower_[i];
      }
      begin++;
      typedef typename std::iterator_traits<PointPointerIter>::value_type P;
      std::for_each(begin, end,set_bounds_from_pointer<Construct_cartesian_const_iterator_d,P,T>(dim, lower_, upper_,construct_it));
      set_max_span();
    }
    
    template <class Construct_cartesian_const_iterator_d,class PointPointerIter> // was PointIter
    Kd_tree_rectangle(int d,  PointPointerIter begin,  PointPointerIter end,const Construct_cartesian_const_iterator_d& construct_it)
      : dim(d), lower_(new FT[d]), upper_(new FT[d])
    {
      update_from_point_pointers<Construct_cartesian_const_iterator_d>(begin,end,construct_it);
    }

    inline int 
    max_span_coord() const 
    { 
      return max_span_coord_; 
    }

    inline FT 
    max_span() const 
    {
      return upper_[max_span_coord_] - lower_[max_span_coord_];
    }

    inline FT 
    min_coord(int i) const 
    {
      CGAL_assertion(lower_ != NULL);
      return lower_[i];
    }
    
    inline FT 
    max_coord(int i) const 
    {
      return upper_[i];
    }
    
    std::ostream& 
    print(std::ostream& s) const 
    {
      s << "Rectangle dimension = " << dim;
      s << "\n lower: ";
      for (int i=0; i < dim; ++i)
	s << lower_[i] << " ";
      // std::copy(lower_, lower_ + dim,
      // 	      std::ostream_iterator<FT>(s," "));
      s << "\n upper: ";
      for (int j=0; j < dim; ++j)
	s << upper_[j] << " ";
      // std::copy(upper_, upper_ + dim,
      //	      std::ostream_iterator<FT>(s," "));
      s << "\n maximum span " << max_span() <<
	" at coordinate " << max_span_coord() << std::endl;
      return s;
    }
    
    // Splits rectangle by modifying itself to lower half 
    // and returns upper half
    //    Kd_tree_rectangle* 
    void
    split(Kd_tree_rectangle& r, int d, FT value)
    {
      CGAL_assertion(d >= 0 && d < dim);
      CGAL_assertion(lower_[d] <= value && value <= upper_[d]);
      
      //Kd_tree_rectangle* r = new Kd_tree_rectangle(*this);
      upper_[d]=value;
      r.lower_[d]=value;
      //return r;
    }
    
    
    ~Kd_tree_rectangle() 
    {
      if (dim) {
	if (lower_) delete [] lower_;
	if (upper_) delete [] upper_;
      }
    }
    
    int 
    dimension() const 
    {
      return dim;
    }

    const T* lower() const {return lower_;}
    const T* upper() const {return upper_;}    
 
    Kd_tree_rectangle<FT,Dynamic_dimension_tag>& 
    operator=(const Kd_tree_rectangle<FT,Dynamic_dimension_tag>& r) 
    {
      CGAL_assertion(dimension() == r.dimension());
      if (this != &r) {
        std::copy(r.lower_, r.lower_+dim, lower_);
	std::copy(r.upper_, r.upper_+dim, upper_);
	set_max_span();
      }
      return *this;
    }

  

  }; // of partial specialization of class Kd_tree_rectangle<FT,0>

  template <class FT, typename D>
  std::ostream& 
  operator<<(std::ostream& s, const Kd_tree_rectangle<FT,D>& r) 
  {
    return r.print(s);
  }

  
} // namespace CGAL

#endif // CGAL_KD_TREE_RECTANGLE_H
