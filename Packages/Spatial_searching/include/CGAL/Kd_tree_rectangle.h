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
// file          : include/CGAL/Kd_tree_rectangle.h
// package       : ASPAS (3.12)
// maintainer    : Hans Tangelder <hanst@cs.uu.nl>
// revision      : 3.0
// revision_date : 2003/07/10 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================

#ifndef CGAL_KD_TREE_RECTANGLE_H
#define CGAL_KD_TREE_RECTANGLE_H
#include <functional>
#include <algorithm>
#include <new>
#include <cassert>

namespace CGAL {

  template <class SearchTraits, class Point_d, class T>
  struct set_bounds : public std::unary_function<Point_d&, void> {
    int dim;
    T *lower;
    T *upper;
    set_bounds(int d, T *l, T *u) : dim(d), lower(l), upper(u) {}
    void operator() (Point_d& p) {
		T h;
		typename SearchTraits::Construct_cartesian_const_iterator_d construct_it;
                typename SearchTraits::Cartesian_const_iterator_d pit = construct_it(p);
		for (int i = 0; i < dim; ++i, ++pit) {
			h=(*pit); 
			if (h < lower[i]) lower[i] = h;
 			if (h > upper[i]) upper[i] = h;
		}
    }
  };

  template <class SearchTraits, class P, class T>
  struct set_bounds_from_pointer : public std::unary_function<P, void> {
    int dim;
    T *lower;
    T *upper;
    set_bounds_from_pointer(int d, T *l, T *u) :
	dim(d), lower(l), upper(u) {}
    void operator() (P p) {
		T h;
		typename SearchTraits::Construct_cartesian_const_iterator_d construct_it;
                typename SearchTraits::Cartesian_const_iterator_d pit = construct_it(*p);
		for (int i = 0; i < dim; ++i, ++pit) {
			h=(*pit);
			if (h < lower[i]) lower[i] = h;
			if (h > upper[i]) upper[i] = h;
		}
    }
  };


  template <class SearchTraits> class Kd_tree_rectangle {
  public:
    typedef typename SearchTraits::FT FT;
    typedef FT T;

  private:

    int dim;
    T *lower_;
    T *upper_;
    int max_span_coord_;

  public:

    inline void set_upper_bound(int i, const FT& x) {
      assert(i >= 0 && i < dim);
      assert(x >= lower_[i]);
      upper_[i] = x;
      set_max_span();
    }

    inline void set_lower_bound(int i, const FT& x) {
      assert(i >= 0 && i < dim);
      assert(x <= upper_[i]);
      lower_[i] = x;
      set_max_span();
    }

    inline void set_max_span() {
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
    
    Kd_tree_rectangle(int d) : 
    dim(d), lower_(new FT[d]), upper_(new FT[d])
    {
      std::fill(lower_, lower_ + dim, 0);
      std::fill(upper_, upper_ + dim, 0);
      set_max_span();
    }

    Kd_tree_rectangle() : dim(0), lower_(0), upper_(0) {}

    
    explicit Kd_tree_rectangle(const Kd_tree_rectangle<SearchTraits>& r) : dim(r.dim),
      lower_(new FT[dim]), upper_(new FT[dim]) {
        std::copy(r.lower_, r.lower_+dim, lower_);
	std::copy(r.upper_, r.upper_+dim, upper_);
	set_max_span();
    }

    template <class PointIter>
    Kd_tree_rectangle(int d,  PointIter begin,  PointIter end)
      : dim(d), lower_(new FT[d]), upper_(new FT[d]) 
    {
      // initialize with values of first point
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it;
      typename SearchTraits::Cartesian_const_iterator_d bit = construct_it(*begin);
      //      typename SearchTraits::Cartesian_const_iterator_d be = construct_it(*begin,1);

	  for (int i=0; i < dim; ++bit, ++i)
	  {
	    lower_[i]=(*bit); upper_[i]=lower_[i];
	  }
	  begin++;
      typedef typename std::iterator_traits<PointIter>::value_type P;
      std::for_each(begin, end, set_bounds<SearchTraits, P,T>(dim, lower_, upper_));
      set_max_span();
    }

    template <class PointPointerIter>
    void update_from_point_pointers(PointPointerIter begin, 
                                    PointPointerIter end, bool empty) {
		if (empty) { // no points
		  for (int i=0; i < dim; ++i)
		  {
			lower_[i]= FT(1); upper_[i]= FT(-1);
		  }
		} else {
          // initialize with values of first point
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it;
      typename SearchTraits::Cartesian_const_iterator_d bit = construct_it(**begin);
	      for (int i=0; i < dim; ++i, ++bit)
		  {
	        lower_[i]= *bit; upper_[i]=lower_[i];
		  }
	      begin++;
          typedef typename 
	  std::iterator_traits<PointPointerIter>::value_type P;
          std::for_each(begin, end,
		    set_bounds_from_pointer<SearchTraits,P,T>(dim, lower_, upper_));
		}
        set_max_span();
    }

    inline int max_span_coord() const { return max_span_coord_; }

    inline FT max_span() const {
      return upper_[max_span_coord_] - lower_[max_span_coord_];
    }

    inline FT min_coord(int i) const {
      return lower_[i];
    }

    inline FT max_coord(int i) const {
      return upper_[i];
    }

    std::ostream& print(std::ostream& s) {
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
split(Kd_tree_rectangle& r, int d, FT value) {
		// assert(d >= 0 && d < dim);
		// assert(lower_[d] <= value && value <= upper_[d]);
                
  //Kd_tree_rectangle* r = new Kd_tree_rectangle(*this);
		upper_[d]=value;
                r.lower_[d]=value;
		//return r;
    }
                      

    ~Kd_tree_rectangle() {
      if (dim) {
	if (lower_) delete [] lower_;
	if (upper_) delete [] upper_;
      }
    }

    int dimension() const {return dim;}

 
    Kd_tree_rectangle<SearchTraits>& operator= (const Kd_tree_rectangle<SearchTraits>& r) {
      
      if (this != &r) {
        std::copy(r.lower_, r.lower_+dim, lower_);
	std::copy(r.upper_, r.upper_+dim, upper_);
	set_max_span();
      }
      return *this;
    }

  

}; // of class Kd_tree_rectangle

  template <class SearchTraits>
    std::ostream& operator<< (std::ostream& s, Kd_tree_rectangle<SearchTraits>& r) {
    return r.print(s);
  }

  
} // namespace CGAL
#endif // CGAL_KD_TREE_RECTANGLE_H

