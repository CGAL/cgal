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
// file          : include/CGAL/Box.h
// package       : ASPAS
// revision      : 1.4 
// revision_date : 2002/16/08 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================

#ifndef CGAL_BOX_H
#define CGAL_BOX_H
#include <functional>
#include <algorithm>
#include <new>

namespace CGAL {

  template <class T> class Points_container;

  // Borland does not like moving this code into the class Box,
  // it gives the error message
  // [C++ Error] E2402 Illegal base class type: formal type
  // 'std::unary_function<Point &,void>' resolves to
  // 'std::unary_function<Point &,void>'
  template <class Point, class T>
  struct set_bounds : public std::unary_function<Point&, void> {
    int dim;
    T *lower;
    T *upper;
    set_bounds(int d, T *l, T *u) : dim(d), lower(l), upper(u) {}
    void operator() (Point& p) {
		for (int i = 0; i < dim; ++i) {
			if (p[i] < lower[i]) lower[i] = p[i];
 			if (p[i] > upper[i]) upper[i] = p[i];
		}
    }
  };

  // Borland does not like moving this code into the class Box,
  // it gives the error message
  // [C++ Error] E2402 Illegal base class type: formal type
  // 'std::unary_function<Point &,void>' resolves to
  // 'std::unary_function<Point &,void>'
  template <class P, class T>
  struct set_bounds_from_pointer : public std::unary_function<P, void> {
    int dim;
    T *lower;
    T *upper;
    set_bounds_from_pointer(int d, T *l, T *u) :
	dim(d), lower(l), upper(u) {}
    void operator() (P p) {
		for (int i = 0; i < dim; ++i) {
			if ((*p)[i] < lower[i]) lower[i] = (*p)[i];
			if ((*p)[i] > upper[i]) upper[i] = (*p)[i];
		}
    }
  };


  template <class T> class Box {
  public:
    typedef T NT;

  private:

    int dim;
    T *lower_;
    T *upper_;
    int max_span_coord_;

  public:

    void set_upper(const int i, const NT& x) {
      assert(i >= 0 && i < dim);
      assert(x >= lower_[i]);
      upper_[i] = x;
      set_max_span();
    }

    void set_lower(const int i, const NT& x) {
      assert(i >= 0 && i < dim);
      assert(x <= upper_[i]);
      lower_[i] = x;
      set_max_span();
    }
  protected:
    void set_max_span() {
      NT span = upper_[0]-lower_[0];
      max_span_coord_ = 0;
      for (int i = 1; i < dim; ++i) {
	NT tmp = upper_[i] - lower_[i];
	if (span < tmp) {
	  span = tmp;
	  max_span_coord_ = i;
	}
      }
    }

  public:
    //Box(const int d) : dim(d) { lower_ = new T[d]; upper_ = new T[d];}
    Box(const int d) : dim(d), lower_(new NT[d]), upper_(new NT[d])
    {
      std::fill(lower_, lower_ + dim, 0);
      std::fill(upper_, upper_ + dim, 0);
      set_max_span();
    }

    Box() : dim(0), lower_(0), upper_(0) {}

    template <class Iter>
    Box(const int d, Iter begin_lower, Iter end_lower,
	Iter begin_upper, Iter end_upper)
      : dim(d) {

      lower_ = new NT[d];
      upper_ = new NT[d];
      std::copy(begin_lower, end_lower, lower_);
      std::copy(begin_upper, end_upper, upper_);
      set_max_span();
    }

    explicit Box(const Box<NT>& b) : dim(b.dim),
      lower_(new NT[dim]), upper_(new NT[dim]) {
        std::copy(b.lower_, b.lower_+dim, lower_);
	std::copy(b.upper_, b.upper_+dim, upper_);
	set_max_span();
    }

    template <class PointIter>
    Box(const int d,  PointIter begin,  PointIter end)
      : dim(d), lower_(new NT[d]), upper_(new NT[d]) {
	  // initialize with values of first point
	  for (int i=0; i < dim; ++i)
	  {
	    lower_[i]=(*begin)[i]; upper_[i]=(*begin)[i];
	  }
	  begin++;
      typedef typename std::iterator_traits<PointIter>::value_type P;
      std::for_each(begin, end, set_bounds<P,T>(dim, lower_, upper_));
      set_max_span();
    }

    template <class PointPointerIter>
    void update_from_point_pointers(PointPointerIter begin, 
                                    PointPointerIter end, bool empty) {
		if (empty) { // no points
		  for (int i=0; i < dim; ++i)
		  {
			lower_[i]= 1.0; upper_[i]=-1.0;
		  }
		} else {
          // initialize with values of first point
	      for (int i=0; i < dim; ++i)
		  {
	        lower_[i]= (*(*begin))[i]; upper_[i]=(*(*begin))[i];
		  }
	      begin++;
          typedef typename 
	  std::iterator_traits<PointPointerIter>::value_type P;
          std::for_each(begin, end,
		    set_bounds_from_pointer<P,T>(dim, lower_, upper_));
		}
        set_max_span();
    }

    inline const int max_span_coord() const { return max_span_coord_; }

    inline const NT max_span() const {
      return upper_[max_span_coord_] - lower_[max_span_coord_];
    }

    inline const NT lower(int i) const {
      assert (i >= 0 && i < dim);
      return lower_[i];
    }

    inline const NT upper(int i) const {
      assert (i >= 0 && i < dim);
      return upper_[i];
    }

    std::ostream& print(std::ostream& s) {
      s << "Box dimension = " << dim;
      s << "\n lower: ";
      std::copy(lower_, lower_ + dim,
       	      std::ostream_iterator<NT>(s, " "));
      s << "\n upper: ";
      std::copy(upper_, upper_ + dim,
	      std::ostream_iterator<NT>(s, " "));
      s << "\n maximum span " << max_span() <<
      " at coordinate " << max_span_coord() << std::endl;
      return s;
    }

    // Splits box by modifying itself to lower half and returns upper half
    Box* split(int d, NT value) {
		assert(d >= 0 && d < dim);
		assert(lower_[d] <= value && value <= upper_[d]);
		Box* b = new Box(*this);
		upper_[d]=value;
                b->lower_[d]=value;
		return b;
    }
                      

    ~Box() {
      if (dim) {
	if (lower_) delete [] lower_;
	if (upper_) delete [] upper_;
      }
    }

    const int dimension() const {return dim;}

    friend class Points_container<T>;

    Box<NT>& operator= (const Box<NT>& b) {
      assert(dim == b.dim);
      if (this != &b) {
        std::copy(b.lower_, b.lower_+dim, lower_);
	std::copy(b.upper_, b.upper_+dim, upper_);
	set_max_span();
      }
      return *this;
    }
  }; // of class Box

  template <class NT>
    std::ostream& operator<< (std::ostream& s, Box<NT>& box) {
    return box.print(s);
  }

  template <class NT, class Point> bool belongs(const Point& p, 
					      const Box<NT>& b) {
    for (int i = 0; i < b.dimension(); ++i)
      if (p[i] < b.lower(i) || p[i] > b.upper(i)) return 0;
    return 1;
  }

  template <class NT, class Point> 
  NT Min_squared_distance_l2_to_box(const Point& p,
					      const Box<NT>& b) {
	NT distance(0.0);
    for (int i = 0; i < b.dimension(); ++i) {
      if (p[i] < b.lower(i)) distance += (b.lower(i)-p[i])*(b.lower(i)-p[i]);
	  if (p[i] > b.upper(i)) distance += (p[i]-b.upper(i))*(p[i]-b.upper(i));
	}
    return distance;
  }

  template <class NT, class Point> 
  NT Max_squared_distance_l2_to_box(const Point& p,
					      const Box<NT>& b) {
	NT distance(0.0);
    for (int i = 0; i < b.dimension(); ++i) {
      if (p[i] >= (b.lower(i)+b.upper(i))/2.0) 
		  distance += (p[i]-b.lower(i))*(p[i]-b.lower(i)); 
	  else
		  distance += (b.upper(i)-p[i])*(b.upper(i)-p[i]);
	}
    return distance;
  }

  template <class NT, class Point> 
  NT Min_distance_linf_to_box(const Point& p,
					      const Box<NT>& b) {
	NT distance(0.0);
    for (int i = 0; i < b.dimension(); ++i) {
      if (b.lower(i) - p[i] > distance)  distance = b.lower(i)-p[i];
	  if (p[i] - b.upper(i) > distance)  distance = p[i]-b.upper(i);
	}
    return distance;
  }

  template <class NT, class Point> 
  NT Max_distance_linf_to_box(const Point& p,
					      const Box<NT>& b) {
	NT distance(0.0);
    for (int i = 0; i < b.dimension(); ++i) {
      if (p[i] >= (b.lower(i)+b.upper(i))/2.0) 
		  if (p[i] - b.lower(i) > distance)  distance = p[i]-b.lower(i); 
	  else
		  if (b.upper(i) - p[i] > distance)  distance = b.upper(i)-p[i];
	}
    return distance;
  }

} // namespace CGAL
#endif // CGAL_BOX_H
