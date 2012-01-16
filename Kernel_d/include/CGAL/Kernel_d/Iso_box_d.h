// Copyright (c) 2002  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Authors       : Hans Tangelder <hanst@cs.uu.nl>, Michael Hoffmann

#ifndef CGAL_ISO_BOX_D_H
#define CGAL_ISO_BOX_D_H

#include <CGAL/basic.h>
#include <CGAL/Handle_for.h>
#include <CGAL/representation_tags.h> 
#include <CGAL/Dimension.h>
#include <functional>
#include <algorithm>
#include <numeric>
#include <cstddef>

namespace CGAL {
  
  namespace Kernel_d {
    struct Begin {};
    struct End {};
    struct Cartesian_end {};
  } // namespace Kernel_d
  
  template < typename Point_, typename Functor_ >
  struct Cartesian_iterator {
    typedef Point_                                    Point;
    typedef Functor_                                  Functor;
    typedef typename Point::Cartesian_const_iterator  Iterator;
    typedef Cartesian_iterator<Point,Functor>         Self;

    typedef typename Functor::result_type  value_type;
    typedef value_type&                    reference;
    typedef value_type*                    pointer;
    typedef std::ptrdiff_t                 difference_type;
    typedef std::input_iterator_tag        iterator_category;
      
  protected:
      
    Iterator pb, qb;
    Functor f;
      
  public:
      
    Cartesian_iterator(const Point& p, const Point& q, Kernel_d::Begin)
      : pb(p.cartesian_begin()), qb(q.cartesian_begin())
    {}
    
    Cartesian_iterator(const Point& p, const Point& q, Kernel_d::End)
      : pb(p.cartesian_end()), qb(q.cartesian_end())
    {}
    
    Cartesian_iterator(const Point& p, const Point& q, Kernel_d::Cartesian_end)
      : pb(p.cartesian_end()), qb(q.cartesian_end())
    {}
    
    Cartesian_iterator(const Point& p, const Point& q, const Functor& f_,
		       Kernel_d::Begin) 
      : pb(p.cartesian_begin()), qb(q.cartesian_begin()), f(f_) 
    {}
      
    Cartesian_iterator(const Point& p, const Point& q, const Functor& f_,
		       Kernel_d::End) 
      : pb(p.cartesian_end()), qb(q.cartesian_end()), f(f_) 
    {}
      
    Cartesian_iterator(const Point& p, const Point& q, const Functor& f_,
		       Kernel_d::Cartesian_end) 
      : pb(p.cartesian_end()), qb(q.cartesian_end()), f(f_) 
    {}
      
    Self& operator++() { ++pb; ++qb; return *this; }
      
    Self operator++(int) {
      Self tmp(*this);
      ++(*this);
      return tmp;
    }
    
    value_type operator*()  const { return f(*pb, *qb); }
    pointer    operator->() const { return &(**this); }
      
    const Functor& functor() const { return f; }
    Iterator base_p() const { return pb; }
    Iterator base_q() const { return qb; }
  };
  
  template < typename Iterator, typename Functor > inline
  bool operator==(const Cartesian_iterator<Iterator,Functor>& x, 
		  const Cartesian_iterator<Iterator,Functor>& y)
  {
    return x.base_p() == y.base_p() && x.base_q() == y.base_q(); 
  }
  
  template < typename Iterator, typename Functor > inline
  bool operator!=(const Cartesian_iterator<Iterator,Functor>& x, 
		  const Cartesian_iterator<Iterator,Functor>& y)
  {
    return ! (x == y);
  }
  
  template < typename Point_, typename Functor_ >
  struct Homogeneous_iterator {
    typedef Point_                                      Point;
    typedef Functor_                                    Functor;
    typedef typename Point::Homogeneous_const_iterator  Iterator;
    typedef Homogeneous_iterator<Point,Functor>         Self;

    typedef typename Functor::result_type  value_type;
    typedef value_type&                    reference;
    typedef value_type*                    pointer;
    typedef std::ptrdiff_t                 difference_type;
    typedef std::input_iterator_tag        iterator_category;
      
  protected:
      
    Iterator pb, qb;
    Functor f;
    typedef typename Kernel_traits<Point>::Kernel::RT RT;
    RT hp, hq; // homogenizing coordinates
    
  public:
      
    Homogeneous_iterator(const Point& p, const Point& q, Kernel_d::Begin)
      : pb(p.homogeneous_begin()), qb(q.homogeneous_begin()),
	hp(p.homogeneous(p.dimension())), hq(q.homogeneous(q.dimension()))
    {}
    
    Homogeneous_iterator(const Point& p, const Point& q, Kernel_d::End)
      : pb(p.homogeneous_end()), qb(q.homogeneous_end()),
	hp(p.homogeneous(p.dimension())), hq(q.homogeneous(q.dimension()))
    {}
    
    Homogeneous_iterator(const Point& p, const Point& q, 
			 Kernel_d::Cartesian_end)
      : pb(p.homogeneous_end()), qb(q.homogeneous_end()),
	hp(p.homogeneous(p.dimension())), hq(q.homogeneous(q.dimension()))
    {
      --pb; --qb;
    }
    
    Homogeneous_iterator(const Point& p, const Point& q, const Functor& f_,
			 Kernel_d::Begin) 
      : pb(p.homogeneous_begin()), qb(q.homogeneous_begin()), f(f_),
	hp(p.homogeneous(p.dimension())), hq(q.homogeneous(q.dimension()))
    {}
    
    Homogeneous_iterator(const Point& p, const Point& q, const Functor& f_,
			 Kernel_d::End) 
      : pb(p.homogeneous_end()), qb(q.homogeneous_end()), f(f_),
	hp(p.homogeneous(p.dimension())), hq(q.homogeneous(q.dimension()))
    {}
      
    Homogeneous_iterator(const Point& p, const Point& q, const Functor& f_,
			 Kernel_d::Cartesian_end) 
      : pb(p.homogeneous_end()), qb(q.homogeneous_end()), f(f_),
	hp(p.homogeneous(p.dimension())), hq(q.homogeneous(q.dimension()))
    {
      --pb; --qb;
    }
      
    Self& operator++() { ++pb; ++qb; return *this; }
      
    Self operator++(int) {
      Self tmp(*this);
      ++(*this);
      return tmp;
    }
    
    value_type operator*()  const { return f(*pb * hq, *qb * hp); }
    pointer    operator->() const { return &(**this); }
      
    const Functor& functor() const { return f; }
    Iterator base_p() const { return pb; }
    Iterator base_q() const { return qb; }
  };
  
  template < typename Iterator, typename Functor > inline
  bool operator==(const Homogeneous_iterator<Iterator,Functor>& x, 
		  const Homogeneous_iterator<Iterator,Functor>& y)
  {
    return x.base_p() == y.base_p() && x.base_q() == y.base_q(); 
  }
  
  template < typename Iterator, typename Functor > inline
  bool operator!=(const Homogeneous_iterator<Iterator,Functor>& x, 
		  const Homogeneous_iterator<Iterator,Functor>& y)
  {
    return ! (x == y);
  }
  
  template < typename Kernel_ > class Iso_box_d;
  
  namespace Kernel_d {

    template < typename RepTag > struct Coordinate_iterator;
    
    template <> struct Coordinate_iterator<Cartesian_tag> {
      template < typename Point, typename Functor > 
      struct Iterator {
	typedef Cartesian_iterator<Point,Functor>  type;
      };
    };
    
    template <> struct Coordinate_iterator<Homogeneous_tag> {
      template < typename Point, typename Functor > 
      struct Iterator {
	typedef Homogeneous_iterator<Point,Functor>  type;
      };
    };

    template < typename Kernel_ > 
    struct Iso_box_d_rep {
      typedef Kernel_                   Kernel;
      friend class Iso_box_d<Kernel>;
      
    protected:

      typedef typename Kernel::Point_d  Point_d;
      Point_d lower, upper;

    public:

      Iso_box_d_rep() {}
      
      template < typename InputIteratorI, typename InputIteratorII >
      Iso_box_d_rep(int dim,
		    InputIteratorI b1, InputIteratorI e1,
		    InputIteratorII b2, InputIteratorII e2)
	: lower(dim, b1, e1), upper(dim, b2, e2)
      {}
      
    };

  } // namespace Kernel_d

  template < typename Kernel_ > 
  class Iso_box_d : public Handle_for< Kernel_d::Iso_box_d_rep<Kernel_> > 
  { 

  public:
    typedef Kernel_                   Kernel;
    typedef Kernel_                   R;
    
  protected:

    typedef Kernel_d::Iso_box_d_rep<Kernel>  Rep;
    typedef Handle_for<Rep>                  Base;
    typedef Iso_box_d<Kernel>                Self;

    using Base::ptr;

    typedef typename Kernel::RT       RT;
    typedef typename Kernel::FT       FT;
    typedef typename Kernel::Point_d  Point_d;
    typedef typename Kernel::Rep_tag  Rep_tag;
   
    typedef CGAL::Kernel_d::Coordinate_iterator<Rep_tag>           CIRT;
    typedef typename CIRT::template Iterator<Point_d,Min<RT> >::type  MinIter;
    typedef typename CIRT::template Iterator<Point_d,Max<RT> >::type  MaxIter;

    typedef Kernel_d::Begin            Begin; 
    typedef Kernel_d::End              End; 
    typedef Kernel_d::Cartesian_end    Cartesian_end; 

    RT volume_nominator() const
    {
      typedef typename CIRT::template Iterator<Point_d,std::minus<RT> >::type
	      Iter;
      Iter b(ptr()->upper, ptr()->lower, Begin());
      Iter e(ptr()->upper, ptr()->lower, Cartesian_end());
      return std::accumulate(b, e, RT(1), std::multiplies<RT>());
    }
    
    RT volume_denominator() const
    {
      RT den = 
	ptr()->lower.homogeneous(dimension()) * 
	ptr()->upper.homogeneous(dimension());
      RT prod = den;
      for (int i = 1; i < dimension(); ++i)
	prod *= den;
      return prod;
    }
    
    FT volume(Homogeneous_tag) const
    { 
      return FT(volume_nominator(), volume_denominator());
    }
    
    FT volume(Cartesian_tag) const
    { 
      return volume_nominator();
    }
    
public:

    typedef CGAL::Dynamic_dimension_tag Ambient_dimension;
    typedef CGAL::Dynamic_dimension_tag Feature_dimension;

    Iso_box_d() {}
    
    Iso_box_d(const Point_d& p, const Point_d& q)
      : Base(Rep(p.dimension(), 
		 MinIter(p, q, Begin()), MinIter(p, q, End()),
		 MaxIter(p, q, Begin()), MaxIter(p, q, End())))
    { 
      CGAL_precondition(p.dimension() == q.dimension());
    }
    
    Bounded_side bounded_side(const Point_d& p) const
    { 
      CGAL_precondition(p.dimension() == dimension());
      typedef typename CIRT::template Iterator<Point_d,Compare<RT> >::type
	      Iter;
      
      Iter il(p, ptr()->lower, Begin());
      Iter ilend(p, ptr()->lower, Cartesian_end());
      Iter iu(p, ptr()->upper, Begin());
      CGAL_assertion_code(Iter iuend(p, ptr()->upper, Cartesian_end()));
      
      for (; il != ilend; ++il, ++iu) {
	CGAL_assertion(iu != iuend);
	Comparison_result low = *il;
	Comparison_result upp = *iu;
	if (low == LARGER && upp == SMALLER) continue;
	if (low == SMALLER || upp == LARGER) return ON_UNBOUNDED_SIDE;
	return ON_BOUNDARY;
      }
      return ON_BOUNDED_SIDE;
    }

    bool has_on_bounded_side(const Point_d& p) const
    { 
      return (bounded_side(p)==ON_BOUNDED_SIDE);
    } 

    bool has_on_unbounded_side(const Point_d& p) const
    {
      return (bounded_side(p)==ON_UNBOUNDED_SIDE); 
    } 

    bool has_on_boundary(const Point_d& p) const
    {
      return (bounded_side(p)==ON_BOUNDARY); 
    } 

    int dimension() const { return ptr()->lower.dimension();}
    
    // FIXME!
    FT min_coord(int i) const { return ptr()->lower[i]; }

    FT max_coord(int i) const { return ptr()->upper[i]; }

    const Point_d& min BOOST_PREVENT_MACRO_SUBSTITUTION () const { return ptr()->lower; }

    const Point_d& max BOOST_PREVENT_MACRO_SUBSTITUTION () const { return ptr()->upper; }

    FT volume() const { return volume(Rep_tag()); }
    
    bool is_degenerate() const
    {
      typedef typename CIRT::
	      template Iterator<Point_d,std::equal_to<RT> >::type Iter;
      // omit homogenizing coordinates
      Iter e(ptr()->lower, ptr()->upper, Cartesian_end());
      return 
	e != std::find(Iter(ptr()->lower, ptr()->upper, Begin()), e, true);
    }

}; // end of class

template < typename Kernel >
inline bool
operator==(const Iso_box_d<Kernel>& b1, Iso_box_d<Kernel>& b2)
{
  CGAL_precondition(b1.dimension() == b2.dimension());
  return (b1.min)() == (b2.min)() && (b1.max)() == (b2.max)();
}

template < typename Kernel >
inline bool
operator!=(const Iso_box_d<Kernel>& b1, Iso_box_d<Kernel>& b2)
{
  return ! (b1 == b2); 
}

} // namespace CGAL

#endif // CGAL_ISO_BOX_D_H
