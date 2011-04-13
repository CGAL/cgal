// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-0.9-I-06 $
// release_date  : $CGAL_Date: 1998/03/11 $
//
// file          : include/CGAL/Polygon_2.h
// source        : 
// revision      : 1.8a
// revision_date : 13 Mar 1998
// author(s)     : Geert-Jan Giezeman <geert@cs.uu.nl>
//                 Wieger Wesselink
//
// coordinator   : Utrecht University
//
// ============================================================================

#ifndef CGAL_POLYGON_2_H
#define CGAL_POLYGON_2_H

#include <CGAL/basic.h>

#ifdef CGAL_CFG_NO_MEMBER_TEMPLATES
#include <vector>
#include <list>
#endif

#include <iterator>

#include <CGAL/circulator.h>
#include <CGAL/enum.h>

#ifdef CGAL_REP_CLASS_DEFINED
#include <CGAL/Polygon_traits_2.h>
#include <CGAL/Aff_transformation_2.h>
#endif // CGAL_REP_CLASS_DEFINED

#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_2_vertex_circulator.h>
#include <CGAL/Polygon_2_edge_iterator.h>
#include <CGAL/Polygon_2_edge_circulator.h>

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------//
//                          Polygon_2
//-----------------------------------------------------------------------//

#if defined(CGAL_POLYGON_2_CACHED) && !defined(CGAL_POLYGON_2_MOD_ITER)

template <class It>
inline
typename std::iterator_traits<It>::pointer address(It it)
{
    return it.operator->();
}     

template <class Struct>
inline
Struct* address(Struct *p)
{
    return p;
}

template <class Struct>
inline
Struct const* address(Struct const*p)
{
    return p;
}

template <class It>
struct Polygon_vertex_iterator_2
{
    typedef typename std::iterator_traits<It>::value_type value_type;
    typedef typename std::iterator_traits<It>::iterator_category
            iterator_category;
    typedef typename std::iterator_traits<It>::difference_type difference_type;
    typedef value_type const * pointer;
    typedef value_type const & reference;
    // more typedefs
    Polygon_vertex_iterator_2() {}
    Polygon_vertex_iterator_2(It it) : m_it(it) {}
    Polygon_vertex_iterator_2 & operator=(Polygon_vertex_iterator_2 const & it)
        { m_it = it.m_it; return *this;}
    Polygon_vertex_iterator_2 operator=(It const & it)
        { m_it = it; return *this;}
    bool operator==(Polygon_vertex_iterator_2 const &o) const
        { return m_it == o.m_it; }
    bool operator!=(Polygon_vertex_iterator_2 const &o) const
        { return m_it != o.m_it; }
    reference operator*() const {return *m_it;}
    pointer operator->() const { return address(m_it);}
    Polygon_vertex_iterator_2 & operator++()
    	{ ++m_it; return *this;}
    Polygon_vertex_iterator_2 operator++(int)
    	{ Polygon_vertex_iterator_2 result = *this; ++m_it; return result;}
    // Bidirectional Iterator methods
    Polygon_vertex_iterator_2 & operator--()
    	{ --m_it; return *this;}
    Polygon_vertex_iterator_2 operator--(int)
    	{ Polygon_vertex_iterator_2 result = *this; --m_it; return result;}
    // Random Access Iterator methods
    Polygon_vertex_iterator_2 & operator+=(difference_type n)
        { m_it += n; return *this;}
    Polygon_vertex_iterator_2 & operator-=(difference_type n)
        { m_it -= n; return *this;}
    difference_type operator-(Polygon_vertex_iterator_2 o)
        { return m_it - o.m_it;}
    bool operator<(Polygon_vertex_iterator_2 o)
        { return m_it < o.m_it;}
    reference operator[](difference_type n)
        { return m_it[n];}

    // Access to internals
    It implementation_it() const {return m_it;}
    // should be private and friend of Polygon_2
private:
    It m_it;
};

template <class It>
inline Polygon_vertex_iterator_2<It>
operator+(Polygon_vertex_iterator_2<It> it,
          typename Polygon_vertex_iterator_2<It>::difference_type n)
{
    Polygon_vertex_iterator_2<It> tmp = it;
    return tmp += n;
}

template <class It>
inline Polygon_vertex_iterator_2<It>
operator+(typename Polygon_vertex_iterator_2<It>::difference_type n,
          Polygon_vertex_iterator_2<It> it)
{
    Polygon_vertex_iterator_2<It> tmp = it;
    return tmp += n;
}

template <class It>
inline Polygon_vertex_iterator_2<It>
operator-(Polygon_vertex_iterator_2<It> it,
          typename Polygon_vertex_iterator_2<It>::difference_type n)
{
    Polygon_vertex_iterator_2<It> tmp = it;
    return tmp -= n;
}

template <class It>
inline typename Polygon_vertex_iterator_2<It>::difference_type
operator-(Polygon_vertex_iterator_2<It> it1, Polygon_vertex_iterator_2<It> it2)
{
    return it1.implementation_it() - it2.implementation_it();
}

#endif // defined(...CACHED)

template <class Traits_P, class Container_P
        = std::vector<CGAL_TYPENAME_MSVC_NULL Traits_P::Point_2> >
class Polygon_2 {

  public:
    //--------------------------------------------------------
    //             Types
    //--------------------------------------------------------

    typedef Traits_P Traits;
    typedef Container_P Container;

    typedef typename Traits_P::FT FT;
    typedef typename Traits_P::Point_2 Point_2;
    typedef typename Traits_P::Segment_2 Segment_2;

    typedef typename Container_P::difference_type difference_type;
    typedef typename Container_P::value_type value_type;
    typedef typename Container_P::pointer pointer;
    typedef typename Container_P::reference reference;
    typedef typename Container_P::const_reference const_reference;


    //-------------------------------------------------------//
    // this intermediary step is required by Sun C++ 4.1
    typedef typename Container_P::iterator iterator;
    typedef typename Container_P::const_iterator const_iterator;
    //-------------------------------------------------------//

#if defined(CGAL_POLYGON_2_CACHED) && !defined(CGAL_POLYGON_2_MOD_ITER)
    typedef Polygon_vertex_iterator_2<typename Container::iterator>
             Vertex_iterator;
    typename Container::iterator get_container_iterator(Vertex_iterator vit)
    {return vit.implementation_it();}
#else // defined(...CACHED)
    typedef typename Container::iterator Vertex_iterator;
    typename Container::iterator get_container_iterator(Vertex_iterator vit)
    {return vit;}
#endif // defined(...CACHED)
    typedef Vertex_iterator Vertex_const_iterator;

    typedef Polygon_circulator<Container_P>
            Vertex_const_circulator;
    typedef Vertex_const_circulator Vertex_circulator;

    typedef Polygon_2_edge_iterator<Traits_P,Container_P>
            Edge_const_iterator;

    typedef Polygon_2_const_edge_circulator<Traits_P,Container_P>
            Edge_const_circulator;

    //--------------------------------------------------------
    //             Creation
    //--------------------------------------------------------

    Polygon_2(Traits p_traits = Traits()) : traits(p_traits), m_flags(CF_EMPTY)
      { }

    Polygon_2(const Polygon_2<Traits_P,Container_P>& polygon)
      : d_container(polygon.d_container), traits(polygon.traits),
        m_flags(CF_EMPTY) { }

    Polygon_2<Traits_P, Container_P>&
    operator=(const Polygon_2<Traits_P,Container_P>& polygon)
    {
      d_container = polygon.d_container;
      invalidate_cache();
      return *this;
    }

    ~Polygon_2()
      { }

   template <class InputIterator>
    Polygon_2(InputIterator first, InputIterator last,
            Traits p_traits = Traits())
        : d_container(first,last), traits(p_traits), m_flags(CF_EMPTY) {}

    //--------------------------------------------------------
    //             Operations
    //--------------------------------------------------------

    void set(Vertex_iterator pos, const Point_2& x)
     { invalidate_cache(); *get_container_iterator(pos) = x; }

#if defined(CGAL_POLYGON_2_CACHED) && !defined(CGAL_POLYGON_2_MOD_ITER)
    void set(Polygon_circulator<Container>const &pos, const Point_2& x)
     { invalidate_cache();
       *pos.mod_iterator() = x;
     }
#endif

    Vertex_iterator insert(Vertex_iterator pos, const Point_2& x)
      { invalidate_cache();
        return d_container.insert(get_container_iterator(pos),x);
      }

#if defined(CGAL_POLYGON_2_CACHED) && !defined(CGAL_POLYGON_2_MOD_ITER)
    Vertex_iterator insert(Vertex_circulator pos, const Point_2& x)
      { invalidate_cache();
        return d_container.insert(pos.mod_iterator(),x);
      }
#endif

#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
    template <class InputIterator>
    void insert(Vertex_iterator pos,
                InputIterator first,
                InputIterator last)
      { d_container.insert(get_container_iterator(pos), first, last); }

#  if defined(CGAL_POLYGON_2_CACHED) && !defined(CGAL_POLYGON_2_MOD_ITER)
    template <class InputIterator>
    void insert(Vertex_circulator pos,
                InputIterator first,
                InputIterator last)
      { d_container.insert(pos.mod_iterator(), first, last); }
#  endif
#endif

    void push_back(const Point_2& x)
      { invalidate_cache(); d_container.insert(d_container.end(), x); }

    void erase(Vertex_iterator pos)
      { invalidate_cache(); d_container.erase(get_container_iterator(pos)); }

#if defined(CGAL_POLYGON_2_CACHED) && !defined(CGAL_POLYGON_2_MOD_ITER)
    void erase(Vertex_circulator pos)
      { invalidate_cache(); d_container.erase(pos.mod_iterator()); }
#endif

    void erase(Vertex_iterator first, Vertex_iterator last)
      { invalidate_cache();
      d_container.erase(get_container_iterator(first),
                        get_container_iterator(last)); }

    void reverse_orientation()
    {
      if (size() <= 1)
        return;
      invalidate_cache(); 
      typename Container_P::iterator i = d_container.begin();
      std::reverse(++i, d_container.end());
    }

    //--------------------------------------------------------
    //             Traversal of a polygon
    //--------------------------------------------------------

//    Vertex_iterator vertices_begin()
//      { return d_container.begin(); }

//    Vertex_iterator vertices_end()
//      { return d_container.end(); }

    Vertex_const_iterator vertices_begin() const
      { return const_cast<Polygon_2&>(*this).d_container.begin(); }

    Vertex_const_iterator vertices_end() const
      { return const_cast<Polygon_2&>(*this).d_container.end(); }

//    Vertex_const_circulator vertices_circulator() const
//      { return Vertex_const_circulator(&d_container, d_container.begin()); }

    Vertex_const_circulator vertices_circulator() const
      { 
        Polygon_2& self = const_cast<Polygon_2&>(*this);
        return Vertex_const_circulator(&self.d_container,
               self.d_container.begin());
      }

    Edge_const_iterator edges_begin() const
      { return Edge_const_iterator(&d_container, d_container.begin()); }

    Edge_const_iterator edges_end() const
      { return Edge_const_iterator(&d_container, d_container.end()); }

    Edge_const_circulator edges_circulator() const
      { return Edge_const_circulator(vertices_circulator()); }

    //--------------------------------------------------------
    //             Predicates
    //--------------------------------------------------------

    bool is_simple() const
    { if (!is_cached(CF_SIMPLE)) {
         m_simple = is_simple_2(d_container.begin(),d_container.end(), traits);
	 mark_cached(CF_SIMPLE);
      }
      return m_simple;
    }

    bool is_convex() const
    { if (!is_cached(CF_CONVEX)) {
         m_convex = is_convex_2(d_container.begin(),d_container.end(), traits);
	 mark_cached(CF_CONVEX);
      }
      return m_convex;
    }

    Orientation orientation() const
    { if (!is_cached(CF_ORIENTATION)) {
         m_orientation = orientation_2(d_container.begin(),
                                       d_container.end(), traits);
	 mark_cached(CF_ORIENTATION);
      }
      return m_orientation;
    }

    Oriented_side oriented_side(const Point_2& value) const
    {
      return oriented_side_2(d_container.begin(), d_container.end(),
                                  value, traits);
    }

    Bounded_side bounded_side(const Point_2& value) const
    {
      CGAL_polygon_precondition(is_simple());
      return bounded_side_2(d_container.begin(), d_container.end(),
                                 value, traits);
    }

    Bbox_2 bbox() const
    { if (!is_cached(CF_BBOX)) {
         m_bbox = bbox_2(d_container.begin(), d_container.end()); 
	 mark_cached(CF_BBOX);
      }
      return m_bbox;
    }

    FT area() const
    { if (!is_cached(CF_AREA)) {
         area_2(d_container.begin(), d_container.end(), m_area, traits);
	 mark_cached(CF_AREA);
      }
      return m_area;
    }

    Vertex_const_iterator left_vertex() const
    { if (!is_cached(CF_LEFT)) {
         Polygon_2 &self = const_cast<Polygon_2&>(*this);
         m_left = left_vertex_2(self.d_container.begin(),
                                self.d_container.end(), traits);
	 mark_cached(CF_LEFT);
      }
      return m_left;
    }

    //Vertex_iterator left_vertex()
    //{
      //return left_vertex_2(d_container.begin(), d_container.end(), traits);
    //}

    Vertex_const_iterator right_vertex() const
    { if (!is_cached(CF_RIGHT)) {
         Polygon_2 &self = const_cast<Polygon_2&>(*this);
         m_right = right_vertex_2(self.d_container.begin(),
                                  self.d_container.end(), traits);
	 mark_cached(CF_RIGHT);
      }
      return m_right;
    }

//    Vertex_iterator right_vertex()
//    {
//      return right_vertex_2(d_container.begin(), d_container.end(), traits);
//    }

    Vertex_const_iterator top_vertex() const
    { if (!is_cached(CF_TOP)) {
         Polygon_2 &self = const_cast<Polygon_2&>(*this);
         m_top = top_vertex_2(self.d_container.begin(),
                              self.d_container.end(), traits);
	 mark_cached(CF_TOP);
      }
      return m_top;
    }

//    Vertex_iterator top_vertex()
//    {
//      return top_vertex_2(d_container.begin(), d_container.end(), traits);
//    }

    Vertex_const_iterator bottom_vertex() const
    { if (!is_cached(CF_BOTTOM)) {
         Polygon_2 &self = const_cast<Polygon_2&>(*this);
         m_bottom = bottom_vertex_2(self.d_container.begin(),
                                    self.d_container.end(), traits);
	 mark_cached(CF_BOTTOM);
      }
      return m_bottom;
    }

//    Vertex_iterator bottom_vertex()
//    {
//      return bottom_vertex_2(d_container.begin(), d_container.end(), traits);
//    }

    bool is_counterclockwise_oriented() const
      { return orientation() == COUNTERCLOCKWISE; }

    bool is_clockwise_oriented() const
      { return orientation() == CLOCKWISE; }

    bool is_collinear_oriented() const
      { return orientation() == COLLINEAR; }

    bool has_on_positive_side(const Point_2& q) const
      { return oriented_side(q) == ON_POSITIVE_SIDE; }

    bool has_on_negative_side(const Point_2& q) const
      { return oriented_side(q) == ON_NEGATIVE_SIDE; }

    bool has_on_boundary(const Point_2& q) const
      { return bounded_side(q) == ON_BOUNDARY; }

    bool has_on_bounded_side(const Point_2& q) const
      { return bounded_side(q) == ON_BOUNDED_SIDE; }

    bool has_on_unbounded_side(const Point_2& q) const
      { return bounded_side(q) == ON_UNBOUNDED_SIDE; }

    //--------------------------------------------------------
    //             Random access methods
    //--------------------------------------------------------

#ifndef CGAL_CFG_NO_LAZY_INSTANTIATION
    const Point_2& vertex(int i) const
      { return *(d_container.begin() + i); }

//    Point_2& vertex(int i)
//      { return *(d_container.begin() + i); }

    const Point_2& operator[](int i) const
      { return vertex(i); }

//    Point_2& operator[](int i)
//      { return vertex(i); }

    Segment_2 edge(int i) const
      { return *(edges_begin() + i); }
#endif

    //--------------------------------------------------------
    //             Miscellaneous
    //--------------------------------------------------------

    int size() const
      { return d_container.size(); }

    bool is_empty() const
      { return d_container.empty(); }

    const Container_P& container() const
      { return d_container; }

    bool identical(const Polygon_2<Traits_P,Container_P> &q) const
      { return this == &q; }


    Traits_P const &traits_member() const { return traits;}

  private:
    enum Cache_flags {CF_EMPTY=0,
            CF_SIMPLE=1<<0, CF_CONVEX=1<<1, CF_ORIENTATION=1<<2,
            CF_BBOX=1<<3, CF_AREA=1<<4, CF_LEFT=1<<5,
	    CF_RIGHT=1<<6, CF_BOTTOM=1<<7, CF_TOP=1<<8};
    Container_P d_container;
    Traits_P traits;
    // cache
    mutable Cache_flags m_flags;
    mutable Bbox_2 m_bbox;
    mutable FT m_area;
    mutable Vertex_iterator m_left, m_right, m_bottom, m_top;
    mutable bool m_simple :1;
    mutable bool m_convex:1;
    mutable Orientation m_orientation:2;
    void invalidate_cache() { m_flags = CF_EMPTY;}
    bool is_cached(Cache_flags f) const { return m_flags & f;}
#if defined(CGAL_POLYGON_2_CACHED)
    void mark_cached(Cache_flags f) const
       { m_flags = Cache_flags((m_flags & ~f) | f); }
#else
    void mark_cached(Cache_flags ) const
       {}
#endif
};

//-----------------------------------------------------------------------//
//               Globally defined operators
//-----------------------------------------------------------------------//

template <class Traits_P, class Container1_P, class Container2_P>
bool operator==( const Polygon_2<Traits_P,Container1_P> &x,
                 const Polygon_2<Traits_P,Container2_P> &y );

template <class Traits_P, class Container1_P, class Container2_P>
inline
bool
operator!=(const Polygon_2<Traits_P,Container1_P> &x,
           const Polygon_2<Traits_P,Container2_P> &y);

#ifdef CGAL_REP_CLASS_DEFINED

CGAL_END_NAMESPACE

#    include <CGAL/Polygon_traits_2.h>

CGAL_BEGIN_NAMESPACE

template <class Transformation, class Traits_P, class Container_P>
Polygon_2<Traits_P,Container_P>
transform(const Transformation& t, const Polygon_2<Traits_P,Container_P>& p);

#endif // CGAL_REP_CLASS_DEFINED

//-----------------------------------------------------------------------//
//               I/O
//-----------------------------------------------------------------------//

template <class Traits_P, class Container_P>
std::istream &operator>>(std::istream &is, Polygon_2<Traits_P,Container_P>& p);

template <class Traits_P, class Container_P>
std::ostream
&operator<<(std::ostream &os, const Polygon_2<Traits_P,Container_P>& p);

//-----------------------------------------------------------------------//
//                         implementation
//-----------------------------------------------------------------------//

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/Polygon_2.C>
#endif

CGAL_BEGIN_NAMESPACE

template <class Traits_P, class Container1_P, class Container2_P>
inline
bool
operator!=(const Polygon_2<Traits_P,Container1_P> &x,
           const Polygon_2<Traits_P,Container2_P> &y)
{
  return !(x==y);
}

CGAL_END_NAMESPACE

#endif

