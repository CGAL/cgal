// ============================================================================
//
// Copyright (c) 1998, 1999, 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : Rectangular_p_center_traits_2.h
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : pcenter.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
// coordinator   : ETH
//
// 2-4-Center Computation for Axis-Parallel 2D-Rectangles
// ============================================================================

#if ! (CGAL_RECTANGULAR_P_CENTER_TRAITS_2_H)
#define CGAL_RECTANGULAR_P_CENTER_TRAITS_2_H 1

#include <CGAL/Point_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/basic_constructions_2.h>



CGAL_BEGIN_NAMESPACE

template < class R >
struct Less_x_2
: public CGAL_STD::binary_function< Point_2< R >, Point_2< R >, bool >
{

  bool operator()(const Point_2< R >& p,
                  const Point_2< R >& q) const
  { return p.x() < q.x(); }
};

template < class R >
struct Less_y_2
: public CGAL_STD::binary_function< Point_2< R >, Point_2< R >, bool >
{
  bool operator()(const Point_2< R >& p,
                  const Point_2< R >& q) const
  { return p.y() < q.y(); }
};

/*
template < class R >
struct Greater_x_2
: public CGAL_STD::binary_function< Point_2< R >, Point_2< R >, bool >
{
  bool operator()(const Point_2< R >& p,
                  const Point_2< R >& q) const
  { return p.x() > q.x(); }
};

template < class R >
struct Greater_y_2
: public CGAL_STD::binary_function< Point_2< R >, Point_2< R >, bool >
{
  bool operator()(const Point_2< R >& p,
                  const Point_2< R >& q) const
  { return p.y() > q.y(); }
};
template < class R >
struct Construct_min_point_2
: public CGAL_STD::unary_function< Iso_rectangle_2< R >, Point_2< R > >
{
  Point_2< R > operator()(const Iso_rectangle_2< R >& r) const
  { return r.min(); }
};

template < class R >
struct Construct_max_point_2
: public CGAL_STD::unary_function< Iso_rectangle_2< R >, Point_2< R > >
{
  Point_2< R > operator()(const Iso_rectangle_2< R >& r) const
  { return r.max(); }
};
*/
template < class A, class S >
struct Select : public CGAL_STD::binary_function< A, A, A > {
  Select() {}
  Select(const S& s) : s_(s) {}
  A operator()(const A& a, const A& b) const
  { return s_(a, b) ? a : b; }
  A operator()(const A& a, const A& b)
  { return s_(a, b) ? a : b; }
private:
  S s_;
};

/*
template < class R >
struct Construct_iso_rectangle_2 {
  Iso_rectangle_2< R >
  operator()(const Point_2< R >& min, const Point_2< R >& max) const
  { return Iso_rectangle_2< R >(min, max); }
};
*/
template < class R >
struct Signed_x_distance_2
: public CGAL_STD::binary_function<
  Point_2< R >, Point_2< R >, typename R::FT >
{
  typename R::FT
  operator()(const Point_2< R >& q1, const Point_2< R >& q2) const
  { return q1.x() - q2.x(); }
};
template < class R >
struct Signed_y_distance_2
: public CGAL_STD::binary_function<
  Point_2< R >, Point_2< R >, typename R::FT >
{
  typename R::FT
  operator()(const Point_2< R >& q1, const Point_2< R >& q2) const
  { return q1.y() - q2.y(); }
};
template < class R >
struct Infinity_distance_2
: public CGAL_STD::binary_function<
  Point_2< R >, Point_2< R >, typename R::FT >
{
  typename R::FT
  operator()(const Point_2< R >& q1, const Point_2< R >& q2) const {
    return max(CGAL_NTS abs(q1.x() - q2.x()),
               CGAL_NTS abs(q1.y() - q2.y()));
  }
};
template < class R >
struct Signed_infinity_distance_2
: public CGAL_STD::binary_function<
  Point_2< R >, Point_2< R >, typename R::FT >
{
  typename R::FT
  operator()(const Point_2< R >& q1, const Point_2< R >& q2) const
  { return max(q1.x() - q2.x(), q1.y() - q2.y()); }
};
template < class R >
struct Construct_corner_2
: public CGAL_STD::binary_function<
  Iso_rectangle_2< R >, unsigned int, Point_2< R > >
{
  Point_2< R >
  operator()(const Iso_rectangle_2< R >& q, unsigned int i) const
  { return q[i]; }
};
template < class R >
struct Construct_point_2_above_right_implicit_point_2 {
  // (p, q, r) |--> (p.x() + r, q.y() + r)
  typedef typename R::FT                     FT;
  typedef Point_2< Cartesian< FT > >         P;

  inline P
  operator()(const P& p, const P& q, FT r) const
  { return P(p.x() + r, q.y() + r); }
};

template < class R >
struct Construct_point_2_above_left_implicit_point_2 {
  // (p, q, r) |--> (p.x() - r, q.y() + r)
  typedef typename R::FT                     FT;
  typedef Point_2< Cartesian< FT > >         P;

  inline P
  operator()(const P& p, const P& q, FT r) const
  { return P(p.x() - r, q.y() + r); }
};

template < class R >
struct Construct_point_2_below_left_implicit_point_2 {
  // (p, q, r) |--> (p.x() - r, q.y() - r)
  typedef typename R::FT                     FT;
  typedef Point_2< Cartesian< FT > >         P;

  inline P
  operator()(const P& p, const P& q, FT r) const
  { return P(p.x() - r, q.y() - r); }
};

template < class R >
struct Construct_point_2_below_right_implicit_point_2 {
  // (p, q, r) |--> (p.x() + r, q.y() - r)
  typedef typename R::FT                     FT;
  typedef Point_2< Cartesian< FT > >         P;

  inline P
  operator()(const P& p, const P& q, FT r) const
  { return P(p.x() + r, q.y() - r); }
};
// Point_2 x Point_2 --> Point_2
// (p, q) |-> projection of p onto the horizontal line through q
template < class R >
struct Construct_projection_onto_horizontal_implicit_line_2
: public CGAL_STD::binary_function< Point_2< R >,
                                    Point_2< R >,
                                    Point_2< R > >
{
  typedef CGAL::Point_2< R > Point_2;
  Point_2 operator()(const Point_2& p, const Point_2& q) const
  { return Point_2(p.x(), q.y()); }
};


template < class R >
struct Rectangular_p_center_default_traits_2 : public R
{
  // -----------------------------------------------------------------
  // types:
  //
  typedef typename R::FT               FT;
  typedef CGAL::Point_2< R >           Point_2;
  typedef CGAL::Iso_rectangle_2< R >   Iso_rectangle_2;

  // -----------------------------------------------------------------
  // predicates:
  //

  // from the kernel
  //typedef typename R::Less_x_2         Less_x_2;
  //typedef typename R::Less_y_2         Less_y_2;
  typedef CGAL::Less_x_2< R > Less_x_2;
  typedef CGAL::Less_y_2< R > Less_y_2;
  Less_x_2 less_x_2_object() const { return Less_x_2(); }
  Less_y_2 less_y_2_object() const { return Less_y_2(); }

  // additions
  struct Greater_x_2 : public std::binary_function< Point_2, Point_2, bool >
  {
    Greater_x_2(const Less_x_2& l) : lessx2(l) {}
  
    bool operator()(const Point_2& p, const Point_2& q) const
    { return lessx2(q, p); }
  private:
    Less_x_2 lessx2;
  };
  struct Greater_y_2 : public std::binary_function< Point_2, Point_2, bool >
  {
    Greater_y_2(const Less_y_2& l) : lessy2(l) {}
  
    bool operator()(const Point_2& p, const Point_2& q) const
    { return lessy2(q, p); }
  private:
    Less_y_2 lessy2;
  };
  
  // object methods:
  Greater_x_2 greater_x_2_object() const
  { return Greater_x_2(less_x_2_object()); }
  Greater_y_2 greater_y_2_object() const
  { return Greater_y_2(less_y_2_object()); }

  // -----------------------------------------------------------------
  // constructions:
  //

  // from the kernel
  typedef typename R::Construct_min_point_2     Construct_min_point_2;
  typedef typename R::Construct_max_point_2     Construct_max_point_2;
  typedef typename R::Construct_iso_rectangle_2 Construct_iso_rectangle_2;

  // additions
  /*
  Now from the Kernel (hopefully :)
  typedef Construct_iso_rectangle_2< R > Construct_iso_rectangle_2;
  Construct_iso_rectangle_2
  construct_iso_rectangle_2_object() const
  { return Construct_iso_rectangle_2(); }
  */

  typedef Signed_x_distance_2< R >  Signed_x_distance_2;
  typedef Signed_y_distance_2< R >  Signed_y_distance_2;

  typedef Infinity_distance_2< R >           Infinity_distance_2;
  typedef Signed_infinity_distance_2< R >    Signed_infinity_distance_2;
  typedef Construct_corner_2< R >            Construct_corner_2;

  typedef Construct_projection_onto_horizontal_implicit_line_2< R >
    Construct_projection_onto_horizontal_implicit_line_2;

  typedef Construct_point_2_above_right_implicit_point_2< R >
    Construct_point_2_above_right_implicit_point_2;
  typedef Construct_point_2_above_left_implicit_point_2< R >
    Construct_point_2_above_left_implicit_point_2;
  typedef Construct_point_2_below_right_implicit_point_2< R >
    Construct_point_2_below_right_implicit_point_2;
  typedef Construct_point_2_below_left_implicit_point_2< R >
    Construct_point_2_below_left_implicit_point_2;

  // get object methods:
  Signed_x_distance_2
  signed_x_distance_2_object() const
  { return Signed_x_distance_2(); }
  Signed_y_distance_2
  signed_y_distance_2_object() const
  { return Signed_y_distance_2(); }
  Infinity_distance_2
  infinity_distance_2_object() const
  { return Infinity_distance_2(); }
  Signed_infinity_distance_2
  signed_infinity_distance_2_object() const
  { return Signed_infinity_distance_2(); }
  Construct_corner_2
  construct_corner_2_object() const
  { return Construct_corner_2(); }
  Construct_projection_onto_horizontal_implicit_line_2
  construct_projection_onto_horizontal_implicit_line_2_object() const
  { return Construct_projection_onto_horizontal_implicit_line_2(); }
  Construct_point_2_above_right_implicit_point_2
  construct_point_2_above_right_implicit_point_2_object() const
  { return Construct_point_2_above_right_implicit_point_2(); }
  Construct_point_2_above_left_implicit_point_2
  construct_point_2_above_left_implicit_point_2_object() const
  { return Construct_point_2_above_left_implicit_point_2(); }
  Construct_point_2_below_left_implicit_point_2
  construct_point_2_below_left_implicit_point_2_object() const
  { return Construct_point_2_below_left_implicit_point_2(); }
  Construct_point_2_below_right_implicit_point_2
  construct_point_2_below_right_implicit_point_2_object() const
  { return Construct_point_2_below_right_implicit_point_2(); }

  //!!! this shouldn't be here as it can be written in terms
  // of known stuff
  struct Min_x_2 : public Select< Point_2, Less_x_2 > {
    Min_x_2(const Less_x_2& l) : Select< Point_2, Less_x_2 >(l) {}
  };
  
  struct Min_y_2 : public Select< Point_2, Less_y_2 > {
    Min_y_2(const Less_y_2& l) : Select< Point_2, Less_y_2 >(l) {}
  };
  
  struct Max_x_2 : public Select< Point_2, Greater_x_2 > {
    Max_x_2(const Greater_x_2& l) : Select< Point_2, Greater_x_2 >(l) {}
  };
  
  struct Max_y_2 : public Select< Point_2, Greater_y_2 > {
    Max_y_2(const Greater_y_2& l) : Select< Point_2, Greater_y_2 >(l) {}
  };
  
  Min_x_2 min_x_2_object() const
  { return Min_x_2(less_x_2_object()); }
  Max_x_2 max_x_2_object() const
  { return Max_x_2(greater_x_2_object()); }
  Min_y_2 min_y_2_object() const
  { return Min_y_2(less_y_2_object()); }
  Max_y_2 max_y_2_object() const
  { return Max_y_2(greater_y_2_object()); }
};

template < class Traits_, class PiercingFunction_ >
struct Rectangular_p_center_matrix_search_traits_2 {
  typedef Traits_                        Traits;
  typedef typename Traits::FT            FT;
  typedef typename Traits::Point_2       Point_2;
  typedef PiercingFunction_              PiercingFunction;
  typedef Staircases< Traits >           LD;
  typedef typename LD::size_type         size_type;

  template < class InputIC >
  Rectangular_p_center_matrix_search_traits_2(
    InputIC f, InputIC l, Traits t, const PiercingFunction& p)
  : ld(f, l, t), pf(p)
#if (!defined(CGAL_OPTIMISATION_NO_ASSERTIONS) && \
  !defined(CGAL_NO_ASSERTIONS) && !defined(NDEBUG))
    , ld_size(ld.size())
#endif // optimisation_assertion_code
  {}

  size_type number_of_points() const { return ld.pts.size(); }

  bool operator()(FT v)
  {
    CGAL_optimisation_assertion(ld.size() == ld_size);
    ld.r = v / FT(2);
    bool ok;
    pf(ld, Wastebasket< Point_2 >(), ok);
    CGAL_optimisation_assertion(ld.size() == ld_size);
    return ok;
  }

  template < class OutputIterator >
  OutputIterator operator()(FT v, OutputIterator o, bool& ok)
  {
    CGAL_optimisation_assertion(ld.size() == ld_size);
    ld.r = v / FT(2);
    OutputIterator n = pf(ld, o, ok);
    CGAL_optimisation_assertion(ld.size() == ld_size);
    return n; //pf(ld, o, ok);
  }

private:
  // data members:
  LD                 ld;
  PiercingFunction   pf;
  CGAL_optimisation_assertion_code(typename LD::size_type ld_size;)

  // copying this would be too inefficient
  Rectangular_p_center_matrix_search_traits_2(
    const Rectangular_p_center_matrix_search_traits_2&)
  {}

}; // Rectangular_p_center_matrix_search_traits_2< ... >

template < class ForwardIterator, class Traits >
Iso_rectangle_2<
  typename std::iterator_traits< ForwardIterator >::value_type::R >
bounding_box_2(ForwardIterator f, ForwardIterator l, const Traits& t)
// PRE: f != l.
{
  CGAL_precondition(f != l);
  //typedef typename R::Less_x_2 Less_x_2;
  //typedef typename R::Less_y_2 Less_y_2;
  typedef typename Traits::Less_x_2                  Less_x_2;
  typedef typename Traits::Less_y_2                  Less_y_2;
  typedef typename Traits::Construct_iso_rectangle_2 Rect;
  typedef typename
    Traits::Construct_projection_onto_horizontal_implicit_line_2
    Pohil;

  Less_x_2 lessx = t.less_x_2_object();
  Less_y_2 lessy = t.less_y_2_object();
  Rect rect = t.construct_iso_rectangle_2_object();
  Pohil pohil =
    t.construct_projection_onto_horizontal_implicit_line_2_object();

  ForwardIterator xmin = f;
  ForwardIterator xmax = f;
  ForwardIterator ymin = f;
  ForwardIterator ymax = f;

  while (++f != l) {
    if (lessx(*f, *xmin)) xmin = f;
    if (lessx(*xmax, *f)) xmax = f;
    if (lessy(*f, *ymin)) ymin = f;
    if (lessy(*ymax, *f)) ymax = f;
  }

  return rect(pohil(*xmin, *ymin), pohil(*xmax, *ymax));
} // bounding_box_2(f, l, t)
template < class ForwardIterator >
inline typename
std::iterator_traits< ForwardIterator >::value_type::R::Iso_rectangle_2
bounding_box_2(ForwardIterator f, ForwardIterator l)
// PRE: f != l.
{
  CGAL_precondition(f != l);
  // that is how it is supposed to be ...
  //typedef typename std::iterator_traits< ForwardIterator >::value_type::R
  //  Traits;
  typedef typename std::iterator_traits< ForwardIterator >::value_type::R R;
  typedef Rectangular_p_center_default_traits_2< R > Traits;
  Traits t;
  return bounding_box_2(f, l, t);
} // bounding_box_2(f, l)
template < class Rectangle, class Traits >
inline Rectangle
construct_bounding_box_union_2(const Rectangle& r1,
                               const Rectangle& r2,
                               const Traits& t)
{
  //typedef typename R::Construct_iso_rectangle_2  Rect;
  //typedef typename R::Less_x_2 Less_x_2;
  //typedef typename R::Less_y_2 Less_y_2;
  typedef typename Traits::Construct_iso_rectangle_2  Rect;
  typedef typename Traits::Construct_min_point_2      Minpt;
  typedef typename Traits::Construct_max_point_2      Maxpt;
  typedef typename Traits::Less_x_2                   Less_x_2;
  typedef typename Traits::Less_y_2                   Less_y_2;
  typedef typename
    Traits::Construct_projection_onto_horizontal_implicit_line_2
    Pohil;

  Less_x_2 lessx = t.less_x_2_object();
  Less_y_2 lessy = t.less_y_2_object();
  Rect     rect  = t.construct_iso_rectangle_2_object();
  Minpt    minpt = t.construct_min_point_2_object();
  Maxpt    maxpt = t.construct_max_point_2_object();
  Pohil    pohil =
    t.construct_projection_onto_horizontal_implicit_line_2_object();

#ifdef __BORLANDC__
  typedef typename Traits::Point_2 Point_2;
  Point_2 bpt1 = lessx(minpt(r1), r2.min()) ? minpt(r1) : minpt(r2);
  Point_2 bpt2 = lessy(minpt(r1), r2.min()) ? minpt(r1) : minpt(r2);
  Point_2 bpt3 = lessx(maxpt(r2), maxpt(r1)) ? maxpt(r1) : maxpt(r2);
  Point_2 bpt4 = lessy(maxpt(r2), maxpt(r1)) ? maxpt(r1) : maxpt(r2);
  return rect(pohil(bpt1, bpt2), pohil(bpt3, bpt4));
#else
  return rect(
    pohil(lessx(minpt(r1), r2.min()) ? minpt(r1) : minpt(r2),
          lessy(minpt(r1), r2.min()) ? minpt(r1) : minpt(r2)),
    pohil(lessx(maxpt(r2), maxpt(r1)) ? maxpt(r1) : maxpt(r2),
          lessy(maxpt(r2), maxpt(r1)) ? maxpt(r1) : maxpt(r2)));
#endif
} // construct_bounding_box_union_2(r1, r2, t)
template < class Rectangle >
inline Rectangle
construct_bounding_box_union_2(const Rectangle& r1, const Rectangle& r2)
{
  typename Rectangle::R t;
  return construct_bounding_box_union_2(r1, r2, t);
} // construct_bounding_box_union_2(r1, r2)


CGAL_END_NAMESPACE


#endif // ! (CGAL_RECTANGULAR_P_CENTER_TRAITS_2_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

