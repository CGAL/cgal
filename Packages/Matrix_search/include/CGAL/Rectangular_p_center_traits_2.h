// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
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
struct Construct_min_2
: public CGAL_STD::unary_function< Iso_rectangle_2< R >, Point_2< R > >
{
  Point_2< R > operator()(const Iso_rectangle_2< R >& r) const
  { return r.min(); }
};

template < class R >
struct Construct_max_2
: public CGAL_STD::unary_function< Iso_rectangle_2< R >, Point_2< R > >
{
  Point_2< R > operator()(const Iso_rectangle_2< R >& r) const
  { return r.max(); }
};
template < class A, class S >
struct Select : public CGAL_STD::binary_function< A, A, A > {
  Select() {}
  Select(S& s) : s_(s) {}
  A operator()(const A& a, const A& b) const
  { return s_(a, b) ? a : b; }
  A operator()(const A& a, const A& b)
  { return s_(a, b) ? a : b; }
private:
  S s_;
};

template < class R >
struct Min_x_2 : public Select< Point_2< R >, Less_x_2< R > > {};

template < class R >
struct Min_y_2 : public Select< Point_2< R >, Less_y_2< R > > {};

template < class R >
struct Max_x_2 : public Select< Point_2< R >, Greater_x_2< R > > {};

template < class R >
struct Max_y_2 : public Select< Point_2< R >, Greater_y_2< R > > {};
template < class R >
struct Construct_iso_rectangle_2 {
  Iso_rectangle_2< R >
  operator()(const Point_2< R >& xmin,
             const Point_2< R >& ymin,
             const Point_2< R >& xmax,
             const Point_2< R >& ymax) const
  {
    typedef Point_2< R >         P;
    typedef Iso_rectangle_2< R > Rect;
    return Rect(P(xmin.x(), ymin.y()), P(xmax.x(), ymax.y()));
  }
};
template < class ForwardIterator >
Iso_rectangle_2<
  typename std::iterator_traits< ForwardIterator >::value_type::R >
bounding_box_2(ForwardIterator f, ForwardIterator l)
// PRE: f != l.
{
  CGAL_precondition(f != l);
  typedef
    typename std::iterator_traits< ForwardIterator >::value_type
  Point;
  typedef typename Point::R R;
  typedef Iso_rectangle_2< R > Return_type;

  ForwardIterator xmin = f;
  ForwardIterator xmax = f;
  ForwardIterator ymin = f;
  ForwardIterator ymax = f;

  while (++f != l) {
    if (Less_x_2< R >()(*f, *xmin))
      xmin = f;
    if (Less_x_2< R >()(*xmax, *f))
      xmax = f;
    if (Less_y_2< R >()(*f, *ymin))
      ymin = f;
    if (Less_y_2< R >()(*ymax, *f))
      ymax = f;
  }

  return Construct_iso_rectangle_2< R >()(*xmin, *ymin, *xmax, *ymax);
} //  bounding_box_2( ... )
template < class Rectangle >
inline Rectangle
construct_bounding_box_union(const Rectangle& r1, const Rectangle& r2)
{
  typedef typename Rectangle::R R;
  return Construct_iso_rectangle_2< R >()(
    Less_x_2< R >()(r1.min(), r2.min()) ? r1.min() : r2.min(),
    Less_y_2< R >()(r1.min(), r2.min()) ? r1.min() : r2.min(),
    Less_x_2< R >()(r2.max(), r1.max()) ? r1.max() : r2.max(),
    Less_y_2< R >()(r2.max(), r1.max()) ? r1.max() : r2.max());
}
template < class R >
struct X_distance_2
: public CGAL_STD::binary_function<
  Point_2< R >, Point_2< R >, typename R::FT >
{
  typename R::FT
  operator()(const Point_2< R >& q1, const Point_2< R >& q2) const
  { return CGAL::abs(q1.x() - q2.x()); }
};
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
struct Y_distance_2
: public CGAL_STD::binary_function<
  Point_2< R >, Point_2< R >, typename R::FT >
{
  typename R::FT
  operator()(const Point_2< R >& q1, const Point_2< R >& q2) const
  { return CGAL::abs(q1.y() - q2.y()); }
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
    return std::max(CGAL::abs(q1.x() - q2.x()),
                    CGAL::abs(q1.y() - q2.y()));
  }
};
template < class R >
struct Signed_infinity_distance_2
: public CGAL_STD::binary_function<
  Point_2< R >, Point_2< R >, typename R::FT >
{
  typename R::FT
  operator()(const Point_2< R >& q1, const Point_2< R >& q2) const
  { return std::max(q1.x() - q2.x(), q1.y() - q2.y()); }
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
struct Construct_iso_rectangle_2_above_right_point_2 {
  typedef typename R::FT                     FT;
  typedef Point_2< Cartesian< FT > >         P;
  typedef Iso_rectangle_2< Cartesian< FT > > Rect;

  inline Rect
  operator()(const P& p, FT r) const
  { return Rect(p, P(p.x() + r, p.y() + r)); }
};

template < class R >
struct Construct_iso_rectangle_2_above_left_point_2 {
  typedef typename R::FT                     FT;
  typedef Point_2< Cartesian< FT > >         P;
  typedef Iso_rectangle_2< Cartesian< FT > > Rect;

  inline Rect
  operator()(const P& p, FT r) const
  { return Rect(P(p.x() - r, p.y()), P(p.x(), p.y() + r)); }
};

template < class R >
struct Construct_iso_rectangle_2_below_right_point_2 {
  typedef typename R::FT                     FT;
  typedef Point_2< Cartesian< FT > >         P;
  typedef Iso_rectangle_2< Cartesian< FT > > Rect;

  inline Rect
  operator()(const P& p, FT r) const
  { return Rect(P(p.x(), p.y() - r), P(p.x() + r, p.y())); }
};

template < class R >
struct Construct_iso_rectangle_2_below_left_point_2 {
  typedef typename R::FT                     FT;
  typedef Point_2< Cartesian< FT > >         P;
  typedef Iso_rectangle_2< Cartesian< FT > > Rect;

  inline Rect
  operator()(const P& p, FT r) const
  { return Rect(P(p.x() - r, p.y() - r), p); }
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
struct Rectangular_p_center_default_traits_2 {
  // -----------------------------------------------------------------
  // types:
  //
  typedef typename R::FT               FT;
  typedef CGAL::Point_2< R >           Point_2;
  typedef CGAL::Iso_rectangle_2< R >   Iso_rectangle_2;

  // -----------------------------------------------------------------
  // predicates:
  //
  typedef Less_x_2< R >                Less_x_2;
  typedef Less_y_2< R >                Less_y_2;
  typedef Greater_x_2< R >             Greater_x_2;
  typedef Greater_y_2< R >             Greater_y_2;

  // get object methods:
  Less_x_2    get_less_x_2() const    { return Less_x_2(); }
  Less_y_2    get_less_y_2() const    { return Less_y_2(); }
  Greater_x_2 get_greater_x_2() const { return Greater_x_2(); }
  Greater_y_2 get_greater_y_2() const { return Greater_y_2(); }

  // -----------------------------------------------------------------
  // constructions:
  //
  typedef X_distance_2< R >         X_distance_2;
  typedef Y_distance_2< R >         Y_distance_2;
  typedef Signed_x_distance_2< R >  Signed_x_distance_2;
  typedef Signed_y_distance_2< R >  Signed_y_distance_2;

  typedef Infinity_distance_2< R >          Infinity_distance_2;
  typedef Signed_infinity_distance_2< R >   Signed_infinity_distance_2;

  typedef Construct_min_2< R >         Construct_min_2;
  typedef Construct_max_2< R >         Construct_max_2;
  typedef Construct_corner_2< R >      Construct_corner_2;
  typedef Construct_projection_onto_horizontal_implicit_line_2< R >
    Construct_projection_onto_horizontal_implicit_line_2;
  typedef Construct_iso_rectangle_2< R >  Construct_iso_rectangle_2;

  typedef Construct_iso_rectangle_2_below_left_point_2< R >
    Construct_iso_rectangle_2_below_left_point_2;
  typedef Construct_iso_rectangle_2_above_left_point_2< R >
    Construct_iso_rectangle_2_above_left_point_2;
  typedef Construct_iso_rectangle_2_below_right_point_2< R >
    Construct_iso_rectangle_2_below_right_point_2;
  typedef Construct_iso_rectangle_2_above_right_point_2< R >
    Construct_iso_rectangle_2_above_right_point_2;

  typedef Construct_point_2_above_right_implicit_point_2< R >
    Construct_point_2_above_right_implicit_point_2;
  typedef Construct_point_2_above_left_implicit_point_2< R >
    Construct_point_2_above_left_implicit_point_2;
  typedef Construct_point_2_below_right_implicit_point_2< R >
    Construct_point_2_below_right_implicit_point_2;
  typedef Construct_point_2_below_left_implicit_point_2< R >
    Construct_point_2_below_left_implicit_point_2;

  // get object methods:
  X_distance_2 get_x_distance_2() const { return X_distance_2(); }
  Y_distance_2 get_y_distance_2() const { return Y_distance_2(); }
  Signed_x_distance_2
  get_signed_x_distance_2() const
  { return Signed_x_distance_2(); }
  Signed_y_distance_2
  get_signed_y_distance_2() const
  { return Signed_y_distance_2(); }
  Infinity_distance_2
  get_infinity_distance_2() const
  { return Infinity_distance_2(); }
  Signed_infinity_distance_2
  get_signed_infinity_distance_2() const
  { return Signed_infinity_distance_2(); }
  Construct_min_2
  get_construct_min_2() const
  { return Construct_min_2(); }
  Construct_max_2
  get_construct_max_2() const
  { return Construct_max_2(); }
  Construct_corner_2
  get_construct_corner_2() const
  { return Construct_corner_2(); }
  Construct_iso_rectangle_2
  get_construct_iso_rectangle_2() const
  { return Construct_iso_rectangle_2(); }
  Construct_projection_onto_horizontal_implicit_line_2
  get_construct_projection_onto_horizontal_implicit_line_2() const
  { return Construct_projection_onto_horizontal_implicit_line_2(); }
  Construct_iso_rectangle_2_below_left_point_2
  get_construct_iso_rectangle_2_below_left_point_2() const
  { return Construct_iso_rectangle_2_below_left_point_2(); }
  Construct_iso_rectangle_2_above_left_point_2
  get_construct_iso_rectangle_2_above_left_point_2() const
  { return Construct_iso_rectangle_2_above_left_point_2(); }
  Construct_iso_rectangle_2_below_right_point_2
  get_construct_iso_rectangle_2_below_right_point_2() const
  { return Construct_iso_rectangle_2_below_right_point_2(); }
  Construct_iso_rectangle_2_above_right_point_2
  get_construct_iso_rectangle_2_above_right_point_2() const
  { return Construct_iso_rectangle_2_above_right_point_2(); }
  Construct_point_2_above_right_implicit_point_2
  get_construct_point_2_above_right_implicit_point_2() const
  { return Construct_point_2_above_right_implicit_point_2(); }
  Construct_point_2_above_left_implicit_point_2
  get_construct_point_2_above_left_implicit_point_2() const
  { return Construct_point_2_above_left_implicit_point_2(); }
  Construct_point_2_below_left_implicit_point_2
  get_construct_point_2_below_left_implicit_point_2() const
  { return Construct_point_2_below_left_implicit_point_2(); }
  Construct_point_2_below_right_implicit_point_2
  get_construct_point_2_below_right_implicit_point_2() const
  { return Construct_point_2_below_right_implicit_point_2(); }

  //!!! this shouldn't be here as it can be written in terms
  // of known stuff
  typedef Min_x_2< R >  Min_x_2;
  typedef Max_x_2< R >  Max_x_2;
  typedef Min_y_2< R >  Min_y_2;
  typedef Max_y_2< R >  Max_y_2;
  Min_x_2 get_min_x_2() const { return Min_x_2(); }
  Max_x_2 get_max_x_2() const { return Max_x_2(); }
  Min_y_2 get_min_y_2() const { return Min_y_2(); }
  Max_y_2 get_max_y_2() const { return Max_y_2(); }
};

template < class _Traits, class _PiercingFunction >
struct Rectangular_p_center_matrix_search_traits_2 {
  typedef _Traits                        Traits;
  typedef typename Traits::FT            FT;
  typedef typename Traits::Point_2       Point_2;
  typedef _PiercingFunction              PiercingFunction;
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

CGAL_END_NAMESPACE


#endif // ! (CGAL_RECTANGULAR_P_CENTER_TRAITS_2_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

