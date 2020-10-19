// Copyright (c) 2012, 2020 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>
//            Ahmed Essam <theartful.ae@gmail.com>

#ifndef CGAL_ARRANGEMENTS_DEMO_UTILS_H
#define CGAL_ARRANGEMENTS_DEMO_UTILS_H

#include "ForwardDeclarations.h"
#include "GraphicsSceneMixin.h"
#include "ArrTraitsAdaptor.h"

#include <CGAL/tags.h>
#include <CGAL/Arr_enums.h>
#include <type_traits>

class QGraphicsScene;

template <typename ArrTraits >
class Arr_compute_y_at_x_2 : public GraphicsSceneMixin
{
public:
  typedef ArrTraits Traits;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef typename ArrTraitsAdaptor< Traits >::CoordinateType CoordinateType;
  // typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Kernel::Line_2                       Line_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Multiplicity                 Multiplicity;
  typedef typename Traits::Intersect_2                  Intersect_2;
  typedef std::pair< typename Traits::Point_2, Multiplicity >
                                                        IntersectionResult;

  /*! Constructor */
  Arr_compute_y_at_x_2( const Traits* );

  CoordinateType
  operator()(const X_monotone_curve_2& curve, const CoordinateType& x);

  double approx(const X_monotone_curve_2& curve, const CoordinateType& x);

protected:
  template <typename TTraits>
  CoordinateType operator()(
    const X_monotone_curve_2& curve, const CoordinateType& x,
    const TTraits* traits_, CGAL::Arr_oblivious_side_tag);

  template <typename TTraits>
  CoordinateType operator()(
    const X_monotone_curve_2& curve, const CoordinateType& x,
    const TTraits* traits_, CGAL::Arr_open_side_tag);

protected:
  const Traits* traits;
  Intersect_2 intersectCurves;
};

template <typename Coefficient_>
class Arr_compute_y_at_x_2< CGAL::Arr_algebraic_segment_traits_2<
                              Coefficient_ > > : public GraphicsSceneMixin
{
public:
  typedef Coefficient_ Coefficient;
  typedef CGAL::Arr_algebraic_segment_traits_2< Coefficient > Traits;
  typedef typename Traits::Algebraic_real_1             CoordinateType;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::Intersect_2                  Intersect_2;
  typedef typename Traits::Multiplicity                 Multiplicity;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

  Arr_compute_y_at_x_2(const Traits* traits_) : traits(traits_) { }

  CoordinateType
  operator()(const X_monotone_curve_2& curve, const CoordinateType& x);

  double approx(const X_monotone_curve_2& curve, const CoordinateType& x);

protected:
  X_monotone_curve_2 makeVerticalLine(const CoordinateType& x);
  const Traits* traits;
};

template <
  typename RatKernel, typename AlgKernel, typename NtTraits,
  typename BoundingTraits>
class Arr_compute_y_at_x_2<CGAL::Arr_Bezier_curve_traits_2<
  RatKernel, AlgKernel, NtTraits, BoundingTraits>> : public GraphicsSceneMixin
{
public:
  typedef CGAL::Arr_Bezier_curve_traits_2<
    RatKernel, AlgKernel, NtTraits, BoundingTraits> Traits;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Traits::Rational Rational;
  typedef typename Traits::Algebraic Algebraic;
  typedef typename Traits::Point_2 Point_2;

  Arr_compute_y_at_x_2(const Traits*) { }

  Algebraic operator()(const X_monotone_curve_2& curve, const Rational& x);
  double approx(const X_monotone_curve_2& curve, const Rational& x);
  Algebraic get_t(const X_monotone_curve_2& curve, const Rational& x);
};

template <typename AlgebraicKernel_d_1>
class Arr_compute_y_at_x_2<
  CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>> :
    public GraphicsSceneMixin
{
public:
  typedef CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1> Traits;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Traits::Algebraic_real_1 Algebraic_real_1;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::Rational Rational;

  Arr_compute_y_at_x_2(const Traits*) { }

  Algebraic_real_1
  operator()(const X_monotone_curve_2& curve, const Algebraic_real_1& x);
  Rational
  operator()(const X_monotone_curve_2& curve, const Rational& x);
  double approx(const X_monotone_curve_2& curve, const Rational& x);
};

template <typename ArrTraits>
class Compute_squared_distance_2_base : public GraphicsSceneMixin
{
};

template <typename ArrTraits >
class Compute_squared_distance_2 :
  public Compute_squared_distance_2_base< ArrTraits >
{ };

template <typename Kernel_ >
class Compute_squared_distance_2< CGAL::Arr_segment_traits_2< Kernel_ > > :
  public Compute_squared_distance_2_base<CGAL::Arr_segment_traits_2<Kernel_> >
{
public:
  typedef Kernel_ Kernel;
  typedef CGAL::Arr_segment_traits_2< Kernel > Traits;
  typedef Compute_squared_distance_2_base< Traits > Superclass;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

  double operator() ( const Point_2& p, const X_monotone_curve_2& c ) const;
};

template <typename Kernel_ >
class Compute_squared_distance_2< CGAL::Arr_linear_traits_2< Kernel_ > > :
  public Compute_squared_distance_2_base<CGAL::Arr_linear_traits_2<Kernel_> >
{
public:
  typedef Kernel_ Kernel;
  typedef CGAL::Arr_linear_traits_2< Kernel > Traits;
  typedef Compute_squared_distance_2_base< Traits > Superclass;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Kernel::Ray_2 Ray_2;
  typedef typename Kernel::Line_2 Line_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

  double operator() ( const Point_2& p, const X_monotone_curve_2& c ) const;
};

template <typename Kernel_ >
class Compute_squared_distance_2< CGAL::Arr_polyline_traits_2< Kernel_ > > :
  public Compute_squared_distance_2_base<CGAL::Arr_polyline_traits_2<Kernel_> >
{
public:
  typedef Kernel_ Kernel;
  typedef CGAL::Arr_polyline_traits_2< Kernel > Traits;
  typedef Compute_squared_distance_2_base< Traits > Superclass;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Traits::Curve_2 Curve_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Curve_2::Subcurve_const_iterator Seg_const_it;

  double operator() ( const Point_2& p, const X_monotone_curve_2& c ) const;
};

template <typename RatKernel, typename AlgKernel, typename NtTraits>
class Compute_squared_distance_2< CGAL::Arr_conic_traits_2< RatKernel,
                                                            AlgKernel,
                                                            NtTraits > > :
  public Compute_squared_distance_2_base< CGAL::Arr_conic_traits_2< RatKernel,
                                                                    AlgKernel,
                                                                    NtTraits > >
{
public:
  typedef AlgKernel                                     Kernel;
  typedef CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > Traits;
  typedef Compute_squared_distance_2_base< Traits >     Superclass;
  // _Conic_point_2< AlgKernel > : public AlgKernel::Point_2
  typedef typename Traits::Point_2                      Conic_point_2;
  typedef typename Kernel::FT                           FT;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Kernel::Segment_2                    Segment_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

public: // methods
  double operator() ( const Point_2& p, const X_monotone_curve_2& c ) const;
};

template <
  typename RatKernel, typename AlgKernel, typename NtTraits,
  typename BoundingTraits>
class Compute_squared_distance_2<CGAL::Arr_Bezier_curve_traits_2<
  RatKernel, AlgKernel, NtTraits, BoundingTraits>> :
    public Compute_squared_distance_2_base<CGAL::Arr_Bezier_curve_traits_2<
      RatKernel, AlgKernel, NtTraits, BoundingTraits>>
{
public:
  typedef CGAL::Arr_Bezier_curve_traits_2<
    RatKernel, AlgKernel, NtTraits, BoundingTraits>
    Traits;
  typedef typename ArrTraitsAdaptor<Traits>::Kernel Kernel;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

public: // methods
  double operator()(const Point_2& p, const X_monotone_curve_2& c) const;
};

template <typename AlgebraicKernel_d_1>
class Compute_squared_distance_2<
  CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>> :
    public Compute_squared_distance_2_base<
      CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>>
{
public:
  typedef CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1> Traits;
  typedef typename ArrTraitsAdaptor<Traits>::Kernel Kernel;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

public:
  double operator()(const Point_2& p, const X_monotone_curve_2& c) const;
};

template <typename Coefficient_ >
class Compute_squared_distance_2< CGAL::Arr_algebraic_segment_traits_2<
                                    Coefficient_ > > :
  public Compute_squared_distance_2_base< CGAL::Arr_algebraic_segment_traits_2<
                                            Coefficient_ > >
{
public:
  typedef Coefficient_                                  Coefficient;
  typedef CGAL::Arr_algebraic_segment_traits_2<Coefficient>
                                                        Traits;
  typedef typename Traits::Bound                        FT;
  typedef typename ArrTraitsAdaptor<Traits>::Kernel     Kernel;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::CKvA_2                       CKvA_2;
  typedef std::pair< double, double >                   Coord_2;
  typedef std::vector< Coord_2 >                        Coord_vec_2;
  typedef typename X_monotone_curve_2::Curve_analysis_2 Curve;
  typedef typename Curve::Polynomial_traits_2           Polynomial_traits_2;
  typedef typename Curve::Polynomial_2                  Polynomial_2;
  typedef typename Polynomial_traits_2::Multivariate_content
                                                        Multivariate_content;
  typedef typename Polynomial_traits_2::Substitute      Substitute;
  typedef typename Polynomial_traits_2::
    Construct_innermost_coefficient_const_iterator_range
      ConstructInnerCoeffIter;

public:
  double operator()(const Point_2& p, const X_monotone_curve_2& c) const;
};

// chcek if arrangement is a model of the concept ArrangementOpenBoundaryTraits_2
template <typename ArrTraits>
struct IsOpenBoundaryArrangement :
    public CGAL::Boolean_tag<
      std::is_convertible<
        typename ArrTraits::Left_side_category,
        CGAL::Arr_open_side_tag>::value &&
      std::is_convertible<
        typename ArrTraits::Bottom_side_category,
        CGAL::Arr_open_side_tag>::value &&
      std::is_convertible<
        typename ArrTraits::Top_side_category,
        CGAL::Arr_open_side_tag>::value &&
      std::is_convertible<
        typename ArrTraits::Right_side_category,
        CGAL::Arr_open_side_tag>::value>
{
};

template <typename ArrTraits, typename=void>
class Param_space_in_x_2
{
public:
  typedef typename ArrTraits::X_monotone_curve_2        X_monotone_curve_2;

  Param_space_in_x_2(const ArrTraits*) {}
  CGAL::Arr_parameter_space
  operator()(const X_monotone_curve_2&, CGAL::Arr_curve_end)
  {
    return CGAL::INTERIOR;
  }
};

template <typename ArrTraits>
class Param_space_in_x_2<
  ArrTraits, std::enable_if_t<IsOpenBoundaryArrangement<ArrTraits>::value>>
{
public:
  typedef typename ArrTraits::X_monotone_curve_2        X_monotone_curve_2;
  typedef typename ArrTraits::Parameter_space_in_x_2    Parameter_space_in_x_2;

  Param_space_in_x_2(const ArrTraits* traits) :
      parameter_space_in_x_2(traits->parameter_space_in_x_2_object())
  {
  }

  CGAL::Arr_parameter_space
  operator()(const X_monotone_curve_2& curve, CGAL::Arr_curve_end curve_end)
  {
    return this->parameter_space_in_x_2(curve, curve_end);
  }

private:
  Parameter_space_in_x_2 parameter_space_in_x_2;
};

template <typename ArrTraits>
class Construct_x_monotone_subcurve_2
{
public:
  typedef typename ArrTraits::X_monotone_curve_2        X_monotone_curve_2;
  typedef typename ArrTraits::Split_2                   Split_2;
  typedef typename ArrTraits::Construct_min_vertex_2    Construct_min_vertex_2;
  typedef typename ArrTraits::Construct_max_vertex_2    Construct_max_vertex_2;
  typedef typename ArrTraits::Compare_x_2               Compare_x_2;
  typedef Param_space_in_x_2<ArrTraits>                 Parameter_space_in_x_2;
  typedef typename ArrTraits::Point_2                   Point_2;

  Construct_x_monotone_subcurve_2( const ArrTraits* traits_ );

  /*
    Return the subcurve of curve bracketed by pLeft and pRight.

    We assume pLeft and pRight don't lie on the curve and always do a vertical
    projection.
  */
  X_monotone_curve_2 operator() ( const X_monotone_curve_2& curve,
                                  const boost::optional<Point_2>& pLeft,
                                  const boost::optional<Point_2>& pRight );

protected:
  const ArrTraits* traits;
  Split_2 split_2;
  Compare_x_2 compare_x_2;
  Arr_compute_y_at_x_2< ArrTraits > compute_y_at_x;
  Construct_min_vertex_2 construct_min_vertex_2;
  Construct_max_vertex_2 construct_max_vertex_2;
  Parameter_space_in_x_2 parameter_space_in_x_2;
}; // class Construct_x_monotone_subcurve_2


/*
 * This specialization for conic traits makes use of X_monotone_curve_2::trim,
 * which is not necessarily available.
 */
template <typename RatKernel, typename AlgKernel, typename NtTraits>
class Construct_x_monotone_subcurve_2< CGAL::Arr_conic_traits_2< RatKernel,
                                                                 AlgKernel,
                                                                 NtTraits > >
{
public:
  typedef CGAL::Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>
                                                        ArrTraits;
  typedef typename ArrTraitsAdaptor<ArrTraits>::Kernel  Kernel;
  typedef typename ArrTraits::X_monotone_curve_2        X_monotone_curve_2;
  typedef typename ArrTraits::Split_2                   Split_2;
  typedef typename ArrTraits::Intersect_2               Intersect_2;
  typedef typename ArrTraits::Multiplicity              Multiplicity;
  typedef typename ArrTraits::Construct_min_vertex_2    Construct_min_vertex_2;
  typedef typename ArrTraits::Construct_max_vertex_2    Construct_max_vertex_2;
  typedef typename ArrTraits::Compare_x_2               Compare_x_2;
  typedef typename Kernel::FT                           FT;
  typedef typename ArrTraitsAdaptor< ArrTraits >::CoordinateType
                                                        CoordinateType;
  typedef typename ArrTraits::Point_2                   Point_2;
  typedef typename Kernel::Point_2                      Kernel_point_2;

  Construct_x_monotone_subcurve_2( const ArrTraits* )
  {
  }

  /*
    Return the subcurve of curve bracketed by pLeft and pRight.
  */
  X_monotone_curve_2 operator() ( const X_monotone_curve_2& curve,
                                  const boost::optional<Point_2>& pLeft,
                                  const boost::optional<Point_2>& pRight );

}; // class Construct_x_monotone_subcurve_2 for Arr_conic_traits_2

template <
  typename RatKernel, typename AlgKernel, typename NtTraits,
  typename BoundingTraits>
class Construct_x_monotone_subcurve_2<CGAL::Arr_Bezier_curve_traits_2<
  RatKernel, AlgKernel, NtTraits, BoundingTraits>>
{
public:
  typedef CGAL::Arr_Bezier_curve_traits_2<
    RatKernel, AlgKernel, NtTraits, BoundingTraits>
                                                        ArrTraits;
  typedef typename ArrTraits::X_monotone_curve_2        X_monotone_curve_2;
  typedef typename ArrTraits::Split_2                   Split_2;
  typedef typename ArrTraits::Intersect_2               Intersect_2;
  typedef typename ArrTraits::Multiplicity              Multiplicity;
  typedef typename ArrTraits::Construct_min_vertex_2    Construct_min_vertex_2;
  typedef typename ArrTraits::Construct_max_vertex_2    Construct_max_vertex_2;
  typedef typename ArrTraits::Compare_x_2               Compare_x_2;
  typedef typename ArrTraits::Point_2                   Point_2;

  Construct_x_monotone_subcurve_2(const ArrTraits* traits_);

  /*
    Return the subcurve of curve bracketed by pLeft and pRight.
  */
  X_monotone_curve_2 operator()(
    const X_monotone_curve_2& curve, const boost::optional<Point_2>& pLeft,
    const boost::optional<Point_2>& pRight);

protected:
  const ArrTraits* traits;
  Split_2 split_2;
  Compare_x_2 compare_x_2;
  Arr_compute_y_at_x_2< ArrTraits > compute_y_at_x;
  Construct_min_vertex_2 construct_min_vertex_2;
  Construct_max_vertex_2 construct_max_vertex_2;
}; // class Construct_x_monotone_subcurve_2 for Arr_conic_traits_2

template <typename AlgebraicKernel_d_1>
class Construct_x_monotone_subcurve_2<
  CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>>
{
public:
  typedef CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>
                                                        Traits;
  typedef typename ArrTraitsAdaptor<Traits>::Kernel     Kernel;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Split_2                      Split_2;
  typedef typename Traits::Intersect_2                  Intersect_2;
  typedef typename Traits::Multiplicity                 Multiplicity;
  typedef typename Traits::Construct_min_vertex_2       Construct_min_vertex_2;
  typedef typename Traits::Construct_max_vertex_2       Construct_max_vertex_2;
  typedef typename Traits::Compare_x_2                  Compare_x_2;
  typedef typename Kernel::FT                           FT;
  typedef typename ArrTraitsAdaptor< Traits >::CoordinateType
                                                        CoordinateType;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Kernel::Point_2                      Kernel_point_2;

  Construct_x_monotone_subcurve_2( const Traits* traits_ );

  /*
    Return the subcurve of curve bracketed by pLeft and pRight.

    We assume pLeft and pRight don't lie on the curve and always do a vertical
    projection.
  */
  X_monotone_curve_2 operator() ( const X_monotone_curve_2& curve,
                                  const boost::optional<Point_2>& pLeft,
                                  const boost::optional<Point_2>& pRight );

protected:
  const Traits* traits;
  Split_2 split_2;
  Compare_x_2 compare_x_2;
  Arr_compute_y_at_x_2< Traits > compute_y_at_x;
  Construct_min_vertex_2 construct_min_vertex_2;
  Construct_max_vertex_2 construct_max_vertex_2;
}; // class Construct_x_monotone_subcurve_2

/**
   Converts between Kernel points and Arrangement points.

   The conversion is not necessarily exact.
*/
template <typename ArrTraits>
class Arr_construct_point_2
{
  typedef ArrTraits                                              Traits;
  typedef typename ArrTraits::Point_2                            Point_2;
  typedef typename ArrTraitsAdaptor< ArrTraits >::CoordinateType CoordinateType;
  typedef typename ArrTraitsAdaptor< ArrTraits >::Kernel         Kernel;
  typedef typename Kernel::Point_2                               Kernel_point_2;
  typedef typename Kernel::FT                                    FT;

public:
  Arr_construct_point_2(const Traits* traits_) : traits(traits_) { }

  template <typename P>
  Point_2 operator()(const P& p)
  {
    return this->operator()(p.x(), p.y());
  }

  template <typename T, typename U>
  Point_2 operator()(const T& x, const U& y)
  {
    return this->operator()(FT{x}, FT{y});
  }

  Point_2 operator()(const Kernel_point_2& pt);
  Point_2 operator()(const FT& x, const FT& y);

protected:
  template <typename TTraits >
  Point_2 operator()(const FT& x, const FT& y, const TTraits*);

  template <typename AlgebraicKernel_d_1>
  Point_2 operator()(
    const FT& x, const FT& y,
    const CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>*);

  const Traits* traits;
};

class Find_nearest_edge_base : public GraphicsSceneMixin
{
public:
  /*! Destructor (virtual) */
  virtual ~Find_nearest_edge_base() {}
};

template <typename Arr_>
class Find_nearest_edge : public Find_nearest_edge_base
{
public: // typedefs
  typedef Arr_ Arrangement;
  typedef typename Arrangement::Geometry_traits_2      ArrTraits;
  typedef Compute_squared_distance_2< ArrTraits > Point_curve_distance;
  typedef typename ArrTraits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename ArrTraitsAdaptor<ArrTraits>::Kernel  Kernel;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Arrangement::Face_const_handle       Face_const_handle;
  typedef typename Arrangement::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Arrangement::Vertex_const_handle     Vertex_const_handle;
  typedef typename Arrangement::Hole_const_iterator     Hole_const_iterator;
  typedef typename Arrangement::Ccb_halfedge_const_circulator
    Ccb_halfedge_const_circulator;

public:
  /*! constructor */
  Find_nearest_edge( Arrangement* arr_ ) :
    Find_nearest_edge_base( ),
    arr( arr_ )
  { }

  /*! Destructor (virtual) */
  virtual ~Find_nearest_edge() {}

public: // member methods
  Halfedge_const_handle operator()( const Point_2& queryPt );

  virtual void setScene( QGraphicsScene* scene_ )
  {
    this->pointCurveDistance.setScene( scene_ );
    Find_nearest_edge_base::setScene( scene_ );
  }

protected: // member methods
  Face_const_handle getFace( const CGAL::Object& obj );

protected: // member fields
  Arrangement* arr;
  Point_curve_distance pointCurveDistance;

}; // class Find_nearest_edge

template <typename Arr_>
class Insert_curve
{
public:
  typedef Arr_                                                   Arrangement;
  typedef typename Arrangement::Geometry_traits_2                ArrTraits;
  typedef typename ArrTraits::Curve_2                            Curve_2;

  void operator()(Arrangement*, const Curve_2&);
};

// free functions gathered here to speed up compilation of other files
// specializing once in one file is better than in multiple files
CGAL::Object createArrangement(demo_types::TraitsType);
void deleteArrangement(demo_types::TraitsType, const CGAL::Object&);
CGAL::Object makeOverlayArrangement(const std::vector<CGAL::Object>&);
void insertCurve(
  demo_types::TraitsType, const CGAL::Object& arr, const CGAL::Object& curve);

#endif // CGAL_ARRANGEMENTS_DEMO_UTILS_H
