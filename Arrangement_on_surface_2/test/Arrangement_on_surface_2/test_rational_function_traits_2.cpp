#include <iostream>

#include <CGAL/config.h>

#if !defined(CGAL_USE_CORE)
int main()
{
//  bool   UNTESTED_TRAITS_AS_CORE_IS_NOT_INSTALLED;
  std::cout << std::endl
            << "NOTE: Core is not installed, "
            << "skipping the test ..."
            << std::endl;
  return 0;
}
#else

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d_1.h>             //Algebraic Kernel
#include <CGAL/Arr_rational_function_traits_2.h>   //Traits
#include <CGAL/Arrangement_2.h>                    //Arrangement
#include <CGAL/Surface_sweep_2_algorithms.h>

typedef CGAL::CORE_arithmetic_kernel::Integer      Number_type;
typedef CGAL::Algebraic_kernel_d_1<Number_type>    AK1;
typedef CGAL::Arr_rational_function_traits_2<AK1>  Traits_2;
typedef CGAL::Arrangement_2<Traits_2>              Arrangement_2;
typedef Traits_2::Point_2                          Point_2;

template <typename Cmp_object>
struct Cmp {
  Cmp<Cmp_object>& operator=(const  Cmp<Cmp_object>&);
  const Cmp_object& m_cmp_object;
  Cmp(const Cmp_object& cmp_object) : m_cmp_object(cmp_object) {}
  bool operator()(const Point_2& p1, const Point_2& p2) const
  { return (m_cmp_object(p1, p2) == CGAL::LARGER); }
};

int main()
{
  // testing that all types are present:

  // Traits
  typedef Traits_2::Polynomial_1         Polynomial_1;
  typedef Traits_2::Algebraic_real_1     Algebraic_real_1;
  typedef Traits_2::Coefficient          Coefficient;
  typedef Traits_2::Bound                Bound;
  typedef Traits_2::Algebraic_kernel_d_1 Algebraic_kernel_d_1;
  typedef Traits_2::Construct_curve_2 Construct_curve_2;
  typedef Traits_2::Construct_x_monotone_curve_2 Construct_x_monotone_curve_2;

  // typedef induced by concept

  typedef Traits_2::Curve_2 Curve_2;
  typedef Traits_2::X_monotone_curve_2 X_monotone_curve_2;
  typedef Traits_2::Point_2 Point_2;

  typedef Traits_2::Compare_x_2 Compare_x_2;
  typedef Traits_2::Compare_xy_2 Compare_xy_2;
  typedef Traits_2::Construct_min_vertex_2 Construct_min_vertex_2;
  typedef Traits_2::Construct_max_vertex_2 Construct_max_vertex_2;
  typedef Traits_2::Is_vertical_2 Is_vertical_2;
  typedef Traits_2::Compare_y_at_x_2 Compare_y_at_x_2;
  typedef Traits_2::Compare_y_at_x_left_2 Compare_y_at_x_left_2;
  typedef Traits_2::Compare_y_at_x_right_2 Compare_y_at_x_right_2;
  typedef Traits_2::Equal_2 Equal_2;
  typedef Traits_2::Parameter_space_in_x_2 Parameter_space_in_x_2;
  typedef Traits_2::Parameter_space_in_y_2 Parameter_space_in_y_2;
  typedef Traits_2::Compare_x_at_limit_2 Compare_x_at_limit_2;
  typedef Traits_2::Compare_x_near_limit_2 Compare_x_near_limit_2;
  typedef Traits_2::Compare_y_near_boundary_2 Compare_y_near_boundary_2;
  typedef Traits_2::Intersect_2 Intersect_2;
  typedef Traits_2::Split_2 Split_2;
  typedef Traits_2::Are_mergeable_2 Are_mergeable_2;
  typedef Traits_2::Merge_2 Merge_2;
  typedef Traits_2::Make_x_monotone_2 Make_x_monotone_2;
  typedef Traits_2::Approximate_2 Approximate_2;

  CGAL_USE_TYPE(Coefficient);
  CGAL_USE_TYPE(Bound);

  CGAL_USE_TYPE(Compare_x_2);
  CGAL_USE_TYPE(Compare_xy_2);
  CGAL_USE_TYPE(Construct_min_vertex_2);
  CGAL_USE_TYPE(Construct_max_vertex_2);
  CGAL_USE_TYPE(Is_vertical_2);
  CGAL_USE_TYPE(Compare_y_at_x_2);
  CGAL_USE_TYPE(Compare_y_at_x_left_2);
  CGAL_USE_TYPE(Compare_y_at_x_right_2);
  CGAL_USE_TYPE(Equal_2);
  CGAL_USE_TYPE(Parameter_space_in_x_2);
  CGAL_USE_TYPE(Parameter_space_in_y_2);
  CGAL_USE_TYPE(Compare_x_at_limit_2);
  CGAL_USE_TYPE(Compare_x_near_limit_2);
  CGAL_USE_TYPE(Compare_y_near_boundary_2);
  CGAL_USE_TYPE(Intersect_2);
  CGAL_USE_TYPE(Split_2);
  CGAL_USE_TYPE(Are_mergeable_2);
  CGAL_USE_TYPE(Merge_2);
  CGAL_USE_TYPE(Make_x_monotone_2);
  CGAL_USE_TYPE(Approximate_2);

  // construction traits
  // default construction
  {
    Traits_2 traits;
    const Algebraic_kernel_d_1* ak (traits.algebraic_kernel_d_1()); (void) ak;
  }

  // construction from pointer to AK
  Algebraic_kernel_d_1 ak;
  Traits_2 traits(&ak);

  traits.cleanup_cache();

  Construct_curve_2 construct_curve_2 = traits.construct_curve_2_object();

  Polynomial_1 P(0);
  Polynomial_1 Q(1);
  Algebraic_real_1 one(1);
  Algebraic_real_1 two(2);
  {
    typedef Construct_curve_2::Polynomial_1 Polynomial_1;
    typedef Construct_curve_2::Algebraic_real_1 Algebraic_real_1;
    typedef Construct_curve_2::Curve_2 Curve_2;
    typedef Construct_curve_2::argument_type argument_type;
    typedef Construct_curve_2::first_argument_type first_argument_type;
    typedef Construct_curve_2::second_argument_type second_argument_type;

    CGAL_USE_TYPE(Polynomial_1);
    CGAL_USE_TYPE(Algebraic_real_1);
    CGAL_USE_TYPE(Curve_2);
    CGAL_USE_TYPE(argument_type);
    CGAL_USE_TYPE(first_argument_type);
    CGAL_USE_TYPE(second_argument_type);
  }
  {Curve_2 curve = construct_curve_2(P);}
  {Curve_2 curve = construct_curve_2(P,one,true);}
  {Curve_2 curve = construct_curve_2(P,one,false);}
  {Curve_2 curve = construct_curve_2(P,one,two);}
  {Curve_2 curve = construct_curve_2(P,Q);}
  {Curve_2 curve = construct_curve_2(P,Q,one,true);}
  {Curve_2 curve = construct_curve_2(P,Q,one,false);}
  //{Curve_2 curve = construct_curve_2(P,Q,one,two);}
  {Curve_2 curve = construct_curve_2(P.begin(),P.end());}
  {Curve_2 curve = construct_curve_2(P.begin(),P.end(),one,true);}
  {Curve_2 curve = construct_curve_2(P.begin(),P.end(),one,false);}
  {Curve_2 curve = construct_curve_2(P.begin(),P.end(),one,two);}
  {Curve_2 curve = construct_curve_2(P.begin(),P.end());}
  {Curve_2 curve = construct_curve_2(P.begin(),P.end(),Q.begin(),Q.end(),one,true);}
  {Curve_2 curve = construct_curve_2(P.begin(),P.end(),Q.begin(),Q.end(),one,false);}
  {Curve_2 curve = construct_curve_2(P.begin(),P.end(),Q.begin(),Q.end(),one,two);}

  Construct_x_monotone_curve_2
    construct_x_monotone_curve_2 = traits.construct_x_monotone_curve_2_object();
  {
    typedef Construct_x_monotone_curve_2::Polynomial_1 Polynomial_1;
    typedef Construct_x_monotone_curve_2::Algebraic_real_1 Algebraic_real_1;
    typedef Construct_x_monotone_curve_2::X_monotone_curve_2
      X_monotone_curve_2;
    typedef Construct_x_monotone_curve_2::argument_type argument_type;
    typedef Construct_x_monotone_curve_2::first_argument_type
      first_argument_type;
    typedef Construct_x_monotone_curve_2::second_argument_type
      second_argument_type;

    CGAL_USE_TYPE(Polynomial_1);
    CGAL_USE_TYPE(Algebraic_real_1);
    CGAL_USE_TYPE(X_monotone_curve_2);
    CGAL_USE_TYPE(argument_type);
    CGAL_USE_TYPE(first_argument_type);
    CGAL_USE_TYPE(second_argument_type);
  }

  {X_monotone_curve_2 x_monotone_curve = construct_x_monotone_curve_2(P);}
  {X_monotone_curve_2 x_monotone_curve = construct_x_monotone_curve_2(P,one,true);}
  {X_monotone_curve_2 x_monotone_curve = construct_x_monotone_curve_2(P,one,false);}
  {X_monotone_curve_2 x_monotone_curve = construct_x_monotone_curve_2(P,one,two);}
  {X_monotone_curve_2 x_monotone_curve = construct_x_monotone_curve_2(P,Q);}
  {X_monotone_curve_2 x_monotone_curve = construct_x_monotone_curve_2(P,Q,one,true);}
  {X_monotone_curve_2 x_monotone_curve = construct_x_monotone_curve_2(P,Q,one,false);}
  {X_monotone_curve_2 x_monotone_curve = construct_x_monotone_curve_2(P,Q,one,two);}
  {X_monotone_curve_2 x_monotone_curve = construct_x_monotone_curve_2(P.begin(),P.end());}
  {X_monotone_curve_2 x_monotone_curve = construct_x_monotone_curve_2(P.begin(),P.end(),one,true);}
  {X_monotone_curve_2 x_monotone_curve = construct_x_monotone_curve_2(P.begin(),P.end(),one,false);}
  {X_monotone_curve_2 x_monotone_curve = construct_x_monotone_curve_2(P.begin(),P.end(),one,two);}
  {X_monotone_curve_2 x_monotone_curve = construct_x_monotone_curve_2(P.begin(),P.end());}
  {X_monotone_curve_2 x_monotone_curve = construct_x_monotone_curve_2(P.begin(),P.end(),Q.begin(),Q.end(),one,true);}
  {X_monotone_curve_2 x_monotone_curve = construct_x_monotone_curve_2(P.begin(),P.end(),Q.begin(),Q.end(),one,false);}
  {X_monotone_curve_2 x_monotone_curve = construct_x_monotone_curve_2(P.begin(),P.end(),Q.begin(),Q.end(),one,two);}

  {
    typedef Curve_2::Polynomial_1 Polynomial_1;
    typedef Curve_2::Algebraic_real_1 Algebraic_real_1;

    CGAL_USE_TYPE(Polynomial_1);
    CGAL_USE_TYPE(Algebraic_real_1);
  }

  {
    Curve_2 curve = construct_curve_2(P,Q);
    {Curve_2 dummy;}
    {Curve_2 dummy(curve);}
    {Curve_2 dummy = curve;}

    assert(CGAL::degree(curve.numerator())  >=0);
    assert(CGAL::degree(curve.denominator())>=0);
    assert(curve.is_continuous());
    assert(curve.left_parameter_space_in_x()==CGAL::ARR_LEFT_BOUNDARY);
    assert(curve.right_parameter_space_in_x()==CGAL::ARR_RIGHT_BOUNDARY);
  }
  {
    Curve_2 curve = construct_curve_2(P,one,two);
    assert(CGAL::degree(curve.numerator())  >=0);
    assert(CGAL::degree(curve.denominator())>=0);
    assert(curve.is_continuous());
    assert(curve.left_parameter_space_in_x()==CGAL::ARR_INTERIOR);
    assert(curve.right_parameter_space_in_x()==CGAL::ARR_INTERIOR);
    assert(one == curve.left_x());
    assert(two == curve.right_x());
  }

  // X_monotone_curve_2
  {
    typedef X_monotone_curve_2::Polynomial_1 Polynomial_1;
    typedef X_monotone_curve_2::Algebraic_real_1 Algebraic_real_1;
    typedef X_monotone_curve_2::Point_2 Point_2;

    CGAL_USE_TYPE(Polynomial_1);
    CGAL_USE_TYPE(Algebraic_real_1);
    CGAL_USE_TYPE(Point_2);

    X_monotone_curve_2 xcurve= construct_x_monotone_curve_2(P,Q,one,two);

    {X_monotone_curve_2 dummy;}
    {X_monotone_curve_2 dummy(xcurve);}
    {X_monotone_curve_2 dummy = xcurve;}

    assert(CGAL::degree(xcurve.numerator()) >= 0);
    assert(CGAL::degree(xcurve.denominator()) >=0);
    assert(xcurve.source_parameter_space_in_x() == CGAL::ARR_INTERIOR);
    assert(xcurve.source_parameter_space_in_y() == CGAL::ARR_INTERIOR);
    xcurve.source();
    xcurve.target();
    xcurve.left();
    xcurve.right();
    xcurve.source_x();
    xcurve.target_x();
    xcurve.left_x();
    xcurve.right_x();
    assert(xcurve.is_left_to_right());
  }

  // Point_2
  {
    typedef Point_2::Polynomial_1 Polynomial_1;
    typedef Point_2::Algebraic_real_1 Algebraic_real_1;
    typedef Point_2::Bound Bound;

    CGAL_USE_TYPE(Polynomial_1);
    CGAL_USE_TYPE(Algebraic_real_1);
    CGAL_USE_TYPE(Bound);

    X_monotone_curve_2 xcurve= construct_x_monotone_curve_2(P, Q, one, two);
    Point_2 p = xcurve.left();
    Point_2 q = xcurve.right();

    p.numerator();
    p.denominator();
    p.to_double();
    p.x();
    p.y();
    p.approximate_absolute_x(1);
    p.approximate_absolute_y(1);
    p.approximate_relative_x(1);
    p.approximate_relative_y(1);

    {Point_2 dummy = p;}
    {Point_2 dummy(p);}
  }

  {
    Arrangement_2 arr;
    const Traits_2* traits = arr.traits();(void) traits;
  }

  {
    const Traits_2 traits;
    traits.cleanup_cache();
    assert(traits.cache().rat_func_map().size()==0);
    assert(traits.cache().rat_pair_map().size()==0);

    Construct_curve_2 construct_curve_2 = traits.construct_curve_2_object();

    Polynomial_1 x = CGAL::shift(Polynomial_1(1), 1);

    {
      std::vector<Curve_2> curves;
      std::vector<X_monotone_curve_2> xcurves;
      std::vector<Point_2> points;
      curves.push_back(construct_curve_2(Polynomial_1(1)));
      curves.push_back(construct_curve_2(x*x-2));
      curves.push_back(construct_curve_2(x*x*x));
      curves.push_back(construct_curve_2(x*x*x, x*x-2));

      for(const Curve_2& curve : curves){
        assert(CGAL::degree(curve.numerator()) >= 0);
      }
      CGAL::compute_subcurves(curves.begin(),curves.end(),
                              std::back_inserter(xcurves),false,traits);
      for(const X_monotone_curve_2& xcurve : xcurves) {
        assert(CGAL::degree(xcurve.numerator()) >= 0);
      }

      CGAL::compute_intersection_points(curves.begin(),curves.end(),
                                        std::back_inserter(points), false,
                                        traits);
      for(const Point_2& point : points) {
        assert(CGAL::degree(point.numerator()) >= 0);
      }
    }
    traits.cleanup_cache();
    assert(traits.cache().rat_func_map().size() == 0);
    assert(traits.cache().rat_pair_map().size() == 0);
    {
      std::vector<Curve_2> curves;
      std::vector<X_monotone_curve_2> xcurves;
      std::vector<Point_2> points;

      curves.push_back(construct_curve_2(Polynomial_1(1)));
      curves.push_back(construct_curve_2(x*x-2));
      curves.push_back(construct_curve_2(x*x*x));
      curves.push_back(construct_curve_2(x*x*x, x*x-2));

      traits.cleanup_cache();
      for(const Curve_2& curve : curves){
        assert(CGAL::degree(curve.numerator()) >= 0);
      }
      CGAL::compute_subcurves(curves.begin(), curves.end(),
                              std::back_inserter(xcurves), false, traits);
      for(const X_monotone_curve_2& xcurve : xcurves) {
        assert(CGAL::degree(xcurve.numerator()) >= 0);
      }
      CGAL::compute_intersection_points(curves.begin(), curves.end(),
                                        std::back_inserter(points), false,
                                        traits);
      for(const Point_2& point : points) {
        assert(CGAL::degree(point.numerator()) >= 0);
      }

      // For some reason some versions of MSVC don't like the auto
      // conversion from Comparison_result to bool.
      //std::sort(points.begin(), points.end(), traits.compare_xy_2_object());
      //std::sort(points.begin(), points.end(), traits.compare_x_2_object());
      std::sort(points.begin(), points.end(),
                Cmp<Traits_2::Compare_x_2>(traits.compare_x_2_object()));
      std::sort(points.begin(), points.end(),
                Cmp<Traits_2::Compare_xy_2>(traits.compare_xy_2_object()));
    }
  }
  return 0;
}

#endif
