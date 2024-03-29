#ifndef CGAL_FILTERED_RATIONAL_H
#define CGAL_FILTERED_RATIONAL_H

#include <CGAL/Exact_rational.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Algebraic_structure_traits.h>

namespace CGAL {

struct Filtered_rational;

template <>
class Algebraic_structure_traits<Filtered_rational>
    : public Algebraic_structure_traits_base<Filtered_rational, Field_tag>
{
  public:
    using Is_exact = Tag_true;
    using Is_numerical_sensitive = Tag_false;
};

struct Filtered_rational : boost::totally_ordered1<Filtered_rational
                                                     //#ifdef _MSC_VER
                 , boost::ordered_field_operators2<Filtered_rational, int
                 , boost::ordered_field_operators2<Filtered_rational, double
                 > >
                                                     //#endif
                 >
{
  Interval_nt<> i;
  Exact_rational rat;

  Filtered_rational()
  {}

  Filtered_rational(int d)
    : i(d), rat(d)
  {}

  Filtered_rational(double d)
    : i(d), rat(d)
  {}

  explicit Filtered_rational(const Exact_rational& rat)
    : i(to_interval(rat)), rat(rat)
  {}


  Filtered_rational(const Interval_nt<>& i, const Exact_rational& rat)
      : i(i), rat(rat)
  {}

  operator Exact_rational() const
  {
    return rat;
  }

 Filtered_rational operator-() const
 {
   return Filtered_rational(-i, -rat);
 }

 Filtered_rational operator+(const Filtered_rational& b) const
 {
     return Filtered_rational(i + b.i, rat + b.rat);
 }


 Filtered_rational operator-(const Filtered_rational& b) const
 {
     return Filtered_rational(i - b.i, rat - b.rat);
 }

 Filtered_rational operator*(const Filtered_rational& b) const
 {
     return Filtered_rational(i * b.i, rat * b.rat);
 }

 Filtered_rational operator/(const Filtered_rational& b) const
 {
     return Filtered_rational(i / b.i, rat / b.rat);
 }



  Filtered_rational& operator+=(const Filtered_rational &a) { return *this = *this + a; }
  Filtered_rational& operator-=(const Filtered_rational &a) { return *this = *this - a; }
  Filtered_rational& operator*=(const Filtered_rational &a) { return *this = *this * a; }
  Filtered_rational& operator/=(const Filtered_rational &a) { return *this = *this / a; }

};


template <> class Real_embeddable_traits< Filtered_rational >
  : public INTERN_RET::Real_embeddable_traits_base< Filtered_rational, CGAL::Tag_true > {
  public:
    class Sgn : public CGAL::cpp98::unary_function<Type, ::CGAL::Sign>
    {
    public:
      ::CGAL::Sign operator()(const Type& x) const
      {
        Uncertain<Sign> u = sign(x.i);
        if(is_indeterminate(u)) {
          return CGAL::sign(x.rat);
        }
        return make_certain(u);
      }
    };

    class To_double : public CGAL::cpp98::unary_function<Type, double>
    {
    public:
      double operator()(const Type& x) const { return to_double(x.i); }
    };

    class To_interval : public CGAL::cpp98::unary_function<Type, std::pair<double, double>>
    {
    public:
      std::pair<double, double> operator()(const Type& x) const { return { x.i.inf(), x.i.sup() }; }
    };

    class Compare : public CGAL::cpp98::binary_function<Type, Type, ::CGAL::Comparison_result>
    {
    public:
      ::CGAL::Comparison_result operator()(const Type& x, const Type& y) const
      {
        Uncertain<Comparison_result> u = compare(x.i, y.i);
        if(is_indeterminate(u)) {
          return CGAL::compare(x.rat, y.rat);
        }
        return make_certain(u);
      }
    };

    class Square : public CGAL::cpp98::unary_function<Type, Type>
    {
    public:
      Type operator()(const Type& x) const
      {
        return Filtered_rational(square(x.i), square(x.rat));
      }
    };
};

bool operator==(const Filtered_rational& a,
                const Filtered_rational& b)
  {
    Uncertain<bool> u = a.i == b.i;
    if(is_indeterminate(u)){
      return a.rat == b.rat;
    }
    return make_certain(u);
  }



bool operator<(const Filtered_rational& a,
                const Filtered_rational& b)
  {
    Uncertain<bool> u = a.i < b.i;
    if(is_indeterminate(u)){
      return a.rat < b.rat;
    }
    return make_certain(u);
  }

  std::ostream& operator<<(std::ostream& os, const Filtered_rational& )
  {
    return os;
  }

  std::istream& operator>>(std::istream& is, const Filtered_rational& )
  {
    return is;
  }

  namespace internal {

      inline bool fit_in_double(const Filtered_rational& fr, double& d)
      {
          return fit_in_double(fr.i, d);
      }
  }
CGAL_DEFINE_COERCION_TRAITS_FOR_SELF(Filtered_rational)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int      ,Filtered_rational)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double   ,Filtered_rational)


} // namespace CGAL

#endif
