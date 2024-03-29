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


  Filtered_rational(const Interval_nt<>& i, const Exact_rational& rat)
      : i(i), rat(rat)
  {}

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

    class Sgn
      : public CGAL::cpp98::unary_function< Type, ::CGAL::Sign > {
      public:
        ::CGAL::Sign operator()( const Type& x ) const {
          return sign(x);
        }
    };

    class To_double
      : public CGAL::cpp98::unary_function< Type, double > {
      public:
        double operator()( const Type& x ) const {
          return to_double(x.i);
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

  Sign sign(const Filtered_rational& a)
  {
    Uncertain<Sign> u = sign(a.i);
     if(is_indeterminate(u)){
       return CGAL::sign(a.rat);
    }
    return make_certain(u);
  }

  Sign compare(const Filtered_rational& a, const Filtered_rational& b){
    Uncertain<Sign> u = compare(a.i, b.i);
     if(is_indeterminate(u)){
       return compare(a.rat, b.rat);
    }
    return make_certain(u);
  }

  Filtered_rational square(const Filtered_rational& a)
  {
    return Filtered_rational(square(a.i), square(a.rat));
  }

  std::pair<double,double> to_interval(const Filtered_rational& b)
{
  return std::make_pair(b.i.inf(), b.i.sup());
}

  double to_double(const Filtered_rational& b)
  {
    return to_double(b.i);
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
