#include <cassert>

namespace CGAL {
template <class NT1, class NT2>
struct checked_NT {
  NT1 x1;
  NT2 x2;
  void verify()const{assert(NT2(x1)==x2);}
  struct pieces{};
  checked_NT(pieces,NT1 y1,NT2 y2):x1(y1),x2(y2){verify();}
  checked_NT():x1(),x2(){verify();}
  checked_NT(checked_NT const&t):x1(t.x1),x2(t.x2){verify();}
  checked_NT& operator=(checked_NT const&t){x1=t.x1;x2=t.x2;verify();return *this;}
  template<class T> checked_NT(T const&t,typename boost::enable_if<is_implicit_convertible<T,NT1>,int>::type=0):x1(t),x2(t){verify();}
  template<class T> explicit checked_NT(T const&t,typename boost::disable_if<is_implicit_convertible<T,NT1>,int>::type=0):x1(t),x2(t){verify();}
  /*TODO: enable_if to restrict the types*/
  template<class T> checked_NT& operator=(T const&t){x1=t;x2=t;verify();return *this;}
  checked_NT operator-()const{return checked_NT(pieces(),-x1,-x2);}
  checked_NT& operator+=(checked_NT const&a){x1+=a.x1;x2+=a.x2;verify();return *this;}
  checked_NT& operator-=(checked_NT const&a){x1-=a.x1;x2-=a.x2;verify();return *this;}
  checked_NT& operator*=(checked_NT const&a){x1*=a.x1;x2*=a.x2;verify();return *this;}
  checked_NT& operator/=(checked_NT const&a){x1/=a.x1;x2/=a.x2;verify();return *this;}
  checked_NT& operator%=(checked_NT const&a){x1%=a.x1;x2%=a.x2;verify();return *this;}
  friend checked_NT operator+(checked_NT const&a,checked_NT const&b){return checked_NT(pieces(),a.x1+b.x1,a.x2+b.x2);}
  friend checked_NT operator-(checked_NT const&a,checked_NT const&b){return checked_NT(pieces(),a.x1-b.x1,a.x2-b.x2);}
  friend checked_NT operator*(checked_NT const&a,checked_NT const&b){return checked_NT(pieces(),a.x1*b.x1,a.x2*b.x2);}
  friend checked_NT operator/(checked_NT const&a,checked_NT const&b){return checked_NT(pieces(),a.x1/b.x1,a.x2/b.x2);}
  friend checked_NT operator%(checked_NT const&a,checked_NT const&b){return checked_NT(pieces(),a.x1%b.x1,a.x2%b.x2);}
  double to_double()const{double a=CGAL::to_double(x1);double b=CGAL::to_double(x2);assert(a==b)/*too strong*/;return a;}
  std::pair<double,double> to_interval()const{std::pair<double,double> a=CGAL::to_interval(x1);std::pair<double,double> b=CGAL::to_interval(x2);assert(a.first<=b.second&&a.second>=b.first)/*overlap*/;return a;}
  friend bool operator<(checked_NT const&a,checked_NT const&b){bool res=a.x1<b.x1;bool other=a.x2<b.x2;assert(res==other);return res;}
  friend bool operator>(checked_NT const&a,checked_NT const&b){bool res=a.x1>b.x1;bool other=a.x2>b.x2;assert(res==other);return res;}
  friend bool operator==(checked_NT const&a,checked_NT const&b){bool res=a.x1==b.x1;bool other=a.x2==b.x2;assert(res==other);return res;}
  friend bool operator!=(checked_NT const&a,checked_NT const&b){bool res=a.x1!=b.x1;bool other=a.x2!=b.x2;assert(res==other);return res;}
  friend bool operator<=(checked_NT const&a,checked_NT const&b){bool res=a.x1<=b.x1;bool other=a.x2<=b.x2;assert(res==other);return res;}
  friend bool operator>=(checked_NT const&a,checked_NT const&b){bool res=a.x1>=b.x1;bool other=a.x2>=b.x2;assert(res==other);return res;}
  friend std::ostream& operator<<(std::ostream&o,const checked_NT&a){return o<<a.x1;}
  /*TODO: check input (output may be different).*/
  friend std::istream& operator>>(std::istream&i,checked_NT&a){i>>a.x1;a.x2=NT2(a.x1);return i;}
  /*TODO: conversion operators*/
};

template<class NT1,class NT2> struct Algebraic_structure_traits<checked_NT<NT1,NT2> >
: Algebraic_structure_traits_base<checked_NT<NT1,NT2>,typename Algebraic_structure_traits<NT1>::Algebraic_category>
{
  typedef Algebraic_structure_traits<NT1> AST1;
  typedef typename AST1::Is_exact Is_exact;
  typedef typename AST1::Is_numerical_sensitive Is_numerical_sensitive;
  typedef checked_NT<NT1,NT2> Type;
  struct Is_zero
    : public std::unary_function< Type, bool > {
      bool operator()( const Type& x ) const {
	bool a=CGAL::is_zero(x.x1);
	bool b=CGAL::is_zero(x.x2);
	assert(a==b);
	return a;
      }
    };
  struct Is_one
    : public std::unary_function< Type, bool > {
      bool operator()( const Type& x ) const {
	bool a=CGAL::is_one(x.x1);
	bool b=CGAL::is_one(x.x2);
	assert(a==b);
	return a;
      }
    };
  struct Is_square
    : public std::unary_function< Type, bool > {
      bool operator()( const Type& x ) const {
	bool a=CGAL::is_square(x.x1);
	bool b=CGAL::is_square(x.x2);
	assert(a==b);
	return a;
      }
    };
  struct Square
    : public std::unary_function< Type, Type > {
      Type operator()( const Type& x ) const {
	return Type(typename Type::pieces(),CGAL::square(x.x1),CGAL::square(x.x2));
      }
    };
  struct Div
    : public std::binary_function< Type, Type, Type > {
      Type operator()( const Type& x, const Type& y ) const {
	return Type(typename Type::pieces(),
	    CGAL::div(x.x1,y.x1),
	    CGAL::div(x.x2,y.x2));
      }
    };
  struct Mod
    : public std::binary_function< Type, Type, Type > {
      Type operator()( const Type& x, const Type& y ) const {
	return Type(typename Type::pieces(),
	    CGAL::mod(x.x1,y.x1),
	    CGAL::mod(x.x2,y.x2));
      }
    };
  struct Integral_division
    : public std::binary_function< Type, Type, Type > {
      Type operator()( const Type& x, const Type& y ) const {
	return Type(typename Type::pieces(),
	    CGAL::integral_division(x.x1,y.x1),
	    CGAL::integral_division(x.x2,y.x2));
      }
    };
};

template<class NT1,class NT2>struct Real_embeddable_traits<checked_NT<NT1,NT2> >
: INTERN_RET::Real_embeddable_traits_base<checked_NT<NT1,NT2>, typename Real_embeddable_traits<NT1>::Is_real_embeddable> {
  typedef checked_NT<NT1,NT2> Type;
  struct Sgn
    : public std::unary_function< Type, ::CGAL::Sign > {
      ::CGAL::Sign operator()( const Type& x ) const {
	::CGAL::Sign a=CGAL::sign(x.x1);
	::CGAL::Sign b=CGAL::sign(x.x2);
	assert(a==b);
	return a;
      }
    };
  struct To_double
    : public std::unary_function< Type, double > {
      double operator()( const Type& x ) const {
	return x.to_double();
      }
    };
  struct To_interval
    : public std::unary_function< Type, std::pair< double, double > > {
      std::pair<double, double> operator()( const Type& x ) const {
	return x.to_interval();
      }
    };
  struct Compare
    : public std::binary_function< Type, Type, Comparison_result > {
      Comparison_result operator()(
	  const Type& x,
	  const Type& y ) const {
	Comparison_result a=CGAL::compare(x.x1,y.x1);
	Comparison_result b=CGAL::compare(x.x2,y.x2);
	assert(a==b);
	return a;
      }
    };

};
}
