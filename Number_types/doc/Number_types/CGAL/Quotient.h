
namespace CGAL {

/*!
\ingroup nt_cgal

\anchor Quotient 

An object of the class `Quotient<NT>` is an element of the 
field of quotients of the integral domain type `NT`. 
If `NT` behaves like an integer, `Quotient<NT>` 
behaves like a rational number. 
\leda's class `rational` (see Section \ref ledant) 
has been the basis for `Quotient<NT>`. 
A `Quotient<NT>` `q` is represented as a pair of 
`NT`s, representing numerator and denominator. 

`NT` must be at least model of concept `IntegralDomainWithoutDivision`. 

`NT` must be a model of concept `RealEmbeddable`. 

\models ::Field 
\models ::RealEmbeddable 
\models ::Fraction 

Operations 
-------------- 

There are two access functions, namely to the 
numerator and the denominator of a quotient. 
Note that these values are not uniquely defined. 
It is guaranteed that `q.numerator()` and 
`q.denominator()` return values `nt_num` and 
`nt_den` such that `q = nt_num/nt_den`, only 
if `q.numerator()` and `q.denominator()` are called 
consecutively wrt `q`, i.e. `q` is not involved in 
any other operation between these calls. 

The stream operations are available as well. 
They assume that corresponding stream operators for type `NT` exist. 

The following functions are added to fulfill the \cgal requirements 
on number types. 

*/
template< typename NT >
class Quotient {
public:

/// \name Creation 
/// @{

/*! 
introduces an uninitialized variable `q`. 
*/ 
Quotient(); 

/*! 
introduces the quotient `t/1`. NT needs to have a constructor from T. 
*/ 
template <class T> Quotient<NT>(const T& t); 

/*! 
introduces the quotient `NT(t.numerator())/NT(t.denominator())`. 
NT needs to have a constructor from T. 
*/ 
template <class T> Quotient<NT>(const Quotient<T>& t); 

/*! 
introduces the quotient `n/d`. 

\pre \f$ d \neq0\f$. 
*/ 
Quotient(const NT& n, const NT& d); 

/// @} 

/// \name Operations 
/// @{

/*! 
returns a numerator of `q`. 
*/ 
NT numerator() const; 

/*! 
returns a denominator of `q`. 
*/ 
NT denominator() const; 

/// @}

}; /* end Quotient */

/*! 
writes `q` to ostream `out` in format "<TT>n/d</TT>", where 
<TT>n</TT>\f$ ==\f$`q.numerator()` and <TT>d</TT>\f$ ==\f$`q.denominator()`. 
\relates Quotient 
*/ 
std::ostream& operator<<(std::ostream& out, const Quotient<NT>& q); 

/*! 
reads `q` from istream `in`. Expected format is 
"<TT>n/d</TT>", where <TT>n</TT> and <TT>d</TT> are of type `NT`. 
A single <TT>n</TT> which is not followed by a <TT>/</TT> is also 
accepted and interpreted as <TT>n/1</TT>. 
\relates Quotient 
*/ 
std::istream& operator>>(std::istream& in, Quotient<NT>& q); 

/*! 
returns some double approximation to `q`. 
\relates Quotient 
*/ 
double to_double(const Quotient<NT>& q); 

/*! 
returns true, if numerator and denominator are valid. 
\relates Quotient 
*/ 
bool is_valid(const Quotient<NT>& q); 

/*! 
returns true, if numerator and denominator are finite. 
\relates Quotient 
*/ 
bool is_finite(const Quotient<NT>& q); 

/*! 
returns the square root of `q`. This is supported if and only if 
`NT` supports the square root as well. 
\relates Quotient 
*/ 
Quotient<NT> sqrt(const Quotient<NT>& q); 

} /* end namespace CGAL */
