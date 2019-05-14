
namespace CGAL {

/*!
\ingroup PkgStlExtensionUtilities

An object of the class `Uncertain` represents an uncertainty on the 
value of type `T`. This uncertainty is represented by a non-empty range of 
values of type `T`. 

The idea is that sometimes you are not sure of the result of a function, and 
you would like to communicate that to the caller. `Uncertain<T>` allows 
just that. 

`Uncertain<T>` is also meant to be used as a drop-in replacement for 
`T` in some template contexts, as much as possible. This is why it 
provides overloaded operators and functions to naturally extend the Boolean 
operations for `Uncertain<bool>` for example, or the operations on 
enumeration types. 

`Uncertain<T>` is used in \cgal as the return type of geometric predicates 
when the number type used is interval arithmetic like `Interval_nt`. End 
users typically do not see it, as it is hidden in the implementation of the 
filtered predicates provided by the various filtered kernels, but it is 
important that providers of predicates that are meant to be filtered by 
`Filtered_predicate`, know about it. 

Note concerning \cgal assertions: assertions checking an expression of type 
`Uncertain<bool>` will trigger an assertion failure only if the assertion 
is certainly `false`. In case of an indeterminate value, the assertion is not 
triggered. This means that we assume, in case of doubt, that there is no 
error. 

It can also be used in other contexts as well, as it is a general tool. 
This can be seen as support for non-deterministic programming. 
Finally, note that this class has some common points with `boost::tribool`. 

\cgalHeading{Parameters}

The parameter `T` can either be `bool` or one of the three-valued 
(-1, 0, 1) enumeration types: `Sign`, `Comparison_result`, 
`Orientation`, `Oriented_side`, `Bounded_side` or `Angle`. 

Some functions are defined only when `T` is `bool` or alternatively 
when it is one of the enumeration types listed previously. 

\sa `CGAL::Interval_nt<bool>` 
*/
template< typename T >
class Uncertain {
public:

/// \name Types 
/// @{ 
/*!
The type `T`. 
*/ 
typedef unspecified_type value_type; 
/// @} 


/// \name Types 
/// @{ 
/*!
The type of the exception 
thrown for uncertain conversions. It is a typedef to the type 
`CGAL::Uncertain_conversion_exception` which derives from `std::range_error`. 
*/ 
typedef unspecified_type Uncertain_conversion_exception; 
/// @} 

/// \name Creation 
/// @{ 
/*!
introduces a certain object with value `T()`. 
*/ 
Uncertain(); 



/// @} 


/// \name Creation 
/// @{ 
/*!
introduces a certain object with value `t`. 
*/ 
Uncertain(T t); 



/// @} 


/// \name Creation 
/// @{ 
/*!
assigns the certain value `t` to `u`. 
*/ 
Uncertain& operator=(T t); 



/// @} 


/// \name Creation 
/// @{ 
/*!
introduces an object representing the range with lower bound `i` and 
upper bound `s`. \pre \f$ i<= s\f$. 
*/ 
Uncertain(T i, T s); 



/// @} 


/// \name Access Functions 
/// The following functions are meant to be used very rarely, they provide ways to inspect 
/// the content of an `Uncertain<T>` object. 
/// @{ 
/*!
returns the lower bound of the range represented by `u`. 
*/ 
T inf() const; 

/*!
returns the upper bound of the range represented by `u`. 
*/ 
T sup() const; 

/*!
returns true whether `u` and `u` are the same range (equality as sets). 
*/ 
bool is_same(Uncertain u) const; 

/// @} 


/// \name Uncertainty Testing and Conversion 
/// There are several ways to extract the content of an `Uncertain` object. The simplest way is to 
/// rely on the implicit conversion from `Uncertain<T>` to `T`. In this case, no special 
/// code has to be written, apart from an exception handler (anywhere higher in the call stack) to manage 
/// the uncertain case. The more involved one is more efficient, but requires manual treatment 
/// of the uncertain case, such as: 
/// 
/// \code{.cpp} 
/// 
/// Uncertain<bool> b = ...; 
/// if (is_certain(b)) 
/// bool cert_b = get_certain(b); // Extract the certain bool it contains 
/// ... 
/// else 
/// ... // b is indeterminate 
/// 
/// \endcode 
/// 
/// Another option is : 
/// 
/// \code{.cpp} 
/// 
/// Uncertain<bool> b = ...; 
/// if (certainly(b)) 
/// ... // b is certainly true 
/// else if (certainly_not(b)) 
/// ... // b is certainly false 
/// else 
/// ... // b is indeterminate 
/// 
/// \endcode 
/// 
/// There are many other handy functions which can be used for easier usage depending 
/// on the context. They are listed in the sequel. 
/// @{ 

/*!
returns `true` iff the value is certain, that is, it is unique, the range 
is a singleton, that is `u.inf() == u.sup()`. 
*/ 
bool is_certain() const; 

/*!
if `u.is_certain()`, then returns the certain value which is represented. 
Otherwise, throws an exception of type `Uncertain_conversion_exception`. 
A profile counter of the number of such exceptions thrown during the execution of 
the program is available with `CGAL_PROFILE`. 
*/ 
T make_certain() const; 

/*!
conversion operator to `T`. It does and returns the same thing as 
`u`. `make_certain()`. Note that relying on the automatic conversion 
can throw exceptions, which defeats the purpose of propagating uncertainty. 
Nevertheless, in many cases, it is hard to avoid it, for example for the 
`&&` and \f$ ||\f$ operators for `bool` (see below). 
*/ 
operator T() const; 

///@}

/*!
returns an indeterminate range. 
*/ 
static Uncertain<T> indeterminate(); 

/*! \name Overloaded Operators
The overloaded operators and functions are defined as preserving the set-inclusion property. 
Similarly to interval arithmetic, the returned range is guaranteed to contain 
the result of the operation over all values of the input range(s). 
In the following documentation we express this as the extension of the corresponding function 
over the type `T`. 
*/
/// @{

/*!
returns the extension of the equality operator over `u` and `v`. 
*/ 
template <class T> 
Uncertain<bool> operator==(Uncertain<T> u, Uncertain<T> v); 

/*!
returns `u == make_uncertain(v)`. 
*/ 
template <class T> 
Uncertain<bool> operator==(Uncertain<T> u, T v); 

/*!
returns `v == u`. 
*/ 
template <class T> 
Uncertain<bool> operator==(T u, Uncertain<T> v); 

/*!
returns the extension of the inequality operator over `u` and `v`. 
*/ 
template <class T> 
Uncertain<bool> operator!=(Uncertain<T> u, Uncertain<T> v); 

/*!
returns `u != make_uncertain(v)`. 
*/ 
template <class T> 
Uncertain<bool> operator!=(Uncertain<T> u, T v); 

/*!
returns `v != u`. 
*/ 
template <class T> 
Uncertain<bool> operator!=(T u, Uncertain<T> v); 

/*!
\name Overloaded Operators for Uncertain<bool>

The overloaded operators and functions are defined as preserving the set-inclusion property. 
Similarly to interval arithmetic, the returned range is guaranteed to contain 
the result of the operation over all values of the input range(s). 
In the following documentation we express this as the extension of the corresponding function 
over the type `T`. 

\note The logical operators \f$ \&\&\f$ and \f$ ||\f$ are not overloaded on purpose. The reason 
is that, when `f() && g()` is evaluated and they return `bool`, then `g()` 
is only evaluated when `f()` returns `true`. One could have a dependency so 
that `g()` has an internal precondition that required that `f()` had returned `true`. 
The overloaded operators for user-defined types can not provide this short-circuiting 
property, and so, if the overloaded operators where provided, then `g()` would 
be evaluated, no matter the result of `f()`, which could lead to an unwanted 
situation, or a performance loss. The \f$ \&\f$ and \f$ |\f$ operators do not have 
this short-circuiting property, and are therefore overloaded safely. 

When translating normal code to use and propagate uncertainty, such as : 

\code{.cpp} 

// Logical AND 
if ( (p.x() == 0) && (p.y() == 0) ) 
... 
else 
... 

// Logical OR 
if ( (q.x() == 0) || (q.y() == 0) ) 
... 
else 
... 

\endcode 

One can do, for example : 

\code{.cpp} 

// Logical AND 
Uncertain<bool> tmp = (p.x() == 0); 
Uncertain<bool> res = certainly_not(tmp) ? make_uncertain(false) : tmp & (p.y() == 0); 

... // Use res 

// Logical OR 
Uncertain<bool> tmp = (q.x() == 0); 
Uncertain<bool> res = certainly(tmp) ? make_uncertain(true) : tmp | (q.y() == 0); 

... // Use res 

\endcode 

This ensures that the first expression is not evaluated twice, and that the second is 
evaluated only if needed. 

This behavior can also be emulated through the use of macros, but only using 
non-standard features ("statement expressions", such as provided by GCC). The 
macros `CGAL_AND` and `CGAL_OR` are provided that perform the lazy 
evaluation of these logical operations. On compilers that do not support 
statement expressions, the macros simply expand to the \f$ \&\&\f$ and \f$ ||\f$ 
operators (which will throw an exception instead of propagating the uncertainty). 

\code{.cpp} 

// Logical AND 
Uncertain<bool> res = CGAL_AND( p.x() == 0 , p.y() == 0 ); 

... // Use res 

// Logical OR 
Uncertain<bool> res = CGAL_OR( q.x() == 0 , q.y() == 0 ); 

... // Use res 

\endcode 

*/
/// @{
/*!

returns the range containing the negated values of `u`. 
*/ 
Uncertain<bool> operator!(Uncertain<bool> u); 

/*!

returns the range containing the values computed as logical or from `u` and `v`. 
*/ 
Uncertain<bool> operator|(Uncertain<bool> u, Uncertain<bool> v); 

/*!

returns `u | make_uncertain(v)`. 
*/ 
Uncertain<bool> operator|(Uncertain<bool> u, bool v); 

/*!

returns `v | u`. 
*/ 
Uncertain<bool> operator|(bool u, Uncertain<bool> v); 

/*!

returns the range containing the values computed as logical and from `u` and `v`. 
*/ 
Uncertain<bool> operator&(Uncertain<bool> u, Uncertain<bool> v); 

/*!

returns `u & make_uncertain(v)`. 
*/ 
Uncertain<bool> operator&(Uncertain<bool> u, bool v); 

/*!

returns `v & u`. 
*/ 
Uncertain<bool> operator&(bool u, Uncertain<bool> v); 


/// @}

/// \name Overloaded Operators and Functions for Uncertain<enum T> Only 
/// @{

/*!

returns the extension of the less-than operator over `u` and `v`. 
*/ 
template <class T> 
Uncertain<bool> operator<(Uncertain<T> u, Uncertain<T> v); 

/*!

returns `u < make_uncertain(v)`. 
*/ 
template <class T> 
Uncertain<bool> operator<(Uncertain<T> u, T v); 

/*!

returns `make_uncertain(u) < v`. 
*/ 
template <class T> 
Uncertain<bool> operator<(T u, Uncertain<T> v); 

/*!

returns the extension of the greater-than operator over `u` and `v`. 
*/ 
template <class T> 
Uncertain<bool> operator>(Uncertain<T> u, Uncertain<T> v); 

/*!

returns `u > make_uncertain(v)`. 
*/ 
template <class T> 
Uncertain<bool> operator>(Uncertain<T> u, T v); 

/*!

returns `make_uncertain(u) > v`. 
*/ 
template <class T> 
Uncertain<bool> operator>(T u, Uncertain<T> v); 

/*!

returns the extension of the less-than or equal operator over `u` and `v`. 
*/ 
template <class T> 
Uncertain<bool> operator<=(Uncertain<T> u, Uncertain<T> v); 

/*!

returns `u <= make_uncertain(v)`. 
*/ 
template <class T> 
Uncertain<bool> operator<=(Uncertain<T> u, T v); 

/*!

returns `make_uncertain(u) <= v`. 
*/ 
template <class T> 
Uncertain<bool> operator<=(T u, Uncertain<T> v); 

/*!

returns the extension of the greater-than or equal operator over `u` and `v`. 
*/ 
template <class T> 
Uncertain<bool> operator>=(Uncertain<T> u, Uncertain<T> v); 

/*!

returns `u > make_uncertain(v)`. 
*/ 
template <class T> 
Uncertain<bool> operator>=(Uncertain<T> u, T v); 

/*!

returns `make_uncertain(u) >= v`. 
*/ 
template <class T> 
Uncertain<bool> operator>=(T u, Uncertain<T> v); 

/*!

returns the extension of the multiplication operator over `u` and `v`. 
This requires `T` to have a multiplication operator as well. 
*/ 
template <class T> 
Uncertain<T> operator*(Uncertain<T> u, Uncertain<T> v); 

/*!

returns `u * make_uncertain(v)`. 
*/ 
template <class T> 
Uncertain<T> operator*(Uncertain<T> u, T v); 

/*!

returns `make_uncertain(u) * v`. 
*/ 
template <class T> 
Uncertain<T> operator<(T u, Uncertain<T> v); 

/*!

returns the extension of the unary minus operator over `u`. 
*/ 
template <class T> 
Uncertain<T> operator-(Uncertain<T> u); 

/*!

returns the extension of the `enum_cast<T>` function over `u`. 
*/ 
template <class T, class U> 
Uncertain<T> enum_cast(Uncertain<U> u); 

/// @}

}; /* end Uncertain */

/*!
returns `u.inf()`. 
\relates Uncertain 
*/ 
template <class T> T inf(Uncertain<T> u); 

/*!
returns `u.sup()`. 
\relates Uncertain 
*/ 
template <class T> T sup(Uncertain<T> u); 

/*!
returns `true`. 
\relates Uncertain 
*/ 
template <class T> bool is_certain(T t); 

/*!
returns `u.is_certain`(). 
\relates Uncertain 
*/ 
template <class T> bool is_certain(Uncertain<T> u); 

/*!
returns `U::indeterminate()` if `U` is `Uncertain<T>`, 
and `U()` otherwise. 

\relates Uncertain 
*/ 
template <class U> U indeterminate(); 

/*!
returns `false`. 
\relates Uncertain 
*/ 
template <class T> bool is_indeterminate(T u); 

/*!
returns `!is_certain(u)`. 
\relates Uncertain 
*/ 
template <class T> bool is_indeterminate(Uncertain<T> u); 

/*!
returns `t`. 
\relates Uncertain 
*/ 
template <class T> T get_certain(T t); 

/*!
returns `u.make_certain`(). \pre `u.is_certain`(). 
\relates Uncertain 
*/ 
template <class T> T get_certain(Uncertain<T> u); 

/*!
returns `t`. 
\relates Uncertain 
*/ 
template <class T> T make_certain(T t); 

/*!
returns `u.make_certain`(). 
\relates Uncertain 
*/ 
template <class T> T make_certain(Uncertain<T> u); 

/*!
returns `Uncertain<T>(u)`. 
\relates Uncertain 
*/ 
template <class T> Uncertain<T> make_uncertain(T t); 

/*!
returns `u`. 
\relates Uncertain 
*/ 
template <class T> Uncertain<T> make_uncertain(Uncertain<T> u); 

/*!
  Boolean operation with 3 arguments.
  \relates CGAL::Uncertain
 */
#define CGAL_AND_3

/*!
  Boolean operation with 3 arguments.
  \relates CGAL::Uncertain
 */
#define CGAL_OR_3

/*!
returns `true` iff `u.is_certain()`, and the `u.make_certain`() 
returns `true`. 
\relates Uncertain 
*/ 
bool certainly(Uncertain<bool> u); 

/*!
returns `u`. 
\relates Uncertain 
*/ 
bool certainly(bool u); 

/*!
returns `true` iff `u.is_certain()` returns `false`, or if 
`u.make_certain`() returns `true`. 
\relates Uncertain 
*/ 
bool possibly(Uncertain<bool> u); 

/*!
returns `u`. 
\relates Uncertain 
*/ 
bool possibly(bool u); 

/*!
returns `true` iff `u.is_certain()`, and the `u.make_certain`() 
returns `false`. 
\relates Uncertain 
*/ 
bool certainly_not(Uncertain<bool> u); 

/*!
returns `!u`. 
\relates Uncertain 
*/ 
bool certainly_not(bool u); 

/*!
returns `true` iff `u.is_certain()` returns `false`, or if 
`u.make_certain`() returns `false`. 
\relates Uncertain 
*/ 
bool possibly_not(Uncertain<bool> u); 

/*!
returns `!u`. 
\relates Uncertain 
*/ 
bool possibly_not(bool u); 

} /* end namespace CGAL */
