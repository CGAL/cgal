
namespace CGAL {

/*!
\ingroup PkgSTLExtensionRef

Some functions can return different types of objects. A typical
\cpp solution to this problem is to derive all possible return
types from a common base class, to return a pointer to this
class and to perform a dynamic cast on this pointer. The class
`Object` provides an abstraction.
An object `obj` of the class `Object` can
represent an arbitrary class. The only operations it provides is
to make copies and assignments, so that you can put them in lists
or arrays. Note that `Object` is NOT a common base class for the
elementary classes. Therefore, there is no
automatic conversion from these classes to `Object`. Rather
this is done with the global function `make_object()`. This
encapsulation mechanism requires the use of `assign` or
`object_cast` to use the functionality of the encapsulated class.

This class is similar in spirit to `std::any`.

\cgalHeading{Example}

In the following example, the object class is used as return value for the
intersection computation, as there are possibly different return values.

\code{.cpp}
{
typedef Cartesian<double> K;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;

Point_2 point;
Segment_2 segment, segment_1, segment_2;

std::cin >> segment_1 >> segment_2;

Object obj = intersection(segment_1, segment_2);

if (assign(point, obj)) {
  //do something with point
} else if (assign(segment, obj)) {
  // do something with segment
}
  //there was no intersection
}
\endcode

A more efficient way to access the object is to use `object_cast`,
which allows to skip a default construction and assignment:

\code{.cpp}
{
typedef Cartesian<double> K;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;

Segment_2 segment_1, segment_2;

std::cin >> segment_1 >> segment_2;

Object obj = intersection(segment_1, segment_2);

if (const Point_2 * point = object_cast<Point_2>(&obj)) {
  // do something with *point
} else if (const Segment_2 * segment = object_cast<Segment_2>(&obj)) {
  // do something with *segment
}
  // there was no intersection
}
\endcode

The intersection routine itself looks roughly as follows:

\code{.cpp}
template < class Kernel >
Object intersection(Segment_2<Kernel> s1, Segment_2<Kernel> s2)
{
  if ( intersection_is_a_point ) {
    Point_2<Kernel> p = ... ;
    return make_object(p);
  } else if ( intersection_is_a_segment ) {
    Segment_2<Kernel> s = ... ;
    return make_object(s);
  }

  // empty intersection
  return Object();
}
\endcode

*/
class Object {
public:

/// \name Creation
/// Objects of type `Object` are normally created using the global function `make_object()`.
/// @{
/*!
introduces an empty object.
*/
Object();

/*!
Copy constructor.
*/
Object(const Object &o);

/*!
Implicit converting constructor for compatibility with
`std::variant`.
*/
Object(std::variant<T...>);


/*!
Implicit converting constructor for compatibility with
`std::optional` and `std::variant`.
 */
Object(std::optional< std::variant<T...> >);

/// @}


/// \name Operations
/// @{
/*!
Assignment.
*/
Object &operator=(const Object &o);

/// @}


/// \name Operations
/// @{
/*!
returns true, if `obj` does not
contain an object of type `T`.
*/
bool empty();

/// @}

/// \name Operations
/// @{
/*!
returns true, iff `obj` contains an object of type `T`.
*/
template<class T> bool is();



/// @}


/// \name Operations
/// @{
/*!
returns the type information of the contained type,
or `typeid(void)` if empty.
*/
const std::type_info & type() const;

/// @}

}; /* end Object */



/*!
Creates an object that contains `t`.

\relates Object
*/
template <class T> Object make_object(const T &t);

/*!

assigns `o` to `c` if `o` was constructed from an object of type `T`.
Returns `true`, if the assignment was possible.  For efficiency
reasons, we recommend using `object_cast` instead.

\relates Object
*/

template <class T> bool assign(T& c, const Object& o);

/*!
Returns a pointer to the object of type `T` stored by `o`,
if any, otherwise returns `NULL`.
\relates Object
*/
template <class T> const T * object_cast(const Object * o);

/*!
Returns a copy of the object of type `T` stored by `o`,
if any, otherwise throws an exception of type `Bad_object_cast`.
\relates Object
*/
template <class T> T object_cast(const Object & o);


} /* end namespace CGAL */
