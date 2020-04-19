/// \defgroup projectionobjects Projection Function Objects
/// \ingroup PkgSTLExtensionRef

namespace CGAL {

/*!
\ingroup projectionobjects



The class `Cast_function_object` applies a C-style type cast to
its argument.



\cgalModels `ProjectionObject`


*/
template< typename Arg, typename Result >
struct Cast_function_object {



/*!
*/
typedef Arg argument_type;




/*!
*/
typedef Result result_type;


/// \name Creation
/// @{
/*!
default constructor.
*/
Cast_function_object();



/// @}


/// \name Operations
/// @{
/*!
returns
`(Result)x`.
*/
result_type& operator()(argument_type& x) const;



/// @}


/// \name Operations
/// @{
/*!
returns `(Result)x`.
*/
const result_type& operator()(const argument_type&
x) const;



/// @}



}; /* end Cast_function_object */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup projectionobjects



The class `Compare_to_less` is used to convert a functor which
returns a `Comparison_result` to a predicate (returning bool) : it
will return true iff the return value of `F` is `SMALLER`.
The class is used in conjunction with the `compare_to_less`
function; see there for an explanation on how exactly the functors
are combined.



\sa `CGAL::compare_to_less()`


*/
template< typename F >
class Compare_to_less {
public:


/// \name Types
/// @{
/*!
type of the composed functor.
*/
typedef unspecified_type Type;
/// @}


}; /* end Compare_to_less */
} /* end namespace CGAL */
namespace CGAL {

/*!
\ingroup projectionobjects


Changes a functor
returning a `Comparison_result` to one which returns a bool.
The returned functor will return `true` iff the original one
returns `SMALLER`.



\sa `CGAL::Compare_to_less<F>`

returns a functor equivalent to
`f`, but which returns a bool instead of a
`Comparison_result`.
*/
template < class F > Compare_to_less< F >
compare_to_less(const F& f);

} /* namespace CGAL */


/// \defgroup STLCreators Creator Function Objects
/// \ingroup PkgSTLExtensionRef

namespace CGAL {

/*!
\ingroup STLCreators


The class `Creator_1` defines types and operations
for creating objects from one argument.

\tparam Arg must be convertible to `Result`.

*/
template< typename Arg, typename Result >
class Creator_1 {
public:


/// \name Requirements
/// @{
/*!
type of argument.
*/
typedef Arg argument_type;
/// @}


/// \name Requirements
/// @{
/*!
type of object to create.
*/
typedef Result result_type;
/// @}

/// \name Requirements
/// @{
/*!
returns
`result_type(a)`.
*/
result_type operator()(argument_type a) const;



/// @}



}; /* end Creator_1 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup STLCreators



The class `Creator_2` defines types and operations
for creating objects from two arguments.


\tparam Result must define a corresponding constructor.

*/
template< typename Arg1, typename Arg2, typename Result >
class Creator_2 {
public:


/// \name Requirements
/// @{
/*!
type of first argument.
*/
typedef Arg1 argument1_type;
/// @}


/// \name Requirements
/// @{
/*!
type of second argument.
*/
typedef Arg2 argument2_type;
/// @}


/// \name Requirements
/// @{
/*!
type of object to create.
*/
typedef Result result_type;
/// @}

/// \name Requirements
/// @{
/*!
returns `result_type(a1, a2)`.
*/
result_type operator()(argument_type1 a1, argument_type2
a2) const;



/// @}



}; /* end Creator_2 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup STLCreators



The class `Creator_3` defines types and operations
for creating objects from three arguments.



\tparam Result must define a corresponding constructor.


*/
template< typename Arg1, typename Arg2, typename Arg3, typename Result >
class Creator_3 {
public:


/// \name Requirements
/// @{
/*!
type of first argument.
*/
typedef Arg1 argument1_type;
/// @}


/// \name Requirements
/// @{
/*!
type of second argument.
*/
typedef Arg2 argument2_type;
/// @}


/// \name Requirements
/// @{
/*!
type of third argument.
*/
typedef Arg3 argument3_type;
/// @}


/// \name Requirements
/// @{
/*!
type of object to create.
*/
typedef Result result_type;
/// @}

/// \name Requirements
/// @{
/*!
returns `result_type(a1, a2,
a3)`.
*/
result_type operator()(argument_type1 a1, argument_type2
a2, argument_type3 a3) const;



/// @}



}; /* end Creator_3 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup STLCreators



The class `Creator_4` defines types and operations
for creating objects from four arguments.



\tparam Result must define a corresponding constructor.


*/
template< typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Result >
class Creator_4 {
public:


/// \name Requirements
/// @{
/*!
type of first argument.
*/
typedef Arg1 argument1_type;
/// @}


/// \name Requirements
/// @{
/*!
type of second argument.
*/
typedef Arg2 argument2_type;
/// @}


/// \name Requirements
/// @{
/*!
type of third argument.
*/
typedef Arg3 argument3_type;
/// @}


/// \name Requirements
/// @{
/*!
type of 4th argument.
*/
typedef Arg4 argument4_type;
/// @}


/// \name Requirements
/// @{
/*!
type of object to create.
*/
typedef Result result_type;
/// @}

/// \name Requirements
/// @{
/*!
returns
`result_type(a1, a2, a3, a4)`.
*/
result_type operator()(argument_type1 a1, argument_type2
a2, argument_type3 a3, argument_type4 a4) const;



/// @}



}; /* end Creator_4 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup STLCreators



The class `Creator_5` defines types and operations
for creating objects from five arguments.



\tparam Result must define a corresponding constructor.


*/
template< typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename Result >
class Creator_5 {
public:


/// \name Requirements
/// @{
/*!
type of first argument.
*/
typedef Arg1 argument1_type;
/// @}


/// \name Requirements
/// @{
/*!
type of second argument.
*/
typedef Arg2 argument2_type;
/// @}


/// \name Requirements
/// @{
/*!
type of third argument.
*/
typedef Arg3 argument3_type;
/// @}


/// \name Requirements
/// @{
/*!
type of 4th argument.
*/
typedef Arg4 argument4_type;
/// @}


/// \name Requirements
/// @{
/*!
type of 5th argument.
*/
typedef Arg5 argument5_type;
/// @}


/// \name Requirements
/// @{
/*!
type of object to create.
*/
typedef Result result_type;
/// @}

/// \name Requirements
/// @{
/*!
returns `result_type(a1, a2, a3, a4, a5)`.
*/
result_type operator()(argument_type1 a1, argument_type2
a2, argument_type3 a3, argument_type4 a4, argument_type5 a5)
const;



/// @}



}; /* end Creator_5 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup STLCreators



The class `Creator_uniform_2` defines types and operations
for creating objects from two arguments of the same type.



\tparam Result must define a constructor from two `Arg` arguments.


*/
template< typename Arg, typename Result >
class Creator_uniform_2 {
public:


/// \name Requirements
/// @{
/*!
type of arguments.
*/
typedef Arge argument_type;
/// @}


/// \name Requirements
/// @{
/*!
type of object to create.
*/
typedef Result result_type;
/// @}

/// \name Requirements
/// @{
/*!
returns `result_type(a1, a2)`.
*/
result_type operator()(argument_type a1, argument_type a2)
const;



/// @}



}; /* end Creator_uniform_2 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup STLCreators



The class `Creator_uniform_3` defines types and operations
for creating objects from three arguments of the same type.



\tparam Result must define a constructor from three `Arg` arguments.


*/
template< typename Arg, typename Result >
class Creator_uniform_3 {
public:


/// \name Requirements
/// @{
/*!
type of arguments.
*/
typedef Arg argument_type;
/// @}


/// \name Requirements
/// @{
/*!
type of object to create.
*/
typedef Result result_type;
/// @}

/// \name Requirements
/// @{
/*!
returns `result_type(a1, a2, a3)`.
*/
result_type operator()(argument_type a1, argument_type a2,
argument_type a3) const;



/// @}



}; /* end Creator_uniform_3 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup STLCreators



The class `Creator_uniform_4` defines types and operations
for creating objects from four arguments of the same type.



\tparam Result must define a constructor from four `Arg` arguments.


*/
template< typename Arg, typename Result >
class Creator_uniform_4 {
public:


/// \name Requirements
/// @{
/*!
type of arguments.
*/
typedef Arg argument_type;
/// @}


/// \name Requirements
/// @{
/*!
type of object to create.
*/
typedef Result result_type;
/// @}

/// \name Requirements
/// @{
/*!
returns
`result_type(a1, a2, a3, a4)`.
*/
result_type operator()(argument_type a1, argument_type a2,
argument_type a3, argument_type a4) const;



/// @}



}; /* end Creator_uniform_4 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup STLCreators



The class `Creator_uniform_5` defines types and operations
for creating objects from five arguments of the same type.



\tparam Result must define a constructor from five `Arg` arguments.


*/
template< typename Arg, typename Result >
class Creator_uniform_5 {
public:


/// \name Requirements
/// @{
/*!
type of arguments.
*/
typedef Arg argument_type;
/// @}


/// \name Requirements
/// @{
/*!
type of object to create.
*/
typedef Result result_type;
/// @}

/// \name Requirements
/// @{
/*!
returns `result_type(a1, a2, a3, a4, a5)`.
*/
result_type operator()(argument_type a1, argument_type a2,
argument_type a3, argument_type a4, argument_type a5)
const;



/// @}



}; /* end Creator_uniform_5 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup STLCreators



The class `Creator_uniform_6` defines types and operations
for creating objects from six arguments of the same type.



\tparam Result must define a constructor from six `Arg` arguments.


*/
template< typename Arg, typename Result >
class Creator_uniform_6 {
public:


/// \name Requirements
/// @{
/*!
type of arguments.
*/
typedef Arg argument_type;
/// @}


/// \name Requirements
/// @{
/*!
type of object to create.
*/
typedef Result result_type;
/// @}

/// \name Requirements
/// @{
/*!
returns `result_type(a1, a2, a3, a4,
a5, a6)`.
*/
result_type operator()(argument_type a1, argument_type a2,
argument_type a3, argument_type a4, argument_type a5,
argument_type a6) const;



/// @}



}; /* end Creator_uniform_6 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup STLCreators



The class `Creator_uniform_7` defines types and operations
for creating objects from seven arguments of the same type.



\tparam Result must define a constructor from seven `Arg` arguments.


*/
template< typename Arg, typename Result >
class Creator_uniform_7 {
public:


/// \name Requirements
/// @{
/*!
type of arguments.
*/
typedef Arg argument_type;
/// @}


/// \name Requirements
/// @{
/*!
type of object to create.
*/
typedef Result result_type;
/// @}

/// \name Requirements
/// @{
/*!
returns
`result_type(a1, a2, a3, a4, a5, a6, a7)`.
*/
result_type operator()(argument_type a1, argument_type a2,
argument_type a3, argument_type a4, argument_type a5,
argument_type a6, argument_type a7) const;



/// @}



}; /* end Creator_uniform_7 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup STLCreators



The class `Creator_uniform_8` defines types and operations
for creating objects from eight arguments of the same type.



\tparam Result must define a constructor from eight `Arg` arguments.


*/
template< typename Arg, typename Result >
class Creator_uniform_8 {
public:


/// \name Requirements
/// @{
/*!
type of arguments.
*/
typedef Arg argument_type;
/// @}


/// \name Requirements
/// @{
/*!
type of object to create.
*/
typedef Result result_type;
/// @}

/// \name Requirements
/// @{
/*!
returns `result_type(a1, a2, a3, a4, a5, a6, a7,
a8)`.
*/
result_type operator()(argument_type a1, argument_type a2,
argument_type a3, argument_type a4, argument_type a5,
argument_type a6, argument_type a7, argument_type a8)
const;



/// @}



}; /* end Creator_uniform_8 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup STLCreators



The class `Creator_uniform_9` defines types and operations
for creating objects from nine arguments of the same type.



\tparam Result must define a constructor from nine `Arg` arguments.


*/
template< typename Arg, typename Result >
class Creator_uniform_9 {
public:


/// \name Requirements
/// @{
/*!
type of arguments.
*/
typedef Arg argument_type;
/// @}


/// \name Requirements
/// @{
/*!
type of object to create.
*/
typedef Result result_type;
/// @}

/// \name Requirements
/// @{
/*!
returns `result_type(a1, a2, a3, a4,
a5, a6, a7, a8, a9)`.
*/
result_type operator()(argument_type a1, argument_type a2,
argument_type a3, argument_type a4, argument_type a5,
argument_type a6, argument_type a7, argument_type a8,
argument_type a9) const;



/// @}



}; /* end Creator_uniform_9 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup STLCreators



The class `Creator_uniform_d` defines types and operations
for creating objects from two arguments of the same type.



\tparam Result must define a constructor from three arguments: one `d` dimension and two `Arg` arguments.


*/
template< typename Arg, typename Result >
class Creator_uniform_d {
public:


/// \name Requirements
/// @{
/*!
type of arguments.
*/
typedef Arg argument_type;
/// @}


/// \name Requirements
/// @{
/*!
type of object to create.
*/
typedef Result result_type;
/// @}

/// \name Requirements
/// @{
/*!
returns `result_type(d, a1, a2)`.
*/
result_type operator()(argument_type a1, argument_type a2)
const;



/// @}



}; /* end Creator_uniform_d */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup projectionobjects



The class `Dereference` dereferences a pointer
(`operator*`).



\cgalModels `ProjectionObject`

*/
template< typename Value >
struct Dereference {


/*!
*/
typedef Value*  argument_type;


/*!
*/
typedef Value result_type;

/// \name Creation
/// @{
/*!
default constructor.
*/
Dereference();



/// @}


/// \name Operations
/// @{
/*!
returns
`*x`.
*/
result_type& operator()(argument_type& x) const;



/// @}


/// \name Operations
/// @{
/*!
returns `*x`.
*/
const result_type& operator()(const argument_type&
x) const;



/// @}



}; /* end Dereference */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup projectionobjects



The class `Get_address` gets the address of an lvalue
(`operator&`).



\cgalModels `ProjectionObject`


*/
template< typename Value >
struct Get_address {


/*!
*/
typedef Value argument_type;

/*!
*/
typedef  Value* result_type;

/// \name Creation
/// @{
/*!
default constructor.
*/
Get_address();



/// @}


/// \name Operations
/// @{
/*!
returns
`&x`.
*/
result_type& operator()(argument_type& x) const;



/// @}


/// \name Operations
/// @{
/*!
returns `&x`.
*/
const result_type& operator()(const argument_type&
x) const;



/// @}



}; /* end Get_address */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup projectionobjects



The class `Identity` represents the identity function
on `Value`.



\cgalModels `ProjectionObject`

*/
template< typename Value >
struct Identity {


/*!
*/
typedef Value argument_type;

/*!
*/
typedef Value result_type;

/// \name Creation
/// @{
/*!
default constructor.
*/
Identity();



/// @}


/// \name Operations
/// @{
/*!
returns
`x`.
*/
result_type& operator()(argument_type& x) const;



/// @}


/// \name Operations
/// @{
/*!
returns `x`.
*/
const result_type& operator()(const argument_type&
x) const;



/// @}



}; /* end Identity */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup projectionobjects



The class `Project_facet` calls the member function
`facet()` on an instance of type `Node`.



\cgalModels `ProjectionObject`

*/
template< typename Node >
struct Project_facet {


/*!
*/
typedef Node argument_type;

/*!
*/
typedef Node::Facet result_type;


/// \name Creation
/// @{
/*!
default constructor.
*/
Project_facet();



/// @}


/// \name Operations
/// @{
/*!
returns `n.facet()`.
*/
result_type& operator()(argument_type& n) const;



/// @}


/// \name Operations
/// @{
/*!
returns `n.facet()`.
*/
const result_type& operator()(const argument_type&
n) const;



/// @}



}; /* end Project_facet */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup projectionobjects



The class `Project_next` calls the member function
`next()` on an instance of type `Node`.



\cgalModels `ProjectionObject`


*/
template< typename Node >
struct Project_next {


/*!
*/
typedef Node* argument_type;

/*!
*/
typedef Node* result_type;

/// \name Creation
/// @{
/*!
default constructor.
*/
Project_next();



/// @}


/// \name Operations
/// @{
/*!
returns  `n->next()`.
*/
result_type& operator()(argument_type& n) const;



/// @}


/// \name Operations
/// @{
/*!
returns `n->next()`.
*/
const result_type& operator()(const argument_type&
n) const;



/// @}



}; /* end Project_next */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup projectionobjects



The class `Project_next_opposite` calls the member functions
`next()->opposite()` on an instance of type `Node`.



\cgalModels `ProjectionObject`

*/
template< typename Node >
struct Project_next_opposite {


/*!
*/
typedef Node* argument_type;

/*!
*/
typedef Node* result_type;


/// \name Creation
/// @{
/*!
default constructor.
*/
Project_next_opposite();



/// @}


/// \name Operations
/// @{
/*!
returns
`n->next()->opposite()`.
*/
result_type& operator()(argument_type& n) const;



/// @}


/// \name Operations
/// @{
/*!
returns `n->next()->opposite()`.
*/
const result_type& operator()(const argument_type&
n) const;



/// @}



}; /* end Project_next_opposite */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup projectionobjects



The class `Project_normal` calls the member function
`normal()` on an instance of type `Node`.



\cgalModels `ProjectionObject`


*/
template< typename Node >
struct Project_normal {


/*!
*/
typedef Node argument_type;


/*!
*/
  typedef Node::Normal result_type;

/// \name Creation
/// @{
/*!
default constructor.
*/
Project_normal();



/// @}


/// \name Operations
/// @{
/*!
returns
`n.normal()`.
*/
result_type& operator()(argument_type& n) const;



/// @}


/// \name Operations
/// @{
/*!
returns `n.normal()`.
*/
const result_type& operator()(const argument_type&
n) const;



/// @}



}; /* end Project_normal */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup projectionobjects



The class `Project_opposite_prev` calls the member functions
`opposite()->prev()` on an instance of type `Node`.



\cgalModels `ProjectionObject`


*/
template< typename Node >
struct Project_opposite_prev {


/*!
*/
typedef Node* argument_type;

/*!
*/
typedef Node* result_type;


/// \name Creation
/// @{
/*!
default constructor.
*/
Project_opposite_prev();



/// @}


/// \name Operations
/// @{
/*!
returns
`n->opposite()->prev()`.
*/
result_type& operator()(argument_type& n) const;



/// @}


/// \name Operations
/// @{
/*!
returns `n->opposite()->prev()`.
*/
const result_type& operator()(const argument_type&
n) const;



/// @}



}; /* end Project_opposite_prev */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup projectionobjects



The class `Project_plane` calls the member function
`plane()` on an instance of type `Node`.



\cgalModels `ProjectionObject`


*/
template< typename Node >
struct Project_plane {


/*!
*/
typedef Node argument_type;

/*!
*/
typedef Node::Plane result_type;


/// \name Creation
/// @{
/*!
default constructor.
*/
Project_plane();



/// @}


/// \name Operations
/// @{
/*!
returns
`n.plane()`.
*/
result_type& operator()(argument_type& n) const;



/// @}


/// \name Operations
/// @{
/*!
returns `n.plane()`.
*/
const result_type& operator()(const argument_type&
n) const;



/// @}



}; /* end Project_plane */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup projectionobjects



The class `Project_point` calls the member function
`point()` on an instance of type `Node`.



\cgalModels `ProjectionObject`

*/
template< typename Node >
struct Project_point {


/*!
*/
typedef Node argument_type;

/*!
*/
typedef Node::Point result_type;


/// \name Creation
/// @{
/*!
default constructor.
*/
Project_point();



/// @}


/// \name Operations
/// @{
/*!
returns
`n.point()`.
*/
result_type& operator()(argument_type& n) const;



/// @}


/// \name Operations
/// @{
/*!
returns `n.point()`.
*/
const result_type& operator()(const argument_type&
n) const;



/// @}



}; /* end Project_point */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup projectionobjects



The class `Project_prev` calls the member function
`prev()` on an instance of type `Node`.



\cgalModels `ProjectionObject`

*/
template< typename Node >
struct Project_prev {


/*!
*/
typedef Node* argument_type;

/*!

*/
typedef Node* result_type;


/// \name Creation
/// @{
/*!
default constructor.
*/
Project_prev();



/// @}


/// \name Operations
/// @{
/*!
returns  `n->prev()`.
*/
result_type& operator()(argument_type& n) const;



/// @}


/// \name Operations
/// @{
/*!
returns `n->prev()`.
*/
const result_type& operator()(const argument_type&
n) const;



/// @}



}; /* end Project_prev */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup projectionobjects



The class `Project_vertex` calls the member function
`vertex()` on an instance of type `Node`.



\cgalModels `ProjectionObject`

*/
template< typename Node >
struct Project_vertex {


/*!
*/
typedef Node argument_type;

/*!
*/
typedef Node::Vertex result_type;


/// \name Creation
/// @{
/*!
default constructor.
*/
Project_vertex();



/// @}


/// \name Operations
/// @{
/*!
returns
`n.vertex()`.
*/
result_type& operator()(argument_type& n) const;



/// @}


/// \name Operations
/// @{
/*!
returns `n.vertex()`.
*/
const result_type& operator()(const argument_type&
n) const;



/// @}



}; /* end Project_vertex */
} /* end namespace CGAL */
