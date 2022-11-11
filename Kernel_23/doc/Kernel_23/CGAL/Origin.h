
namespace CGAL {

/*!
\ingroup PkgKernel23Ref

\cgal defines a symbolic constant
`NULL_VECTOR` to construct zero length vectors.
`Null_vector` is the type of this constant.

\sa `CGAL::Vector_2<Kernel>`
\sa `CGAL::Vector_3<Kernel>`

*/

class Null_vector {
public:
}; /* end Null_vector */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgKernel23Ref

\cgal defines a symbolic constant
`ORIGIN` which denotes the point at the origin.
`Origin` is the type of this constant.
It is used in the conversion between points and vectors.

\sa `CGAL::Point_2<Kernel>`
\sa `CGAL::Point_3<Kernel>`
\sa `CGAL::Vector_2<Kernel>`
\sa `CGAL::Vector_3<Kernel>`
\sa `CGAL::operator+`
\sa `CGAL::operator-`

*/

class Origin {
public:
}; /* end Origin */

/*!
A symbolic constant which denotes the point at the origin.  This
constant is used in the conversion between points and vectors.

\cgalHeading{Example}

\code
Point_2< Cartesian<Exact_NT> > p(1.0, 1.0), q;
Vector2< Cartesian<Exact_NT> > v;
v = p - ORIGIN;
q = ORIGIN + v;
assert( p == q );
\endcode

\sa `CGAL::Point_2<Kernel>`
\sa `CGAL::Point_3<Kernel>`
*/
const Origin ORIGIN;

} /* end namespace CGAL */
