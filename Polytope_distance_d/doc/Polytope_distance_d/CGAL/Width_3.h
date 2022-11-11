
namespace CGAL {

/*!
\ingroup PkgPolytopeDistanceDRef

Given a set of points \f$ \mathcal{S}=\left\{p_1,\ldots , p_n\right\}\f$ in
\f$ \mathbb{R}^3\f$. The width of \f$ \mathcal{S}\f$, denoted as \f$ \mathcal{W(S)}\f$, is defined
as the minimum distance between two parallel planes of support of
\f$ \mathit{conv(\mathcal{S})}\f$; where \f$ \mathit{conv(\mathcal{S})}\f$ denotes
the convex hull of \f$ \mathcal{S}\f$. The width in direction \f$ \mathbf{d}\f$,
denoted as \f$ \mathcal{W}_d\mathcal{(S)}\f$, is the distance between two
parallel planes of support of \f$ \mathit{conv(\mathcal{S})}\f$, which are
orthogonal to \f$ \mathbf{d}\f$.

Subject to the applications of the width algorithm, several objects
might be interesting:
<OL>
<LI>The two parallel planes of support such that the distance
between them is as small as possible. These planes are called
width-planes in further considerations.
<LI>The width \f$ \mathcal{W(S)}\f$, i.e., the distance between the
width-planes.
<LI>The direction \f$ \mathbf{d}_{opt}\f$ such that
\f$ \mathcal{W(S)}=\mathcal{W}_{d_{opt}}\mathcal{(S)}\f$
</OL>

<I>Note:</I> There might be several optimal build directions. Hence
neither the width-planes nor the direction \f$ \mathbf{d}_{opt}\f$ are
unique - only the width is.

\tparam Traits must be a model for `WidthTraits_3`.

We provide the model `Width_default_traits_3<Kernel>` based on a
three-dimensional \cgal kernel.

\sa `CGAL::Width_default_traits_3<K>`
\sa `WidthTraits_3`

\cgalHeading{Implementation}

Since the width of the point set \f$ \mathcal{S}\f$ and the width of the convex
hull of \f$ \mathcal{S}\f$ (\f$ \mathit{conv(\mathcal{S})}\f$) is the same, the
algorithm uses the 3D convex hull algorithm \cgal provides.

The width-algorithm is not incremental and therefore inserting and
erasing points cause not an 'automatic' update of the width. Instead
you have to run the width-algorithm again even if the point set is
extended by only one new point.

<b>Large Numbers.</b>

Because there is no need for dividing values during the algorithm, the
numbers can get really huge (all the computations are made using a lot
of multiplications). Therefore it is strongly recommended to use a
number type that can handle numbers of arbitrary length (e.g.,
`leda_integer` in combination with the homogeneous representation
of the points). But these large numbers have a disadvantage:
Operations on them are slower as greater the number gets. Therefore it
is possible to shorten the numbers by using the compiler flag
<span class="textsc">-Dsimplify</span>. For using this option it is required that
the underlying number type provides the 'modulo' operation.

<b>Information Output during the Computations.</b>

If during the algorithm the program should output some information
(e.g., during the debugging phase) you can turn on the output
information by giving the compiler flag <span class="textsc">debug</span>. In the file
<TT>width_assertions.h</TT> you can turn on/off the output of some
functions and additional informations by changing the defined values
from 0 (no output) to 1 (output available). But then it is required
that the `operator<<()` has to been overloaded for `Point_3`,
`Plane_3`, `Vector_3` and `RT`.

\cgalHeading{Example}

\cgalExample{Polytope_distance_d/width_simplex.cpp}

*/
template< typename Traits >
class Width_3 {
public:

/// \name Types
/// @{

/*!
traits class.
*/
typedef unspecified_type Traits;

/*!
point type.
*/
typedef typename Traits::Point_3 Point_3;

/*!
plane type.
*/
typedef typename Traits::Plane_3 Plane_3;

/*!
vector type.
*/
typedef typename Traits::Vector_3 Vector_3;

/*!
algebraic ring type.
*/
typedef typename Traits::RT RT;

/*!
traits
class for the 3D convex hull algorithm.
*/
typedef typename Traits::ChullTraits ChullTraits;

/// @}

/// \name Creation
/// @{

/*!

creates a variable `width` initialized to the width of \f$ \mathcal{S}\f$ -
with \f$ \mathcal{S}\f$ being the set of points in the range
[`first`,`beyond`).

\tparam InputIterator has `Point_3` as value type.
*/
template < class InputIterator >
Width_3( InputIterator first, InputIterator beyond);

/*!
creates a variable `width` initialized to
the width of the polyhedron \f$ P\f$. Note that the vertex point coordinates
are altered!
\pre \f$ P\f$ is a convex polyhedron.

\tparam Polyhedron is a
`CGAL::Polyhedron_3` with facets supporting plane equations
where `Polyhedron::Point_3` \f$ \equiv\f$ `Point_3` and
`Polyhedron::Plane_3` \f$ \equiv\f$ `Plane_3`.
*/
template < class Polyhedron >
Width_3( Polyhedron& P);

/// @}

/// \name Access Functions
/// @{

/*!
returns the squared width. For the reason of exact
computation not the width itself is stored, but the <I>squared</I>
width as a fraction: The numerator in `width_num` and the
denominator in `width_denom`. The width of the point set
\f$ \mathcal{S}\f$ is
\f$ \sqrt{\frac{width\_num}{width\_denom}}\f$.
*/
void get_squared_width ( RT& width_num, RT&
width_denom );

/*!

The planes `e1` and `e2` are the two parallel supporting
planes, which distance is minimal (among all such planes).
*/
void get_width_planes ( Plane_3& e1, Plane_3& e2 );

/*!
The returned coefficients `A,B,C,D,K` have the
property that width-plane `e1` is given by the equation
\f$ Ax+By+Cz+D=0\f$ and width-plane `e2` by \f$ Ax+By+Cz+K=0\f$.
*/
void get_width_coefficients ( RT& A,RT& B,RT& C,
RT& D, RT& K );

/*!
returns a
direction \f$ \mathbf{d}_{opt}\f$ such that the width-planes `e1` and
`e2` are perpendicular to \f$ \mathbf{d}_{opt}\f$. The width of the
point set is minimal in this direction.
*/
Vector_3 get_build_direction ( );

/*!
All the build directions are stored in the vector
`dir`. It might happen that a certain body has several
different build directions, but it is also possible to have only one
build direction.
*/
void get_all_build_directions ( std::vector<Vector_3>&
dir );

/*!
returns
the number of optimal solutions, i.e., the number of optimal build
directions.
*/
int get_number_of_optimal_solutions( );

/// @}

}; /* end Width_3 */
} /* end namespace CGAL */
