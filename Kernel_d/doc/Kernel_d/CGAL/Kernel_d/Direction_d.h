namespace CGAL {

/*!
\ingroup PkgKernelDKernelObjs

A `Direction_d` is a vector in the \f$ d\f$-dimensional vector space
where we forget about its length. We represent directions in
\f$ d\f$-dimensional space as a tuple \f$ (h_0,\ldots,h_d)\f$ of variables of
type `RT` which we call the homogeneous coordinates of the
direction. The coordinate \f$ h_d\f$ must be positive. The %Cartesian
coordinates of a direction are \f$ c_i = h_i/h_d\f$ for \f$ 0 \le i < d\f$,
which are of type `FT`. Two directions are equal if their
%Cartesian coordinates are positive multiples of each other. Directions
are in one-to-one correspondence to points on the unit sphere.

\cgalHeading{Downward compatibility}

We provide the operations of the lower dimensional interface `dx()`,
`dy()`, `dz()`.

\cgalHeading{Implementation}

Directions are implemented by arrays of integers as an item type. All
operations like creation, initialization, tests, inversion, input and
output on a direction \f$ d\f$ take time \cgalBigO{d.\mathit{dimension}()}.
`dimension()`, coordinate access and conversion take constant
time. The space requirement is \cgalBigO{d.\mathit{dimension}()}.

*/
template< typename Kernel >
class Direction_d {
public:

/// \name Types
/// @{

/*!
the linear algebra layer.
*/
typedef unspecified_type LA;

/*!

a read-only iterator for the deltas of `dir`.
*/
typedef unspecified_type Delta_const_iterator;

/*!
construction tag.
*/
typedef unspecified_type Base_direction;

/// @}

/// \name Creation
/// @{

/*!
introduces a variable `dir` of
type `Direction_d<Kernel>`.
*/
Direction_d<Kernel>();

/*!

introduces a variable `dir` of type `Direction_d<Kernel>`
initialized to the direction of `v`.
*/
Direction_d<Kernel>(Vector_d<Kernel> v);

/*!

introduces a variable `dir` of type `Direction_d<Kernel>` in
dimension `d` with representation tuple `set [first,last)`.
\pre `d` is nonnegative, `[first,last)` has `d` elements.
\tparam InputIterator has `RT` as value type.
*/
template <class InputIterator>
Direction_d<Kernel>(int d, InputIterator first, InputIterator last);

/*!
returns
a variable `dir` of type `Direction_d<Kernel>` initialized to the
direction of the \f$ i\f$-th base vector of dimension \f$ d\f$.

\pre \f$ 0 \leq i < d\f$.
*/
Direction_d<Kernel>(int d, Base_direction, int i);

/*!

introduces a variable `dir` of type `Direction_d<Kernel>` in
\f$ 2\f$-dimensional space.
*/
Direction_d<Kernel>(RT x, RT y);

/*!

introduces a variable `dir` of type `Direction_d<Kernel>` in
\f$ 3\f$-dimensional space.
*/
Direction_d<Kernel>(RT x, RT y, RT z);

/// @}

/// \name Operations
/// @{

/*!

returns the dimension of `dir`.
*/
int dimension() ;

/*!

returns the \f$ i\f$-th component of `dir`.
\pre \f$ 0 \leq i < d\f$.
*/
RT delta(int i) ;

/*!

returns the \f$ i\f$-th delta of `dir`.
\pre \f$ 0 \leq i < d\f$.
*/
RT operator[](int i) ;

/*!

returns an iterator pointing to the first delta of `dir`.
*/
Delta_const_iterator deltas_begin() ;

/*!

returns an iterator pointing beyond the last delta of `dir`.
*/
Delta_const_iterator deltas_end() ;

/*!

returns a vector pointing in direction `dir`.
*/
Vector_d<Kernel> vector() ;

/*!

returns true iff `dir.delta(i)==0` for all \f$ 0\leq i < d\f$.
*/
bool is_degenerate() ;

/*!

returns \f$ t(p)\f$.
*/
Direction_d<Kernel> transform(const Aff_transformation_d<Kernel>& t) ;

/*!

returns the direction opposite to `dir`.
*/
Direction_d<Kernel> opposite() ;

/*!

returns the direction opposite to `dir`.
*/
Direction_d<Kernel> operator- () ;

/// @}

}; /* end Direction_d */
} /* end namespace CGAL */
