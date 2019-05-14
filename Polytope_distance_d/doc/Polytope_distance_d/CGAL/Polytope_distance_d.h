
namespace CGAL {

/*!
\ingroup PkgOptimalDistances

An object of the class `Polytope_distance_d` represents the (squared) distance 
between two convex polytopes, given as the convex hulls of two finite point 
sets in \f$ d\f$-dimensional Euclidean space \f$ \E^d\f$. For point sets \f$ P\f$ and \f$ Q\f$ 
we denote by \f$ pd(P,Q)\f$ the distance between the convex hulls of \f$ P\f$ and 
\f$ Q\f$. Note that \f$ pd(P,Q)\f$ can be 
degenerate, 
i.e.\ \f$ pd(P,Q)=\infty\f$ if \f$ P\f$ or \f$ Q\f$ is empty. 

Two inclusion-minimal subsets \f$ S_P\f$ of \f$ P\f$ and \f$ S_Q\f$ of \f$ Q\f$ with 
\f$ pd(S_P,S_Q)=pd(P,Q)\f$ are called <I>pair of support 
sets</I>, the 
points in \f$ S_P\f$ and \f$ S_Q\f$ are the <I>support points</I>. A pair of support 
sets has size at most \f$ d+2\f$ (by size we mean \f$ |S_P|+|S_Q|\f$). The distance 
between the two polytopes is <I>realized</I> by a pair of points \f$ p\f$ and 
\f$ q\f$ lying in the convex hull of \f$ S_P\f$ and \f$ S_Q\f$, repectively, 
i.e.\ \f$ \sqrt{||p-q||}=pd(P,Q)\f$. In general, neither the support sets nor the 
realizing points are necessarily unique. 

The underlying algorithm can cope with all kinds of input, e.g. \f$ P\f$ and \f$ Q\f$ 
may be in non-convex position or points may occur more than once. The 
algorithm computes a pair of support sets \f$ S_P\f$ and \f$ S_Q\f$ with realizing 
points \f$ p\f$ and \f$ q\f$ which remain fixed until the next set, insert, or clear 
operation. 


\tparam Traits must be a model for `PolytopeDistanceDTraits`. 

We provide the models `Polytope_distance_d_traits_2`, 
`Polytope_distance_d_traits_3`, and `Polytope_distance_d_traits_d` using the 
two-, three-, and \f$ d\f$-dimensional \cgal kernel, respectively. 

\sa `CGAL::Polytope_distance_d_traits_2<K,ET,NT>` 
\sa `CGAL::Polytope_distance_d_traits_3<K,ET,NT>` 
\sa `CGAL::Polytope_distance_d_traits_d<K,ET,NT>` 
\sa `PolytopeDistanceDTraits` 

\cgalHeading{Implementation}

The problem of finding the distance between two convex polytopes given as 
the convex hulls of two finite point sets can be formulated as an 
optimization problem with linear constraints and a convex quadratic 
objective function. The solution is obtained using our exact solver 
for quadratic programs \cgalCite{gs-eegqp-00}. 

The creation time is almost always linear in the number of points. Access 
functions and predicates take constant time, inserting a point might take 
up to linear time. The clear operation and the check for validity each 
take linear time. 

*/
template< typename Traits >
class Polytope_distance_d {
public:

/// \name Types 
/// @{

/*!
typedef to `Traits::Point_d`. 
Point type used to represent the input points. 
*/ 
typedef unspecified_type Point; 

/*!
typedef to `Traits::FT`. 
Number type used to return the squared distance 
between the two polytopes. 
*/ 
typedef unspecified_type FT; 

/*!
typedef to `Traits::ET`. 
Number type used to do the exact computations in the underlying 
solver for quadratic programs (cf. <B>Implementation</B>). 
*/ 
typedef unspecified_type ET; 

/*!

non-mutable model of the \stl concept <I>RandomAccessIterator</I> 
with value type `Point`. Used to access the points 
of the two polytopes. 
*/ 
typedef unspecified_type Point_iterator; 

/*!

non-mutable model of the \stl concept <I>RandomAccessIterator</I> 
with value type `Point`. Used to access the support points. 
*/ 
typedef unspecified_type Support_point_iterator; 

/*!

non-mutable model of the \stl concept <I>RandomAccessIterator</I> 
with value type `int`. Used to access the indices of the 
support points in the provided input order (starting from 0 
in both point sets). 
*/ 
typedef unspecified_type Support_point_index_iterator; 

/*!

non-mutable model of the \stl concept <I>RandomAccessIterator</I> 
with value type `ET`. Used to access the coordinates of 
the realizing points. 
*/ 
typedef unspecified_type Coordinate_iterator; 

/// @} 

/// \name Creation 
/// @{

/*!

initializes `poly_dist` to \f$ pd(\emptyset, \emptyset)\f$. 
*/ 
Polytope_distance_d( const Traits& traits = Traits(), 
int verbose = 0, 
std::ostream& stream = std::cout); 

/*!

initializes `poly_dist` to \f$ pd(P,Q)\f$ with \f$ P\f$ and \f$ Q\f$ being the 
sets of points in the range [`p_first`,`p_last`) and 
[`q_first`,`q_last`), respectively. 
\tparam InputIterator1 has `Point` as value type.
\tparam InputIterator2 has `Point` as value type.
\pre All points have the same dimension. 

\attention If `verbose` is set to \f$ 1\f$, \f$ 2\f$, or
\f$ 3\f$ then some, more, or full verbose output of the underlying
solver for quadratic programs is written to `stream`, resp.
*/ 
template < class InputIterator1, class InputIterator2 > 
Polytope_distance_d( InputIterator1 p_first, 
InputIterator1 p_last, 
InputIterator2 q_first, 
InputIterator2 q_last, 
const Traits& traits = Traits(), 
int verbose = 0, 
std::ostream& stream = std::cout); 

/// @} 

/// \name Access Functions 
/// @{

/*!

returns the dimension of the points in \f$ P\f$ and \f$ Q\f$. 
If `poly_dist` is \f$ pd(\emptyset, \emptyset)\f$, the ambient dimension is \f$ -1\f$. 
*/ 
int ambient_dimension( ) const; 

/*!

returns the number of all points of `poly_dist`, i.e.\ \f$ |P|+|Q|\f$. 
*/ 
int number_of_points( ) const; 

/*!

returns the number of points in \f$ P\f$. 
*/ 
int number_of_points_p( ) const; 

/*!

returns the number of points in \f$ Q\f$. 
*/ 
int number_of_points_q( ) const; 

/*!

returns the number of support points of `poly_dist`, i.e.\ \f$ |S_P|+|S_Q|\f$. 
*/ 
int number_of_support_points( ) const; 

/*!

returns the number of support points in \f$ S_P\f$. 
*/ 
int number_of_support_points_p( ) const; 

/*!

returns the number of support points in \f$ S_Q\f$. 
*/ 
int number_of_support_points_q( ) const; 

/*!

returns an iterator referring to the first point in \f$ P\f$. 
*/ 
Point_iterator points_p_begin( ) const; 

/*!

returns the corresponding past-the-end iterator. 
*/ 
Point_iterator points_p_end( ) const; 

/*!

returns an iterator referring to the first point in \f$ Q\f$. 
*/ 
Point_iterator points_q_begin( ) const; 

/*!

returns the corresponding past-the-end iterator. 
*/ 
Point_iterator points_q_end( ) const; 

/*!

returns an iterator referring to the first support point in \f$ S_P\f$. 
*/ 
Support_point_iterator support_points_p_begin( ) const; 

/*!

returns the corresponding past-the-end iterator. 
*/ 
Support_point_iterator support_points_p_end( ) const; 

/*!

returns an iterator referring to the first support point in \f$ S_Q\f$. 
*/ 
Support_point_iterator support_points_q_begin( ) const; 

/*!

returns the corresponding past-the-end iterator. 
*/ 
Support_point_iterator support_points_q_end( ) const; 

/*!

returns an iterator referring to the index of the first support point in \f$ P\f$. 
*/ 
Support_point_index_iterator support_points_p_indices_begin( ) const; 

/*!

returns the corresponding past-the-end iterator. 
*/ 
Support_point_index_iterator support_points_p_indices_end( ) const; 

/*!

returns an iterator referring to the index of the first support point in \f$ Q\f$. 
*/ 
Support_point_index_iterator support_points_q_indices_begin( ) const; 

/*!

returns the corresponding past-the-end iterator. 
*/ 
Support_point_index_iterator support_points_q_indices_end( ) const; 

/*!

returns the realizing point of \f$ P\f$. 
An implicit conversion from `ET` to `RT` must be available.
\pre \f$ pd(P,Q)\f$ is finite. 
*/ 
Point realizing_point_p( ) const; 

/*!

returns the realizing point of \f$ Q\f$. 
An implicit conversion from `ET` to `RT` must be available.
\pre \f$ pd(P,Q)\f$ is finite. 
*/ 
Point realizing_point_q( ) const; 

/*!

returns the squared distance of `poly_dist`, i.e.\ \f$ (pd(P,Q))^2\f$. 
An implicit conversion from `ET` to `RT` must be available.
\pre \f$ pd(P,Q)\f$ is finite. 
*/ 
FT squared_distance( ) const; 

/*!

returns an iterator referring to the first coordinate of the 
realizing point of \f$ P\f$. 

\note
The coordinates have a rational 
representation, i.e.\ the first \f$ d\f$ elements of the iterator 
range are the numerators and the \f$ (d\!+\!1)\f$-st element is the 
common denominator. 
*/ 
Coordinate_iterator 
realizing_point_p_coordinates_begin() const; 

/*!

returns the corresponding past-the-end iterator. 
*/ 
Coordinate_iterator 
realizing_point_p_coordinates_end() const; 

/*!

returns an iterator referring to the first coordinate of the 
realizing point of \f$ Q\f$. 

\note
The coordinates have a rational 
representation, i.e.\ the first \f$ d\f$ elements of the iterator 
range are the numerators and the \f$ (d\!+\!1)\f$-st element is the 
common denominator. 
*/ 
Coordinate_iterator 
realizing_point_q_coordinates_begin() const; 

/*!

returns the corresponding past-the-end iterator. 
*/ 
Coordinate_iterator 
realizing_point_q_coordinates_end() const; 

/*!

returns the numerator of the squared distance of `poly_dist`. 
*/ 
ET squared_distance_numerator( ) const; 

/*!

returns the denominator of the squared distance of `poly_dist`. 
*/ 
ET squared_distance_denominator( ) const; 

/// @} 

/// \name Predicates 
/// @{

/*!

returns `true`, if \f$ pd(P,Q)\f$ is finite, 
i.e.\ none of the two polytopes is empty. 
*/ 
bool is_finite( ) const; 

/*!

returns `true`, if \f$ pd(P,Q)\f$ is zero, 
i.e.\ the two polytopes intersect (this implies degeneracy). 
*/ 
bool is_zero( ) const; 

/*!

returns `true`, iff \f$ pd(P,Q)\f$ is degenerate, 
i.e.\  \f$ pd(P,Q)\f$ is not finite. 
*/ 
bool is_degenerate( ) const; 

/// @} 

/// \name Modifiers 
/// @{

/*!

resets `poly_dist` to \f$ pd(\emptyset, \emptyset)\f$. 
*/ 
void clear( ); 

/*!

sets `poly_dist` to \f$ pd(P,Q)\f$ with \f$ P\f$ and \f$ Q\f$ being the sets of 
points in the ranges [`p_first`,`p_last`) and 
[`q_first`,`q_last`), respectively. 
\tparam InputIterator1 has `Point` as value type.
\tparam InputIterator2 has `Point` as value type.
\pre All points have the same dimension. 
*/ 
template < class InputIterator1, class InputIterator2 > 
void set( InputIterator1 p_first, 
InputIterator1 p_last, 
InputIterator2 q_first, 
InputIterator2 q_last ); 

/*!

sets `poly_dist` to \f$ pd(P,Q)\f$ with \f$ P\f$ being the set of points 
in the range [`p_first`,`p_last`) (\f$ Q\f$ remains unchanged). 
\tparam InputIterator has `Point` as value type.
\pre All points in \f$ P\f$ have dimension `poly_dist``.ambient_dimension()` if \f$ Q\f$ is not empty. 
*/ 
template < class InputIterator > 
void set_p( InputIterator p_first, 
InputIterator p_last ); 

/*!

sets `poly_dist` to \f$ pd(P,Q)\f$ with \f$ Q\f$ being the set of points 
in the range [`q_first`,`q_last`) (\f$ P\f$ remains unchanged). 
\tparam InputIterator has `Point` as value type.
\pre All points in \f$ Q\f$ have dimension `poly_dist``.ambient_dimension()` if \f$ P\f$ is not empty. 
*/ 
template < class InputIterator > 
void set_q( InputIterator q_first, 
InputIterator q_last ); 

/*!

inserts `p` into \f$ P\f$. 
\pre The dimension of `p` equals `poly_dist``.ambient_dimension()` if `poly_dist` is not \f$ pd(\emptyset, \emptyset)\f$. 
*/ 
void insert_p( const Point& p); 

/*!

inserts `q` into \f$ Q\f$. 
\pre The dimension of `q` equals `poly_dist``.ambient_dimension()` if `poly_dist` is not \f$ pd(\emptyset, \emptyset)\f$. 
*/ 
void insert_q( const Point& q); 

/*!

inserts the points in the range [`p_first`,`p_last`) 
and [`q_first`,`q_last`) into \f$ P\f$ and \f$ Q\f$, respectively, 
and recomputes the (squared) distance. 
\tparam InputIterator1 has `Point` as value type.
\tparam InputIterator2 has `Point` as value type.
\pre All points have the same dimension. If `poly_dist` is not \f$ pd(\emptyset, \emptyset)\f$, this dimension must be equal to `poly_dist``.ambient_dimension()`. 
*/ 
template < class InputIterator1, class InputIterator2 > 
void insert( InputIterator1 p_first, 
InputIterator1 p_last, 
InputIterator2 q_first, 
InputIterator2 q_last ); 

/*!

inserts the points in the range [`p_first`,`p_last`) into 
\f$ P\f$ and recomputes the (squared) distance (\f$ Q\f$ remains unchanged). 
\tparam InputIterator has `Point` as value type.
\pre All points have the same dimension. If `poly_dist` is not empty, this dimension must be equal to `poly_dist``.ambient_dimension()`. 
*/ 
template < class InputIterator > 
void insert_p( InputIterator p_first, 
InputIterator p_last ); 

/*!

inserts the points in the range [`q_first`,`q_last`) into 
\f$ Q\f$ and recomputes the (squared) distance (\f$ P\f$ remains unchanged). 
\tparam InputIterator has `Point` as value type.
\pre All points have the same dimension. If `poly_dist` is not empty, this dimension must be equal to `poly_dist``.ambient_dimension()`. 
*/ 
template < class InputIterator > 
void insert_q( InputIterator q_first, 
InputIterator q_last ); 

/// @} 

/// \name Validity Check 
/// @{

/*!

returns `true`, iff `poly_dist` is valid. If `verbose` is 
`true`, some messages concerning the performed checks are 
written to standard error stream. The second parameter 
`level` is not used, we provide it only for consistency 
with interfaces of other classes. 

An object `poly_dist` is valid, iff \f$ ldots\f$ 
- `poly_dist` contains all points of its defining set \f$ P\f$, 
- `poly_dist` is the smallest sphere containing its support set \f$ S\f$, 
- and \f$ S\f$ is minimal, i.e.\ no support point is redundant.

*/ 
bool is_valid( bool verbose = false, 
int level = 0 ) const; 

/// @} 

/// \name Miscellaneous 
/// @{

/*!

returns a const reference to the traits class object. 
*/ 
const Traits& traits( ) const; 


/// @}

}; /* end Polytope_distance_d */

/*!

writes `poly_dist` to output stream `os`. 
An overload of `operator<<` must be defined for `Point_d`.
\relates Polytope_distance_d 
*/ 
std::ostream& 
operator << ( std::ostream& os, const Polytope_distance_d<Traits>& poly_dist); 

/*!

reads `poly_dist` from input stream `is`. 
An overload of `operator>>` must be defined for `Point_d`.
\relates Polytope_distance_d 
*/ 
std::istream& 
operator >> ( std::istream& is, Polytope_distance_d<Traits> poly_dist&); 

} /* end namespace CGAL */
