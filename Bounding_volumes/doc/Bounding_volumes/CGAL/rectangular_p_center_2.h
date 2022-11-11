namespace CGAL {

/*!
\ingroup PkgBoundingVolumesRef

The class `Rectangular_p_center_default_traits_2` defines types and operations
needed to compute rectilinear \f$ p\f$-centers of a planar point set
using the function `rectangular_p_center_2()`.


\tparam K must be a model for `Kernel`.

\cgalModels `RectangularPCenterTraits_2`

\sa `CGAL::rectangular_p_center_2()`

*/
template< typename K >
class Rectangular_p_center_default_traits_2 {
public:

/// \name Types
/// @{

/*!
typedef to `K::FT`.
*/
typedef unspecified_type FT;

/*!
typedef to `K::Point_2`.
*/
typedef unspecified_type Point_2;

/*!
typedef to `K::Iso_rectangle_2`.
*/
typedef unspecified_type Iso_rectangle_2;

/*!
typedef to `K::Less_x_2`.
*/
typedef unspecified_type Less_x_2;

/*!
typedef to `K::Less_y_2`.
*/
typedef unspecified_type Less_y_2;

/*!
typedef to
`K::Construct_vertex_2`.
*/
typedef unspecified_type Construct_vertex_2;

/*!
typedef to
`K::Construct_iso_rectangle_2`.
*/
typedef unspecified_type Construct_iso_rectangle_2;

/*!
adaptable binary function
class: `Point_2` \f$ \times\f$ `Point_2` \f$ \rightarrow\f$
`FT` returns the signed distance of two points'
\f$ x\f$-coordinates.
*/
typedef unspecified_type Signed_x_distance_2;

/*!
adaptable binary function
class: `Point_2` \f$ \times\f$ `Point_2` \f$ \rightarrow\f$
`FT` returns the signed distance of two points'
\f$ y\f$-coordinates.
*/
typedef unspecified_type Signed_y_distance_2;

/*!
adaptable binary function
class: `Point_2` \f$ \times\f$ `Point_2` \f$ \rightarrow\f$
`FT` returns the \f$ ||\cdot||_{\infty}\f$ distance of two
points.
*/
typedef unspecified_type Infinity_distance_2;

/*!
adaptable binary
function class: `Point_2` \f$ \times\f$ `Point_2`
\f$ \rightarrow\f$ `FT` returns the signed \f$ ||\cdot||_{\infty}\f$
distance of two points.
*/
typedef unspecified_type Signed_infinity_distance_2;

/*!

3-argument function class: `Point_2` \f$ \times\f$ `Point_2`
\f$ \times\f$ `FT` \f$ \rightarrow\f$ `Point_2`. For arguments
\f$ (p,\,q,\,r)\f$ it returns the lower-left corner of the iso-oriented
square with sidelength \f$ r\f$ and upper-right corner at the
intersection of the vertical line through \f$ p\f$ and the horizontal
line through \f$ q\f$.
*/
typedef unspecified_type Construct_point_2_below_left_implicit_point_2;

/*!

3-argument function class: `Point_2` \f$ \times\f$ `Point_2`
\f$ \times\f$ `FT` \f$ \rightarrow\f$ `Point_2`. For arguments
\f$ (p,\,q,\,r)\f$ it returns the lower-right corner of the
iso-oriented square with sidelength \f$ r\f$ and upper-left corner at
the intersection of the vertical line through \f$ p\f$ and the
horizontal line through \f$ q\f$.
*/
typedef unspecified_type Construct_point_2_below_right_implicit_point_2;

/*!

3-argument function class: `Point_2` \f$ \times\f$ `Point_2`
\f$ \times\f$ `FT` \f$ \rightarrow\f$ `Point_2`. For arguments
\f$ (p,\,q,\,r)\f$ it returns the upper-right corner of the
iso-oriented square with sidelength \f$ r\f$ and lower-left corner at
the intersection of the vertical line through \f$ p\f$ and the
horizontal line through \f$ q\f$.
*/
typedef unspecified_type Construct_point_2_above_right_implicit_point_2;

/*!

3-argument function class: `Point_2` \f$ \times\f$ `Point_2`
\f$ \times\f$ `FT` \f$ \rightarrow\f$ `Point_2`. For arguments
\f$ (p,\,q,\,r)\f$ it returns the upper-left corner of the iso-oriented
square with sidelength \f$ r\f$ and lower-right corner at the
intersection of the vertical line through \f$ p\f$ and the horizontal
line through \f$ q\f$.
*/
typedef unspecified_type Construct_point_2_above_left_implicit_point_2;

/// @}

/// \name Operations
/// For every function class listed above there is a member function
/// to fetch the corresponding function object.
/// @{

/*!

*/
Inf_distance_2 inf_distance_2_object() const;

/*!

*/
Signed_inf_distance_2
signed_inf_distance_2_object() const;

/*!

*/
Construct_vertex_2
construct_vertex_2_object() const;

/*!

*/
Construct_iso_rectangle_2
construct_iso_rectangle_2_object() const;

/*!

*/
Construct_iso_rectangle_2_below_left_point_2
construct_iso_rectangle_2_below_left_point_2_object() const;

/*!

*/
Construct_iso_rectangle_2_above_left_point_2
construct_iso_rectangle_2_above_left_point_2_object() const;

/*!

*/
Construct_iso_rectangle_2_below_right_point_2
construct_iso_rectangle_2_below_right_point_2_object() const;

/*!

*/
Construct_iso_rectangle_2_above_right_point_2
construct_iso_rectangle_2_above_right_point_2_object() const;

/// @}

};

/*!
\ingroup PkgBoundingVolumesRef

Computes rectilinear
\f$ p\f$-centers of a planar point set, i.e.\ a set of \f$ p\f$ points such
that the maximum minimal \f$ L_{\infty}\f$-distance between both sets is
minimized.

More formally the problem can be defined as follows.

Given a finite set \f$ \mathcal{P}\f$ of points, compute a
point set \f$ \mathcal{C}\f$ with \f$ |\mathcal{C}| \le p\f$ such that the
\f$ p\f$-radius of \f$ \mathcal{P}\f$,
\f[
rad_p(\mathcal{P}) := \max_{P \in \mathcal{P}} \min_{Q \in
\mathcal{C}} || P - Q ||_\infty
\f]
is minimized. We can interpret \f$ \mathcal{C}\f$ as the best
approximation (with respect to the given metric) for \f$ \mathcal{P}\f$
with at most \f$ p\f$ points.

computes rectilinear `p`-centers for the point set described by
the range [`f`, `l`), sets `r` to the corresponding
\f$ p\f$-radius, writes the at most `p` center points to `o` and
returns the past-the-end iterator of this sequence.

\pre 2 \f$ \le\f$ `p` \f$ \le\f$ 4.

The geometric types and operations to be used for the computation
are specified by the traits class parameter `t`. This parameter
can be omitted if `ForwardIterator` refers to a point type from
the 2D-Kernel. In this case, a default traits class
(`Rectangular_p_center_default_traits_2<K>`) is used.

<OL>
<LI><I>Either: (if no traits parameter is given)</I> Value type
of `ForwardIterator` must be `CGAL::Point_2<K>` for some
representation class `K` and `FT` must be equivalent to
`K::FT`,
<LI><I>Or: (if a traits parameter is specified)</I> `Traits`
must be a model for `RectangularPCenterTraits_2`.
<LI>`OutputIterator` must accept the value type of
`ForwardIterator` as value type.
</OL>

\sa `RectangularPCenterTraits_2`
\sa `CGAL::Rectangular_p_center_default_traits_2<K>`
\sa `CGAL::sorted_matrix_search()`

\cgalHeading{Implementation}

The runtime is linear for \f$ p \in \{2,\,3\}\f$ and
\f$ \mathcal{O}(n \cdot \log n)\f$ for \f$ p = 4\f$ where \f$ n\f$ is the number of
input points. These runtimes are worst case optimal. The \f$ 3\f$-center
algorithm uses a prune-and-search technique described in
\cgalCite{cgal:h-slacr-99}. The \f$ 4\f$-center implementation uses sorted matrix
search \cgalCite{fj-fkppc-83}, \cgalCite{fj-gsrsm-84} and fast algorithms for
piercing rectangles \cgalCite{sw-rpppp-96}.

\cgalHeading{Example}

The following code generates a random set of ten points
and computes its two-centers.

\cgalExample{Rectangular_p_center_2/rectangular_p_center_2.cpp}

*/
template < class ForwardIterator, class
OutputIterator, class FT, class Traits > OutputIterator
rectangular_p_center_2(ForwardIterator f, ForwardIterator l, OutputIterator o,
                       FT& r, int p, const Traits& t = Default_traits);
} /* namespace CGAL */

