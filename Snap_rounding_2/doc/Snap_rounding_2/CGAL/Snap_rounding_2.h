namespace CGAL {

/*!
\ingroup PkgSnapRounding2Ref

\tparam Traits must be a model of `SnapRoundingTraits_2`.
\tparam InputIterator must be an iterator with value type `Traits::Segment_2`.
\tparam OutputContainer must be a container with a method `push_back(const OutputContainer::value_type& c)`,
where `OutputContainer::value_type` must be a container with a method `push_back(const Traits::Point_2& p)`

\param begin,end The first two parameters denote the iterator range
of the input segments.

\param output_container is a reference to a
container of the output polylines. Since a polyline is composed of a
sequence of points, a polyline is a container itself.

\param do_isr The fifth parameter determines whether to apply ISR or SR.

\param pixel_size The fourth parameter denotes the pixel size `w`. The plane will be
tiled with square pixels of width `w` such that the origin is the
center of a pixel. `pixel_size` must have a positive value.

\param int_output The sixth parameter denotes the output
representation. If the value of the sixth parameter is `true`
then the centers of pixels constitute the integer grid, and hence the
vertices of the output polylines will be integers. For example, the
coordinates of the center of the pixel to the right of the pixel
containing the origin will be `(1,0)` regardless of the pixel width.
If the value of the sixth parameter is `false`, then the centers
of hot pixels (and hence the vertices of the output polylines) will
bear their original coordinates, which may not necessarily be
integers. In the latter case, the coordinates of the center of the
pixel to the right of the pixel containing the origin, for example,
will be `(w,0)`.

\param number_of_kd_trees The seventh parameter is briefly described later on this page; for a
detailed description see \cgalCite{cgal:hp-isr-02}.

Snap Rounding (SR, for short) is a well known method for converting
arbitrary-precision arrangements of segments into a fixed-precision
representation \cgalCite{gght-srlse-97}, \cgalCite{gm-rad-98}, \cgalCite{h-psifp-99}. In
the study of robust geometric computing, it can be classified as a
finite precision approximation technique. Iterated Snap Rounding (ISR,
for short) is a modification of SR in which each vertex is at least
half-the-width-of-a-pixel away from any non-incident edge
\cgalCite{cgal:hp-isr-02}. This package supports both methods. Algorithmic
details and experimental results are given in \cgalCite{cgal:hp-isr-02}.

Given a finite collection \f$ \Sc\f$ of segments in the plane, the
arrangement of \f$ \Sc\f$ denoted \f$ \Ac(\Sc)\f$ is the subdivision of the plane
into vertices, edges, and faces induced by \f$ \Sc\f$. A <I>vertex</I> of the arrangement is either a segment endpoint or
the intersection of two segments. Given an arrangement of segments
whose vertices are represented with arbitrary-precision coordinates,
SR proceeds as follows. We tile the plane
with a grid of unit squares, <I>pixels</I>, each centered at a point
with integer coordinates. A pixel is <I>hot</I> if it contains a
vertex of the arrangement. Each vertex of the arrangement is replaced
by the center of the hot pixel containing it and each edge \f$ e\f$ is
replaced by the polygonal chain through the centers of the hot pixels
met by \f$ e\f$, in the same order as they are met by \f$ e\f$.

In a snap-rounded arrangement, the distance between a vertex and
a non-incident edge can be extremely small compared with the width of a
pixel in the grid used for rounding. ISR is a modification of SR which
makes a vertex and a non-incident edge well separated (the distance
between each is at least half-the-width-of-a-pixel). However, the
guaranteed quality of the approximation in ISR degrades. For more
details on ISR see \cgalCite{cgal:hp-isr-02}.

The traits used here must support (arbitrary-precision) rational number type as
this is a basic requirement of SR.



\cgalHeading{About the Number of kd-Trees}

A basic query used in the algorithm is to report the hot pixels of
size \f$ w\f$ that a certain segment \f$ s\f$ intersects. An alternative way to
do the same is to query the hot pixels' centers contained in a
Minkowski sum of \f$ s\f$ with a pixel of width \f$ w\f$ centered at the origin;
we denote this Minkowski sum by \f$ M(s)\f$. Since efficiently implementing
this kind of query is difficult, we use an orthogonal range-search
structure instead. We query with the bounding box \f$ B(M(s))\f$ of \f$ M(s)\f$
in a two-dimensional kd-tree which stores the centers of hot
pixels. Since \f$ B(M(s))\f$ in general is larger than \f$ M(s)\f$, we still
need to filter out the hot pixels which do not intersect \f$ s\f$.

While this approach is easy to implement with CGAL, it may incur
considerable overhead since the area of \f$ B(M(s))\f$ may be much larger
than the area of \f$ M(s)\f$, possibly resulting in many redundant hot pixels
to filter out. Our heuristic solution, which we describe next, is to
use a cluster of kd-trees rather than just one. The cluster includes
several kd-trees, each has the plane, and hence the centers of hot
pixels, rotated by a different angle in the first quadrant of the
plane; for our purpose, a rotation by angles outside this quadrant
is symmetric to a rotation by an angle in the first quadrant.

Given a parameter \f$ c\f$, the angles of rotation are \f$ (i - 1)
\frac{\pi}{2c}, i=1,\ldots,c\f$, and we construct a kd-tree
corresponding to each of these angles. Then for a query segment \f$ s\f$,
we choose the kd-tree for which the area of \f$ B(M(s))\f$ is the smallest,
in order to (potentially) get less hot pixels to filter out. Since
constructing many kd-trees may be costly, our algorithm avoids
building a kd-tree which it expects to be queried a relatively small
number of times (we estimate this number in advance). How many
kd-trees should be used? It is difficult to provide a simple
answer for that. There are inputs for which the time to build more
than one kd-tree is far greater than the time saved by having to
filter out less hot pixels (sparse arrangements demonstrate this
behavior), and there are inputs which benefit from using several
kd-trees. Thus, the user can control the number of kd-trees
with the parameter  `number_of_kd_trees`. Typically, but not
always, one kd-tree (the default) is sufficient.
*/

template < class Traits, class InputIterator, class OutputContainer >
void
snap_rounding_2(
InputIterator begin,
InputIterator end,
OutputContainer& output_container,
typename Traits::FT pixel_size,
bool do_isr = true,
bool int_output = true,
unsigned int number_of_kd_trees = 1);

} /* namespace CGAL */

