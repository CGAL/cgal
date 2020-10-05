
namespace CGAL {

/*!
\ingroup PkgInscribedAreasRef

Given a set of points in the plane, the class `Largest_empty_iso_rectangle_2` is a data
structure that maintains an iso-rectangle with the largest area among
all iso-rectangles that are inside a given bounding box( iso-rectangle), and
that do not contain any point of the point set.

\tparam T must be a model of the concept `LargestEmptyIsoRectangleTraits_2`.

\cgalHeading{Implementation}

The algorithm is an implementation of \cgalCite{o-naler-90}. The runtime of an
insertion or a removal is \f$ O(\log n)\f$. A query takes \f$ O(n^2)\f$ worst
case time and \f$ O(n \log n)\f$ expected time. The working storage is \f$
O(n)\f$.

*/
template< typename T >
class Largest_empty_iso_rectangle_2 {
public:

/// \name Types
/// @{

/*!

*/
typedef T Traits;

/*!

*/
typedef Traits::Point_2 Point_2;

/*!

*/
typedef Traits::Iso_rectangle_2 Iso_rectangle_2;

/*!
Iterator over the points.

This iterator allows to enumerate the points. It is non mutable,
bidirectional and its value type is `Point_2`. It is invalidated by
any insertion or removal of a point.
*/
typedef unspecified_type const_iterator;

/// @}

/// \name Creation
/// @{

/*!
Constructor. The iso-rectangle `b` is the bounding rectangle.
*/
Largest_empty_iso_rectangle_2<Traits>
(const Iso_rectangle_2 &b);

/*!
Constructor. The iso-rectangle whose lower left and upper right points are `p` and
`q` respectively is the bounding rectangle.
*/
Largest_empty_iso_rectangle_2<Traits>
(const Point_2 p,const Point_2 q);

/*!
Constructor. The iso-rectangle whose lower left point and upper right points are (0,0)
and (1,1) respectively is the bounding rectangle.
*/
Largest_empty_iso_rectangle_2<Traits>
();

/*!
Copy constructor.
*/
Largest_empty_iso_rectangle_2<Traits>
(const Largest_empty_iso_rectangle_2<Traits> tr);

/// @}

/// \name Assignment
/// @{

/*!

*/
Largest_empty_iso_rectangle_2<T>
operator=(const Largest_empty_iso_rectangle_2<T> & tr);

/// @}

/// \name Access Functions
/// @{

/*!
Returns a const reference to the traits object.
*/
const Traits & traits() const;

/*!
Returns an iterator to the beginning of the point set.
*/
const_iterator begin() const;

/*!
Returns a past-the-end iterator for the point set.
*/
const_iterator end() const;

/// @}

/// \name Queries
/// @{

/*!
Returns the four points that define the largest empty iso-rectangle.
(Note that these points are not necessarily on a corner of an iso-rectangle.)
*/
Quadruple<Point_2, Point_2, Point_2, Point_2>
get_left_bottom_right_top();

/*!
Returns the largest empty iso-rectangle. (Note that the two
points defining the iso-rectangle are not necessarily part of
the point set.)
*/
Iso_rectangle_2 get_largest_empty_iso_rectangle();

/*!
Returns the iso-rectangle passed in the constructor.
*/
Iso_rectangle_2 get_bounding_box();

/// @}

/// \name Insertion
/// @{

/*!
Inserts point `p` in the point set, if it is not already in the set.
*/
void
insert(const Point_2& p);

/*!
Inserts point `p` in the point set, if it is not already in the set.
*/
void
push_back(const Point_2& p);

/*!
Inserts the points in the range `[first, last)`, and returns the number of inserted points.

\tparam InputIterator must be an iterator with value type `Point_2`.
*/
template < class InputIterator >
int
insert(InputIterator first, InputIterator last);

/// @}

/// \name Removal
/// @{

/*!
Removes point `p`.
Returns false iff `p` is not in the point set.
*/
bool remove(const Point_2& p);

/*!
Removes all points.
*/
void clear();

/// @}

}; /* end Largest_empty_iso_rectangle_2 */
} /* end namespace CGAL */
