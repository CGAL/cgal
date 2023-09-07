namespace CGAL {
namespace Box_intersection_d {

/*!
\ingroup PkgBoxIntersectionDClasses

`Box_with_handle_d` is a generic iso-oriented bounding box in dimension \f$ D\f$
that stores additionally a handle to some underlying geometric object.
It provides in each dimension an interval with lower and upper
endpoints represented with the number type `NT`. This class is
designed to work smoothly with the algorithms for intersecting
sequences of iso-oriented boxes. For degeneracy handling, the boxes
need to provide a unique `id`-number. The policy parameter
`IdPolicy` offers several choices.


\tparam NT number type for the box boundaries, needs to be a model
of the `Assignable` and the `LessThanComparable` concept.
\tparam D the dimension of the box.
\tparam Handle Handle concept, e.g., a pointer, an iterator, or a circulator.
\tparam IdPolicy specifies how the `id`-number will be
provided and can be one of the following types, where
`ID_FROM_HANDLE` is the default for this parameter:
  - `ID_NONE`: no `id`-number is provided. This can be useful
    to have this class as a base class for different
    implementations of `id`-numbers than the ones provided
    here.
  - `ID_EXPLICIT`: the `id`-number is stored explicitly in
    the box and automatically created and assigned at construction
    time of the box. Note that copying a box (copy-constructor and
    assignment) does not create a new `id`-number but keeps
    the old one, which is the behavior needed by the
    `CGAL::box_self_intersection_d()` algorithm. This is therefore
the safe default implementation.
  - `ID_FROM_BOX_ADDRESS`: casts the address of the box into a
    `std::ptrdiff_t` to create the `id`-number. This works fine
    if the intersection algorithms work effectively with pointers
    to boxes, but not in the case where the algorithms work with
    box values, because the algorithms modify the order of the
    boxes, and the `CGAL::box_self_intersection_d()` algorithm
    creates copies of the boxes that would not have identical
    `id`-numbers.
  - `ID_FROM_HANDLE`: casts the address of the value of the
    handle into a `std::ptrdiff_t` to create the
    `id`-number. Works in many conceivable settings, e.g.,
    it works with boxes copied by value or by pointer, and
    the self intersection test. It will not work if there
    is no one-to-one mapping between boxes and the geometry that
    is referred to with the handles, i.e., this `id`-number
    scheme fails if a geometric object creates several boxes with
    the same handle value. Note that this option is not
    available for the `CGAL::Box_intersection_d::Box_d` type
    that does not store a handle.

\cgalModels{BoxIntersectionBox_d}

\sa \link PkgBoxIntersectionD_box_intersection_d `CGAL::box_intersection_d()` \endlink
\sa \link PkgBoxIntersectionD_box_self_intersection_d `CGAL::box_self_intersection_d()` \endlink
\sa \link PkgBoxIntersectionD_box_intersection_all_pairs_d `CGAL::box_intersection_all_pairs_d()` \endlink
\sa \link PkgBoxIntersectionD_box_self_intersection_all_pairs_d `CGAL::box_self_intersection_all_pairs_d()` \endlink
\sa `CGAL::Box_intersection_d::Box_traits_d<BoxHandle>`
\sa `BoxIntersectionTraits_d`

*/
template< typename NT, typename int D, typename Handle, typename IdPolicy >
class Box_with_handle_d {
public:

/// \name Types
/// @{

/*!
number type to represent the box
boundaries. Allowed are the built-in types `int`, `unsigned
int`, `float`, and `double`.
*/
typedef unspecified_type NT;

/*!
type for the box `id`-number.
*/
typedef std::size_t ID;

/// @}

/// \name Creation
/// @{

/*!
%Default constructor. No
particular initialization.
*/
Box_with_handle_d();

/*!
initializes to the
complete or the empty space. If empty, all interval starting (end)
points will be set to positive (negative) infinity, sets handle to \f$ h\f$.
*/
Box_with_handle_d(bool complete, Handle h);

/*!
initializes
the box intervals to [`lo[i]`,`hi[i]`], \f$ 0 \leq i < D\f$ and
sets the handle to \f$ h\f$.
\pre `lo[i]` \f$ <\f$ `hi[i]` for \f$ 0 \leq i < D\f$.
*/
Box_with_handle_d(NT lo[D], NT hi[D], Handle h);

/*!
constructs
from bbox and sets the handle to \f$ h\f$, exists iff \f$ D=2\f$ and `NT`\f$
\equiv\f$`double`.
*/
Box_with_handle_d( const Bbox_2& bbox, Handle h);

/*!
constructs
from bbox and sets the handle to \f$ h\f$, exists iff \f$ D=3\f$ and `NT`\f$
\equiv\f$`double`.
*/
Box_with_handle_d( const Bbox_3& bbox, Handle h);

/// @}

/// \name Modifiers
/// @{

/*!
initializes to the complete or
the empty space. If empty, all interval starting(end) points will be
set to positive(negative) infinity.
*/
void init( bool complete = false);

/*!
extend `box` to contain the
old `box` and `point`.
*/
void extend(NT point[D]);

/// @}

/// \name Access Functions
/// @{

/*!
returns the handle stored in `box`.
*/
Handle handle() const;

/*!
returns \f$ D\f$, the dimension of the box.
*/
static int dimension();

/*!
returns a unique box id, see the
`IdPolicy` template parameter above for the different
choices. Does not exist if `ID_NONE` has been chosen for the
`IdPolicy`.
*/
std::size_t id();

/*!
returns the lower boundary in dimension `d`, \f$ 0 \leq\f$`d`\f$ < D\f$.
*/
NT min_coord( int d) const;

/*!
returns the upper boundary in dimension `d`, \f$ 0 \leq\f$`d`\f$ < D\f$.
*/
NT max_coord( int d) const;

/*!
returns the bounding box iff
\f$ D=2\f$ and `NT`\f$ \equiv\f$`double`.
*/
const Bbox_2& bbox() const;

/*!
returns the bounding box iff
\f$ D=3\f$ and `NT`\f$ \equiv\f$`double`.
*/
const Bbox_3& bbox() const;

/// @}

}; /* end Box_with_handle_d */
} /* end namespace Box_intersection_d */
} /* end namespace CGAL */
