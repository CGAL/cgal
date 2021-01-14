namespace CGAL {
namespace Box_intersection_d {
/*!
\ingroup PkgBoxIntersectionDClasses

`Box_d` is a generic iso-oriented bounding box in dimension \f$ D\f$.
It provides in each dimension an interval with lower and upper
endpoints represented with the number type `NT`. This class is
designed to work smoothly with the algorithms for intersecting
sequences of iso-oriented boxes. For degeneracy handling, the boxes
need to provide a unique `id`-number. The policy parameter
`IdPolicy` offers several choices. The template parameters have to
comply with the following requirements:

\tparam NT is the number type for the box boundaries. It must meet the requierements
of the concepts `Assignable` and `LessThanComparable`.
\tparam D is an integer and the dimension of the box.
\tparam IdPolicy specifies how the `id`-number will be
provided and can be one of the following types, where
`ID_EXPLICIT` is the default for this parameter:
<UL>
<LI>`ID_NONE`: no `id`-number is provided. This can be useful
if `Box_d` is used as a base class for a different
implementation of `id`-numbers than the ones provided
here.
<LI>`ID_EXPLICIT`: the `id`-number is stored explicitly in
the box and automatically created and assigned at construction
time of the box. Note that copying a box (copy-constructor and
assignment) does not create a new `id`-number but keeps
the old one, which is the behavior needed by the
`box_self_intersection_d()` algorithm. This is therefore
the safe default implementation.
<LI>`ID_FROM_BOX_ADDRESS`: casts the address of the box into a
`std::ptrdiff_t` to create the `id`-number. This works fine
if the intersection algorithms work effectively with pointers
to boxes, but not in the case where the algorithms work with
box values, because the algorithms modify the order of the
boxes, and the `box_self_intersection_d()` algorithm
creates copies of the boxes that would not have identical
`id`-numbers.
</UL>

\cgalModels `BoxIntersectionBox_d`

\sa \link PkgBoxIntersectionD_box_intersection_d `CGAL::box_intersection_d()` \endlink
\sa \link PkgBoxIntersectionD_box_self_intersection_d `CGAL::box_self_intersection_d()` \endlink
\sa \link PkgBoxIntersectionD_box_intersection_all_pairs_d `CGAL::box_intersection_all_pairs_d()` \endlink
\sa \link PkgBoxIntersectionD_box_self_intersection_all_pairs_d `CGAL::box_self_intersection_all_pairs_d()` \endlink
\sa `CGAL::Box_intersection_d::Box_with_handle_d<NT, int D, Handle, IdPolicy>`
\sa `CGAL::Box_intersection_d::Box_traits_d<BoxHandle>`
\sa `BoxIntersectionTraits_d`

*/
template< typename NT, typename int D, typename IdPolicy >
class Box_d {
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
%Default constructor. No particular initialization.
*/
Box_d();

/*!
Constructor initialized to the complete or the empty space.
If empty, all interval starting(end) points will be
set to positive(negative) infinity.
*/
Box_d(bool complete);

/*!
initializes the box
intervals to [`lo[i]`,`hi[i]`], \f$ 0 \leq i < D\f$.
\pre `lo[i]` \f$ <\f$ `hi[i]` for \f$ 0 \leq i < D\f$.
*/
Box_d(NT lo[D], NT hi[D]);

/*!
constructs from `bbox`.
Requirements: \f$ D=2\f$ and `NT`\f$ \equiv\f$`double`.
*/
Box_d( const Bbox_2& bbox);

/*!
constructs from `bbox`.
Requirements: \f$ D=3\f$ and `NT`\f$ \equiv\f$`double`.
*/
Box_d( const Bbox_3& bbox);

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

/*!
returns \f$ D\f$, the dimension of the box.
*/
static int dimension();

/// @}

/// \name Access Functions
/// @{

/*!
returns a unique box id, see the
`IdPolicy` template parameter above for the different
choices.
\tparam IdPolicy \f$ \neq\f$`ID_NONE`
*/
std::size_t id();

/*!
returns the lower boundary in dimension `d`
\pre \f$ 0 \leq\f$`d`\f$ < D\f$.
*/
NT min_coord(int d) const;

/*!
returns the upper boundary in dimension `d`
\pre \f$ 0 \leq\f$`d`\f$ < D\f$.
*/
NT max_coord(int d) const;

/*!
returns the bounding box
Requirements: \f$ D=2\f$ and `NT`\f$ \equiv\f$`double`
*/
const Bbox_2& bbox() const;

/*!
returns the bounding box
Requirements: \f$ D=3\f$ and `NT`\f$ \equiv\f$`double`
*/
const Bbox_3& bbox() const;

/*!
extends `box` to the smallest
box that additionally contains the point represented by coordinates in `p`.
*/
void extend(NT p[N]);

/*!
extends `box` to the smallest
box that additionally contains the point represented by coordinate intervals in `p`.
*/
void extend(std::pair<NT,NT> p[N]);

/// @}

}; /* end Box_d */
} /* Box_intersection_d */
} /* end namespace CGAL */
