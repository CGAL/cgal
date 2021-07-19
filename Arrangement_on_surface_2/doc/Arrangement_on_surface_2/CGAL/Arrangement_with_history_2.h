
namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2Ref

\anchor arr_refarr_with_hist

An object `arr` of the class `Arrangement_with_history_2` represents the
planar subdivision induced by a set of input curves \f$ \cal C\f$.
The arrangement is represented as a doubly-connected edge-list (\dcel).
As is the case for the `Arrangement_2<Traits,Dcel>`, each \dcel
vertex is associated with a point and each edge is associated with an
\f$ x\f$-monotone curve whose interior is disjoint from all other edges and
vertices. Each such \f$ x\f$-monotone curve is a subcurve of some
\f$ C \in \cal C\f$ - or may represent an overlap among several curves
in \f$ \cal C\f$.

The `Arrangement_with_history_2` class-template extends the `Arrangement_2`
class-template by keeping an additional container of input curves
representing \f$ \cal C\f$, and by maintaining a cross-mapping between these
curves and the arrangement edges they induce. This way it is possible
to determine the inducing curve(s) of each arrangement edge. This mapping
also allows the traversal of input curves, and the traversal of edges
induced by each curve.

The `Arrangement_with_history_2` template has two parameters:
<UL>
<LI>The `Traits` template-parameter should be instantiated with
a model of the `ArrangementTraits_2` concept. The traits
class defines the `Curve_2` type, which represents an input curve.
It also defines the types of \f$ x\f$-monotone curves and two-dimensional
points, namely `ArrangementTraits_2::X_monotone_curve_2` and `ArrangementTraits_2::Point_2`,
respectively, and supports basic geometric predicates on them.
<LI>The `Dcel` template-parameter should be instantiated with
a class that is a model of the `ArrangementDcelWithRebind` concept. The
value of this parameter is by default
`Arr_default_dcel<Traits>`.
</UL>

\sa `ArrangementDcel`
\sa `Arr_default_dcel<Traits>`
\sa `ArrangementTraits_2`
\sa `Arrangement_2<Traits,Dcel>`
\sa `insertion functions`
\sa `removal functions`
\sa `overlaying arrangements`

*/
template< typename Traits, typename Dcel >
class Arrangement_with_history_2 : public Arrangement_2<Traits,Dcel> {
public:

/// \name Types

/// @{

/*!
a private type used as an abbreviation of the `Arrangement_with_history_2` type hereafter.
*/
typedef Arrangement_with_history_2<Traits_2,Dcel> Self;

/*!
the traits class in use.
*/
typedef unspecified_type Traits_2;

/*!
the \dcel representation of the arrangement.
*/
typedef unspecified_type Dcel;

/*!
the point type, as defined by the traits class.
*/
typedef typename Traits_2::Point_2 Point_2;

/*!
the \f$ x\f$-monotone curve type, as defined by the traits class.
*/
typedef typename Traits_2::X_monotone_curve_2 X_monotone_curve_2;

/*!
the curve type, as defined by the traits class.
*/
typedef typename Traits_2::Curve_2 Curve_2;


/// @}

/*! \name

In addition, the nested types `Vertex`, `Halfedge` and `Face` are defined, as well as all handle, iterator and circulator types, as defined by the `Arrangement_2` class-template.

*/
/// @{


/*!
a handle for an input curve.
*/
typedef unspecified_type Curve_handle;

/*!
a bidirectional iterator over the
curves that induce the arrangement. Its value-type is
`Curve_2`.
*/
typedef unspecified_type Curve_iterator;

/*!
an iterator over the edges induced by an input curve.
Its value type is `Halfedge_handle`.
*/
typedef unspecified_type Induced_edge_iterator;

/*!
an iterator for the curves that originate a given arrangement edge.
Its value type is `Curve_handle`.
*/
typedef unspecified_type Originating_curve_iterator;

/// @}

/// \name Creation
/// @{

/*!
constructs an empty arrangement containing one unbounded face,
which corresponds to the
whole plane.
*/
Arrangement_with_history_2<Traits, Dcel>();

/*!
copy constructor.
*/
Arrangement_with_history_2<Traits, Dcel>(const Self& other);

/*!
constructs an empty arrangement that uses the given `traits`
instance for performing the geometric predicates.
*/
Arrangement_with_history_2<Traits, Dcel>(Traits_2 *traits);

/// @}

/// \name Assignment Methods
/// @{

/*!
assignment operator.
*/
Self& operator= (other);

/*!
assigns the contents of another arrangement.
*/
void assign (const Self& other);

/*!
clears the arrangement.
*/
void clear ();

/// @}

/*! \name Access Functions for Input Curves

See the `Arrangement_2` reference pages for the full list.
*/

/// @{

/*!
returns the number of input curves that induce the arrangement.
*/
Size number_of_curves() const;

/*!
returns the begin-iterator of the curves inducing the arrangement.
*/
Curve_iterator curves_begin();

/*!
returns the past-the-end iterator of the curves inducing the arrangement.
*/
Curve_iterator curves_end();

/*!
returns the number of arrangement edges induced by the curve `ch`.
*/
Size number_of_induced_edges (Curve_handle ch) const;

/*!
returns the begin-iterator of the edges induced by the curve `ch`.
*/
Induced_edge_iterator
induced_edges_begin (Curve_handle ch) const;

/*!
returns the past-the-end iterator of the edges induced by the curve `ch`.
*/
Induced_edge_iterator
induced_edges_end (Curve_handle ch) const;

/*!
returns the number of input curves that originate the edge `e`.
*/
Size number_of_originating_curves (Halfedge_handle e) const;

/*!
returns the begin-iterator of the curves originating the edge `e`.
*/
Originating_curve_iterator
originating_curves_begin (Halfedge_handle e) const;

/*!
returns the past-the-end iterator of the curves originating the edge
`e`.
*/
Originating_curve_iterator
originating_curves_end (Halfedge_handle e) const;

/// @}

/*! \name Modifying Arrangement Edges

The following functions override their counterparts in the
`Arrangement_2` class, as they also maintain the cross-relationships
between the input curves and the edges they induce.

See the `Arrangement_2` reference pages for the full list of functions
for modifying arrangement vertices
*/

/// @{

/*!
splits the edge `e` into two edges (more precisely, into two halfedge
pairs), at a given split point `p`.
The function returns a handle for the halfedge whose source is the same
as `e->source()` and whose target vertex is the split point.
\pre `p` lies in the interior of the curve associated with `e`.
*/
Halfedge_handle split_edge (Halfedge_handle e,
const Point_2& p);

/*!
merges the edges represented by `e1` and `e2` into
a single edge.
The function returns a handle for one of the merged halfedges.
\pre `e1` and `e2` share a common end-vertex, of degree \f$ 2\f$, and the \f$ x\f$-monotone curves associated with `e1` and `e2` are mergeable into a single \f$ x\f$-monotone curves.
*/
Halfedge_handle merge_edge (Halfedge_handle e1,
Halfedge_handle e2);

/*!
removes the edge `e` from the arrangement. Since the `e` may
be the only edge incident to its source vertex (or its target vertex),
this vertex can be removed as well. The flags `remove_source` and
`remove_target` indicate whether the endpoints of `e` should be
removed, or whether they should be left as isolated vertices in the
arrangement.
If the operation causes two faces to merge, the merged face is returned.
Otherwise, the face to which the edge was incident is returned.
*/
Face_handle remove_edge(Halfedge_handle e,
bool remove_source = true,
bool remove_target = true);

/// @}

}; /* end Arrangement_with_history_2 */

/*!
\ingroup PkgArrangementOnSurface2Insert

Inserts the given curve `c` into the arrangement with history `arr`,
and returns a handle to the inserted curve. `c` is subdivided into
\f$ x\f$-monotone subcurves (and perhaps isolated points). Each subcurve is in
turn inserted into the arrangement by locating its left endpoint and
computing its zone until reaching the right endpoint.

The given point-location object `pl` is used to locate the left
endpoints of the \f$ x\f$-monotone curves. By default, the function uses the
"walk along line" point-location strategy  -  namely an instance of
the class `Arr_walk_along_line_point_location<Arrangement_2<Traits,Dcel> >`.

\pre If provided, `pl` is attached to the given arrangement `arr`.
*/
template<typename Traits, typename Dcel,
         typename PointLocation>
typename Arrangement_with_history_2<Traits,Dcel>::Curve_handle
insert (Arrangement_with_history_2<Traits,Dcel>& arr,
        const typename Traits::Curve_2& c,
        const PointLocation& pl = walk_pl);

/*!
\ingroup PkgArrangementOnSurface2Insert
Aggregately inserts the curves in the range `[first,last)` into the
arrangement with history `arr` using the sweep-line framework.

*/
template<class Traits, class Dcel, InputIterator>
void insert(Arrangement_with_history_2<Traits,Dcel>& arr,
            InputIterator first, InputIterator last);


/*!
\ingroup PkgArrangementOnSurface2Funcs

Removes a given curve from a given arrangement.

The curve is specified by its handle `ch`, from
the arrangement `arr`, by deleting all the edges it induces. The
function returns the number of deleted edges.

*/
template <class Traits, class Dcel>
Size remove_curve (Arrangement_with_history_2<Traits,Dcel>& arr,
                   typename Arrangement_with_history_2<Traits,Dcel>::Curve_handle ch);


/*!
\addtogroup PkgArrangementOnSurface2Overlay

Computes the overlay of two arrangements with history `arr1` and
`arr2`, and sets the output arrangement with history `res` to
represent the overlaid arrangement. The function also constructs a
consolidated set of curves that induce `res`.
\pre `res` does not refer to either `arr1` or `arr2` (that is, "self overlay" is not supported).
*/
template<typename Traits, typename Dcel1, typename Dcel2,
         typename ResDcel, typename OverlayTraits>
void overlay (const Arrangement_with_history_2<Traits,Dcel1>& arr1,
              const Arrangement_with_history_2<Traits,Dcel2>& arr2,
              Arrangement_with_history_2<Traits,ResDcel>& res,
              OverlayTraits& ovl_tr);


/*!
\addtogroup PkgArrangementOnSurface2Overlay

Computes the (simple) overlay of two arrangements with history `arr1`
and `arr2`, and sets the output arrangement with history `res` to
represent the overlaid arrangement. The function also constructs a
consolidated set of curves that induce `res`. It employs the default
overlay-traits, which practically does nothing.
\pre `res` does not refer to either `arr1` or `arr2` (that is, "self overlay" is not supported).
*/
template<typename Traits, typename Dcel1, typename Dcel2,
                           typename ResDcel>
           void overlay (const Arrangement_with_history_2<Traits,Dcel1>& arr1,
                         const Arrangement_with_history_2<Traits,Dcel2>& arr2,
                         Arrangement_with_history_2<Traits,ResDcel>& res);

} /* namespace CGAL */

