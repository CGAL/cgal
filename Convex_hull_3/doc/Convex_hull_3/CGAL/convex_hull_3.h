namespace CGAL {

/*!
\ingroup PkgConvexHull3Functions

\brief computes the convex hull of the set of points in the range
[`first`, `last`). The polygon mesh `pm` is cleared, then
the convex hull is stored in `pm`. Note that the convex hull will be triangulated,
that is `pm` will contain only triangular facets.
if the convex hull is a point or a segment, endpoints will be added in `pm` as isolated vertices.

\tparam InputIterator must be an input iterator with a value type  equivalent to `Traits::Point_3`.
\tparam PolygonMesh must be a model of `MutableFaceGraph`.
\tparam Traits must be a model of the concept `ConvexHullTraits_3`.
For the purposes of checking the postcondition that the convex hull
is valid, `Traits` must also be a model of the concept
`IsStronglyConvexTraits_3`.

If the kernel `R` of the points determined by the value type of `InputIterator`
is a kernel with exact predicates but inexact constructions
(in practice we check `R::Has_filtered_predicates_tag` is `Tag_true` and `R::FT` is a floating point type),
then the default traits class of `::convex_hull_3()` is `Convex_hull_traits_3<R>`, and `R` otherwise.

\attention The user must include the header file of the `PolygonMesh` type.

\cgalHeading{Implementation}

The algorithm implemented by these functions is the quickhull algorithm of
Barber <I>et al.</I> \cgalCite{bdh-qach-96}.


*/

template <class InputIterator, class PolygonMesh, class Traits>
void convex_hull_3(InputIterator first, InputIterator last, PolygonMesh& pm, const Traits& ch_traits = Default_traits);


/*!
\ingroup PkgConvexHull3Functions
 * \brief computes the convex hull of the points associated to the vertices of `g`.
 * The polygon mesh `pm` is cleared, then
 * the convex hull is stored in `pm`. Note that the convex hull will be triangulated,
 * that is `pm` will contain only triangular faces.
 * if the convex hull is a point or a segment, endpoints will be added in `pm` as isolated vertices.
 *
 * \tparam VertexListGraph a model of `VertexListGraph`.
 * \tparam PolygonMesh must be a model of `MutableFaceGraph`. an internal property map for `CGAL::vertex_point_t` must be available.
 * \tparam NamedParameters a sequence of named parameters
 *
 * \param g the graph
 * \param pm the `PolygonMesh` that will contain the convex hull
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `g`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<VertexListGraph>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `VertexListGraph`.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 * \attention The user must include the header file of the `PolygonMesh` and `VertexListGraph` types.
 */
template <class VertexListGraph, class PolygonMesh, class NamedParameters = parameters::Default_named_parameters>
void convex_hull_3(const VertexListGraph& g,
                   PolygonMesh& pm,
                   const NamedParameters& np = parameters::default_values());

/*!
\ingroup PkgConvexHull3Functions

\brief computes the convex hull of the set of points in the range
[`first`, `last`). The result, which may be a point, a segment,
a triangle, or a polygon mesh, is stored in `ch_object`.
In the case the result is a polygon mesh, the convex hull will be triangulated,
that is the polygon mesh will contain only triangular facets.

\tparam InputIterator must be an input iterator with a value type  equivalent to `Traits::Point_3`.
\tparam Traits must be model of the concept `ConvexHullTraits_3`.
For the purposes of checking the postcondition that the convex hull
is valid, `Traits` must also be a model of the concept
`IsStronglyConvexTraits_3`.   Furthermore, `Traits` must define a type `PolygonMesh` that is a model of
`MutableFaceGraph`.

If the kernel `R` of the points determined by the value type  of `InputIterator`
is a kernel with exact predicates but inexact constructions
(in practice we check `R::Has_filtered_predicates_tag` is `Tag_true` and `R::FT` is a floating point type),
then the default traits class of `convex_hull_3()` is `Convex_hull_traits_3<R>`, and `R` otherwise.

\attention The user must include the header file of the `PolygonMesh` type.
*/

template <class InputIterator, class Traits>
void convex_hull_3(InputIterator first, InputIterator last,
Object& ch_object,
const Traits& ch_traits = Default_traits);

/*!
\ingroup PkgConvexHull3Functions

\brief computes the convex hull of the set of points in the range
[`first`, `last`). The result, which may be a point, a segment,
a triangle, or a triangle mesh, is stored as an indexed triangle soup
in a vector of points and a vector of index triples.


\tparam InputIterator must be an input iterator with a value type  equivalent to `Traits::Point_3`.
\tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type.
\tparam PolygonRange a model of the concepts `RandomAccessContainer` and `BackInsertionSequence`
                     whose `value_type` is itself a model of the concepts `RandomAccessContainer`
                     and `BackInsertionSequence` whose `value_type` is an unsigned integer type
                     convertible to `std::size_t`
\tparam Traits must be model of the concept `ConvexHullTraits_3`.
        For the purposes of checking the postcondition that the convex hull
        is valid, `Traits` must also be a model of the concept
        `IsStronglyConvexTraits_3`.

If the kernel `R` of the points determined by the value type  of `InputIterator`
is a kernel with exact predicates but inexact constructions
(in practice we check `R::Has_filtered_predicates_tag` is `Tag_true` and `R::FT` is a floating point type),
then the default traits class of `convex_hull_3()` is `Convex_hull_traits_3<R>`, and `R` otherwise.


*/

  template <class InputIterator, class PointRange, class PolygonRange class Traits>
void convex_hull_3(InputIterator first, InputIterator last,
                   PointRange& vertices,
                   PolygonRange& faces,
                   const Traits& ch_traits = Default_traits);

/*!
\ingroup PkgConvexHull3Functions

\brief copies in `out` the points on the convex hull of the points in `range`.

\tparam InputRange a range of `Traits::Point_3`, model of `ConstRange`.
        Its iterator type is `InputIterator`.
\tparam OutputIterator must be an output iterator where points of type `Traits::Point_3` can be put.
\tparam Traits must be model of the concept `ConvexHullTraits_3`.

If the kernel `R` of the points from `InputRange`
is a kernel with exact predicates but inexact constructions
(in practice we check `R::Has_filtered_predicates_tag` is `Tag_true` and `R::FT` is a floating point type),
then the default traits class used is `Convex_hull_traits_3<R>`, and `R` otherwise.

@param range the range of input points.
@param out an output iterator where the extreme points will be put.
@param traits an instance of `Traits`.

\sa `CGAL::Extreme_points_traits_adapter_3`
*/
template <class InputRange, class OutputIterator, class Traits>
OutputIterator
extreme_points_3(InputRange range,
              OutputIterator out,
              const Traits& traits);



} /* namespace CGAL */
