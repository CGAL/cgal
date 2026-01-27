namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2Funcs
 *
 * produces the symbolic vertical decomposition of a given arrangement. More
 * precisely, this function performs a batched vertical ray-shooting query from
 * every arrangement vertex, and pairs each vertex with a pair of polymorphic
 * objects, one corresponds to the arrangement feature that lies below it, and
 * the other corresponds to the feature that lies above it.
 *
 * The finite arrangement vertices and the features they "see", if exist,
 * that are, the query results, are inserted in ascending \f$xy\f$-lexicographic
 * order (of the query vertex) into an output container given through an output
 * iterator. If the vertex is the top end-vertex of a vertical edge, we say that
 * there is no feature below it; similarly, if it is the bottom end-vertex of a
 * vertical edge, we say that there is no feature above it.  Each feature, if
 * exists, is represented by a discriminated union container that holds an
 * object of one of the following types:
 *
 * <UL>
 * <LI> `Arrangement_on_surface_2::Halfedge_const_handle`, if the vertex is
 *      located above (or below) an edge. The given halfedge is always directed
 *      from right to left.  In case there is no concrete edge below (or above)
 *      the vertex, and the arrangement is unbounded, then the object returned
 *      is a <I>fictitious</I> halfedge.
 * <LI> `Arrangement_on_surface_2::Face_const_handle`, in case there is no edge
 *      below (or above) the vertex, and the arrangement is bounded.
 * <LI> `Arrangement_on_surface_2::Vertex_const_handle`, in case the vertex is
 *      located vertically above (or below) another arrangement vertex.
 * </UL>
 *
 * The output of this function can be readily used for inserting vertical walls
 * and physically decomposing the arrangement into pseudo-trapezoids.
 *
 * \param arr The arrangement.
 * \param oi The output iterator that points at the output container.
 * \return The past-the-end iterator of the output container.
 *
 * \cgalHeading{Requirements}
 *
 * \pre Dereferencing `oi` must yield an object of type
 * `std::pair<Arrangement_on_surface_2::Vertex_const_handle,
 *            std::pair<std::optional<Type,std::optional<Type>>>`,
 * where `Type` is
 * `std::variant<Arrangement_on_surface_2::Vertex_const_handle, Arrangement_on_surface_2::Halfedge_const_handle, Arrangement_on_surface_2::Face_const_handle>`.
 */
template <typename Traits, typename TopologyTraits, typename OutputIterator>
OutputIterator
decompose(const Arrangement_on_surface_2<GeometryTraits,TopologyTraits>& arr,
          OutputIterator oi);

} /* namespace CGAL */
