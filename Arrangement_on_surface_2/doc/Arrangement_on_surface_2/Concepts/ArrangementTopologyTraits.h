/*! \ingroup PkgArrangementOnSurface2ConceptsTopologyTraits
 * \cgalConcept
 *
 * A geometry traits class encapsulates the definitions of the geometric
 * entities and implements the geometric predicates and constructions needed by
 * instances of the `CGAL::Arrangement_on_surface_2` class template and by the
 * peripheral algorithms that operate on objects of such instances. Essentially,
 * it maintains the doubly-connected connected edge list (DCEL) used by the
 * arrangement.
 *
 * The package contains three topology traits, as follows:
 * <UL>
 * <li> `CGAL::Arr_spherical_topology_traits_2`&mdash;can serve as a topology traits
 * for an arrangement of planar bounded curves.
 * <li> `CGAL::Arr_bounded_planar_topology_traits_2`&mdash;can serve as a topology traits
 * for an arrangement of planar unbounded curves.
 * <li> `CGAL::Arr_unb_planar_topology_traits_2`&mdash;can serve as a topology traits
 * for an arrangement of arcs of great circles embedded on a sphere.
 * </ul>
 *
 * At this point we do not expose all the requirements of this concept.
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Arr_bounded_planar_topology_traits_2<GeometryTraits_2,Dcel>}
 * \cgalHasModels{CGAL::Arr_unb_planar_topology_traits_2<GeometryTraits_2,Dcel>}
 * \cgalHasModels{CGAL::Arr_spherical_topology_traits_2<GeometryTraits_2,Dcel>}
 * \cgalHasModelsEnd
 *
 * \sa `Arrangement_on_surface_2<GeometryTraits_2,TopologyTraits>`
 */

class ArrangementTopologyTraits {
public:

  /// \name Types
  /// @{

  /*! */
  typedef unspecified_type Geometry_traits_2;

  /*! */
  typedef unspecified_type Dcel;
  /// @}

  /// \name Creation
  /// @{

  /*! constructs default. */
  Arr_topology_traits();

  /*! constructs from a geometry-traits object. */
  Arr_topology_traits(const Geometry_traits_2* geometry_traits);

  /// @}

  /// \name Access Functions
  /// @{

  /*! obtains the (const) DCEL. */
  const Dcel& dcel() const;

  /*! obtains the (non-const) DCEL. */
  Dcel& dcel();

  /// @}
};
