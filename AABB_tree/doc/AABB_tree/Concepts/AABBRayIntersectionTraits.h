/*!
\ingroup PkgAABBTreeConcepts
\cgalConcept

The concept `AABBRayIntersectionTraits` is a refinement of the concept
`AABBTraits`. In addition to the types and functions required by
`AABBTraits` it also requires function objects to calculate the
distance of an intersection along a ray.

\cgalHasModelsBegin
\cgalHasModels{CGAL::AABB_traits_2<AABBGeomTraits_2,AABBPrimitive>}
\cgalHasModels{CGAL::AABB_traits_3<AABBGeomTraits_3,AABBPrimitive>}
\cgalHasModelsEnd

\sa `CGAL::AABB_tree<AABBTraits>`
\sa `AABBPrimitive`

*/
class AABBRayIntersectionTraits {
public:
  /*!
    Type of a ray.
  */
  typedef unspecified_type Ray;

  /*!
    Type of a vector.
  */
  typedef unspecified_type Vector;

  /*!
    A functor object to construct the source point of a ray. Provides the operator:
    `Point operator()(const Ray&);`
   */
  typedef unspecified_type Construct_source;

  /*!
   */
  Construct_source construct_source_object();

  /*!
    A model of `CartesianConstIterator2` or `CartesianConstIterator3`, must compatible with `Vector`.
   */
  typedef unspecified_type Cartesian_const_iterator;

  /*!
    A model of `ConstructCartesianConstIterator2` or `ConstructCartesianConstIterator3`, must compatible with `Vector`.
   */
  typedef unspecified_type  Construct_cartesian_const_iterator;

  /*!
   */
  Construct_source construct_cartesian_const_iterator_object();

  /*!
    A functor object to construct a vector having the same direction as a ray. Provides the operator:
    `Vector operator()(const Ray&);`
   */
  typedef unspecified_type Construct_vector;

  /*!
   */
  Construct_source construct_vector_object();


  /*!
    A functor object to compute the distance between the source of a ray and its
    closest intersection point between the ray and a primitive or a bounding box.
    An empty `std::optional` is returned, if there is no intersection.
    When there is an intersection, an object of type `FT` is returned such that
    if `i1` and `i2` are two intersection points, then `i1` is closer to the source
    of the ray than `i2` iff `n1 < n2`, `n1` and `n2` being the numbers returned for `i1` and `i2`
    respectively.

    Provides the operators:
    `std::optional<FT> operator()(const Ray& r, const Bounding_box& bbox)`.
    `std::optional<std::pair<FT, Intersection_and_primitive_id<Ray>::%Type > >
     operator()(const Ray& r, const Primitive& primitive)`.

    A common algorithm to compute the intersection between a bounding box and a ray is <A
    HREF="https://education.siggraph.org/static/HyperGraph/raytrace/rtinter3.htm">the
    slab method</A>.
  */
  typedef unspecified_type Intersection_distance;

  /*!
    returns the intersection distance functor.
  */
  Intersection_distance intersection_distance_object() const ;
};
