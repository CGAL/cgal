/*!
\ingroup PkgAABBTreeConcepts
\cgalConcept

The concept `AABBRayIntersectionTraits` is a refinement of the concept
`AABBTraits`. In addition to the types and functions required by
`AABBTraits` it also requires function objects to calculate the
distance of an intersection along a ray.

\cgalRefines{AABBTraits}

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
    Type of bounding box.
  */
  typedef AABBTraits::Bounding_box Bounding_box;

  /*!
    A functor object to construct the source point of a ray. Provides the operator:
    `Point operator()(const Ray&);`
   */
  typedef unspecified_type Construct_source;

  /*!
    returns the `Construct_source` functor.
   */
  Construct_source construct_source_object();

  /*!
    A model of `CartesianConstIterator_2` or `CartesianConstIterator_3`, depending on the dimension of `Vector`.
   */
  typedef unspecified_type Cartesian_const_iterator;

  /*!
    A model of `ConstructCartesianConstIterator_2` or `ConstructCartesianConstIterator_3`, depending on the dimension of  `Vector`.
   */
  typedef unspecified_type  Construct_cartesian_const_iterator;

  /*!
    returns the `Construct_cartesian_const_iterator` functor.
   */
  Construct_cartesian_const_iterator construct_cartesian_const_iterator_object();

  /*!
    A functor object to construct a vector having the same direction as a ray. Provides the operator:
    `Vector operator()(const Ray&);`
   */
  typedef unspecified_type Construct_vector;

  /*!
    returns the `Construct_vector` functor.
   */
  Construct_vector construct_vector_object();


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
