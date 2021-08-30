/*!
\ingroup PkgAABBTreeConcepts
\cgalConcept

The concept `AABBRayIntersectionTraits` is a refinement of the concept
`AABBTraits`. In addition to the types and functions required by
`AABBTraits` it also requires function objects to calculate the
distance of an intersection along a ray.

\cgalHasModel `CGAL::AABB_traits<AABBGeomTraits,AABBPrimitive>`

\sa `CGAL::AABB_traits<AABBGeomTraits,AABBPrimitive>`
\sa `CGAL::AABB_tree<AABBTraits>`
\sa `AABBPrimitive`

*/
class AABBRayIntersectionTraits {
public:

  /*!
    Type of a 3D ray.
  */
  typedef unspecified_type Ray_3;


  /*!
    A functor object to compute the distance between the source of a ray and its
    closest intersection point between the ray and a primitive or a bounding box.
    An empty `boost::optional` is returned, if there is no intersection.
    When there is an intersection, an object of type `FT` is returned such that
    if `i1` and `i2` are two intersection points, then `i1` is closer to the source
    of the ray than `i2` iff `n1 < n2`, `n1` and `n2` being the numbers returned for `i1` and `i2`
    respectively.

    Provides the operators:
    `boost::optional<FT> operator()(const Ray_3& r, const Bounding_box& bbox)`.
    `boost::optional<std::pair<FT, Intersection_and_primitive_id<Ray_3>::%Type > >
     operator()(const Ray_3& r, const Primitive& primitive)`.

    A common algorithm to compute the intersection between a bounding box and a ray is <A
    HREF="http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter3.htm">the
    slab method</A>.
  */
  typedef unspecified_type Intersection_distance;

  /*!
    returns the intersection distance functor.
  */
  Intersection_distance intersection_distance_object() const ;
};
