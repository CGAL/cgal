/*!
\ingroup PkgAABB_treeConcepts
\cgalConcept

The concept `AABBRayIntersectionTraits` is a refinement of the Model
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

    A functor object to compute the smallest ray parameter, if any, at
    which a bounding box and a ray intersect. Provides the operator:
    `boost::optional<FT> operator()(const Ray_3& r, const Bounding_box& bbox)`.

    A common algorithm for this is <A
    HREF="http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter3.htm">the
    slab method</A>.
  */
  typedef unspecified_type Intersection_distance;

  /*!
    Returns the intersection distance functor. 
  */ 
  Intersection_distance intersection_distance_object() const ;
};
