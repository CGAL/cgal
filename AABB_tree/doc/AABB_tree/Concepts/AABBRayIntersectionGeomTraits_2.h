3
/*!
\ingroup PkgAABBTreeConcepts
\cgalConcept

The concept `AABBRayIntersectionGeomTraits_2` is a refinement of the
concept `AABBGeomTraits_2`. In addition to the types required by
`AABBGeomTraits_2` it also requires types and functors necessary to
define the Intersection_distance functor.

\cgalRefines{AABBGeomTraits}

\cgalHasModelsBegin
\cgalHasModelsBare{All models of the concept `Kernel`}
\cgalHasModelsEnd

\sa `CGAL::AABB_traits_2<AABBGeomTraits_2,AABBPrimitive>`
\sa `CGAL::AABB_tree<AABBTraits>`
\sa `AABBPrimitive`

*/
class AABBRayIntersectionGeomTraits_2 {
public:
  /*!
    Type of a 2D ray.
  */
  typedef unspecified_type Ray_2;

  /*!
    Type of a 2D vector.
  */
  typedef unspecified_type Vector_2;

  /*!
    A functor object to construct the source point of a ray. Provides the operator:
    `Point_2 operator()(const Ray_2&);`
   */
  typedef unspecified_type Construct_source_2;

  /*!
   */
  Construct_source_2 construct_source_2_object();

  /*!
    A model of `CartesianConstIterator2`.
   */
  typedef unspecified_type Cartesian_const_iterator_2;

  /*!
    A model of `ConstructCartesianConstIterator2`.
   */
  typedef unspecified_type  Construct_cartesian_const_iterator_2;

  /*!
   */
  Construct_source_2 construct_cartesian_const_iterator_2_object();

  /*!
    A functor object to construct a vector giving the direction of a ray. Provides the operator:
    `Vector_2 operator()(const Ray_2&);`
   */
  typedef unspecified_type Construct_vector_2;

  /*!
   */
  Construct_source_2 construct_vector_2_object();
};