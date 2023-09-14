/*!
\ingroup PkgAABBTreeConcepts
\cgalConcept

The concept `AABBRayIntersectionGeomTraits` is a refinement of the
concept `AABBGeomTraits`. In addition to the types required by
`AABBGeomTraits` it also requires types and functors necessary to
define the Intersection_distance functor.

\cgalRefines{AABBGeomTraits}

\cgalHasModelsBegin
\cgalHasModelsBare{All models of the concept `Kernel`}
\cgalHasModelsEnd

\sa `CGAL::AABB_traits<AABBGeomTraits,AABBPrimitive>`
\sa `CGAL::AABB_tree<AABBTraits>`
\sa `AABBPrimitive`

*/
class AABBRayIntersectionGeomTraits {
public:
  /*!
    Type of a 3D ray.
  */
  typedef unspecified_type Ray_3;

  /*!
    Type of a 3D vector.
  */
  typedef unspecified_type Vector_3;

  /*!
    A functor object to construct the source point of a ray. Provides the operator:
    `Point_3 operator()(const Ray_3&);`
   */
  typedef unspecified_type Construct_source_3;

  /*!
   */
  Construct_source_3 construct_source_3_object();

  /*!
    A model of `CartesianConstIterator3`.
   */
  typedef unspecified_type Cartesian_const_iterator_3;

  /*!
    A model of `ConstructCartesianConstIterator3`.
   */
  typedef unspecified_type  Construct_cartesian_const_iterator_3;

  /*!
   */
  Construct_source_3 construct_cartesian_const_iterator_3_object();

  /*!
    A functor object to construct a vector giving the direction of a ray. Provides the operator:
    `Vector_3 operator()(const Ray_3&);`
   */
  typedef unspecified_type Construct_vector_3;

  /*!
   */
  Construct_source_3 construct_vector_3_object();
};
