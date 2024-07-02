/*!
\ingroup PkgAABBTreeConcepts
\cgalConcept

The concept `AABBRayIntersectionGeomTraits_3` is a refinement of the
concept `AABBGeomTraits_3`. In addition to the types required by
`AABBGeomTraits_3` it also requires types and functors necessary to
define the Intersection_distance functor.

\cgalRefines{AABBGeomTraits_3}

\cgalHasModelsBegin
\cgalHasModelsBare{All models of the concept `Kernel`}
\cgalHasModelsEnd

\sa `CGAL::AABB_traits_3<AABBGeomTraits_3,AABBPrimitive>`
\sa `CGAL::AABB_tree<AABBTraits>`
\sa `AABBPrimitive`

*/
class AABBRayIntersectionGeomTraits_3 {
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
    returns the `Construct_source_3` functor.
   */
  Construct_source_3 construct_source_3_object();

  /*!
    A model of `CartesianConstIterator_3`.
   */
  typedef unspecified_type Cartesian_const_iterator_3;

  /*!
    A model of `ConstructCartesianConstIterator_3`.
   */
  typedef unspecified_type  Construct_cartesian_const_iterator_3;

  /*!
    returns the `Construct_cartesian_const_iterator_3` functor.
   */
  Construct_cartesian_const_iterator_3 construct_cartesian_const_iterator_3_object();

  /*!
    A functor object to construct a vector having the same direction as a ray. Provides the operator:
    `Vector_3 operator()(const Ray_3&);`
   */
  typedef unspecified_type Construct_vector_3;

  /*!
    returns the `Construct_vector_3` functor.
   */
  Construct_vector_3 construct_vector_3_object();
};
