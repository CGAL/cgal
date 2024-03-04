/*!
\ingroup PkgAABBTreeConcepts
\cgalConcept

The concept `AABBRayIntersectionGeomTraits_3` is a refinement of the
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
class AABBRayIntersectionGeomTraits_3 {
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
    @todo update me
    A model of `CartesianConstIterator3`.
   */
  typedef unspecified_type Cartesian_const_iterator;

  /*!
    @todo update me
    A model of `ConstructCartesianConstIterator3`.
   */
  typedef unspecified_type  Construct_cartesian_const_iterator;

  /*!
   */
  Construct_source construct_cartesian_const_iterator_object();

  /*!
    A functor object to construct a vector giving the direction of a ray. Provides the operator:
    `Vector operator()(const Ray&);`
   */
  typedef unspecified_type Construct_vector;

  /*!
   */
  Construct_source construct_vector_object();
};
