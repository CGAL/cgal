namespace CGAL {

namespace Scale_space_reconstruction_3 {

/** \ingroup PkgScaleSpaceReconstruction3Concepts
 * \cgalConcept
 *
 * Concept describing a meshing algorithm used to produce the output
 * mesh of the scale space reconstruction algorithm.
 *
 * A mesher is a functor that can be applied to a range of
 * points and that returns a set of facets.
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Scale_space_reconstruction_3::Alpha_shape_mesher}
 * \cgalHasModels{CGAL::Scale_space_reconstruction_3::Advancing_front_mesher}
 * \cgalHasModelsEnd
 *
 */
class Mesher
{
public:

  /**
   * \tparam InputIterator iterator over input points.
   * \tparam OutputIterator output iterator to which
  `std::array<std::size_t, 3>` can be assigned.
   * \param begin iterator over the first input point.
   * \param end past-the-end iterator over the input points.
   * \param output iterator where facets (triples of point indices)
   * are stored.
   */
  template <typename InputIterator, typename OutputIterator>
  void operator() (InputIterator begin, InputIterator end,
                   OutputIterator output);
};

}


}
