namespace CGAL {

namespace Scale_space_reconstruction_3 {

/** \ingroup PkgScaleSpaceReconstruction3Concepts
 * \cgalConcept
 *
 * Concept describing a smoothing algorithm used to construct the
 * scales of the scale space reconstruction algorithm.
 *
 * A smoother is a functor that can be applied to a range of
 * points. The container is directly modified, points are replaced by
 * their smoothed versions.
 *
 * \note the functor can be applied several times by the scale space
 * reconstruction algorithm.
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Scale_space_reconstruction_3::Weighted_PCA_smoother}
 * \cgalHasModels{CGAL::Scale_space_reconstruction_3::Jet_smoother}
 * \cgalHasModelsEnd
 */
class Smoother
{
public:

  /**
   * \tparam InputIterator iterator over input points.
   * \param begin iterator over the first input point.
   * \param end past-the-end iterator over the input points.
   */
  template <typename InputIterator>
  void operator() (InputIterator begin, InputIterator end);
};

}


}
