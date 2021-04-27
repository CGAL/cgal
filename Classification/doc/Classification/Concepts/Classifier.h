namespace CGAL
{

namespace Classification
{

/*!
\ingroup PkgClassificationConcepts
\cgalConcept

Concept describing a classifier used by classification functions (see
`CGAL::Classification::classify()`, `CGAL::Classification::classify_with_local_smoothing()` and
`CGAL::Classification::classify_with_graphcut()`).

\cgalHasModel `CGAL::Classification::Sum_of_weighted_features_classifier`
\cgalHasModel `CGAL::Classification::ETHZ::Random_forest_classifier`
\cgalHasModel `CGAL::Classification::OpenCV::Random_forest_classifier`

*/
class Classifier
{
public:

  /*!
    \brief Computes for each label indexed from 0 to `out.size()`, the
    probability (between 0 and 1) that the item at `item_index`
    belongs to this label.
   */
  void operator() (std::size_t item_index, std::vector<float>& out) const;


};

} // namespace Classification

} // namespace CGAL
