namespace CGAL
{

namespace Classification
{

/*!
\ingroup PkgClassificationConcepts
\cgalConcept

Concept describing a neighbor query used for classification.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Classification::Point_set_neighborhood::K_neighbor_query}
\cgalHasModels{CGAL::Classification::Point_set_neighborhood::Sphere_neighbor_query}
\cgalHasModels{CGAL::Classification::Mesh_neighborhood::One_ring_neighbor_query}
\cgalHasModels{CGAL::Classification::Mesh_neighborhood::N_ring_neighbor_query}
\cgalHasModelsEnd

*/
class NeighborQuery
{
public:

  /*!
    \brief Type of the data that is classified.
  */
  typedef unspecified_type value_type;

  /*!

    \brief puts in `output` the indices of the neighbors of `query`.

    \tparam OutputIterator An output iterator accepting `std::size_t`
    values.
   */
  template <typename OutputIterator>
  OutputIterator operator() (const value_type& query, OutputIterator output) const;

};

} // namespace Classification

} // namespace CGAL
