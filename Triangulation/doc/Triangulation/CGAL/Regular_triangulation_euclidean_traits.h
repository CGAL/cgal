
namespace CGAL {

  /*!
  \ingroup PkgTriangulationsTriangulationClasses

  The class `Regular_triangulation_euclidean_traits` is designed as a traits 
  class for the class
  `Regular_triangulation<RegularTriangulationTraits, TriangulationDataStructure>`.

  \tparam K must be a model of the `Kernel_d` concept. We recommend to use 
            `Epick_d`.

  \tparam Weight is optional. If is it not provided, `K::RT` will be used.

  \cgalModels `RegularTriangulationTraits`
  */
template < class K, class Weight = typename K::RT >
class Regular_triangulation_euclidean_traits
  : public K
{
};

} //namespace CGAL

