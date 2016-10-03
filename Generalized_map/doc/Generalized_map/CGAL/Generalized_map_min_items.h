
namespace CGAL {

/*!
\ingroup PkgGeneralizedMapsClasses

The class `Generalized_map_min_items` is a model of the `BasicMapItems` concept. It defines the type of darts which is a `GMap_dart<d,GMap>`. The `Generalized_map_min_items` has a template argument for the dimension of the generalized map. In this class, no attribute is enabled.

\cgalModels `BasicMapItems`

\tparam d the dimension of the generalized map.

\cgalHeading{Example}

The following example shows the implementation of the `Generalized_map_min_items` class.

\code{.cpp}
template <unsigned int d>
struct Generalized_map_min_items
{
  template <class GMap>
  struct Dart_wrapper
  {
    typedef CGAL::GMap_dart<d, GMap> Dart;
    typedef CGAL::cpp11::tuple<>     Attributes;
  };
};
\endcode

\sa `Generalized_map<d,Items,Alloc>`
\sa `GMap_dart<d,GMap>`

*/
template< typename d >
class Generalized_map_min_items {
public:

/// @}

}; /* end Generalized_map_min_items */
} /* end namespace CGAL */
