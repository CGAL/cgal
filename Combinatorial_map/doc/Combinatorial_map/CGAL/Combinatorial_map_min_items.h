
namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsClasses

The class `Combinatorial_map_min_items` is a model of the `BasicMapItems` concept. It defines the type of darts which is a `Combinatorial_map_dart<d,CMap>`. The `Combinatorial_map_min_items` has a template argument for the dimension of the combinatorial map. In this class, no attribute is enabled.

\cgalModels `BasicMapItems`

\tparam d the dimension of the combinatorial map.

\cgalHeading{Example}

The following example shows the implementation of the `Combinatorial_map_min_items` class.

\code{.cpp}
template <unsigned int d>
struct Combinatorial_map_min_items
{
  template <class CMap>
  struct Dart_wrapper
  {
    typedef CGAL::Combinatorial_map_dart<d, CMap> Dart;
    typedef CGAL::cpp11::tuple<> Attributes;
  };
};
\endcode

\sa `Combinatorial_map<d,Items,Alloc>`
\sa `Combinatorial_map_dart<d,CMap>`

*/
template< unsigned int d >
class Combinatorial_map_min_items {
public:

/// @}

}; /* end Combinatorial_map_min_items */
} /* end namespace CGAL */
