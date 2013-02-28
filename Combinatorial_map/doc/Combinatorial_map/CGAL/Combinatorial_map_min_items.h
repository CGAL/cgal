
namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsClasses

The class `Combinatorial_map_min_items` is a model of the `CombinatorialMapItems`
concept. It defines the type of darts which is a
`Dart<d,CMap>`. The `Combinatorial_map_min_items` has a
template argument for the dimension of the combinatorial map.
In this class, no attribute is enabled.

\cgalModels `CombinatorialMapItems`

\cgalHeading{Example}

The following example shows the implementation of the
`Combinatorial_map_min_items` class.

\code{.cpp}
template <unsigned int d>
struct Combinatorial_map_min_items
{
  template <class CMap>
  struct Dart_wrapper
  {
    typedef CGAL::Dart<d, CMap> Dart;
    typedef CGAL::cpp11::tuple<> Attributes;
  };
};
\endcode

\sa `Combinatorial_map<d,Items,Alloc>`
\sa `Dart<d,CMap>`

*/
template< typename d >
class Combinatorial_map_min_items {
public:

/// @}

}; /* end Combinatorial_map_min_items */
} /* end namespace CGAL */
