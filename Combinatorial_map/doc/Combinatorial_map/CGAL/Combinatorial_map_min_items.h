
namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsClasses

The class `Combinatorial_map_min_items` defines the type of darts which is a `Dart<d,CMap>`. The `Combinatorial_map_min_items` has a template argument for the dimension of the combinatorial map. In this class, no attribute is enabled.

\tparam d the dimension of the combinatorial map.

\deprecated This class is deprecated since CGAL 4.9. Users are required to use class `Generic_map_min_items` instead, where the `Dart` type is no more defined, but replaced by the `Dart_info` type. `CGAL_CMAP_DART_DEPRECATED` can be defined to keep the old behavior.

\cgalHeading{Example}

The following example shows the implementation of the `Combinatorial_map_min_items` class.

\code{.cpp}
template <unsigned int d>
struct Combinatorial_map_min_items
{
  template <class CMap>
  struct Dart_wrapper
  {
    typedef CGAL::Dart<d, CMap> Dart;
    typedef std::tuple<> Attributes;
  };
};
\endcode

\sa `Generic_map_items`

*/
template< unsigned int d >
struct Combinatorial_map_min_items {

}; /* end Combinatorial_map_min_items */
} /* end namespace CGAL */
