
namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsClasses

The class `Generic_map_min_items` defines `void` as the information associated with darts, and no attribute is enabled.

\cgalModels `GenericMapItems`

\cgalHeading{Example}

The following example shows the implementation of the `Generic_map_min_items` class.

\code{.cpp}
struct Generic_map_min_items
{
  template <class CMap>
  struct Dart_wrapper
  {};
};
\endcode

\sa `Combinatorial_map<d,Items,Alloc>`
\sa `Generalized_map<d,Items,Alloc>`

*/
struct Generic_map_min_items {

}; /* end Generic_map_min_items */
} /* end namespace CGAL */
