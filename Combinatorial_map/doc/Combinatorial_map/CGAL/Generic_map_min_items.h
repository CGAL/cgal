
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
  {
    typedef void Dart_info;
    typedef CGAL::cpp11::tuple<> Attributes;
  };
};
\endcode

\sa `Combinatorial_map<d,Items,Alloc>`
\sa `Generalized_map<d,Items,Alloc>`

*/
class Generic_map_min_items {
public:

/// @}

}; /* end Combinatorial_map_min_items */
} /* end namespace CGAL */
