
namespace CGAL {

/*!
\ingroup PkgCombinatorialMaps

The class `Combinatorial_map_min_items` is a model of the `CombinatorialMapItems` 
concept. It defines the type of darts which is a 
`CGAL::Dart<d,CMap>`. The `Combinatorial_map_min_items` has a 
template argument for the dimension of the combinatorial map. 
In this class, no attribute is enabled. 

\models ::CombinatorialMapItems 

### Example ###

The following example shows the implementation of the 
`CGAL::Combinatorial_map_min_items` class. 

\code{.cpp} 
template <unsigned int d> 
struct Combinatorial_map_min_items 
{ 
  template <class CMap> 
  struct Dart_wrapper 
  { 
    typedef CGAL::Dart<d, CMap> Dart; 
    typedef CGAL::cpp0x::tuple<> Attributes; 
  }; 
}; 
\endcode 

\sa `CGAL::Combinatorial_map<d,Items,Alloc>`
\sa `CGAL::Dart<d,CMap>`

*/
template< typename d >
class Combinatorial_map_min_items {
public:

/// @}

}; /* end Combinatorial_map_min_items */
} /* end namespace CGAL */
