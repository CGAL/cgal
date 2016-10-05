
namespace CGAL {

/*!
\ingroup PkgLinearCellComplexClasses

The class `Linear_cell_complex_for_generalized_map_min_items` defines the type of darts, which is a `Generalized_map_dart<d,LCC>`, and the attributes used. In this class, 0-attributes are enabled and associated with `Cell_attribute_with_point`.

\cgalModels `LinearCellComplexItems`

\tparam d the dimension of the generalized map.

\cgalHeading{Example}

The following example shows one implementation of the `Linear_cell_complex_for_generalized_map_min_items` class.

\code{.cpp}

template <unsigned int d>
struct Linear_cell_complex_for_generalized_map_min_items
{
  template <class LCC>
  struct Dart_wrapper
  {
    typedef CGAL::Generalized_map_dart<d, LCC> Dart;
    typedef CGAL::Cell_attribute_with_point<LCC> Vertex_attrib;
    typedef CGAL::cpp11::tuple<Vertex_attrib> Attributes;
  };
};

\endcode

\sa `CGAL::Linear_cell_complex_for_generalized_map<d,d2,LCCTraits,Items,Alloc>`
\sa `CGAL::Generalized_map_dart<d,CMap>`

*/
template< typename d >
class Linear_cell_complex_for_generalized_map_min_items {
public:

/// @}

}; /* end Linear_cell_complex_for_generalized_map_min_items */
} /* end namespace CGAL */
