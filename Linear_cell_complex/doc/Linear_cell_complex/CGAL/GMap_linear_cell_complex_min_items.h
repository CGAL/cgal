
namespace CGAL {

/*!
\ingroup PkgLinearCellComplexClasses

The class `GMap_linear_cell_complex_min_items` defines the type of darts, which is a `GMap_dart<d,LCC>`, and the traits class used. In this class, 0-attributes are enabled and associated with `Cell_attribute_with_point`.

\cgalModels `LinearCellComplexItems`

\tparam d the dimension of the generalized map.

\cgalHeading{Example}

The following example shows one implementation of the `GMap_linear_cell_complex_min_items` class.

\code{.cpp}

template <unsigned int d>
struct GMap_linear_cell_complex_min_items
{
  template <class LCC>
  struct Dart_wrapper
  {
    typedef CGAL::GMap_dart<d, LCC> Dart;
    typedef CGAL::Cell_attribute_with_point<LCC> Vertex_attrib;
    typedef CGAL::cpp11::tuple<Vertex_attrib> Attributes;
  };
};

\endcode

\sa `CGAL::GMap_linear_cell_complex<d,d2,LCCTraits,Items,Alloc>`
\sa `CGAL::GMap_dart<d,CMap>`

*/
template< typename d >
class GMap_linear_cell_complex_min_items {
public:

/// @}

}; /* end GMap_Linear_cell_complex_min_items */
} /* end namespace CGAL */
