
namespace CGAL {

/*!
\ingroup PkgLinearCellComplexClasses

The class `Linear_cell_complex_min_items` defines the type of darts, which is a
\ref CombinatorialMapItems::Dart_wrapper "Dart_wrapper::Dart<d,LCC>", and the traits class used. In
this class, 0-attributes are enabled and associated with
`Cell_attribute_with_point`.

\cgalModels `LinearCellComplexItems`

\tparam d the dimension of the combinatorial map.

\cgalHeading{Example}

The following example shows one implementation of the
`Linear_cell_complex_min_items` class.

\code{.cpp}

template <unsigned int d>
struct Linear_cell_complex_min_items
{
  template <class LCC>
  struct Dart_wrapper
  {
    typedef CGAL::Dart<d, LCC> Dart;
    typedef CGAL::Cell_attribute_with_point<LCC> Vertex_attrib;
    typedef CGAL::cpp11::tuple<Vertex_attrib> Attributes;
  };
};

\endcode

\sa `CGAL::Linear_cell_complex<d,d2,LCCTraits,Items,Alloc>`
\sa `CGAL::Dart<d,CMap>`

*/
template< typename d >
class Linear_cell_complex_min_items {
public:

/// @}

}; /* end Linear_cell_complex_min_items */
} /* end namespace CGAL */
