
namespace CGAL {

/*!
\ingroup PkgLinearCellComplexClasses

The class `Linear_cell_complex_min_items` defines `void` as the information associated with darts, and the attributes used. In this class, 0-attributes are enabled and associated with `Cell_attribute_with_point`.

\cgalModels `LinearCellComplexItems`

\deprecated Before CGAL 4.9, this class was templated by the dimension of the darts, and users must define the type of darts used (see also deprecated class `Combinatorial_map_min_items`). `CGAL_CMAP_DART_DEPRECATED` can be defined to keep the old behavior (only possible with `Combinatorial_map` and not for `Generalized_map`).

\cgalHeading{Example}

The following example shows one implementation of the `Linear_cell_complex_min_items` class.

\code{.cpp}

struct Linear_cell_complex_min_items
{
  template <class LCC>
  struct Dart_wrapper
  {
    typedef CGAL::Cell_attribute_with_point<LCC> Vertex_attrib;
    typedef std::tuple<Vertex_attrib> Attributes;
  };
};

\endcode

\sa `CGAL::Linear_cell_complex_for_combinatorial_map<d,d2,LCCTraits,Items,Alloc>`
\sa `CGAL::Linear_cell_complex_for_generalized_map<d,d2,LCCTraits,Items,Alloc>`

*/
struct Linear_cell_complex_min_items {

}; /* end Linear_cell_complex_min_items */

} /* end namespace CGAL */
