
namespace CGAL {

/*!
\ingroup BGLGraphExternalIndices

The class `Linear_cell_complex_bgl_min_items` defines `void` as the information associated with darts, darts have ids and 0- and 2-attributes are enabled and have ids.

\cgalModels `LinearCellComplexItems`

\cgalHeading{Example}

The following example shows one implementation of the `Linear_cell_complex_bgl_min_items` class.

\code{.cpp}

struct Linear_cell_complex_bgl_min_items
{
  template <class LCC>
  struct Dart_wrapper
  {
    typedef CGAL::Tag_true Darts_with_id;
    typedef CGAL::Cell_attribute_with_point_and_id<LCC> Vertex_attrib;
    typedef CGAL::Cell_attribute_with_id<LCC> Face_attrib;
    typedef std::tuple<Vertex_attrib, void, Face_attrib> Attributes;
  };
};

\endcode

\sa `CGAL::Linear_cell_complex_min_item`

*/
struct Linear_cell_complex_bgl_min_items {

}; /* end Linear_cell_complex_bgl_min_items */

} /* end namespace CGAL */
