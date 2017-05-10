
namespace CGAL {

/*!
\ingroup PkgBGLHelper

\cgalModifBegin

The class `Linear_cell_complex_bgl_min_items` defines `void` as the information associated with darts, and the attributes used. In this class, 0- and 2-attributes are enabled and have ids.

\cgalModels `LinearCellComplexItems`

\cgalHeading{Example}

The following example shows one implementation of the `Linear_cell_complex_bgl_min_items` class.

\code{.cpp}

struct Linear_cell_complex_bgl_min_items
{
  template <class LCC>
  struct Dart_wrapper
  {
    typedef CGAL::Cell_attribute_with_point_and_id<LCC> Vertex_attrib;
    typedef CGAL::Cell_attribute_with_id<LCC> Face_attrib;
    typedef CGAL::cpp11::tuple<Vertex_attrib, void, Face_attrib> Attributes;
  };
};

\endcode

\sa `CGAL::Linear_cell_complex_for_generalized_map<d,d2,LCCTraits,Items,Alloc>`

\cgalModifEnd

*/
class Linear_cell_complex_bgl_min_items {
public:

/// @}

}; /* end Linear_cell_complex_bgl_min_items */

} /* end namespace CGAL */
