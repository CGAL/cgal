
namespace CGAL {

/*!
\ingroup PkgLinearCellComplexClasses

The auxiliary class `Linear_cell_complex_incremental_builder` supports the incremental
construction of linear cell complexes.

\tparam LCC a linear cell complex
*/

template < class LCC >
class Linear_cell_complex_incremental_builder_3
{
  typedef LCC_ LCC;
  typedef typename LCC::Dart_descriptor             DH;
  typedef typename LCC::Vertex_attribute_descriptor VAH;
  typedef typename LCC::Point                       Point_3;
  typedef typename LCC::size_type                   size_type;

  /// \name Creation
  /// @{

  /*!
   * Constructor
   */
  Linear_cell_complex_incremental_builder_3(LCC & alcc);

/// @}

  /*!
\name Surface Creation

To build a linear cell complex, the following regular expression gives
the correct and allowed order and nesting of method calls from this
section:

\code
begin_surface ( add_vertex  | ( begin_facet  add_vertex_to_facet  end_facet ) ) end_surface
\endcode
*/
/// @{


  /*!
   *
   */
  void begin_surface();


  /*!
   * adds a new vertex for `p` and returns its handle.
   */
  VAH add_vertex(const Point_3& p);

  /*
   * starts a new facet and returns its handle.
   */
  void begin_facet();

  /*!
   *
   */
  void add_vertex_to_facet(size_type i);

  /*!
   *  End of the facet. Returns the first dart of this facet.
   */
  DH end_facet();


  /*!
   * End of the surface construction. Returns one dart of the created surface.
   */
  DH end_surface();

/// @}

/// \name Additional Operations
/// @{

/*!
 * is a synonym for `begin_facet()`, a call to `add_vertex_to_facet()` for each
 * value in the range `[first,beyond)`, and a call to `end_facet()`.
 */
  DH add_facet(std::initializer_list<size_type> l);


/// @}

 };

} // namespace CGAL
