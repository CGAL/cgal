namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsOperations

\deprecated Inserts a 0-cell in the 1-cell containing `dh`. Deprecated. Use `cm.insert_cell_0_in_cell_1()` instead.
*/
template < class CMap >
typename CMap::Dart_handle insert_cell_0_in_cell_1(CMap& cm,
typename CMap::Dart_handle dh);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsOperations

\deprecated Inserts a 0-cell in the 2-cell containing `dh`. Deprecated. Use `cm.insert_cell_0_in_cell_2()` instead.
*/
template <class CMap>
typename CMap::Dart_handle insert_cell_0_in_cell_2(CMap & cm,
typename CMap::Dart_handle dh);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsOperations

\deprecated Inserts a 1-cell in the 2-cell containing `dh1` and `dh2`. Deprecated. Use `cm.insert_cell_1_in_cell_2()` instead.
*/
template < class CMap >
typename CMap::Dart_handle insert_cell_1_in_cell_2(CMap& cm,
typename CMap::Dart_handle dh1,typename CMap::Dart_handle dh2);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsOperations

\deprecated Inserts a 2-cell along the path of 1-cells containing darts given by the range `[afirst,alast)`. Deprecated. Use `cm.insert_cell_2_in_cell_3()` instead.
*/
template <class CMap, class InputIterator>
typename CMap::Dart_handle insert_cell_2_in_cell_3(CMap & cm,
InputIterator afirst, InputIterator alast);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsOperations

\deprecated Inserts a 1-cell in a the 2-cell containing `dh`, the 1-cell
being attached only by one of its extremity to the 0-cell containing `dh`. Deprecated. Use `cm.insert_dangling_cell_1_in_cell_2()` instead.
*/
template < class CMap >
typename CMap::Dart_handle insert_dangling_cell_1_in_cell_2(CMap& cm,
typename CMap::Dart_handle dh);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsOperations

\deprecated Returns true iff it is possible to insert a 1-cell in `cm` between `dh1` and `dh2`. Deprecated. Use `cm.is_insertable_cell_1_in_cell_2()` instead.
*/
template < class CMap >
bool is_insertable_cell_1_in_cell_2(const CMap & cm,
typename CMap::Dart_const_handle dh1,
typename CMap::Dart_const_handle dh2);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsOperations

\deprecated Returns true iff it is possible to insert a 2-cell in `cm` along the path
of darts given by the range `[afirst,alast)`. Deprecated. Use `cm.is_insertable_cell_2_in_cell_3()` instead.
*/
template <class CMap, class InputIterator>
bool is_insertable_cell_2_in_cell_3(const CMap & cm,
InputIterator afirst, InputIterator alast);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsOperations

\deprecated Returns true iff the <I>i</I>-cell containing `dh` can be removed. Deprecated. Use `cm.is_removable()` instead.
*/
template <class CMap, unsigned int i>
bool is_removable(const CMap& cm, typename CMap::Dart_const_handle dh);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsOperations

\deprecated Removes the <I>i</I>-cell containing `dh`. Deprecated. Use `cm.remove_cell()` instead.
*/
template <class CMap, unsigned int i>
typename CMap::size_type remove_cell(CMap& cm, typename CMap::Dart_handle dh);

} /* namespace CGAL */

