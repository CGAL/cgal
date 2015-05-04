namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsOperations

Inserts a 0-cell in the 1-cell containing `dh`.
Returns a handle on one dart belonging to the new 0-cell.
\pre \ref CombinatorialMap::dimension "CMap::dimension"\f$ \geq\f$ 1 and
   `*dh`\f$ \in\f$\ref CombinatorialMap::darts "cm.darts()".

See example in \cgalFigureRef{figinsertvertex}.

If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`,
if 1-attributes are non `void`,
\ref CellAttribute::On_split "Attribute_type<1>::type::On_split"(<I>a</I>,<I>a'</I>) is called, with <I>a</I> the original 1-attribute associated with <I>dh</I> and <I>a'</I> the new 1-attribute created during the operation. If set, the dynamic onsplit function of 1-attributes is also called on <I>a</I> and <I>a'</I>.

\cgalAdvancedBegin
If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are
not updated; thus the combinatorial map can be no more valid after this operation.
\cgalAdvancedEnd

\sa `CGAL::insert_cell_0_in_cell_2<CMap>`
\sa `CGAL::insert_cell_1_in_cell_2<CMap>`
\sa `CGAL::insert_dangling_cell_1_in_cell_2<CMap>`
\sa `CGAL::insert_cell_2_in_cell_3<CMap,InputIterator>`
\sa `CGAL::remove_cell<CMap,i>`

*/
template < class CMap >
typename CMap::Dart_handle insert_cell_0_in_cell_1(CMap& cm,
typename CMap::Dart_handle dh);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsOperations

Inserts a 0-cell in the 2-cell containing `dh`.
The 2-cell is split in triangles, one for each initial edge of the facet.
Returns a handle on one dart belonging to the new 0-cell.
\pre \ref CombinatorialMap::dimension "CMap::dimension"\f$ \geq\f$ 2 and
   `*dh`\f$ \in\f$\ref CombinatorialMap::darts "cm.darts()".

See example in \cgalFigureRef{figtriangulate}.

If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`,
if 2-attributes are non `void`,
\ref CellAttribute::On_split "Attribute_type<2>::type::On_split"(<I>a</I>,<I>a'</I>) is called, with <I>a</I> the original 2-attribute associated with `dh` and <I>a'</I> each new 2-attribute created during the operation. If set, the dynamic onsplit function of 2-attributes is also called on <I>a</I> and <I>a'</I>.

\cgalAdvancedBegin
If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are
not updated; thus the combinatorial map can be no more valid after this operation.
\cgalAdvancedEnd

\sa `CGAL::insert_cell_0_in_cell_2<CMap>`
\sa `CGAL::insert_cell_1_in_cell_2<CMap>`
\sa `CGAL::insert_dangling_cell_1_in_cell_2<CMap>`
\sa `CGAL::insert_cell_2_in_cell_3<CMap,InputIterator>`
\sa `CGAL::remove_cell<CMap,i>`


*/
template <class CMap>
typename CMap::Dart_handle insert_cell_0_in_cell_2(CMap & cm,
typename CMap::Dart_handle dh);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsOperations

Inserts a 1-cell in the 2-cell containing `dh1` and `dh2`.
Returns \f$ \beta_0\f$(`dh1`), a handle on one dart belonging to the new 1-cell.
\pre `::is_insertable_cell_1_in_cell_2<Map>(cm,dh1,dh2)`.

See example in \cgalFigureRef{figinsertedge}.

If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`,
if 2-attributes are non `void`, \ref CellAttribute::On_split "Attribute_type<2>::type::On_split"(<I>a</I>,<I>a'</I>) is called, with <I>a</I> the original 2-attribute associated with `dh` and <I>a'</I> the new 2-attribute created during the operation. If set, the dynamic onsplit function of 2-attributes is also called on <I>a</I> and <I>a'</I>.

\cgalAdvancedBegin
If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are
not updated; thus the combinatorial map can be no more valid after this operation.
\cgalAdvancedEnd

\sa `CGAL::is_insertable_cell_1_in_cell_2<CMap>`
\sa `CGAL::insert_cell_0_in_cell_1<CMap>`
\sa `CGAL::insert_cell_0_in_cell_2<CMap>`
\sa `CGAL::insert_dangling_cell_1_in_cell_2<CMap>`
\sa `CGAL::insert_cell_2_in_cell_3<CMap,InputIterator>`
\sa `CGAL::remove_cell<CMap,i>`


*/
template < class CMap >
typename CMap::Dart_handle insert_cell_1_in_cell_2(CMap& cm,
typename CMap::Dart_handle dh1,typename CMap::Dart_handle dh2);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsOperations

Inserts a 2-cell along the path of 1-cells containing darts given
by the range `[afirst,alast)`.
Returns a handle on one dart belonging to the new 2-cell.
\pre `::is_insertable_cell_2_in_cell_3<Map>(cm,afirst,alast)`.

See example in \cgalFigureRef{figinsertface}.

If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`,
if 3-attributes are non `void`, \ref CellAttribute::On_split "Attribute_type<3>::type::On_split"(<I>a</I>,<I>a'</I>) is called, with <I>a</I> the original 3-attribute associated with `dh` and <I>a'</I> the new 3-attribute created during the operation. If set, the dynamic onsplit function of 3-attributes is also called on <I>a</I> and <I>a'</I>.

\cgalAdvancedBegin
If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are
not updated; thus the combinatorial map can be no more valid after this operation.
\cgalAdvancedEnd

\sa `CGAL::is_insertable_cell_2_in_cell_3<CMap,InputIterator>`
\sa `CGAL::insert_cell_0_in_cell_1<CMap>`
\sa `CGAL::insert_cell_0_in_cell_2<CMap>`
\sa `CGAL::insert_cell_1_in_cell_2<CMap>`
\sa `CGAL::insert_dangling_cell_1_in_cell_2<CMap>`
\sa `CGAL::remove_cell<CMap,i>`


*/
template <class CMap, class InputIterator>
typename CMap::Dart_handle insert_cell_2_in_cell_3(CMap & cm,
InputIterator afirst, InputIterator alast);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsOperations

Inserts a 1-cell in a the 2-cell containing `dh`, the 1-cell
being attached only by one of its extremity to the 0-cell containing `dh`.
Returns a handle on the dart belonging to the new 1-cell and to the new 0-cell.
\pre \ref CombinatorialMap::dimension "CMap::dimension"\f$ \geq\f$ 2 and
   `*dh`\f$ \in\f$\ref CombinatorialMap::darts "cm.darts()".

See example in \cgalFigureRef{figinsertedge}.

\cgalAdvancedBegin
If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are
not updated; thus the combinatorial map can be no more valid after this operation.
\cgalAdvancedEnd

\sa `CGAL::insert_cell_0_in_cell_1<CMap>`
\sa `CGAL::insert_cell_0_in_cell_2<CMap>`
\sa `CGAL::insert_cell_1_in_cell_2<CMap>`
\sa `CGAL::insert_cell_2_in_cell_3<CMap,InputIterator>`
\sa `CGAL::remove_cell<CMap,i>`

*/
template < class CMap >
typename CMap::Dart_handle insert_dangling_cell_1_in_cell_2(CMap& cm,
typename CMap::Dart_handle dh);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsOperations

Returns true iff it is possible to insert a 1-cell in `cm`
between `dh1` and `dh2`.

This is possible if `dh1`\f$ \neq\f$`dh2` and `dh1`\f$ \in\f$\f$ \langle{}\f$\f$ \beta_1\f$\f$ \rangle{}\f$(`dh2`).
\pre \ref CombinatorialMap::dimension "CMap::dimension"\f$ \geq\f$ 2,
   `*dh1`\f$ \in\f$\ref CombinatorialMap::darts "cm.darts()", and
   `*dh2`\f$ \in\f$\ref CombinatorialMap::darts "cm.darts()".

\sa `CGAL::insert_cell_1_in_cell_2<CMap>`
\sa `CGAL::is_insertable_cell_2_in_cell_3<CMap,InputIterator>`

*/
template < class CMap >
bool is_insertable_cell_1_in_cell_2(const CMap & cm,
typename CMap::Dart_const_handle dh1,
typename CMap::Dart_const_handle dh2);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsOperations

Returns true iff it is possible to insert a 2-cell in `cm` along the path
of darts given by the range `[afirst,alast)`. The 2-cell can be inserted
iff each couple of consecutive darts of the path <I>a1</I> and <I>a2</I> belong to the
same vertex and the same volume, and if the path is closed.
\pre \ref CombinatorialMap::dimension "CMap::dimension"\f$ \geq\f$ 3.

\sa `CGAL::insert_cell_2_in_cell_3<CMap,InputIterator>`
\sa `CGAL::is_insertable_cell_1_in_cell_2<CMap>`

*/
template <class CMap, class InputIterator>
bool is_insertable_cell_2_in_cell_3(const CMap & cm,
InputIterator afirst, InputIterator alast);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsOperations

Returns true iff the <I>i</I>-cell containing `dh` can be removed.

An <I>i</I>-cell can be removed if `i`==\ref CombinatorialMap::dimension "CMap::dimension"
or if `i`==\ref CombinatorialMap::dimension "CMap::dimension"-1 or
if `i`\f$ <\f$\ref CombinatorialMap::dimension "CMap::dimension"-1 and the <I>i</I>-cell containing `dh`
is incident to at most two (<I>i+1</I>)-cells.
\pre 0\f$ \leq\f$`i`\f$ \leq\f$\ref CombinatorialMap::dimension "CMap::dimension"
  and `*dh`\f$ \in\f$\ref CombinatorialMap::darts "cm.darts()".

\sa `CGAL::remove_cell<CMap,i>`

*/
template <class CMap, unsigned int i>
bool is_removable(const CMap& cm, typename CMap::Dart_const_handle dh);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsOperations

Removes the <I>i</I>-cell containing `dh`.
Returns the number of darts removed from `cm`.
\pre `::is_removable<CMap,i>(cm,dh)`.

See examples in \cgalFigureRef{figinsertvertex}, \cgalFigureRef{figinsertedge} and \cgalFigureRef{figinsertface}.

If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`,
if `i`\f$ <\f$\ref CombinatorialMap::dimension "CMap::dimension", and <I>i+1</I>-attributes are non `void`, and if there are two distinct (<I>i+1</I>)-cells around dart `dh`, \ref CellAttribute::On_merge "Attribute_type<i+1>::type::On_merge"(<I>a1</I>,<I>a2</I>) is called, with <I>a1</I> the (<I>i+1</I>)-attribute associated to `dh`, and <I>a2</I> the (<I>i+1</I>)-attribute associated to \f$ \beta_{i+1}\f$(<I>dh</I>). If set, the dynamic onmerge function of <I>i+1</I>-attributes is also called on <I>a1</I> and <I>a2</I>.

If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`,
if a <I>j</I>-cell is disconnected in two <I>j</I>-cells during the operation, and if <I>j</I>-attributes are non void, \ref CellAttribute::On_split "Attribute_type<j>::type::On_split"(<I>a</I>,<I>a'</I>) is called with <I>a</I> the original <I>j</I>-attribute and <I>a'</I> the new <I>j</I>-attribute created due to the disconnection. If set, the dynamic onsplit function of <i>j</i>-attributes is also called on <I>a</I> and <I>a'</I>.

\cgalAdvancedBegin
If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are
not updated; thus the combinatorial map can be no more valid after this operation.
\cgalAdvancedEnd

\sa `CGAL::is_removable<CMap,i>`
\sa `CGAL::insert_cell_0_in_cell_1<CMap>`
\sa `CGAL::insert_cell_0_in_cell_2<CMap>`
\sa `CGAL::insert_cell_1_in_cell_2<CMap>`
\sa `CGAL::insert_dangling_cell_1_in_cell_2<CMap>`
\sa `CGAL::insert_cell_2_in_cell_3<CMap,InputIterator>`
*/
template <class CMap, unsigned int i>
typename CMap::size_type remove_cell(CMap& cm, typename CMap::Dart_handle dh);

} /* namespace CGAL */

