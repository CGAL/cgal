
/*!
\ingroup PkgCombinatorialMapsConcepts
\cgalConcept

The concept `CellAttribute` represents a non void attribute associated
with a cell of a combinatorial map. It can keep a handle to one
dart of its associated cell, and can contain any information.

\cgalHasModel \ref CGAL::Cell_attribute "CGAL::Cell_attribute<CMap,Info_,Tag,OnMerge,OnSplit>"

\sa `CombinatorialMap`
\sa `CombinatorialMapItems`

*/

class CellAttribute {
public:

/// \name Types
/// @{

/*!
%Dart handle type.
*/
typedef unspecified_type Dart_handle;

/*!
%Dart const handle type.
*/
typedef unspecified_type Dart_const_handle;

/*!
Type of the information contained in the attribute.
If `void`, the cell attribute does not have any information.
*/
typedef unspecified_type Info;

/*!
Equals to \ref CGAL::Tag_true "Tag_true" to enable the storage of a
`Dart_handle` of the associated cell, \ref CGAL::Tag_false "Tag_false" otherwise.
*/
typedef unspecified_type Supports_cell_dart;

/*!
Functor called before merging two attributes. Must be
a model of the  [Binary Function](http://www.sgi.com/tech/stl/BinaryFunction.html)
concept having two references
to a model of `CellAttribute` as type of both arguments and `void` as
return type.
*/
typedef unspecified_type On_merge;

/*!
Functor called after an attribute was split in two.
Must be a model of the [Binary Function](http://www.sgi.com/tech/stl/BinaryFunction.html)
concept having two references
to a model of `CellAttribute` as type of both arguments and `void` as
return type.
*/
typedef unspecified_type On_split;

/// @}

/// \name Creation
/// @{

/*!

*/
CellAttribute();

/*!
Constructor initializing the information of this attribute by the
copy constructor \ref Info "Info(info)".
Defined only if \ref Info "Info" is different from `void`.
*/
Cell_attribute(const Info& info);

/// @}

/// \name Access Member Functions
/// @{

/*!
Returns one dart of the cell associated to this attribute.
`NULL` if \ref Supports_cell_dart "Supports_cell_dart" is equal to
\ref CGAL::Tag_false "Tag_false".
*/
Dart_handle dart();

/*!
Returns one dart of the cell associated to this attribute,
when it is const.
`NULL` if \ref Supports_cell_dart "Supports_cell_dart" is equal to
\ref CGAL::Tag_false "Tag_false".
*/
Dart_const_handle dart() const;

/*!
Sets the dart of the cell associated to this attribute to `ahandle`,
if \ref Supports_cell_dart "Supports_cell_dart" is equal to
\ref CGAL::Tag_true "Tag_true".
Otherwise, this method does nothing.
\pre `ahandle` belongs to the cell associated to this attribute.

*/
void set_dart(Dart_handle ahandle);

/*!
Returns the information of this attribute.
Defined only if \ref Info "Info" is different from `void`.
*/
Info& info();

/*!
Returns the information of this attribute, when it is const.
Defined only if \ref Info "Info" is different from `void`.
*/
const Info& info() const;

/// @}

}; /* end CellAttribute */

