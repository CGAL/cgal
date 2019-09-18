
/*!
\ingroup PkgCombinatorialMapsConcepts
\cgalConcept

The concept `CellAttribute` represents a non void attribute associated with a cell of a generic map. It can keep a handle to one dart of its associated cell, and can contain any information.

\cgalHasModel \link CGAL::Cell_attribute `CGAL::Cell_attribute<Map,Info_,Tag,OnMerge,OnSplit>`\endlink

\sa `GenericMap`
\sa `GenericMapItems`

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
Type of the information contained in the attribute. If `void`, the cell attribute does not have any information.
*/
typedef unspecified_type Info;

/*!
Equals to \link CGAL::Tag_true `Tag_true`\endlink to enable the storage of a `Dart_handle` of the associated cell, \link CGAL::Tag_false `Tag_false`\endlink otherwise.
*/
typedef unspecified_type Supports_cell_dart;

/*!
Functor called before merging two attributes. It must be a binary functor taking as argument two references to a model of `CellAttribute`.
*/
typedef unspecified_type On_merge;

/*!
Functor called after an attribute was split in two. It must be a binary functor taking as argument two references to a model of `CellAttribute`.
*/
typedef unspecified_type On_split;

/// @}

/// \name Creation
/// @{

/*!

*/
CellAttribute();

/*!
Constructor initializing the information of this attribute by the copy constructor \link Info `Info(info)`\endlink. Defined only if \link Info `Info`\endlink is different from `void`.
*/
Cell_attribute(const Info& info);

/// @}

/// \name Access Member Functions
/// @{

/*!
Returns one dart of the cell associated to this attribute. `nullptr` if \link Supports_cell_dart `Supports_cell_dart`\endlink is equal to \link CGAL::Tag_false `Tag_false`\endlink.
*/
Dart_handle dart();

/*!
Returns one dart of the cell associated to this attribute, when it is const. `nullptr` if \link Supports_cell_dart `Supports_cell_dart`\endlink is equal to \link CGAL::Tag_false `Tag_false`\endlink.
*/
Dart_const_handle dart() const;

/*!
Sets the dart of the cell associated to this attribute to `ahandle`, if \link Supports_cell_dart `Supports_cell_dart`\endlink is equal to \link CGAL::Tag_true `Tag_true`\endlink. Otherwise, this method does nothing. \pre `ahandle` belongs to the cell associated to this attribute.

*/
void set_dart(Dart_handle ahandle);

/*!
Returns the information of this attribute. Defined only if \link Info `Info`\endlink is different from `void`.
*/
Info& info();

/*!
Returns the information of this attribute, when it is const. Defined only if \link Info `Info`\endlink is different from `void`.
*/
const Info& info() const;

/// @}

}; /* end CellAttribute */

