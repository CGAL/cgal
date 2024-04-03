/*!
\ingroup PkgCombinatorialMapsConcepts
\cgalConcept

The concept `GenericMap` defines a <I>d</I>-dimensional generic map. This concept is defined only to factorize the common notions between \link CombinatorialMap `CombinatorialMap`\endlink and \link GeneralizedMap `GeneralizedMap`\endlink concepts.

\cgalRefines{DefaultConstructible}

A generic map has a set of darts <I>D</I>, and functions \f$ f_0\f$,\f$ \ldots\f$,\f$ f_{d}\f$ that link these darts between them.

\cgalHasModelsBegin
\cgalHasModelsBare{\link CGAL::Combinatorial_map `CGAL::Combinatorial_map<d,Items,Alloc>`\endlink}
\cgalHasModelsBare{\link CGAL::Generalized_map `CGAL::Generalized_map<d,Items,Alloc>`\endlink}
\cgalHasModelsEnd

\sa `CombinatorialMap`
\sa `GeneralizedMap`

*/

class GenericMap {
public:

/// \name Creation
/// @{

/*!
Construct a new generic map from another one.
The new generic map is created by copying the darts and the non void attributes of bmap. Map must be a model of `GenericMap` concept, which can be defined with a different dimension and/or different attributes than `*this`. In this case, only dimensions that are common to `bmap` and `*this`, and only non void i-attributes of `bmap` whose info type is the same to the info of non void i-attributes of `*this`, are copied.
*/
template<typename Map>
GenericMap(const Map& bmap);

/// @}

/// \name Types
/// @{

/*!
%Dart type.
*/
typedef unspecified_type Dart;

/*!
%Dart descriptor type.
*/
typedef unspecified_type Dart_descriptor;

/*!
%Dart const descriptor type. When using indices, Dart_const_descriptor is equal to Dart_descriptor.
*/
typedef unspecified_type Dart_const_descriptor;

/*!
Information associated to each dart.
*/
typedef unspecified_type Dart_info;

/*!
Size type (an unsigned integral type).
*/
typedef unspecified_type size_type;

/// @}

/// \name Constants
/// @{

/*!
The dimension <I>d</I> of the generic map.
*/
static const unsigned int dimension;

/*!
The number of available Boolean marks of the generic map.
*/
static const size_type NB_MARKS;

/*!
A null descriptor. This descriptor is not valid and should not be used to access a dart or an attribute.
*/
static unspecified_type null_descriptor;

/// @}

/// \name Types for Attributes
/// @{

/*!
The tuple of cell attributes.
It contains at most \link GenericMap::dimension `dimension`\endlink`+1` types
(one for each possible cell of the generic map). Each type of
the tuple must be either a model of the `CellAttribute` concept or
`void`.  The first type corresponds to 0-attributes, the second to
1-attributes and so on.  If the <i>i <sup>th</sup></i> type in the tuple
is `void`, <I>(i-1)</I>-attributes are disabled. Otherwise,
<I>(i-1)</I>-attributes are enabled and have the given type. If the
size of the tuple is <I>k</I>, with <I>k</I> \f$
< \f$ \link GenericMap::dimension `dimension`\endlink`+1`, \f$ \forall \f$ <I>i</I>: <I>k</I> \f$
\leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink,
<I>i</I>-attributes are disabled.
*/
typedef unspecified_type Attributes;

  /*!
    `Attribute_type<i>::%type` is the type of <I>i</I>-attributes, a model of `CellAttribute` concept.
    \link CellAttribute::Dart_descriptor `Attribute_type<i>::type::Dart_descriptor`\endlink is equal to
    \link GenericMap::Dart_descriptor `Dart_descriptor`\endlink, and
    \link CellAttribute::Dart_const_descriptor `Attribute_type<i>::type::Dart_const_descriptor`\endlink is equal to
    \link GenericMap::Dart_const_descriptor `Dart_const_descriptor`\endlink.
    \pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink and
    <I>i</I>-attributes are non `void`.
    \note It can be implemented using a nested template class.
  */
  template <unsigned int i>
  using Attribute_type = unspecified_type;

  /*!
  `Attribute_descriptor<i>::%type` is a descriptor to <I>i</I>-attributes.
  \pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink
       and <I>i</I>-attributes are non void.
  \note It can be implemented using a nested template class.
  */
  template <unsigned int i>
  using Attribute_descriptor = unspecified_type;

  /*!
  `Attribute_descriptor<i>::%type` is a const descriptor to <I>i</I>-attributes.
  \pre 0\f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink
       and <I>i</I>-attributes are non void.
  \note It can be implemented using a nested template class.
  */
  template <unsigned int i>
  using Attribute_const_descriptor = unspecified_type;

/// @}

/// \name Range Types
/// @{

/*!
%Range of all the darts of the generic map.
This type is a model of `Range` concept, its iterator type is bidirectional and its value type is
  \link GenericMap::Dart `Dart`\endlink.
*/
typedef unspecified_type Dart_range;

/*!
Const range of all the darts of the generic map.
This type is a model of `ConstRange` concept, its iterator type is bidirectional and its value type is
  \link GenericMap::Dart `Dart`\endlink.
*/
typedef unspecified_type Dart_const_range;


/*!
`Attribute_range<i>::%type` is the range of all the <I>i</I>-attributes.
  `Attribute_range<i>::%type` is a model of `Range` concept, its iterator type is bidirectional and its value
  type is \link GenericMap::Attribute_type `Attribute_type<i>::type` \endlink.
  \pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink
       and <I>i</I>-attributes are non void.
  \note It can be implemented using a nested template class.
*/
  template <unsigned int i>
  using Attribute_range = unspecified_type;


/*! `Attribute_const_range<i>::%type` is the const range of all the <I>i</I>-attributes.
   `Attribute_const_range<i>::%type` is a model of `ConstRange` concept, its iterator type is bidirectional and
   its value type is \link GenericMap::Attribute_type `Attribute_type<i>::type`\endlink.
  \pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink
       and <I>i</I>-attributes are non void.
  \note It can be implemented using a nested template class.
*/
template <unsigned int i>
using Attribute_const_range = unspecified_type;

/*!
%Range of all the darts of the `<I...>` orbit.
This type is a model of `Range` concept, its iterator type is forward and its value type is
 \link GenericMap::Dart `Dart`\endlink.
*/
template<unsigned int... I>
using Dart_of_orbit_range = unspecified_type;

/*!
Const range of all the darts of the `<I...>` orbit.
This type is a model of `ConstRange` concept, its iterator type is forward and its value type is
  \link GenericMap::Dart `Dart`\endlink.
*/
template<unsigned int ... I>
using Dart_of_orbit_const_range = unspecified_type;

/*!
%Range of all the darts of an <I>i</I>-cell.
Cells are considered in <I>dim</I> dimension, with 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ <I>dim+1</I> and
 0 \f$ \leq \f$ <I>dim</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink. If <I>i==dim+1</I>,
range of all the darts of a connected component.
This type is a model of `Range` concept, its iterator type is forward and its value type is
  \link GenericMap::Dart `Dart`\endlink.
*/
template<unsigned int i,unsigned int dim=dimension>
using Dart_of_cell_range = unspecified_type;

/*!
Const range of all the darts of the <I>i</I>-cell.
Cells are considered in <I>dim</I> dimension, with 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ <I>dim+1</I> and
 0 \f$ \leq \f$ <I>dim</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink. If <I>i==dim+1</I>,
range of all the darts of a connected component.
This type is a model of `ConstRange` concept, its iterator type is forward and its value type is
  \link GenericMap::Dart `Dart`\endlink.
*/
template<unsigned int i,unsigned int dim=dimension>
using Dart_of_cell_const_range = unspecified_type;

/*!
%Range of one dart of each <I>i</I>-cell incident to one <I>j</I>-cell.
Cells are considered in <I>dim</I> dimension,
with 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ <I>dim+1</I>, 0 \f$ \leq \f$ <I>j</I> \f$ \leq \f$ <I>dim+1</I> and
0 \f$ \leq \f$ <I>dim</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink. If <I>i</I>==<I>dim+1</I>,
consider each connected component instead of each <I>i</I>-cell. If <I>j</I>==<I>dim+1</I>,
consider one connected component instead of one <I>j</I>-cell.
This type is a model of `Range` concept, its iterator type is forward and its value type is
  \link GenericMap::Dart `Dart`\endlink.
*/
template<unsigned int i,unsigned int j,unsigned int dim=dimension>
using One_dart_per_incident_cell_range = unspecified_type;

/*!
Const range of one dart of each <I>i</I>-cell incident to one <I>j</I>-cell.
Cells are considered in <I>dim</I> dimension,
with 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ <I>dim+1</I>, 0\f$ \leq \f$ <I>j</I> \f$ \leq \f$ <I>dim+1</I> and
0\f$ \leq \f$ <I>dim</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink. If <I>i==dim+1</I>,
consider each connected component instead of each <I>i</I>-cell. If <I>j==dim+1</I>,
consider one connected component instead of one <I>j</I>-cell.
This type is a model of `ConstRange` concept, its iterator type is forward and its value type is
  \link GenericMap::Dart `Dart`\endlink.
*/
template<unsigned int i,unsigned int j,unsigned int dim=dimension>
using One_dart_per_incident_cell_const_range = unspecified_type;

/*!
%Range of one dart of each <I>i</I>-cell of the generic map.
Cells are considered in <I>dim</I> dimension,
with 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ <I>dim+1</I> and
0\f$ \leq \f$ <I>dim</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink.
If <I>i==dim+1</I>, consider each connected component instead of each <I>i</I>-cell.
This type is a model of `Range` concept, its iterator type is forward and its value type is
  \link GenericMap::Dart `Dart`\endlink.
*/
template<unsigned int i,unsigned int dim=dimension>
using One_dart_per_cell_range = unspecified_type;

/*!
Const range of one dart of each <I>i</I>-cell of the generic map.
Cells are considered in <I>dim</I> dimension,
with 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ <I>dim+1</I> and
0\f$ \leq \f$ <I>dim</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink.
If <I>i==dim+1</I>, consider each connected component instead of each <I>i</I>-cell.
This type is a model of `ConstRange` concept, its iterator type is forward and its value type is
  \link GenericMap::Dart `Dart`\endlink.
*/
template<unsigned int i,unsigned int dim=dimension>
using One_dart_per_cell_const_range = unspecified_type;

/// @}

/// \name Access Member Functions
/// @{

/*!
Returns `true` iff the generic map is empty, i.e.\ it contains no dart.
*/
bool is_empty() const;

/*!
Returns `true` iff the generic map is without <I>i</I>-boundary.

The map is without <I>i</I>-boundary if there is no `i`-free dart.
\pre 1\f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink.
*/
bool is_without_boundary(unsigned int i) const;

/*!
Returns `true` iff the generic map is without boundary in all dimensions.
*/
bool is_without_boundary() const;

/*!
Returns the number of darts in the generic map.
*/
size_type number_of_darts() const;

/*!
Returns the number of <I>i</I>-attributes in the generic map.
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink,
  and <I>i</I>-attributes are non void.
*/
template <unsigned int i>
size_type number_of_attributes() const;

/*!
Returns an upper bound of the id of dart descriptors if indices are used or 0 otherwise.
*/
size_type upper_bound_on_dart_ids() const;

/*!
Returns an upper bound of the id of <I>i</I>-attributes descriptors if indices are used or 0 otherwise.
*/
template <unsigned int i>
size_type upper_bound_on_attribute_ids() const;

/*! Returns `true` if `d` is a descriptor of a used dart (i.e.\ valid).
 */
bool is_dart_used(Dart_const_descriptor d) const;

/*!
Returns `true` iff dart `d` is <I>i</I>-free.
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink.
*/
bool is_free(Dart_const_descriptor d, unsigned int i) const;

/*!
Returns `true` iff dart `d` is <I>i</I>-free.
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink.
*/
template<unsigned int i>
bool is_free(Dart_const_descriptor d) const;

/*!
Returns the highest dimension <I>i</I> such that dart `d` is not <I>i</I>-free. -1 if `d` is free for any dimension.
*/
int highest_nonfree_dimension(Dart_const_descriptor d) const;

/*!
Returns a descriptor to a dart belonging to the other vertex of the edge containing dart `d` (but not necessarily to the same edge). `null_descriptor` if such a dart does not exist.
*/
Dart_descriptor other_extremity(Dart_descriptor d);

/*!
Returns a const descriptor to a dart belonging to the other vertex of the edge containing dart `d`, when the dart is const (but not necessarily to the same edge). `null_descriptor` if such a dart does not exist.
*/
Dart_const_descriptor other_extremity(Dart_const_descriptor d) const;

/*!
Returns a descriptor to a dart belonging to the next edge after `d`, that does not belong to the same <I>0</I>-cell than `d` and that belongs to the same <I>i</I>-cell than `d`for each <I>i</I>, 2\f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink.
*/
Dart_descriptor next(Dart_descriptor d);

/*!
Returns a const descriptor to a dart belonging to the next edge after `d`, that does not belong to the same <I>0</I>-cell than `d` and that belongs to the same <I>i</I>-cell than `d`for each <I>i</I>, 2\f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink.
*/
Dart_const_descriptor next(Dart_const_descriptor d) const;

/*!
Returns a descriptor to a dart belonging to the previous edge before `d`, that does not belong to the same <I>0</I>-cell than `d` and that belongs to the same <I>i</I>-cell than `d`for each <I>i</I>, 2\f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink.
*/
Dart_descriptor previous(Dart_descriptor d);

/*!
Returns a const descriptor to a dart belonging to the previous edge before `d`, that does not belong to the same <I>0</I>-cell than `d` and that belongs to the same <I>i</I>-cell than `d`for each <I>i</I>, 2\f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink.
*/
Dart_const_descriptor previous(Dart_const_descriptor d) const;

/*!
Returns a descriptor to a dart belonging to the opposite <I>i</I>-cell than `d`. This dart does not belong to the same <I>0</I>-cell than `d`, nor to the same <I>i</I>-cell, but belongs to the same <I>i</I>-cell than `d`for each <I>j</I>, 2\f$ \leq \f$ <I>j</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink, <I>j</I> \f$ \neq \f$ <I>i</I>.
\pre 2\f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink.
*/
template<unsigned int i>
Dart_descriptor opposite(Dart_descriptor d);

/*!
Returns a const descriptor to a dart belonging to the opposite <I>i</I>-cell than `d`. This dart does not belong to the same <I>0</I>-cell than `d`, nor to the same <I>i</I>-cell, but belongs to the same <I>i</I>-cell than `d`for each <I>j</I>, 2\f$ \leq \f$ <I>j</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink, <I>j</I> \f$ \neq \f$ <I>i</I>.
\pre 2\f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink.
*/
template<unsigned int i>
Dart_const_descriptor opposite(Dart_const_descriptor d) const;

/*!
Displays on `os` the number of elements of the generic map.
Its number of darts,
its number of <I>i</I>-cells, for each <I>i</I>,
0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink,
and its number of connected components.

Example of output for a 3D combinatorial map containing two disjoint combinatorial tetrahedra:

<TT>\#Darts=24, \#0-cells=8, \#1-cells=12, \#2-cells=8, \#3-cells=2, \#ccs=2</TT>
*/
std::ostream& display_characteristics(std::ostream & os) const;

/// @}

/// \name Attributes Access Member Functions
/// @{
///

/*!
Returns the information associated to dart `d`.
\pre `Dart_info` is not `void`.
*/
Dart_info& info(Dart_descriptor d);
/*!
Returns the information associated to dart `d`, when the dart is const.
\pre `Dart_info` is not `void`.
*/
const Dart_info& info(Dart_const_descriptor d) const;

/*!
Returns a descriptor to the <I>i</I>-attribute associated to dart `d`.
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink, and <I>i</I>-attributes are non `void`.
*/
template <unsigned int i> Attribute_descriptor<i>::type attribute(Dart_descriptor d);

/*!
Returns a const descriptor to the <I>i</I>-attribute associated to dart `d`, when the dart is const.
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink, and <I>i</I>-attributes are non `void`.
*/
template <unsigned int i>
Attribute_const_descriptor<i>::type attribute(Dart_const_descriptor d) const;

/*!
Returns one dart of the cell associated to the <I>i</I>-attribute `*ah`.
`nullptr` if \link CellAttribute::Supports_cell_dart `Supports_cell_dart`\endlink of <I>i</I>-attributes  is equal to \link CGAL::Tag_false `Tag_false`\endlink.
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink, <I>i</I>-attributes are non `void` and `ah`!=nullptr.
*/
template<unsigned int i>
Dart_descriptor dart_of_attribute(typename Attribute_descriptor<i>::type ah);

/*!
Returns one dart of the cell associated to the const <I>i</I>-attribute `*ah`.
`nullptr` if \link CellAttribute::Supports_cell_dart `Supports_cell_dart`\endlink of <I>i</I>-attributes is equal to \link CGAL::Tag_false `Tag_false`\endlink.
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink, <I>i</I>-attributes are non `void` and `ah`!=nullptr.
*/
template<unsigned int i>
Dart_const_descriptor dart_of_attribute(typename Attribute_const_descriptor<i>::type ah) const;

/*!
Returns the information of the <I>i</I>-attribute `*ah`.
Defined only if \link CellAttribute::Info `Info`\endlink of <I>i</I>-attributes is not `void`.
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink, <I>i</I>-attributes are non `void` and `ah`!=nullptr.
*/
template <unsigned int i>
Attribute_type<i>::type::Info& info_of_attribute(typename Attribute_descriptor<i>::type ah);

/*!
Returns the information of the const <I>i</I>-attribute `*ah`.
Defined only if \link CellAttribute::Info `Info`\endlink of <I>i</I>-attributes is not `void`.
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink, <I>i</I>-attributes are non `void` and `ah`!=nullptr.
*/
template <unsigned int i>
const Attribute_type<i>::type::Info& info_of_attribute(typename Attribute_const_descriptor<i>::type ah) const;

/*!
A shortcut for \link GenericMap::info_of_attribute `info_of_attribute<i>`\endlink`(`\link GenericMap::attribute `attribute<i>`\endlink`(adart))`.
\pre \link GenericMap::attribute `attribute<i>`\endlink`(adart)!=nullptr`.
*/
template<unsigned int i>
typename Attribute_type<i>::type::Info & info(Dart_descriptor adart);

/*!
A shortcut for \link GenericMap::info_of_attribute(typename Attribute_const_descriptor<i>::type)const `info_of_attribute<i>`\endlink`(`\link GenericMap::attribute(Dart_const_descriptor)const `attribute<i>`\endlink`(adart))` for const descriptor.
\pre \link GenericMap::attribute(Dart_const_descriptor)const `attribute<i>`\endlink`(adart)!=nullptr`.
*/
template<unsigned int i>
const typename Attribute_type<i>::type::Info & info(Dart_const_descriptor adart) const;

/*!
A shortcut for \link GenericMap::dart_of_attribute `dart_of_attribute<i>`\endlink`(`\link GenericMap::attribute `attribute<i>`\endlink`(adart))`.
\pre `attribute<i>(adart)!=nullptr`.
*/
template<unsigned int i>
Dart_descriptor & dart(Dart_descriptor adart);

/*!
A shortcut for \link GenericMap::dart_of_attribute(typename Attribute_const_descriptor<i>::type)const `dart_of_attribute<i>`\endlink`(`\link GenericMap::attribute(Dart_const_descriptor)const `attribute<i>`\endlink`(adart))` for const descriptor.
\pre `attribute<i>(adart)!=nullptr`.
*/
template<unsigned int i>
Dart_const_descriptor dart(Dart_const_descriptor adart) const;

/*! Returns `true` if ah points to a used i-attribute (i.e.\ valid).
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink, and <I>i</I>-attributes are non `void`.
 */
template<unsigned int i>
bool is_attribute_used(typename Attribute_const_descriptor<i>::type ah) const;

/// @}

/// \name Transformations Between Descriptors and Instances
/// @{

/*!
Returns the dart descriptors of `d`.
\pre `d`\f$ \in \f$ `darts()`.
*/
Dart_descriptor dart_descriptor(Dart& d);

/*!
Returns the dart const descriptors of `d`.
\pre `d`\f$ \in \f$ `darts()`.
*/
Dart_const_descriptor dart_descriptor(const Dart& d) const;

/// @}

/// \name Range Access Member Functions
/// @{

/*!
Returns a range of all the darts in the generic map.
*/
Dart_range& darts();

/*!
Returns a const range of all the darts in the generic map.
*/
Dart_const_range& darts() const;

/*!
Returns a range of all the <I>i</I>-attributes in the generic map.
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink,
    and <I>i</I>-attributes are non void.
*/
template<unsigned int i> Attribute_range<i>::type & attributes();

/*!
Returns a const range of all the <I>i</I>-attributes in the generic map.
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink,
     and <I>i</I>-attributes are non void.
*/
template<unsigned int i> Attribute_const_range<i>::type & attributes() const;

/*!
Returns a range of all the darts of the orbit `<I...>(d)`.
The first element in the range points onto `d`.
\pre `d`\f$ \in \f$ `darts()` and `I...` is a sequence of integers \f$ i_1\f$,\f$ \ldots\f$,\f$ i_k\f$,
    such that 0\f$ \leq \f$\f$ i_1\f$\f$ < \f$\f$ i_2\f$\f$ < \f$\f$ \ldots\f$\f$ < \f$\f$ i_k\f$\f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink.
*/
template<unsigned int ... I> Dart_of_orbit_range darts_of_orbit(Dart_descriptor d);

/*!
Returns a const range of all the darts of the orbit `<I...>(d)`.
The first element in the range points onto `d`.
\pre Same as for the non const version.
*/
template<unsigned int ... I> Dart_of_orbit_const_range darts_of_orbit(Dart_const_descriptor d) const;

/*!
Returns a range of all the darts of the <I>i</I>-cell containing `d`.
The first element in the range points onto `d`.
<I>i</I>-cells are considered in <I>dim</I> dimension. If <I>i==dim+1</I>,
range of all the darts of the connected component containing `d`.
\pre `d`\f$ \in \f$ `darts()`, 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ <I>dim+1</I> and
    0\f$ \leq \f$ <I>dim</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink.
*/
template<unsigned int i,unsigned int dim=dimension> Dart_of_cell_range darts_of_cell(Dart_descriptor d);

/*!
Returns a const range of all the darts of the <I>i</I>-cell containing `d`.
The first element in the range points onto `d`.
<I>i</I>-cells are considered in <I>dim</I> dimension. If <I>i==dim+1</I>,
const range of all the darts of the connected component containing `d`.
\pre Same as for the non const version.
*/
template<unsigned int i,unsigned int dim=dimension> Dart_of_cell_const_range darts_of_cell(Dart_const_descriptor d) const;

/*!
Returns a range of one dart of each <I>i</I>-cell incident to the <I>j</I>-cell containing `d`.
The first element in the range points onto `d`.
Cells are considered in <I>dim</I> dimension. If <I>i==dim+1</I>,
consider each connected component instead of each <I>i</I>-cell. If <I>j==dim+1</I>,
consider the connected component containing `d` instead of the <I>j</I>-cell.
\pre `d`\f$ \in \f$ `darts()`, 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ <I>dim+1</I>,
      0\f$ \leq \f$ <I>j</I> \f$ \leq \f$ <I>dim+1</I> and
      0\f$ \leq \f$ <I>dim</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink.
*/
template<unsigned int i,unsigned int j,unsigned int dim=dimension> One_dart_per_incident_cell_range one_dart_per_incident_cell(Dart_descriptor d);

/*!
Returns a const range of one dart of each <I>i</I>-cell incident to the <I>j</I>-cell containing `d`.
The first element in the range points onto `d`.
Cells are considered in <I>dim</I> dimension. If <I>i==dim+1</I>,
consider each connected component instead of each <I>i</I>-cell. If <I>j==dim+1</I>,
consider the connected component containing `d` instead of the <I>j</I>-cell.
\pre Same as for the non const version.
*/
template<unsigned int i,unsigned int j,unsigned int dim=dimension> One_dart_per_incident_cell_const_range one_dart_per_incident_cell(Dart_const_descriptor d) const;

/*!
Returns a range of one dart of each <I>i</I>-cell in the generic map.
Cells are considered in <I>dim</I> dimension. If <I>i==dim+1</I>,
range of one dart of each connected component in the generic map.
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ <I>dim+1</I> and
  0\f$ \leq \f$ <I>dim</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink.
*/
template<unsigned int i,unsigned int dim=dimension> One_dart_per_cell_range one_dart_per_cell();

/*!
Returns a const range of one dart of each <I>i</I>-cell in the generic map.
Cells are considered in <I>dim</I> dimension. If <I>i==dim+1</I>,
const range of one dart of each connected component in the generic map.
\pre Same as for the non const version.
*/
template<unsigned int i,unsigned int dim=dimension> One_dart_per_cell_const_range one_dart_per_cell() const;

/// @}

/// \name Modifiers
/// @{

/*!
Creates a new dart in the generic map, and returns the corresponding descriptor.
Calls the constructor of dart having `T1` as parameter.
A new dart is initialized to be <I>i</I>-free,
\f$ \forall \f$ <I>i</I>: 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink, and to have no associated attribute for each non void attribute.
Overloads of this member function are defined that take from zero to nine arguments.
With zero argument, `%create_dart()` creates a new dart by using the default constructor.

*/
template<typename T1>
Dart_descriptor create_dart(T1 t1);

/*!
Removes `d` from the generic map.
\pre `d`\f$ \in \f$ `darts()`.

*/
void erase_dart(Dart_descriptor d);

/*!
Creates a new <I>i</I>-attribute in the generic map, and returns the corresponding descriptor.
Calls the constructor of <I>i</I>-attribute having `T1` as parameter.
Overloads of this member function are defined that take from zero to nine arguments.
With zero argument, `create_attribute<i>()` creates a new <I>i</I>-attribute by using the default constructor.
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink,
    and <I>i</I>-attributes are non void.
*/
template<unsigned int i,typename T1> Attribute_descriptor<i>::type create_attribute(T1 t1);

/*!
Removes the <I>i</I>-attribute `*ah` from the generic map.
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink,
<I>i</I>-attributes are non void, and
`*ah`\f$ \in \f$ \link GenericMap::attributes() `attributes<i>()`\endlink.
*/
template <unsigned int i> void erase_attribute(Attribute_descriptor<i>::type ah);

/*!
Associates the <I>i</I>-attribute of all the darts of the <I>i</I>-cell containing `d` to `a`.
\pre `d`\f$ \in \f$ `darts()`, 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink,
<I>i</I>-attributes are non void, and `*ah`\f$ \in \f$ \link GenericMap::attributes() `attributes<i>()`\endlink.
*/
template <unsigned int i> void set_attribute(Dart_descriptor d, Attribute_descriptor<i>::type a);

/*!
Deletes all the darts and all the attributes of the generic map.
*/
void clear();

/*!
Assignment operator.
All darts and attributes are duplicated, and the former generic map is deleted.
*/
GenericMap& operator= (const GenericMap& bmap);

/*!
Swap the current generic map with `bmap`.
There is no copy of darts and attributes thus this method runs in constant time.
*/
void swap(GenericMap& bmap);

/// @}

/// \name Attributes management
/// @{

/*!
Returns the status of the management of the attributes of the generic map. <code>true</code> if the high level operations update the non void attributes (default value); <code>false</code> otherwise.
*/
bool are_attributes_automatically_managed() const;

/*!
Set the status of the management of the attributes of the generic map.

\cgalAdvancedBegin
After calling `set_automatic_attributes_management(false)`, all high level operations will not update non void attributes, until the call of `set_automatic_attributes_management(true)`. The call of `set_automatic_attributes_management(true)` call the \link GenericMap::correct_invalid_attributes `correct_invalid_attributes()`\endlink function.
\cgalAdvancedEnd

*/
void set_automatic_attributes_management(bool update_attributes);

/*!
Correct the invalid attributes of the generic map.
We can have invalid attribute either if we have called \link GenericMap::set_automatic_attributes_management `set_automatic_attributes_management(false)`\endlink before to use some modification operations or if we have modified the generic map by using low level operations.

\f$ \forall i \f$, 0 \f$ \leq \f$ i \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink such that the i-attributes are non void, \f$ \forall \f$ d \f$ \in \f$ `darts()`:
 - if there exists a dart `d2` in the same i-cell than `d` with a different i-attribute, then the i-attribute of `d2` is set to the i-attribute of `d`;
 - if there exists a dart `d2` in a different i-cell than `d` with the same i-attribute, then the i-attribute of all the darts in i-cell(`d`) is set to a new i-attribute (copy of the original attribute);
 - ensures that \link GenericMap::dart_of_attribute `dart_of_attribute(d)`\endlink \f$ \in \f$ i-cell(`d`).
*/
void correct_invalid_attributes();

/// @}

/// \cond SKIP_IN_MANUAL boost::function \endcond

/// \name Dynamic On-Merge/On-Split functors
/// @{

/*!
  Return the current dynamic on-split function associated with i-attributes.
  This is a boost::function returning void and having two references to \link GenericMap::Attribute_type `Attribute_type<i>::type`\endlink as parameters.
  The on-split function is returned by reference so that we can modify it.
*/
  template<int i>
  boost::function<void(typename Attribute_type< i >::type&,
                       typename Attribute_type< i >::type&)>&
  onsplit_function();

/*!
  Return the current dynamic on-split function associated with i-attributes, when *this is const.
  This is a boost::function returning void and having two references to \link GenericMap::Attribute_type `Attribute_type<i>::type`\endlink as parameters.
*/
  template<int i>
  const boost::function<void(typename Attribute_type< i >::type&,
                             typename Attribute_type< i >::type&)>&
  onsplit_function() const;

/*!
  Return the current dynamic on-merge function associated with i-attributes.
  This is a boost::function returning void and having two references to \link GenericMap::Attribute_type `Attribute_type<i>::type`\endlink as parameters.
  The on-merge function is returned by reference so that we can modify it.
*/
  template<int i>
  boost::function<void(typename Attribute_type< i >::type&,
                       typename Attribute_type< i >::type&)>&
  onmerge_function();

/*!
  Return the current dynamic on-merge function associated with i-attributes, when *this is const.
  This is a boost::function returning void and having two references to \link GenericMap::Attribute_type `Attribute_type<i>::type`\endlink as parameters.
*/
  template<int i>
  const boost::function<void(typename Attribute_type< i >::type&,
                             typename Attribute_type< i >::type&)>&
  onmerge_function() const;

/// @}

/// \name Boolean Marks
/// @{

/*!
Reserves a new mark. Returns its
index. If there is no more available free mark, throw the exception Exception_no_more_available_mark.
*/
size_type get_new_mark() const;

/*!
Returns `true` iff `m` is a reserved mark of the generic map.
\pre 0\f$ \leq \f$ <I>m</I> \f$ < \f$ \link GenericMap::NB_MARKS `NB_MARKS`\endlink.
*/
bool is_reserved(size_type m) const;

/*!
Returns `true` iff dart `d` is marked for `m`.
\pre \link GenericMap::is_reserved `is_reserved(m)`\endlink and `d`\f$ \in \f$ `darts()`.
*/
bool is_marked(Dart_const_descriptor d, size_type m) const;

/*!
Marks dart `d` for `m`.
\pre \link GenericMap::is_reserved `is_reserved(m)`\endlink and `d`\f$ \in \f$ `darts()`.
*/
void mark(Dart_const_descriptor d, size_type m) const;

/*!
Unmarks dart `d` for the mark `m`.
\pre \link GenericMap::is_reserved `is_reserved(m)`\endlink and `d`\f$ \in \f$ `darts()`.
*/
void unmark(Dart_const_descriptor d, size_type m) const;

/*!
Inverse the mark `m` for all the darts of the generic map.
All the marked darts become unmarked and all the unmarked darts
become marked.
\pre \link GenericMap::is_reserved `is_reserved(m)`\endlink.
*/
void negate_mark(size_type m) const;

/*!
Unmarks all the darts of the generic map for `m`.
\pre \link GenericMap::is_reserved `is_reserved(m)`\endlink.
*/
void unmark_all(size_type m) const;

/*!
Returns the number of marked darts for `m`.
\pre \link GenericMap::is_reserved `is_reserved(m)`\endlink.
*/
size_type number_of_marked_darts(size_type m) const;

/*!
Return the number of unmarked darts for `m`.
\pre \link GenericMap::is_reserved `is_reserved(m)`\endlink.
*/
size_type number_of_unmarked_darts(size_type m) const;

/*!
Frees mark `m`.
\pre \link GenericMap::is_reserved `is_reserved(m)`\endlink.
*/
void free_mark(size_type m) const;

/// @}

/// \name Constructions
/// @{

/*!
Creates a combinatorial hexahedron (six combinatorial quadrangles 2-sewn together), and adds it in the generic map. Returns a descriptor on one dart of this combinatorial hexahedron.
\pre `dimension` \f$\geq\f$ 2.

\sa `make_edge()`
\sa `make_combinatorial_polygon()`
\sa `make_combinatorial_tetrahedron()`

*/
Dart_descriptor make_combinatorial_hexahedron();

/*!
Creates a combinatorial polygon of length `lg` (a cycle of `lg` edges), and adds it in the generic map. Returns a descriptor on one dart of this combinatorial polygon.
\pre `dimension`\f$ \geq\f$ 1 and `lg`\f$ >\f$ 0.

\sa `make_edge()`
\sa `make_combinatorial_tetrahedron()`
\sa `make_combinatorial_hexahedron()`
*/
Dart_descriptor make_combinatorial_polygon(unsigned int lg);

/*!
Creates a combinatorial tetrahedron (four combinatorial triangles 2-sewn together), and adds it in the generic map. Returns a descriptor on one dart of this combinatorial tetrahedron.
\pre `dimension`\f$ \geq\f$ 2.

\sa `make_edge()`
\sa `make_combinatorial_polygon()`
\sa `make_combinatorial_hexahedron()`
*/
Dart_descriptor make_combinatorial_tetrahedron();

/*!
Creates an isolated edge (two darts sewn to represent one edge and two vertices) and adds it in the generic map. Returns a descriptor on one dart of this edge.
\pre `dimension`\f$ \geq\f$ 2.

\sa `make_combinatorial_polygon()`
\sa `make_combinatorial_tetrahedron()`
\sa `make_combinatorial_hexahedron()`
*/
Dart_descriptor make_edge();

/// @}

/// \name Operations
/// @{

/*!
Inserts a 0-cell in the 1-cell containing `d`. Returns `next(d)`, a descriptor on one dart belonging to the new 0-cell.
\pre \link GenericMap::dimension `dimension`\endlink \f$ \geq\f$ 1 and `d`\f$ \in \f$ \link GenericMap::darts `darts()`\endlink.

See examples for combinatorial map in \cgalFigureRef{fig_cmap_insert_vertex} and for generalized map in \cgalFigureRef{fig_gmap_insert_vertex}.

If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`, if 1-attributes are non `void`, \link CellAttribute::On_split `Attribute_type<1>::type::On_split`\endlink(<I>a</I>,<I>a'</I>) is called, with <I>a</I> the original 1-attribute associated with <I>d</I> and <I>a'</I> the new 1-attribute created during the operation. If set, the dynamic on-split function of 1-attributes is also called on <I>a</I> and <I>a'</I>.

\cgalAdvancedBegin
If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are not updated; thus the generic map can be no more valid after this operation.
\cgalAdvancedEnd

\sa `insert_cell_0_in_cell_2()`
\sa `insert_cell_1_in_cell_2()`
\sa `insert_cell_1_between_two_cells_2()`
\sa `insert_dangling_cell_1_in_cell_2()`
\sa `insert_cell_2_in_cell_3<InputIterator>()`
\sa `remove_cell<i>()`
*/
Dart_descriptor insert_cell_0_in_cell_1(Dart_descriptor d);

/*!
Inserts a 0-cell in the 2-cell containing `d`. The 2-cell is split in triangles, one for each initial edge of the facet. Returns a descriptor on one dart belonging to the new 0-cell.
\pre \link GenericMap::dimension `dimension`\endlink \f$ \geq\f$ 2 and `d`\f$ \in \f$ \link GenericMap::darts `darts()`\endlink.

See examples for combinatorial map in \cgalFigureRef{fig_cmap_triangulation} and for generalized map in \cgalFigureRef{fig_gmap_triangulation}.

If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`, if 2-attributes are non `void`, \link CellAttribute::On_split `Attribute_type<2>::type::On_split`\endlink(<I>a</I>,<I>a'</I>) is called, with <I>a</I> the original 2-attribute associated with `d` and <I>a'</I> each new 2-attribute created during the operation. If set, the dynamic on-split function of 2-attributes is also called on <I>a</I> and <I>a'</I>.

\cgalAdvancedBegin
If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are not updated; thus the generic map can be no more valid after this operation.
\cgalAdvancedEnd

\sa `insert_cell_0_in_cell_2()`
\sa `insert_cell_1_in_cell_2()`
\sa `insert_cell_1_between_two_cells_2()`
\sa `insert_dangling_cell_1_in_cell_2()`
\sa `insert_cell_2_in_cell_3<InputIterator>()`
\sa `remove_cell<i>()`
*/
Dart_descriptor insert_cell_0_in_cell_2(Dart_descriptor d);

/*!
Inserts a 1-cell in the 2-cell containing `d1` and `d2`. Returns `previous(d1)`, a descriptor on one dart belonging to the new 1-cell.
\pre `is_insertable_cell_1_in_cell_2(d1,d2)`.

See examples for combinatorial map in \cgalFigureRef{fig_cmap_insert_edge} and for generalized map in \cgalFigureRef{fig_gmap_insert_edge}.

If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`, if 2-attributes are non `void`, \link CellAttribute::On_split `Attribute_type<2>::type::On_split`\endlink(<I>a</I>,<I>a'</I>) is called, with <I>a</I> the original 2-attribute associated with `d1` and <I>a'</I> the new 2-attribute created during the operation. If set, the dynamic on-split function of 2-attributes is also called on <I>a</I> and <I>a'</I>.

\cgalAdvancedBegin
If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are not updated; thus the generic map can be no more valid after this operation.
\cgalAdvancedEnd

\sa `is_insertable_cell_1_in_cell_2()`
\sa `insert_cell_0_in_cell_1()`
\sa `insert_cell_0_in_cell_2()`
\sa `insert_cell_1_between_two_cells_2()`
\sa `insert_dangling_cell_1_in_cell_2()`
\sa `insert_cell_2_in_cell_3<InputIterator>()`
\sa `remove_cell<i>()`
*/
Dart_descriptor insert_cell_1_in_cell_2(Dart_descriptor d1, Dart_descriptor d2);

/*!
Inserts a 1-cell between the 2-cell containing `d1` and the one containing `d2`. Returns `previous(d1)`, a descriptor on one dart belonging to the new 1-cell.
\pre `is_insertable_cell_1_between_two_cells_2(d1,d2)`.

If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`, call \link CellAttribute::On_merge `Attribute_type<i>::type::On_merge`\endlink(<I>a</I>,<I>a'</I>) is called for all enabled i-attributes, for i>=2, with <I>a</I> the original 2-attribute associated with `d1` and <I>a'</I> the original 2-attribute associated with `d2`. If set, the dynamic on-merge function of i-attributes is also called on <I>a</I> and <I>a'</I>.

\cgalAdvancedBegin
If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are not updated; thus the generic map can be no more valid after this operation.
\cgalAdvancedEnd

\sa `is_insertable_cell_1_between_two_cells_2()`
\sa `insert_cell_0_in_cell_1()`
\sa `insert_cell_0_in_cell_2()`
\sa `insert_cell_1_in_cell_2()`
\sa `insert_dangling_cell_1_in_cell_2()`
\sa `insert_cell_2_in_cell_3<InputIterator>()`
\sa `remove_cell<i>()`
*/
Dart_descriptor insert_cell_1_between_two_cells_2(Dart_descriptor d1, Dart_descriptor d2);

/*! Call `insert_cell_1_in_cell_2()` if `is_insertable_cell_1_in_cell_2(d1, d2)`, otherwise call `insert_cell_1_between_two_cells_2()`.
\sa `insert_cell_1_in_cell_2()`
\sa `insert_cell_1_between_two_cells_2()`
\sa `is_insertable_cell_1_in_cell_2()`
*/
Dart_descriptor insert_cell_1(Dart_descriptor d1, Dart_descriptor d2);

/*!
Inserts a 2-cell along the path of 1-cells containing darts given by the range `[afirst,alast)`. Returns `opposite<2>(*afirst)`, a descriptor on one dart belonging to the new 2-cell.
\pre `is_insertable_cell_2_in_cell_3(afirst,alast)`.

See examples for combinatorial map in \cgalFigureRef{fig_cmap_insert_facet} and for generalized map in \cgalFigureRef{fig_gmap_insert_facet}.

If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`, if 3-attributes are non `void`, \link CellAttribute::On_split `Attribute_type<3>::type::On_split`\endlink(<I>a</I>,<I>a'</I>) is called, with <I>a</I> the original 3-attribute associated with `d` and <I>a'</I> the new 3-attribute created during the operation. If set, the dynamic on-split function of 3-attributes is also called on <I>a</I> and <I>a'</I>.

\cgalAdvancedBegin
If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are not updated; thus the generic map can be no more valid after this operation.
\cgalAdvancedEnd

\sa `is_insertable_cell_2_in_cell_3<InputIterator>()`
\sa `insert_cell_0_in_cell_1()`
\sa `insert_cell_0_in_cell_2()`
\sa `insert_cell_1_in_cell_2()`
\sa `insert_cell_1_between_two_cells_2()`
\sa `insert_dangling_cell_1_in_cell_2()`
\sa `remove_cell<i>()`
*/
template <class InputIterator>
Dart_descriptor insert_cell_2_in_cell_3(InputIterator afirst, InputIterator alast);

/*!
Inserts a 1-cell in a the 2-cell containing `d`, the 1-cell being attached only by one of its extremity to the 0-cell containing `d`. Returns `previous(d)`, a descriptor on the dart belonging to the new 1-cell and to the new 0-cell.
\pre \link GenericMap::dimension `dimension`\endlink \f$ \geq\f$ 2 and `d`\f$ \in \f$ \link GenericMap::darts `darts()`\endlink.

See examples for combinatorial map in \cgalFigureRef{fig_cmap_insert_edge} and for generalized map in \cgalFigureRef{fig_gmap_insert_edge}.

\cgalAdvancedBegin
If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are not updated; thus the generic map can be no more valid after this operation.
\cgalAdvancedEnd

\sa `insert_cell_0_in_cell_1()`
\sa `insert_cell_0_in_cell_2()`
\sa `insert_cell_1_in_cell_2()`
\sa `insert_cell_1_between_two_cells_2()`
\sa `insert_cell_2_in_cell_3<InputIterator>()`
\sa `remove_cell<i>()`

*/
Dart_descriptor insert_dangling_cell_1_in_cell_2(Dart_descriptor d);

/*!
Returns `true` iff it is possible to insert a 1-cell in the generic map between `d1` and `d2`.

This is possible if `d1`\f$ \neq \f$ `d2` and `d1` can be reached from `d2` by using some `previous` and `next` calls.
\pre \link GenericMap::dimension `dimension`\endlink \f$ \geq\f$ 2, `d1`\f$ \in \f$ \link GenericMap::darts `darts()`\endlink, and `d2`\f$ \in \f$ \link GenericMap::darts `darts()`\endlink.

\sa `insert_cell_1_in_cell_2()`
\sa `is_insertable_cell_2_in_cell_3<InputIterator>()`

*/
bool is_insertable_cell_1_in_cell_2(Dart_const_descriptor d1, Dart_const_descriptor d2);

/*!
Returns `true` iff it is possible to insert a 1-cell in the generic map between `d1` and `d2`.

This is possible if `d1`\f$ \neq \f$ `d2` and `d1` can not be reached from `d2` by using some `previous` and `next` calls.
\pre \link GenericMap::dimension `dimension`\endlink \f$ \geq\f$ 2, `d1`\f$ \in \f$ \link GenericMap::darts `darts()`\endlink, and `d2`\f$ \in \f$ \link GenericMap::darts `darts()`\endlink.

\sa `insert_cell_1_between_two_cells_2()`

*/
bool is_insertable_cell_1_between_two_cells_2(Dart_const_descriptor d1, Dart_const_descriptor d2);

/*!
Returns `true` iff it is possible to insert a 2-cell in the generic map along the path of darts given by the range `[afirst,alast)`. The 2-cell can be inserted iff the ordered list of darts form a closed path of edges inside a same volume.
\pre \link GenericMap::dimension `dimension`\endlink \f$ \geq\f$ 3.

\sa `insert_cell_2_in_cell_3<InputIterator>()`
\sa `is_insertable_cell_1_in_cell_2()`

*/
template <class InputIterator>
bool is_insertable_cell_2_in_cell_3(InputIterator afirst, InputIterator alast);

/*!
Returns `true` iff the <I>i</I>-cell containing `d` can be removed.

An <I>i</I>-cell can be removed if `i`==\link GenericMap::dimension `dimension`\endlink or if `i`==\link GenericMap::dimension `dimension`\endlink-1 or if `i`\f$ < \f$ \link GenericMap::dimension `dimension`\endlink-1 and the <I>i</I>-cell containing `d` is incident to at most two (<I>i+1</I>)-cells.
\pre 0\f$ \leq \f$ `i`\f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink and `d`\f$ \in \f$ \link GenericMap::darts `darts()`\endlink.

\sa `remove_cell<i>()`
*/
template <unsigned int i>
bool is_removable(Dart_const_descriptor d);

/*!
Removes the <I>i</I>-cell containing `d`. Returns the number of darts removed from the generic map.
\pre `is_removable<i>(d)`.

See examples in \cgalFigureRef{fig_cmap_insert_vertex}, \cgalFigureRef{fig_cmap_insert_edge} and \cgalFigureRef{fig_cmap_insert_facet}.

If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`, if `i`\f$ < \f$ \link GenericMap::dimension `dimension`\endlink, and <I>i+1</I>-attributes are non `void`, and if there are two distinct (<I>i+1</I>)-cells around dart `d`, \link CellAttribute::On_merge `Attribute_type<i+1>::type::On_merge`\endlink(<I>a1</I>,<I>a2</I>) is called, with <I>a1</I> the (<I>i+1</I>)-attribute associated to `d`, and <I>a2</I> the (<I>i+1</I>)-attribute associated to \f$ \beta_{i+1}\f$(<I>d</I>). If set, the dynamic on-merge function of <I>i+1</I>-attributes is also called on <I>a1</I> and <I>a2</I>.

If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`, if a <I>j</I>-cell is disconnected in two <I>j</I>-cells during the operation, and if <I>j</I>-attributes are non void, \link CellAttribute::On_split `Attribute_type<j>::type::On_split`\endlink(<I>a</I>,<I>a'</I>) is called with <I>a</I> the original <I>j</I>-attribute and <I>a'</I> the new <I>j</I>-attribute created due to the disconnection. If set, the dynamic on-split function of <i>j</i>-attributes is also called on <I>a</I> and <I>a'</I>.

\cgalAdvancedBegin
If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are not updated; thus the generic map can be no more valid after this operation.
\cgalAdvancedEnd

\sa `is_removable<i>()`
\sa `insert_cell_0_in_cell_1()`
\sa `insert_cell_0_in_cell_2()`
\sa `insert_cell_1_in_cell_2()`
\sa `insert_cell_1_between_two_cells_2()`
\sa `insert_dangling_cell_1_in_cell_2()`
\sa `insert_cell_2_in_cell_3<InputIterator>()`
*/
template <unsigned int i>
size_type remove_cell(Dart_descriptor d);

/// @}

}; /* end GenericMap */


