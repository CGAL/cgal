/*!
\ingroup PkgCombinatorialMapsConcepts
\cgalConcept

The concept `CombinatorialMap` defines a <I>d</I>-dimensional combinatorial map.

\cgalHasModel \ref  CGAL::Combinatorial_map "CGAL::Combinatorial_map<d,Items,Alloc>"

*/
 //@{
class CombinatorialMap {
public:

/// \name Creation
/// @{

/*!
%Default constructor creating an empty combinatorial map.
*/
CombinatorialMap();

/*!
Construct a new combinatorial map from another one.
The new combinatorial map is created by copying the darts and the non void attributes of cmap. CMap must be a model of `CombinatorialMap` concept, which can be defined with a different dimension and/or different attributes than `*this`. In this case, only permutations that are common to `cmap` and `*this`, and only non void i-attributes of `cmap` whose info type is the same to the info of non void i-attributes of `*this`, are copied.
*/
template<typename CMap>
CombinatorialMap(const CMap& cmap);

/// @}

/// \name Types
/// @{

/*!
%Dart type, a model of the `Dart` concept.
*/
typedef unspecified_type Dart;

/*!
%Dart handle type, equal to `Dart::Dart_handle`.
*/
typedef unspecified_type Dart_handle;

/*!
%Dart const handle type, equal to `Dart::Dart_const_handle`.
*/
typedef unspecified_type Dart_const_handle;

/*!
Size type (an unsigned integral type).
*/
typedef unspecified_type size_type;

/// @}

/// \name Constants
/// @{

/*!
The dimension <I>d</I> of the combinatorial map, equal to `Dart::dimension`.
*/
static unsigned int dimension;

/*!
The number of available Boolean marks of the combinatorial map.
*/
static size_type NB_MARKS;

/*!
The null dart handle constant.
A dart `d` is <I>i</I>-free if `beta(d, i)==null_dart_handle`.
Note that `*null_dart_handle`\f$ \notin\f$`darts()`.
*/
Dart_handle null_dart_handle;

/// @}

/// \name Types for Attributes
/// @{

/*!
The tuple of cell attributes.
It contains at most \ref CombinatorialMap::dimension "dimension"`+1` types
(one for each possible cell of the combinatorial map). Each type of
the tuple must be either a model of the `CellAttribute` concept or
`void`.  The first type corresponds to 0-attributes, the second to
1-attributes and so on.  If the <i>i <sup>th</sup></i> type in the tuple
is `void`, <I>(i-1)</I>-attributes are disabled. Otherwise,
<I>(i-1)</I>-attributes are enabled and have the given type. If the
size of the tuple is <I>k</I>, with <I>k</I>\f$
<\f$\ref CombinatorialMap::dimension "dimension"`+1`, \f$ \forall\f$<I>i</I>: <I>k</I>\f$
\leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension",
<I>i</I>-attributes are disabled.
*/
typedef unspecified_type Attributes;

  /*!
    `Attribute_type<i>::%type` is the type of <I>i</I>-attributes, a model of `CellAttribute` concept.
    \ref CellAttribute::Dart_handle "Attribute_type<i>::type::Dart_handle" is equal to
    \ref CombinatorialMap::Dart_handle "Dart_handle", and
    \ref CellAttribute::Dart_const_handle "Attribute_type<i>::type::Dart_const_handle" is equal to
    \ref CombinatorialMap::Dart_const_handle "Dart_const_handle".
    \pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension" and
    <I>i</I>-attributes are non `void`.
    \note It can be implemented using a nested template class.
  */
  template <unsigned int i>
  using Attribute_type = unspecified_type;

  /*!
  `Attribute_handle<i>::%type` is a handle to <I>i</I>-attributes, equal to \link Dart::Attribute_handle `Dart::Attribute_handle<i>::type` \endlink.
  \pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension"
       and <I>i</I>-attributes are non void.
  \note It can be implemented using a nested template class.
  */
  template <unsigned int i>
  using Attribute_handle = unspecified_type;

  /*!
  `Attribute_handle<i>::%type` is a const handle to <I>i</I>-attributes, equal to \link Dart::Attribute_const_handle `Dart::Attribute_const_handle<i>::type` \endlink.
  \pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension"
       and <I>i</I>-attributes are non void.
  \note It can be implemented using a nested template class.
  */
  template <unsigned int i>
  using Attribute_const_handle = unspecified_type;

/// @}

/// \name Range Types
/// @{

/*!
%Range of all the darts of the combinatorial map.
This type is a model of `Range` concept, its iterator type is bidirectional and its value type is
  \ref CombinatorialMap::Dart "Dart".
*/
typedef unspecified_type Dart_range;

/*!
Const range of all the darts of the combinatorial map.
This type is a model of `ConstRange` concept, its iterator type is bidirectional and its value type is
  \ref CombinatorialMap::Dart "Dart".
*/
typedef unspecified_type Dart_const_range;


/*!
`Attribute_range<i>::%type` is the range of all the <I>i</I>-attributes.
  `Attribute_range<i>::%type` is a model of `Range` concept, its iterator type is bidirectional and its value
  type is \link CombinatorialMap::Attribute_type `Attribute_type<i>::type` \endlink.
  \pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension"
       and <I>i</I>-attributes are non void.
  \note It can be implemented using a nested template class.
*/
  template <unsigned int i>
  using Attribute_range = unspecified_type;


/*! `Attribute_const_range<i>::%type` is the const range of all the <I>i</I>-attributes.
   `Attribute_const_range<i>::%type` is a model of `ConstRange` concept, its iterator type is bidirectional and
   its value type is \link CombinatorialMap::Attribute_type `Attribute_type<i>::type`\endlink.
  \pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension"
       and <I>i</I>-attributes are non void.
  \note It can be implemented using a nested template class.
*/
template <unsigned int i>
using Attribute_const_range = unspecified_type;

/*!
%Range of all the darts of the `<Beta...>` orbit.
This type is a model of `Range` concept, its iterator type is forward and its value type is
 \ref CombinatorialMap::Dart "Dart".
*/
template<unsigned int... Beta>
using Dart_of_orbit_range = unspecified_type;

/*!
Const range of all the darts of the `<Beta...>` orbit.
This type is a model of `ConstRange` concept, its iterator type is forward and its value type is
  \ref CombinatorialMap::Dart "Dart".
*/
template<unsigned int ... Beta>
using Dart_of_orbit_const_range = unspecified_type;

/*!
%Range of all the darts of an <I>i</I>-cell.
Cells are considered in <I>dim</I> dimension, with 0\f$ \leq\f$<I>i</I>\f$ \leq\f$<I>dim+1</I> and
 0\f$ \leq\f$<I>dim</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension". If <I>i==dim+1</I>,
range of all the darts of a connected component.
This type is a model of `Range` concept, its iterator type is forward and its value type is
  \ref CombinatorialMap::Dart "Dart".
*/
template<unsigned int i,unsigned int dim=dimension>
using Dart_of_cell_range = unspecified_type;

/*!
Const range of all the darts of the <I>i</I>-cell.
Cells are considered in <I>dim</I> dimension, with 0\f$ \leq\f$<I>i</I>\f$ \leq\f$<I>dim+1</I> and
 0\f$ \leq\f$<I>dim</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension". If <I>i==dim+1</I>,
range of all the darts of a connected component.
This type is a model of `ConstRange` concept, its iterator type is forward and its value type is
  \ref CombinatorialMap::Dart "Dart".
*/
template<unsigned int i,unsigned int dim=dimension>
using Dart_of_cell_const_range = unspecified_type;

/*!
%Range of one dart of each <I>i</I>-cell incident to one <I>j</I>-cell.
Cells are considered in <I>dim</I> dimension,
with 0\f$ \leq\f$<I>i</I>\f$ \leq\f$<I>dim+1</I>, 0\f$ \leq\f$<I>j</I>\f$ \leq\f$<I>dim+1</I> and
0\f$ \leq\f$<I>dim</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension". If <I>i</I>==<I>dim+1</I>,
consider each connected component instead of each <I>i</I>-cell. If <I>j</I>==<I>dim+1</I>,
consider one connected component instead of one <I>j</I>-cell.
This type is a model of `Range` concept, its iterator type is forward and its value type is
  \ref CombinatorialMap::Dart "Dart".
*/
template<unsigned int i,unsigned int j,unsigned int dim=dimension>
using One_dart_per_incident_cell_range = unspecified_type;

/*!
Const range of one dart of each <I>i</I>-cell incident to one <I>j</I>-cell.
Cells are considered in <I>dim</I> dimension,
with 0\f$ \leq\f$<I>i</I>\f$ \leq\f$<I>dim+1</I>, 0\f$ \leq\f$<I>j</I>\f$ \leq\f$<I>dim+1</I> and
0\f$ \leq\f$<I>dim</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension". If <I>i==dim+1</I>,
consider each connected component instead of each <I>i</I>-cell. If <I>j==dim+1</I>,
consider one connected component instead of one <I>j</I>-cell.
This type is a model of `ConstRange` concept, its iterator type is forward and its value type is
  \ref CombinatorialMap::Dart "Dart".
*/
template<unsigned int i,unsigned int j,unsigned int dim=dimension>
using One_dart_per_incident_cell_const_range = unspecified_type;

/*!
%Range of one dart of each <I>i</I>-cell of the combinatorial map.
Cells are considered in <I>dim</I> dimension,
with 0\f$ \leq\f$<I>i</I>\f$ \leq\f$<I>dim+1</I> and
0\f$ \leq\f$<I>dim</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension".
If <I>i==dim+1</I>, consider each connected component instead of each <I>i</I>-cell.
This type is a model of `Range` concept, its iterator type is forward and its value type is
  \ref CombinatorialMap::Dart "Dart".
*/
template<unsigned int i,unsigned int dim=dimension>
using One_dart_per_cell_range = unspecified_type;

/*!
Const range of one dart of each <I>i</I>-cell of the combinatorial map.
Cells are considered in <I>dim</I> dimension,
with 0\f$ \leq\f$<I>i</I>\f$ \leq\f$<I>dim+1</I> and
0\f$ \leq\f$<I>dim</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension".
If <I>i==dim+1</I>, consider each connected component instead of each <I>i</I>-cell.
This type is a model of `ConstRange` concept, its iterator type is forward and its value type is
  \ref CombinatorialMap::Dart "Dart".
*/
template<unsigned int i,unsigned int dim=dimension>
using One_dart_per_cell_const_range = unspecified_type;

/// @}

/// \name Access Member Functions
/// @{

/*!
Returns true iff the combinatorial map is empty, i.e.\ it contains no dart.
*/
bool is_empty() const;

/*!
Returns true iff the combinatorial map is valid.

A combinatorial map is valid (see Sections \ref
sseccombimapanddarts and \ref sseccombimapvalidity) if for all its darts `d`
\f$\in\f$`darts()`:

- `d` is 0-free, or \f$ \beta_1(\beta_0(d))=d\f$;
- `d` is 1-free, or \f$ \beta_0(\beta_1(d))=d\f$;
- \f$ \forall\f$<I>i</I>, 2\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension":
   `d` is <i>i</i>-free, or \f$ \beta_i(\beta_i(d))=d\f$;
- \f$ \forall\f$<I>i</I>, <I>j</I>,
  0\f$ \leq\f$<I>i</I>\f$ <\f$<I>i</I>+2\f$ \leq\f$<I>j</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension"
  such that <I>j</I>\f$ \geq\f$ 3: \f$ \beta_j(\beta_i(d))=\varnothing\f$ or ;
\f$ \beta_j(\beta_i(\beta_j(\beta_i(d))))=d\f$;
- \f$ \forall\f$<I>i</I>, 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension"
  such that <I>i</I>-attributes are non void:
  + \f$ \forall\f$<I>d2</I> in the same <I>i</I>-cell than <I>d</I>: <I>d</I> and <I>d2</I> have the same <I>i</I>-attribute;
  + \f$ \forall\f$<I>d2</I>  in a different <I>i</I>-cell than <I>d</I>: <I>d</I> and <I>d2</I> have different <I>i</I>-attributes.
*/
bool is_valid() const;

/*!
Returns true iff the combinatorial map is without <I>i</I>-boundary.

The map is without <I>i</I>-boundary if there is no `i`-free dart.
\pre 1\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension".
*/
bool is_without_boundary(unsigned int i) const;

/*!
Returns true iff the combinatorial map is without boundary in all dimensions.
*/
bool is_without_boundary() const;

/*!
Returns the number of darts in the combinatorial map.
*/
size_type number_of_darts() const;

/*!
Returns the number of <I>i</I>-attributes in the combinatorial map.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension",
  and <I>i</I>-attributes are non void.
*/
template <unsigned int i>
size_type number_of_attributes() const;

/*! Returns true if dh points to a used dart (i.e.\ valid).
 */
bool is_dart_used(Dart_const_handle dh) const;

/*!
Returns \f$ \beta_j\f$(\f$ \beta_i\f$(`*dh`)).
Overloads of this member function are defined that take from one to nine integer as arguments.
For each function, betas are applied in the same order as their indices are given as parameters.

For example `beta(dh,1)`=\f$ \beta_1\f$(`*dh`),
and `beta(dh,1,2,3,0)`=\f$ \beta_0\f$(\f$ \beta_3\f$(\f$ \beta_2\f$(\f$ \beta_1\f$(`*dh`)))).
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension",
  0\f$ \leq\f$<I>j</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension"
  and `*dh`\f$ \in\f$`darts()`.
*/
Dart_handle beta(Dart_handle dh, int i, int j);

/*!
Returns \f$ \beta_j\f$(\f$ \beta_i\f$(`*dh`)).
Overloads of this member function are defined that take from one to nine integer as arguments.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension",
     0\f$ \leq\f$<I>j</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension"
     and `*dh`\f$ \in\f$`darts()`.

*/
Dart_const_handle beta(Dart_const_handle dh, int i, int j) const;

/*!
Returns \f$ \beta_j\f$(\f$ \beta_i\f$(`*dh`)).
Overloads of this member function are defined that take from one to nine integer as template arguments.
For each function, betas are applied in the same order as their indices are given as template arguments.

For example `beta<1>(dh)`=\f$ \beta_1\f$(`*dh`),
and `beta<1,2,3,0>(dh)`=\f$ \beta_0\f$(\f$ \beta_3\f$(\f$ \beta_2\f$(\f$ \beta_1\f$(`*dh`)))).
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension",
  0\f$ \leq\f$<I>j</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension"
  and `*dh`\f$ \in\f$`darts()`.
*/
template<int i, int j>
Dart_handle beta(Dart_handle dh);

/*!
Returns \f$ \beta_j\f$(\f$ \beta_i\f$(`*dh`)).
Overloads of this member function are defined that take from one to nine integer as template arguments.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension",
     0\f$ \leq\f$<I>j</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension"
     and `*dh`\f$ \in\f$`darts()`.

*/
template<int i, int j>
Dart_const_handle beta(Dart_const_handle dh) const;

/*!
Returns true iff dart `*dh` is <I>i</I>-free.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension".
*/
bool is_free(Dart_const_handle dh, unsigned int i) const;

/*!
Returns true iff dart `*dh` is <I>i</I>-free.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension".
*/
template<unsigned int i>
bool is_free(Dart_const_handle dh) const;

/*!
Returns the highest dimension <I>i</I> such that dart `*dh` is not <I>i</I>-free. -1 if `dh` is free for any dimension.
*/
int highest_nonfree_dimension(Dart_const_handle dh) const;

/*!
Returns a handle to a dart belonging to the same edge than dart `*dh`, and not to the same vertex. `NULL` if such a dart does not exist.
*/
Dart_handle opposite(Dart_handle dh);

/*!
Returns a const handle to a dart belonging to the same edge than dart `*dh`, and not to the same vertex, when the dart is const. `NULL` if such a dart does not exist.
*/
Dart_const_handle opposite(Dart_const_handle dh) const;

/*!
Returns a handle to a dart belonging to the other vertex of the edge containing dart `*dh` (but contrary to `opposite()` not necessarily to the same edge). `NULL` if such a dart does not exist.
*/
Dart_handle other_extremity(Dart_handle dh);

/*!
Returns a const handle to a dart belonging to the other vertex of the edge containing dart `*dh`, when the dart is const (but contrary to `opposite()` not necessarily to the same edge). `NULL` if such a dart does not exist.
*/
Dart_const_handle other_extremity(Dart_const_handle dh) const;

/*!
Displays on `os` the number of elements of the combinatorial map.
Its number of darts,
its number of <I>i</I>-cells, for each <I>i</I>,
0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension",
and its number of connected components.

Example of output for a 3D combinatorial map containing two disjoint
combinatorial tetrahedra:

<TT>\#Darts=24, \#0-cells=8, \#1-cells=12, \#2-cells=8, \#3-cells=2, \#ccs=2</TT>
*/
std::ostream& display_characteristics(std::ostream & os) const;

/// @}

/// \name Attributes Access Member Functions
/// @{
///
/*!
Returns a handle to the <I>i</I>-attribute associated to dart `*dh`.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension", and <I>i</I>-attributes are non `void`.
*/
template <unsigned int i> Attribute_handle<i>::type attribute(Dart_handle dh);

/*!
Returns a const handle to the <I>i</I>-attribute associated to dart `*dh`, when the dart is const.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension", and <I>i</I>-attributes are non `void`.
*/
template <unsigned int i>
Attribute_const_handle<i>::type attribute(Dart_const_handle dh) const;

/*!
Returns one dart of the cell associated to the <I>i</I>-attribute `*ah`.
`NULL` if \ref CellAttribute::Supports_cell_dart "Supports_cell_dart" of <I>i</I>-attributes  is equal to \ref CGAL::Tag_false "Tag_false".
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension", <I>i</I>-attributes are non `void` and `ah`!=NULL.
*/
template<unsigned int i>
Dart_handle dart_of_attribute(typename Attribute_handle<i>::type ah);

/*!
Returns one dart of the cell associated to the const <I>i</I>-attribute `*ah`.
`NULL` if \ref CellAttribute::Supports_cell_dart "Supports_cell_dart" of <I>i</I>-attributes is equal to \ref CGAL::Tag_false "Tag_false".
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension", <I>i</I>-attributes are non `void` and `ah`!=NULL.
*/
template<unsigned int i>
Dart_const_handle dart_of_attribute(typename Attribute_const_handle<i>::type ah) const;

/*!
Returns the information of the <I>i</I>-attribute `*ah`.
Defined only if \ref CellAttribute::Info "Info" of <I>i</I>-attributes is not `void`.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension", <I>i</I>-attributes are non `void` and `ah`!=NULL.
*/
template <unsigned int i>
Attribute_type<i>::type::Info& info_of_attribute(typename Attribute_handle<i>::type ah);

/*!
Returns the information of the const <I>i</I>-attribute `*ah`.
Defined only if \ref CellAttribute::Info "Info" of <I>i</I>-attributes is not `void`.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension", <I>i</I>-attributes are non `void` and `ah`!=NULL.
*/
template <unsigned int i>
const Attribute_type<i>::type::Info& info_of_attribute(typename Attribute_const_handle<i>::type ah) const;

/*!
A shorcut for \link CombinatorialMap::info_of_attribute `info_of_attribute<i>`\endlink`(`\link CombinatorialMap::attribute `attribute<i>`\endlink`(adart))`.
\pre \ref CombinatorialMap::attribute "attribute<i>"`(adart)!=NULL`.
*/
template<unsigned int i>
typename Attribute_type<i>::type::Info & info(Dart_handle adart);

/*!
A shorcut for \link CombinatorialMap::info_of_attribute(typename Attribute_const_handle<i>::type)const `info_of_attribute<i>`\endlink`(`\link CombinatorialMap::attribute(Dart_const_handle)const `attribute<i>`\endlink`(adart))` for const handle.
\pre \link CombinatorialMap::attribute(Dart_const_handle)const `attribute<i>`\endlink`(adart)!=NULL`.
*/
template<unsigned int i>
const typename Attribute_type<i>::type::Info & info(Dart_const_handle adart) const;

/*!
A shorcut for \link CombinatorialMap::dart_of_attribute `dart_of_attribute<i>`\endlink`(`\link CombinatorialMap::attribute `attribute<i>`\endlink`(adart))`.
\pre `attribute<i>(adart)!=NULL`.
*/
template<unsigned int i>
Dart_handle & dart(Dart_handle adart);

/*!
A shorcut for \link CombinatorialMap::dart_of_attribute(typename Attribute_const_handle<i>::type)const `dart_of_attribute<i>`\endlink`(`\link CombinatorialMap::attribute(Dart_const_handle)const `attribute<i>`\endlink`(adart))` for const handle.
\pre `attribute<i>(adart)!=NULL`.
*/
template<unsigned int i>
Dart_const_handle dart(Dart_const_handle adart) const;

/*! Returns true if ah points to a used i-attribute (i.e.\ valid).
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension", and <I>i</I>-attributes are non `void`.
 */
template<unsigned int i>
bool is_attribute_used(typename Attribute_const_handle<i>::type ah) const;

/// @}

/// \name Transformations Between Handles and Instances
/// @{

/*!
Returns the dart handle of `d`.
\pre `d`\f$ \in\f$`darts()`.
*/
Dart_handle dart_handle(Dart& d);

/*!
Returns the dart const handle of `d`.
\pre `d`\f$ \in\f$`darts()`.
*/
Dart_const_handle dart_handle(const Dart& d) const;

/// @}

/// \name Range Access Member Functions
/// @{

/*!
Returns a range of all the darts in the combinatorial map.
*/
Dart_range& darts();

/*!
Returns a const range of all the darts in the combinatorial map.
*/
Dart_const_range& darts() const;

/*!
Returns a range of all the <I>i</I>-attributes in the combinatorial map.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension",
    and <I>i</I>-attributes are non void.
*/
template<unsigned int i> Attribute_range<i>::type & attributes();

/*!
Returns a const range of all the <I>i</I>-attributes in the combinatorial map.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension",
     and <I>i</I>-attributes are non void.
*/
template<unsigned int i> Attribute_const_range<i>::type & attributes() const;

/*!
Returns a range of all the darts of the orbit `<Beta...>(*dh)`.
The first element in the range points onto `*dh`.
\pre `*dh`\f$ \in\f$`darts()` and `Beta...` is a sequence of integers \f$ i_1\f$,\f$ \ldots\f$,\f$ i_k\f$,
    such that 0\f$ \leq\f$\f$ i_1\f$\f$ <\f$\f$ i_2\f$\f$ <\f$\f$ \ldots\f$\f$ <\f$\f$ i_k\f$\f$ \leq\f$\ref CombinatorialMap::dimension "dimension",
     and (\f$ i_1\f$\f$ \neq\f$ 0 or \f$ i_2\f$\f$ \neq\f$ 1).
*/
template<unsigned int ... Beta> Dart_of_orbit_range darts_of_orbit(Dart_handle dh);

/*!
Returns a const range of all the darts of the orbit `<Beta...>(*dh)`.
The first element in the range points onto `*dh`.
\pre Same as for the non const version.
*/
template<unsigned int ... Beta> Dart_of_orbit_const_range darts_of_orbit(Dart_const_handle dh) const;

/*!
Returns a range of all the darts of the <I>i</I>-cell containing `*dh`.
The first element in the range points onto `*dh`.
<I>i</I>-cells are considered in <I>dim</I> dimension. If <I>i==dim+1</I>,
range of all the darts of the connected component containing `dh`.
\pre `*dh`\f$ \in\f$`darts()`, 0\f$ \leq\f$<I>i</I>\f$ \leq\f$<I>dim+1</I> and
    0\f$ \leq\f$<I>dim</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension".
*/
template<unsigned int i,unsigned int dim=dimension> Dart_of_cell_range darts_of_cell(Dart_handle dh);

/*!
Returns a const range of all the darts of the <I>i</I>-cell containing `*dh`.
The first element in the range points onto `*dh`.
<I>i</I>-cells are considered in <I>dim</I> dimension. If <I>i==dim+1</I>,
const range of all the darts of the connected component containing `*dh`.
\pre Same as for the non const version.
*/
template<unsigned int i,unsigned int dim=dimension> Dart_of_cell_const_range darts_of_cell(Dart_const_handle dh) const;

/*!
Returns a range of one dart of each <I>i</I>-cell incident to the <I>j</I>-cell containing `*dh`.
The first element in the range points onto `*dh`.
Cells are considered in <I>dim</I> dimension. If <I>i==dim+1</I>,
consider each connected component instead of each <I>i</I>-cell. If <I>j==dim+1</I>,
consider the connected component containing `*dh` instead of the <I>j</I>-cell.
\pre `*dh`\f$ \in\f$`darts()`, 0\f$ \leq\f$<I>i</I>\f$ \leq\f$<I>dim+1</I>,
      0\f$ \leq\f$<I>j</I>\f$ \leq\f$<I>dim+1</I> and
      0\f$ \leq\f$<I>dim</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension".
*/
template<unsigned int i,unsigned int j,unsigned int dim=dimension> One_dart_per_incident_cell_range one_dart_per_incident_cell(Dart_handle dh);

/*!
Returns a const range of one dart of each <I>i</I>-cell incident to the <I>j</I>-cell containing `*dh`.
The first element in the range points onto `*dh`.
Cells are considered in <I>dim</I> dimension. If <I>i==dim+1</I>,
consider each connected component instead of each <I>i</I>-cell. If <I>j==dim+1</I>,
consider the connected component containing `*dh` instead of the <I>j</I>-cell.
\pre Same as for the non const version.
*/
template<unsigned int i,unsigned int j,unsigned int dim=dimension> One_dart_per_incident_cell_const_range one_dart_per_incident_cell(Dart_const_handle dh) const;

/*!
Returns a range of one dart of each <I>i</I>-cell in the combinatorial map.
Cells are considered in <I>dim</I> dimension. If <I>i==dim+1</I>,
range of one dart of each connected component in the combinatorial map.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$<I>dim+1</I> and
  0\f$ \leq\f$<I>dim</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension".
*/
template<unsigned int i,unsigned int dim=dimension> One_dart_per_cell_range one_dart_per_cell();

/*!
Returns a const range of one dart of each <I>i</I>-cell in the combinatorial map.
Cells are considered in <I>dim</I> dimension. If <I>i==dim+1</I>,
const range of one dart of each connected component in the combinatorial map.
\pre Same as for the non const version.
*/
template<unsigned int i,unsigned int dim=dimension> One_dart_per_cell_const_range one_dart_per_cell() const;

/// @}

/// \name Modifiers
/// @{

/*!
Creates a new dart in the combinatorial map, and returns the corresponding handle.
Calls the constructor of dart having `T1` as parameter.
A new dart is initialized to be <I>i</I>-free,
\f$ \forall\f$<I>i</I>: 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension",
  and to have no associated attribute for each non void attribute.
Overloads of this member function are defined that take from zero to nine arguments.
With zero argument, `%create_dart()` creates a new dart by using the default constructor.

*/
template<typename T1>
Dart_handle create_dart(T1 t1);

/*!
Removes `*dh` from the combinatorial map.
\pre `*dh`\f$ \in\f$`darts()`.

*/
void erase_dart(Dart_handle dh);

/*!
Creates a new <I>i</I>-attribute in the combinatorial map, and returns the corresponding handle.
Calls the constructor of <I>i</I>-attribute having `T1` as parameter.
Overloads of this member function are defined that take from zero to nine arguments.
With zero argument, `create_attribute<i>()` creates a new <I>i</I>-attribute by using the default constructor.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension",
    and <I>i</I>-attributes are non void.
*/
template<unsigned int i,typename T1> Attribute_handle<i>::type create_attribute(T1 t1);

/*!
Removes the <I>i</I>-attribute `*ah` from the combinatorial map.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension",
<I>i</I>-attributes are non void, and
`*ah`\f$ \in\f$\ref CombinatorialMap::attributes() "attributes<i>()".
*/
template <unsigned int i> void erase_attribute(Attribute_handle<i>::type ah);

/*!
Associates the <I>i</I>-attribute of all the darts of the <I>i</I>-cell containing `*dh` to `*ah`.
\pre `*dh`\f$ \in\f$`darts()`, 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension",
<I>i</I>-attributes are non void, and `*ah`\f$ \in\f$\ref CombinatorialMap::attributes() "attributes<i>()".
*/
template <unsigned int i> void set_attribute(Dart_handle dh, Attribute_handle<i>::type ah);

/*!
Deletes all the darts and all the attributes of the combinatorial map.
*/
void clear();

/*!
Assignment operator.
All darts and attributes are duplicated, and the former combinatorial map is deleted.
*/
CombinatorialMap& operator= (const CombinatorialMap& cmap);

/*!
Swap the current combinatorial map with `cmap`.
There is no copy of darts and attributes thus this method runs in constant time.
*/
void swap(CombinatorialMap& cmap);

/// @}

/// \name Attributes management
/// @{

/*!
Returns the status of the management of the attributes of the combinatorial map. <code>true</code> if the high level operations update the non void attributes (default value); <code>false</code> otherwise.
*/
bool are_attributes_automatically_managed() const;

/*!
Set the status of the managment of the attributes of the combinatorial map.

\cgalAdvancedBegin
After calling `set_automatic_attributes_management(false)`, all high level operations will not update non void attributes, until the call of `set_automatic_attributes_management(true)`. The call of `set_automatic_attributes_management(true)` call the \link CombinatorialMap::correct_invalid_attributes `correct_invalid_attributes()`\endlink function.
\cgalAdvancedEnd

*/
void set_automatic_attributes_management(bool update_attributes);

/*!
Correct the invalid attributes of the combinatorial map.
We can have invalid attribute either if we have called \link CombinatorialMap::set_automatic_attributes_management `set_automatic_attributes_management(false)`\endlink before to use some modification operations or if we have modified the combinatorial map by using low level operations.

\f$ \forall i \f$, 0 \f$ \leq \f$ i \f$ \leq \f$ \ref CombinatorialMap::dimension "dimension" such that the i-attributes are non void, \f$ \forall \f$ d \f$ \in\f$`darts()`:
 - if there exists a dart `d2` in the same i-cell than `d` with a different i-attribute, then the i-attribute of `d2` is set to the i-attribute of `d`;
 - if there exists a dart `d2` in a different i-cell than `d` with the same i-attribute, then the i-attribute of all the darts in i-cell(`d`) is set to a new i-attribute (copy of the original attribute);
 - ensures that \link CombinatorialMap::dart_of_attribute `dart_of_attribute(d)`\endlink \f$ \in \f$ i-cell(`d`).
*/
void correct_invalid_attributes();

/// @}

/// \name Operations
/// @{

/*!
Returns true iff `*dh1` can be <I>i</I>-sewn with `*dh2` by keeping the combinatorial map valid.
This is true if there is
a bijection <I>f</I> between all the darts of the orbit
<I>D1</I>=\f$ \langle{}\f$\f$ \beta_1\f$,\f$ \ldots\f$,\f$ \beta_{i-2}\f$,\f$ \beta_{i+2}\f$,\f$ \ldots\f$,\f$ \beta_d\f$\f$ \rangle{}\f$(<I>*dh1</I>) and
<I>D2</I>=\f$ \langle{}\f$\f$ \beta_1\f$,\f$ \ldots\f$,\f$ \beta_{i-2}\f$,\f$ \beta_{i+2}\f$,\f$ \ldots\f$,\f$ \beta_d\f$\f$ \rangle{}\f$(<I>*dh2</I>)
satisfying: <I>f</I>(<I>*dh1</I>)=<I>*dh2</I>, and for all <I>e</I>\f$ \in\f$<I>D1</I>, for all <I>j</I>\f$ \in\f${1,\f$ \ldots\f$,<I>i</I>-2,<I>i</I>+2,\f$ \ldots\f$,<I>d</I>},
<I>f</I>(\f$ \beta_j\f$(<I>e</I>))=\f$ \beta_j^{-1}\f$(<I>f</I>(<I>e</I>)).
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension",
  `*dh1`\f$ \in\f$`darts()`, and `*dh2`\f$ \in\f$`darts()`.
*/
template <unsigned int i> bool is_sewable(Dart_const_handle dh1, Dart_const_handle dh2) const;

/*!
  <I>i</I>-sew darts `*dh1` and `*dh2`, by keeping the combinatorial map valid.
  Links by \f$ \beta_i\f$
two by two all the darts of the orbit
<I>D1</I>=\f$ \langle{}\f$\f$ \beta_1\f$,\f$ \ldots\f$,\f$ \beta_{i-2}\f$,\f$ \beta_{i+2}\f$,\f$ \ldots\f$,\f$ \beta_d\f$\f$ \rangle{}\f$(`*dh1`) and
<I>D2</I>=\f$ \langle{}\f$\f$ \beta_0\f$,\f$ \beta_2\f$,\f$ \ldots\f$,\f$ \beta_{i-2}\f$,\f$ \beta_{i+2}\f$,\f$ \ldots\f$,\f$ \beta_d\f$\f$ \rangle{}\f$(`*dh2`)
such that <I>d2</I>=<I>f</I>(<I>d1</I>).

<I>f</I> is the bijection between <I>D1</I> and <I>D2</I>
satisfying: <I>f</I>(<I>*dh1</I>)=<I>*dh2</I>, and for all <I>e</I>\f$ \in\f$<I>D1</I>, for all
<I>j</I>\f$ \in\f${1,\f$ \ldots\f$,<I>i</I>-2,<I>i</I>+2,\f$ \ldots\f$,<I>d</I>},
<I>f</I>(\f$ \beta_j\f$(<I>e</I>))=\f$ \beta_j^{-1}\f$(<I>f</I>(<I>e</I>)).

If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`,
when necessary, non void
attributes are updated to ensure the validity of the combinatorial map: for each
<I>j</I>-cells <I>c1</I> and <I>c2</I> which are merged into one <I>j</I>-cell during
the sew, the two associated attributes <I>attr1</I> and <I>attr2</I> are
considered. If one attribute is
NULL and the other not, the non NULL attribute is associated to all
the darts of the resulting cell. When the two attributes are non
NULL, functor \ref CellAttribute::On_merge "Attribute_type<i>::type::On_merge"
is called on
the two attributes <I>attr1</I> and <I>attr2</I>. If set, the dynamic onmerge function of <i>i</i>-attributes is also called on <I>attr1</I> and <I>attr2</I>. Then, the attribute
<I>attr1</I> is associated to all darts of the resulting
<I>j</I>-cell. Finally, attribute <I>attr2</I> is removed from the combinatorial map.
\pre \ref CombinatorialMap::is_sewable "is_sewable<i>(dh1,dh2)".

\cgalAdvancedBegin
If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are
not updated; thus the combinatorial map can be no more valid after this operation.
\cgalAdvancedEnd

*/
template <unsigned int i> void sew(Dart_handle dh1,Dart_handle dh2);

/*!
  <I>i</I>-unsew darts `*dh` and \f$ \beta_i\f$`(*dh)`, by keeping the combinatorial map valid.
Unlinks by \f$ \beta_i\f$ all the darts in the
orbit
\f$ \langle{}\f$\f$ \beta_1\f$,\f$ \ldots\f$,\f$ \beta_{i-2}\f$,\f$ \beta_{i+2}\f$,\f$ \ldots\f$,\f$ \beta_d\f$\f$ \rangle{}\f$(`*dh`). If
\link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`,
when necessary, non void
attributes are updated to ensure the validity of the combinatorial map: for each
<I>j</I>-cell <I>c</I> split in two <I>j</I>-cells <I>c1</I> and <I>c2</I> by the
operation, if <I>c</I> is associated to a <I>j</I>-attribute <I>attr1</I>, then
this attribute is duplicated into <I>attr2</I>, and all the darts
belonging to <I>c2</I> are associated with this new attribute. Finally,
the functor \ref CellAttribute::On_split "Attribute_type<i>::type::On_split"
is called on the
two attributes <I>attr1</I> and <I>attr2</I>. If set, the dynamic onsplit function of <i>i</i>-attributes is also called on <I>attr1</I> and <I>attr2</I>.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension",
     `*dh`\f$ \in\f$`darts()` and `*dh` is not <I>i</I>-free.

\cgalAdvancedBegin
If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are
not updated thus the combinatorial map can be no more valid after this operation.
\cgalAdvancedEnd
*/
template <unsigned int i> void unsew(Dart_handle dh);

/*!
Links `*dh1` and `*dh2` by \f$ \beta_i\f$.
The combinatorial map can be no more valid after this operation. If
\link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`,
non void attributes of `*dh1` and
`*dh2` are updated: if one dart has an attribute and the second
dart not, the non null attribute is associated to the dart having a null attribute.
If both darts have an attribute,
the attribute of `*dh1` is associated to `*dh2`.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension",
    `*dh1`\f$ \in\f$`darts()`, `*dh2`\f$ \in\f$`darts()` and (<I>i</I>\f$ <\f$ 2 or `dh1`\f$ \neq\f$`dh2`).
*/
template <unsigned int i> void link_beta(Dart_handle dh1, Dart_handle dh2);

/*!
Unlinks `*dh` and \f$ \beta_i\f$(`*dh`) by \f$ \beta_i\f$.
The combinatorial map can be no more valid after this operation.
Attributes of `*dh` and \f$ \beta_i\f$(`*dh`)
are not modified.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension",
     `*dh`\f$ \in\f$`darts()`, and `*dh` is not <I>i</I>-free.
*/
template <unsigned int i> void unlink_beta(Dart_handle dh);


/*!
  Reverse the orientation (swap \f$ \beta_0\f$ and \f$ \beta_1\f$ links) of the entire map.
*/
void reverse_orientation();

/*!
    Reverse the orientation (swap \f$ \beta_0\f$ and \f$ \beta_1\f$ links) of the connected component containing the given dart.
*/
void reverse_orientation_connected_component(Dart_handle adart);

/// @}

/// \name Dynamic Onmerge/Onsplit functors
/// @{

/*!
  Return the current dynamic onsplit function associated with i-attributes.
  This is a boost:function returning void and having two references to \link CombinatorialMap::Attribute_type `Attribute_type<i>::type`\endlink as parameters.
  The onsplit function is returned by reference so that we can modify it.
*/
  template<int i>
  boost::function<void(typename Attribute_type< i >::type&,
                       typename Attribute_type< i >::type&)>&
  onsplit_function();

/*!
  Return the current dynamic onsplit function associated with i-attributes, when *this is const.
  This is a boost:function returning void and having two references to \link CombinatorialMap::Attribute_type `Attribute_type<i>::type`\endlink as parameters.
*/
  template<int i>
  const boost::function<void(typename Attribute_type< i >::type&,
                             typename Attribute_type< i >::type&)>&
  onsplit_function() const;

/*!
  Return the current dynamic onmerge function associated with i-attributes.
  This is a boost:function returning void and having two references to \link CombinatorialMap::Attribute_type `Attribute_type<i>::type`\endlink as parameters.
  The onmerge function is returned by reference so that we can modify it.
*/
  template<int i>
  boost::function<void(typename Attribute_type< i >::type&,
                       typename Attribute_type< i >::type&)>&
  onmerge_function();

/*!
  Return the current dynamic onmerge function associated with i-attributes, when *this is const.
  This is a boost:function returning void and having two references to \link CombinatorialMap::Attribute_type `Attribute_type<i>::type`\endlink as parameters.
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
Returns true iff `m` is a reserved mark of the combinatorial map.
\pre 0\f$ \leq\f$<I>m</I>\f$ <\f$\ref  CombinatorialMap::NB_MARKS "NB_MARKS".
*/
bool is_reserved(size_type m) const;

/*!
Returns true iff `*dh` is marked for `m`.
\pre \ref CombinatorialMap::is_reserved "is_reserved(m)" and `*dh`\f$ \in\f$`darts()`.
*/
bool is_marked(Dart_const_handle dh, size_type m) const;

/*!
Marks `*dh` for `m`.
\pre \ref CombinatorialMap::is_reserved "is_reserved(m)" and `*dh`\f$ \in\f$`darts()`.
*/
void mark(Dart_const_handle dh, size_type m) const;

/*!
Unmarks `*dh` for the mark `m`.
\pre \ref CombinatorialMap::is_reserved "is_reserved(m)" and `*dh`\f$ \in\f$`darts()`.
*/
void unmark(Dart_const_handle dh, size_type m) const;

/*!
Inverse the mark `m` for all the darts of the combinatorial map.
All the marked darts become unmarked and all the unmarked darts
become marked.
\pre \ref CombinatorialMap::is_reserved "is_reserved(m)".
*/
void negate_mark(size_type m) const;

/*!
Unmarks all the darts of the combinatorial map for `m`.
\pre \ref CombinatorialMap::is_reserved "is_reserved(m)".
*/
void unmark_all(size_type m) const;

/*!
Returns the number of marked darts for `m`.
\pre \ref CombinatorialMap::is_reserved "is_reserved(m)".
*/
size_type number_of_marked_darts(size_type m) const;

/*!
Return the number of unmarked darts for `m`.
\pre \ref CombinatorialMap::is_reserved "is_reserved(m)".
*/
size_type number_of_unmarked_darts(size_type m) const;

/*!
Frees mark `m`.
\pre \ref CombinatorialMap::is_reserved "is_reserved(m)".
*/
void free_mark(size_type m) const;

/// @}

/// \name Constructions
/// @{

/*!
Creates a combinatorial hexahedron (six combinatorial quadrangles linked together by \f$ \beta_2\f$), and adds it in the combinatorial map. Returns a handle on one dart of this combinatorial hexahedron.
\pre `dimension` \f$\geq\f$ 2.

\sa `make_edge`
\sa `make_combinatorial_polygon`
\sa `make_combinatorial_tetrahedron`

*/
Dart_handle make_combinatorial_hexahedron();

/*!
Creates a combinatorial polygon of length `lg` (`lg` darts linked by \f$ \beta_1\f$), and adds it in the combinatorial map. Returns a handle on one dart of this combinatorial polygon.
\pre `dimension`\f$ \geq\f$ 1 and `lg`\f$ >\f$ 0.

\sa `make_edge`
\sa `make_combinatorial_tetrahedron`
\sa `make_combinatorial_hexahedron`
*/
Dart_handle make_combinatorial_polygon(unsigned int lg);

/*!
Creates a combinatorial tetrahedron (four combinatorial triangles linked together by \f$ \beta_2\f$), and adds it in the combinatorial map. Returns a handle on one dart of this combinatorial tetrahedron.
\pre `dimension`\f$ \geq\f$ 2.

\sa `make_edge`
\sa `make_combinatorial_polygon`
\sa `make_combinatorial_hexahedron`
*/
Dart_handle make_combinatorial_tetrahedron();

/*!
Creates an isolated edge (two darts linked by \f$ \beta_2\f$) and adds it in the combinatorial map. Returns a handle on one dart of this edge.
\pre `dimension`\f$ \geq\f$ 2.

\sa `make_combinatorial_polygon`
\sa `make_combinatorial_tetrahedron`
\sa `make_combinatorial_hexahedron`
*/
Dart_handle make_edge();

/// @}

/// \name Operations
/// @{

/*!
Inserts a 0-cell in the 1-cell containing `dh`. Returns a handle on one dart belonging to the new 0-cell.
\pre \ref CombinatorialMap::dimension "dimension"\f$ \geq\f$ 1 and
   `*dh`\f$ \in\f$\ref CombinatorialMap::darts "darts()".

See example in \cgalFigureRef{figinsertvertex}.

If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`, if 1-attributes are non `void`, \ref CellAttribute::On_split "Attribute_type<1>::type::On_split"(<I>a</I>,<I>a'</I>) is called, with <I>a</I> the original 1-attribute associated with <I>dh</I> and <I>a'</I> the new 1-attribute created during the operation. If set, the dynamic onsplit function of 1-attributes is also called on <I>a</I> and <I>a'</I>.

\cgalAdvancedBegin
If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are not updated; thus the combinatorial map can be no more valid after this operation.
\cgalAdvancedEnd

\sa `insert_cell_0_in_cell_2`
\sa `insert_cell_1_in_cell_2`
\sa `insert_dangling_cell_1_in_cell_2`
\sa `insert_cell_2_in_cell_3<InputIterator>`
\sa `remove_cell<i>`
*/
Dart_handle insert_cell_0_in_cell_1(Dart_handle dh);

/*!
Inserts a 0-cell in the 2-cell containing `dh`. The 2-cell is split in triangles, one for each initial edge of the facet. Returns a handle on one dart belonging to the new 0-cell.
\pre \ref CombinatorialMap::dimension "dimension"\f$ \geq\f$ 2 and `*dh`\f$ \in\f$\ref CombinatorialMap::darts "darts()".

See example in \cgalFigureRef{figtriangulate}.

If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`, if 2-attributes are non `void`, \ref CellAttribute::On_split "Attribute_type<2>::type::On_split"(<I>a</I>,<I>a'</I>) is called, with <I>a</I> the original 2-attribute associated with `dh` and <I>a'</I> each new 2-attribute created during the operation. If set, the dynamic onsplit function of 2-attributes is also called on <I>a</I> and <I>a'</I>.

\cgalAdvancedBegin
If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are not updated; thus the combinatorial map can be no more valid after this operation.
\cgalAdvancedEnd

\sa `insert_cell_0_in_cell_2`
\sa `insert_cell_1_in_cell_2`
\sa `insert_dangling_cell_1_in_cell_2`
\sa `insert_cell_2_in_cell_3<InputIterator>`
\sa `remove_cell<i>`
*/
Dart_handle insert_cell_0_in_cell_2(Dart_handle dh); 

/*!
Inserts a 1-cell in the 2-cell containing `dh1` and `dh2`. Returns \f$ \beta_0\f$(`dh1`), a handle on one dart belonging to the new 1-cell.
\pre `is_insertable_cell_1_in_cell_2(dh1,dh2)`.

See example in \cgalFigureRef{figinsertedge}.

If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`, if 2-attributes are non `void`, \ref CellAttribute::On_split "Attribute_type<2>::type::On_split"(<I>a</I>,<I>a'</I>) is called, with <I>a</I> the original 2-attribute associated with `dh` and <I>a'</I> the new 2-attribute created during the operation. If set, the dynamic onsplit function of 2-attributes is also called on <I>a</I> and <I>a'</I>.

\cgalAdvancedBegin
If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are not updated; thus the combinatorial map can be no more valid after this operation.
\cgalAdvancedEnd

\sa `is_insertable_cell_1_in_cell_2`
\sa `insert_cell_0_in_cell_1`
\sa `insert_cell_0_in_cell_2`
\sa `insert_dangling_cell_1_in_cell_2`
\sa `insert_cell_2_in_cell_3<InputIterator>`
\sa `remove_cell<i>`
*/
Dart_handle insert_cell_1_in_cell_2(Dart_handle dh1, Dart_handle dh2);

/*!
Inserts a 2-cell along the path of 1-cells containing darts given by the range `[afirst,alast)`. Returns a handle on one dart belonging to the new 2-cell.
\pre `is_insertable_cell_2_in_cell_3(afirst,alast)`.

See example in \cgalFigureRef{figinsertface}.

If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`, if 3-attributes are non `void`, \ref CellAttribute::On_split "Attribute_type<3>::type::On_split"(<I>a</I>,<I>a'</I>) is called, with <I>a</I> the original 3-attribute associated with `dh` and <I>a'</I> the new 3-attribute created during the operation. If set, the dynamic onsplit function of 3-attributes is also called on <I>a</I> and <I>a'</I>.

\cgalAdvancedBegin
If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are not updated; thus the combinatorial map can be no more valid after this operation.
\cgalAdvancedEnd

\sa `is_insertable_cell_2_in_cell_3<InputIterator>`
\sa `insert_cell_0_in_cell_1`
\sa `insert_cell_0_in_cell_2`
\sa `insert_cell_1_in_cell_2`
\sa `insert_dangling_cell_1_in_cell_2`
\sa `remove_cell<i>`
*/
template <class InputIterator>
Dart_handle insert_cell_2_in_cell_3(InputIterator afirst, InputIterator alast);  

/*!
Inserts a 1-cell in a the 2-cell containing `dh`, the 1-cell being attached only by one of its extremity to the 0-cell containing `dh`. Returns a handle on the dart belonging to the new 1-cell and to the new 0-cell.
\pre \ref CombinatorialMap::dimension "dimension"\f$ \geq\f$ 2 and `*dh`\f$ \in\f$\ref CombinatorialMap::darts "darts()".

See example in \cgalFigureRef{figinsertedge}.

\cgalAdvancedBegin
If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are not updated; thus the combinatorial map can be no more valid after this operation.
\cgalAdvancedEnd

\sa `insert_cell_0_in_cell_1`
\sa `insert_cell_0_in_cell_2`
\sa `insert_cell_1_in_cell_2`
\sa `insert_cell_2_in_cell_3<InputIterator>`
\sa `remove_cell<i>`

*/
Dart_handle insert_dangling_cell_1_in_cell_2(Dart_handle dh);

/*!
Returns true iff it is possible to insert a 1-cell in the combinatorial map between `dh1` and `dh2`.

This is possible if `dh1`\f$ \neq\f$`dh2` and `dh1`\f$ \in\f$\f$ \langle{}\f$\f$ \beta_1\f$\f$ \rangle{}\f$(`dh2`).
\pre \ref CombinatorialMap::dimension "dimension"\f$ \geq\f$ 2, `*dh1`\f$ \in\f$\ref CombinatorialMap::darts "darts()", and `*dh2`\f$ \in\f$\ref CombinatorialMap::darts "darts()".

\sa `insert_cell_1_in_cell_2`
\sa `is_insertable_cell_2_in_cell_3<InputIterator>`

*/
bool is_insertable_cell_1_in_cell_2(Dart_const_handle dh1, Dart_const_handle dh2);

/*!
Returns true iff it is possible to insert a 2-cell in the combinatorial map along the path of darts given by the range `[afirst,alast)`. The 2-cell can be inserted iff each couple of consecutive darts of the path <I>a1</I> and <I>a2</I> belong to the same vertex and the same volume, and if the path is closed.
\pre \ref CombinatorialMap::dimension "dimension"\f$ \geq\f$ 3.

\sa `insert_cell_2_in_cell_3<InputIterator>`
\sa `is_insertable_cell_1_in_cell_2`

*/
template <class InputIterator>
bool is_insertable_cell_2_in_cell_3(InputIterator afirst, InputIterator alast);

/*!
Returns true iff the <I>i</I>-cell containing `dh` can be removed.

An <I>i</I>-cell can be removed if `i`==\ref CombinatorialMap::dimension "dimension" or if `i`==\ref CombinatorialMap::dimension "dimension"-1 or if `i`\f$ <\f$\ref CombinatorialMap::dimension "dimension"-1 and the <I>i</I>-cell containing `dh` is incident to at most two (<I>i+1</I>)-cells.
\pre 0\f$ \leq\f$`i`\f$ \leq\f$\ref CombinatorialMap::dimension "dimension" and `*dh`\f$ \in\f$\ref CombinatorialMap::darts "darts()".

\sa `remove_cell<i>`
*/
template <unsigned int i>
bool is_removable(Dart_const_handle dh);

/*!
Removes the <I>i</I>-cell containing `dh`. Returns the number of darts removed from the combinatorial map.
\pre `is_removable<i>(dh)`.

See examples in \cgalFigureRef{figinsertvertex}, \cgalFigureRef{figinsertedge} and \cgalFigureRef{figinsertface}.

If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`, if `i`\f$ <\f$\ref CombinatorialMap::dimension "dimension", and <I>i+1</I>-attributes are non `void`, and if there are two distinct (<I>i+1</I>)-cells around dart `dh`, \ref CellAttribute::On_merge "Attribute_type<i+1>::type::On_merge"(<I>a1</I>,<I>a2</I>) is called, with <I>a1</I> the (<I>i+1</I>)-attribute associated to `dh`, and <I>a2</I> the (<I>i+1</I>)-attribute associated to \f$ \beta_{i+1}\f$(<I>dh</I>). If set, the dynamic onmerge function of <I>i+1</I>-attributes is also called on <I>a1</I> and <I>a2</I>.

If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`, if a <I>j</I>-cell is disconnected in two <I>j</I>-cells during the operation, and if <I>j</I>-attributes are non void, \ref CellAttribute::On_split "Attribute_type<j>::type::On_split"(<I>a</I>,<I>a'</I>) is called with <I>a</I> the original <I>j</I>-attribute and <I>a'</I> the new <I>j</I>-attribute created due to the disconnection. If set, the dynamic onsplit function of <i>j</i>-attributes is also called on <I>a</I> and <I>a'</I>.

\cgalAdvancedBegin
If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are not updated; thus the combinatorial map can be no more valid after this operation.
\cgalAdvancedEnd

\sa `is_removable<i>`
\sa `insert_cell_0_in_cell_1`
\sa `insert_cell_0_in_cell_2`
\sa `insert_cell_1_in_cell_2`
\sa `insert_dangling_cell_1_in_cell_2`
\sa `insert_cell_2_in_cell_3<InputIterator>`
*/
template <unsigned int i>
size_type remove_cell(Dart_handle dh);
  
} /* namespace CGAL */

/// @}
  
}; /* end CombinatorialMap */
//@}

