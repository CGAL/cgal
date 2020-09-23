
/*!
\ingroup PkgLinearCellComplexConcepts
\cgalConcept

The concept `LinearCellComplex` represents a linear cell complex in dimension `d`, in an ambient space of dimension `d2`. This is a model of the concept of `GenericMap` adding a requirement to ensure that each vertex of the map is associated with a model of `CellAttributeWithPoint`.

\cgalRefines `GenericMap`

\cgalHasModel \link CGAL::Linear_cell_complex_for_combinatorial_map `CGAL::Linear_cell_complex_for_combinatorial_map<d,d2,LCCTraits,Items,Alloc>`\endlink
\cgalHasModel \link CGAL::Linear_cell_complex_for_generalized_map `CGAL::Linear_cell_complex_for_generalized_map<d,d2,LCCTraits,Items,Alloc>`\endlink

\sa `LinearCellComplexItems`
\sa `LinearCellComplexTraits`

*/
class LinearCellComplex
{
public:

/// \name Creation
/// @{

/*!
%Default constructor creating an empty linear cell complex.
*/
LinearCellComplex();

/// @}

/// \name Constants
/// @{

/*!
Ambient dimension, must be > 1.
*/
static unsigned int ambient_dimension;

/// @}

/// \name Types
/// @{

/// Traits class, a model of the `LinearCellComplexTraits` concept.
typedef unspecified_type Traits;

/*!
*/
typedef Traits::FT FT;

/*!
*/
typedef Traits::Point Point;

/*!
*/
typedef Traits::Vector Vector;

/*!
%Type of 0-attributes, a model of `CellAttributeWithPoint` concept
(a shortcut for \link GenericMap::Attribute_type `Attribute_type<0>::type` \endlink).
*/
typedef unspecified_type Vertex_attribute;

/*!
%Handle through 0-attributes
(a shortcut for \link GenericMap::Attribute_handle `Attribute_handle<0>::type` \endlink).
*/
typedef unspecified_type Vertex_attribute_handle;

/*!
Const handle through 0-attributes
(a shortcut for \link GenericMap::Attribute_const_handle `Attribute_const_handle<0>::type` \endlink).
*/
typedef unspecified_type Vertex_attribute_const_handle;

/*!
%Range of all the 0-attributes, a model of the `Range` concept
(a shortcut for \link GenericMap::Attribute_range `Attribute_range<0>::type` \endlink).
%Iterator inner type is bidirectional iterator and value type is \link LinearCellComplex::Vertex_attribute `Vertex_attribute`\endlink.
*/
typedef unspecified_type Vertex_attribute_range;

/*!
%Const range of all the 0-attributes, a model of the `ConstRange` concept
a shortcut for \link GenericMap::Attribute_const_range `Attribute_const_range<0>::type` \endlink).
%Iterator inner type is bidirectional iterator and value type is \link LinearCellComplex::Vertex_attribute `Vertex_attribute`\endlink.
*/
typedef unspecified_type Vertex_attribute_const_range;

/// @}

/// \name %Range Access Member Functions
/// @{

/*!
Returns a range of all the 0-attributes in this linear cell complex
(a shortcut for \link GenericMap::attributes `attributes<0>()`\endlink).
*/
Vertex_attribute_range& vertex_attributes();

/*!
Returns a const range of all the 0-attributes in this linear cell complex
(a shortcut for \link GenericMap::attributes `attributes<0>() const`\endlink).
*/
Vertex_attribute_const_range& vertex_attributes() const;

/// @}

/// \name Access Member Functions
/// @{

/*!
Returns true iff this linear cell complex is valid.

A linear cell complex `lcc` is valid if it is a valid generic map (cf. `GenericMap::is_valid()`), and if for each dart handle <I>dh</I> such that `*dh` \f$ \in \f$ \link GenericMap::darts `darts()`\endlink: \link GenericMap::attribute `dh->attribute<0>()`\endlink`!=NULL`.
*/
bool is_valid() const;

/*!
Returns the number of 0-attributes in this linear cell complex
(a shortcut for \link GenericMap::number_of_attributes `number_of_attributes<0>()`\endlink).
*/
size_type number_of_vertex_attributes() const;

/*!
Returns the 0-attribute associated with `dh`.
*/
Vertex_attribute_handle vertex_attribute(Dart_handle dh);

/*!
Returns the 0-attribute associated with `dh`, when `dh` is const.
*/
Vertex_attribute_const_handle vertex_attribute(Dart_const_handle dh);

/*!
Returns the point in the 0-attribute `vh`.
*/
Point& point_of_vertex_attribute(Vertex_attribute_handle vh);

/*!
Returns the point in the 0-attribute `vh`, when `vh` is const.
*/
const Point& point_of_vertex_attribute(Vertex_attribute_const_handle vh) const;

/*!
Returns the point in the 0-attribute associated with `dh`.
*/
Point& point(Dart_handle dh);

/*!
Returns the point in the 0-attribute associated with `dh`, when `dh` is const.
*/
const Point& point(Dart_const_handle dh) const;

/// @}

/// \name Modifiers
/// @{

/*!
Creates a new dart in this linear cell complex, sets its associated 0-attribute to `vh` and returns the corresponding handle.
\pre `*vh` \f$ \in \f$ `vertex_attributes()`.
*/
Dart_handle create_dart(Vertex_attribute_handle vh);

/*!
Creates a new dart in this linear cell complex, creates a new 0-attribute initialized with `apoint`, sets the associated 0-attribute of the new dart to this new 0-attribute, and returns the corresponding handle.
*/
Dart_handle create_dart(const Point& apoint);

/*!
Creates a new 0-attribute in this linear cell complex, and returns the corresponding handle (a shortcut for \link GenericMap::create_attribute `create_attribute<0>(t1)`\endlink). Calls the constructor of \link LinearCellComplex::Vertex_attribute `Vertex_attribute`\endlink having `T1` as parameter. Overloads of this member function are defined that take from zero to nine arguments. With zero argument, `create_vertex_attribute()` creates a new 0-attribute by using the default constructor.
*/
template<typename T1> Vertex_attribute_handle
create_vertex_attribute(T1 t1);

/*!
Removes the 0-attribute pointed to by `vh` from this linear cell complex (a shortcut for \link GenericMap::erase_attribute `erase_attribute<0>(vh)`\endlink).
\pre `*vh` \f$ \in \f$ `vertex_attributes()`.

*/
void erase_vertex_attribute(Vertex_attribute_handle vh);

/*!
Associates the 0-attribute of all the darts of the 0-cell containing `dh` to `vh` (a shortcut for \link GenericMap::set_attribute `set_attribute<0>(dh,vh)`\endlink).
\pre `*dh` \f$ \in \f$ \link GenericMap::darts `darts()`\endlink and `*vh` \f$ \in \f$ `vertex_attributes()`.

*/
void set_vertex_attribute(Dart_handle dh, Vertex_attribute_handle vh);

/// @}

/// \name Attributes management
/// @{
/*!
Correct the invalid attributes of the linear cell complex. We can have invalid attribute either if we have called \link GenericMap::set_automatic_attributes_management `set_automatic_attributes_management(false)`\endlink before to use some modification operations or if we have modified the combinatorial map by using low level operations.

The validation process of a linear cell complex validates its generic map (cf. \link GenericMap::correct_invalid_attributes `correct_invalid_attributes()`\endlink), and for each dart `d` having no vertex attribute, a new vertex attribute is created, with its Point initialized to `CGAL::Origin`, and all the darts of the 0-cell containing `d` are linked with the new attribute.
*/
void correct_invalid_attributes();

/// @}

/// \name Operations
/// @{

/*!
Returns the barycenter of the <I>i</I>-cell containing `dh`.
\pre 1 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink and `*dh` \f$ \in \f$ \link GenericMap::darts `darts()`\endlink.

*/
template<unsigned int i> Point barycenter(Dart_const_handle dh) const;

/*!
Inserts a point, copy of `p`, in the <I>i</I>-cell containing `dh`.
Returns a handle on one dart of this cell.
\pre <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink \f$ \leq \f$ 2 and `*dh` \f$ \in \f$ \link GenericMap::darts `darts()`\endlink.

If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`, if <I>i</I>-attributes are non void, \link CellAttribute::On_split `Attribute_type<i>::type::On_split`\endlink(<I>a</I>,<I>a'</I>) is called, with <I>a</I> the original <I>i</I>-attribute associated with <I>dh</I> and <I>a'</I> each new <I>i</I>-attribute created during the operation.

\cgalAdvancedBegin
If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are not updated; thus the combinatorial map can be no more valid after this operation.
\cgalAdvancedEnd

*/
template <unsigned int i> Dart_handle insert_point_in_cell(Dart_handle dh, Point p);

/*!
Inserts a point in the barycenter of the <I>i</I>-cell containing `dh`.
Returns a handle on one dart of this cell.
\pre <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink \f$ \leq \f$ 2 and `*dh` \f$ \in \f$ \link GenericMap::darts `darts()`\endlink.

If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`, if <I>i</I>-attributes are non void, \link CellAttribute::On_split `Attribute_type<i>::type::On_split`\endlink(<I>a</I>,<I>a'</I>) is called, with <I>a</I> the original <I>i</I>-attribute associated with <I>dh</I> and <I>a'</I> each new <I>i</I>-attribute created during the operation.

\cgalAdvancedBegin
If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are not updated; thus the combinatorial map can be no more valid after this operation.
\cgalAdvancedEnd

*/
template <unsigned int i> Dart_handle insert_barycenter_in_cell(Dart_handle dh);

/*!
Inserts a 1-cell in the 2-cell containing `dh`, the 1-cell being attached only by one of its vertex to the 0-cell containing `dh`. The second vertex is associated with a new 0-attribute containing a copy of `p` as point. Returns a handle on one dart belonging to the new 0-cell.
\pre 2 \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink and `*dh` \f$ \in \f$ \link GenericMap::darts `darts()`\endlink.

\cgalAdvancedBegin
If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are not updated; thus the combinatorial map can be no more valid after this operation.
\cgalAdvancedEnd

*/
Dart_handle insert_dangling_cell_1_in_cell_2(Dart_handle dh, Point p);

/// @}

/// \name Constructions
/// @{

/*!
Creates an isolated segment in this linear cell complex (two darts linked by \f$ \beta_2\f$) having `p0`, `p1` as points.
Returns a handle on the dart associated with `p0`.
\pre \link GenericMap::dimension `dimension`\endlink \f$ \geq\f$ 2.

\image html lcc_make_segment.svg "Example of r=lcc.make_segment(p0,p1), left for combinatorial map as combinatorial data-structure, right for generalized maps."
\image latex lcc_make_segment.svg "Example of r=lcc.make_segment(p0,p1), left for combinatorial map as combinatorial data-structure, right for generalized maps."
*/
Dart_handle make_segment(const Point& p0, const Point& p1);

/*!
Creates an isolated triangle in this linear cell complex having `p0`, `p1`, `p2` as points.
Returns a handle on the dart associated with `p0` and with edge [`p0`,`p1`].
\pre \link GenericMap::dimension `dimension`\endlink \f$ \geq\f$ 1.

\image html lcc_make_triangle.svg "Example of r=lcc.make_triangle(p0,p1,p2), left for combinatorial map as combinatorial data-structure, right for generalized maps."
\image latex lcc_make_triangle.svg "Example of r=lcc.make_triangle(p0,p1,p2), left for combinatorial map as combinatorial data-structure, right for generalized maps."
*/
Dart_handle make_triangle(const Point& p0, const Point& p1, const Point& p2);

/*!
Creates an isolated quadrangle in this linear cell complex having `p0`, `p1`, `p2`, `p3` as points.
Returns a handle on the dart associated with `p0` and with edge [`p0`,`p1`].
\pre \link GenericMap::dimension `dimension`\endlink \f$ \geq\f$ 1.

\image html lcc_make_quadrilateral.svg "Example of r=lcc.make_quadrangle(p0,p1,p2,p3), left for combinatorial map as combinatorial data-structure, right for generalized maps."
\image latex lcc_make_quadrilateral.svg "Example of r=lcc.make_quadrangle(p0,p1,p2,p3), left for combinatorial map as combinatorial data-structure, right for generalized maps."
*/
Dart_handle make_quadrangle(const Point& p0,const Point& p1,const Point& p2,const Point& p3);

/*!
Creates an isolated tetrahedron in this linear cell complex having `p0`, `p1`,`p2`,`p3` as points.
Returns a handle on the dart associated with `p0`, with edge [`p0`,`p1`] and belonging to the 2-cell having `p0`, `p1`, `p2` as points.
\pre \link GenericMap::dimension `dimension`\endlink \f$ \geq\f$ 2.

\image html lcc_make_tetrahedron.svg "Example of r=lcc.make_tetrahedron(p0,p1,p2,p3), left for combinatorial map as combinatorial data-structure, right for generalized maps."
\image latex lcc_make_tetrahedron.svg "Example of r=lcc.make_tetrahedron(p0,p1,p2,p3), left for combinatorial map as combinatorial data-structure, right for generalized maps."
*/
Dart_handle make_tetrahedron(const Point& p0,const Point& p1,const Point& p2,const Point& p3);

/*!
Creates an isolated hexahedron in this linear cell complex having `p0`, `p1`, `p2`, `p3`, `p4`, `p5`, `p6`, `p7` as points.
Returns a handle on the dart associated with `p0`, with edge [`p0`,`p5`] and belonging to the 2-cell having `p0`, `p5`, `p6`, `p1` as points.
\pre \link GenericMap::dimension `dimension`\endlink \f$ \geq \f$ 2.

\image html lcc_make_hexahedron.svg "Example of r=lcc.make_hexahedron(p0,p1,p2,p3,p4,p5,p6,p7), left for combinatorial map as combinatorial data-structure, right for generalized maps."
\image latex lcc_make_hexahedron.svg "Example of r=lcc.make_hexahedron(p0,p1,p2,p3,p4,p5,p6,p7), left for combinatorial map as combinatorial data-structure, right for generalized maps."
*/
Dart_handle make_hexahedron(const Point& p0,const Point& p1,const Point& p2,
                            const Point& p3,const Point& p4,const Point& p5,
                            const Point& p6,const Point& p7);
/// @}

}; /* end LinearCellComplex */
