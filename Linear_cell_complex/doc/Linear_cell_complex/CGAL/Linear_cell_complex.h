
namespace CGAL {

/*!
\ingroup PkgLinearCellComplexClasses

The class `Linear_cell_complex` represents a linear cell complex in dimension `d`,
in an ambient space of dimension `d2`. This is a model of the concept of
`CombinatorialMap` adding a requirement to ensure that
each vertex of the map is associated with a
model of `CellAttributeWithPoint`.

\cgalModels `CombinatorialMap`

\tparam d an integer for the dimension of the combinatorial map,
\tparam d2 an integer for the dimension of the ambient space,
\tparam LCCTraits must be a model of the `LinearCellComplexTraits` concept,
satisfying \ref LinearCellComplexTraits::ambient_dimension "LCCTraits::ambient_dimension"`==d2`,
\tparam Items must be a model of the `LinearCellComplexItems` concept,
\tparam Alloc has to match the standard allocator requirements.

There are four default template arguments:
`d2` is equal to `d`,
`LCCTraits` is equal to `CGAL::Linear_cell_complex_traits<d2>`,
`Items` is equal to `CGAL::Linear_cell_complex_min_items<d>` and
`Alloc` is `CGAL_ALLOCATOR(int)`.

\cgalAdvancedBegin
Note that there is an additional, and undocumented, template
parameter `CMap` for
`Linear_cell_complex<d,d2,LCCTraits,Items,Alloc,CMap>` allowing to
inherit from any model of the `CombinatorialMap` concept.
\cgalAdvancedEnd

\sa `CombinatorialMap`
\sa `CGAL::Combinatorial_map<d,Items,Alloc>`
\sa `Dart`
\sa `LinearCellComplexItems`
\sa `CGAL::Linear_cell_complex_min_items<d>`
\sa `LinearCellComplexTraits`
\sa `CGAL::Linear_cell_complex_traits<d,K>`

\deprecated Since \cgal 4.4, `vertex_attribute` and `point` methods are no more static. You can define the `CGAL_CMAP_DEPRECATED` macro to keep the old behavior.

*/
template< typename d, typename d2, typename LCCTraits, typename Items, typename Alloc >
class Linear_cell_complex  : public Combinatorial_map<d,Items,Alloc>
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
static unsigned int ambient_dimension = d2;

/// @}

/// \name Types
/// @{

/*!

*/
typedef Linear_cell_complex<d,d2,LCCTraits,Items,Alloc> Self;

/*!
The type of dart, must satisfy \ref Dart::dimension "Dart::dimension"`==d`.
*/
typedef Items::Dart_wrapper<Self>::Dart Dart;

/*!

*/
typedef LCCTraits Traits;

/*!

*/
typedef Items Items;

/*!

*/
typedef Alloc Alloc;

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
(a shortcut for \link CombinatorialMap::Attribute_type `Attribute_type<0>::type` \endlink).
*/
typedef unspecified_type Vertex_attribute;

/*!
%Handle through 0-attributes
(a shortcut for \link CombinatorialMap::Attribute_handle `Attribute_handle<0>::type` \endlink).
*/
typedef unspecified_type Vertex_attribute_handle;

/*!
Const handle through 0-attributes
(a shortcut for \link CombinatorialMap::Attribute_const_handle `Attribute_const_handle<0>::type` \endlink).
*/
typedef unspecified_type Vertex_attribute_const_handle;

/*!
%Range of all the 0-attributes, a model of the `Range` concept
(a shortcut for \link CombinatorialMap::Attribute_range `Attribute_range<0>::type` \endlink).
%Iterator inner type is bidirectional iterator and value type is \ref Linear_cell_complex::Vertex_attribute "Vertex_attribute".
*/
typedef unspecified_type Vertex_attribute_range;

/*!
%Const range of all the 0-attributes, a model of the `ConstRange` concept
a shortcut for \link CombinatorialMap::Attribute_const_range `Attribute_const_range<0>::type` \endlink).
%Iterator inner type is bidirectional iterator and value type is \ref Linear_cell_complex::Vertex_attribute "Vertex_attribute".
*/
typedef unspecified_type Vertex_attribute_const_range;

/// @}

/// \name %Range Access Member Functions
/// @{

/*!
Returns a range of all the 0-attributes in this linear cell complex
(a shortcut for \ref CombinatorialMap::attributes "attributes<0>()").
*/
Vertex_attribute_range& vertex_attributes();

/*!
Returns a const range of all the 0-attributes in this linear cell complex
(a shortcut for \ref CombinatorialMap::attributes "attributes<0>() const").
*/
Vertex_attribute_const_range& vertex_attributes() const;

/// @}

/// \name Access Member Functions
/// @{

/*!
Returns true iff this linear cell complex is valid.

A linear cell complex `lcc` is valid if it is a valid combinatorial
map (cf. `CombinatorialMap::is_valid()`), and if for each dart handle <I>dh</I> such that
`*dh`\f$\in\f$\ref CombinatorialMap::darts "darts()": \ref Dart::attribute "dh->attribute<0>()"`!=NULL`.
*/
bool is_valid() const;

/*!
Returns the number of 0-attributes in this linear cell complex
(a shortcut for \ref CombinatorialMap::number_of_attributes "number_of_attributes<0>()").
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
const Point& point(Dart_const_handle dh);

/// @}

/// \name Modifiers
/// @{

/*!
Creates a new dart in this linear cell complex, sets its associated 0-attribute
to `vh` and returns the corresponding handle.
\pre `*vh`\f$ \in\f$`vertex_attributes()`.
*/
Dart_handle create_dart(Vertex_attribute_handle vh);

/*!
Creates a new dart in this linear cell complex, creates a new 0-attribute
initialized with `apoint`, sets the associated 0-attribute
of the new dart to this new 0-attribute,
and returns the corresponding handle.
*/
Dart_handle create_dart(const Point& apoint);

/*!
Creates a new 0-attribute in this linear cell complex,
and returns the corresponding handle (a shortcut for
\ref CombinatorialMap::create_attribute "create_attribute<0>(t1)").
Calls the constructor of
\ref Linear_cell_complex::Vertex_attribute "Vertex_attribute" having `T1` as parameter. Overloads of this
member function are defined that take from zero to nine arguments.
With zero argument, `create_vertex_attribute()` creates a new
0-attribute by using the default constructor.
*/
template<typename T1> Vertex_attribute_handle
create_vertex_attribute(T1 t1);

/*!
Removes the 0-attribute pointed to by `vh` from this linear cell complex
(a shortcut for \ref CombinatorialMap::erase_attribute "erase_attribute<0>(vh)").
\pre `*vh`\f$ \in\f$`vertex_attributes()`.

*/
void erase_vertex_attribute(Vertex_attribute_handle vh);

/*!
Associates the 0-attribute of all the darts of the 0-cell
containing `dh` to `vh`
(a shortcut for \ref CombinatorialMap::set_attribute "set_attribute<0>(dh,vh)").
\pre `*dh`\f$ \in\f$\ref CombinatorialMap::darts "darts()" and `*vh`\f$ \in\f$`vertex_attributes()`.

*/
void set_vertex_attribute(Dart_handle dh, Vertex_attribute_handle vh);

/// @}

/// \name Attributes management
/// @{
/*!
Correct the invalid attributes of the linear cell complex.
We can have invalid attribute either if we have called \link CombinatorialMap::set_automatic_attributes_management `set_automatic_attributes_management(false)`\endlink before to use some modification operations or if we have modified the combinatorial map by using low level operations.

The validation process of a linear cell complex validates its combinatorial map (cf. \link CombinatorialMap::correct_invalid_attributes `correct_invalid_attributes()`\endlink), and for each dart `d` having no vertex attribute, a new vertex attribute is created, with its Point initialized to `CGAL::Origin`, and all the darts of the 0-cell containing `d` are linked with the new attribute.
*/
void correct_invalid_attributes();

/// @}

/// \name Operations
/// @{

/*!
Returns the barycenter of the <I>i</I>-cell containing `dh`.
\pre 1\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension" and `*dh`\f$ \in\f$\ref CombinatorialMap::darts "darts()".

*/
template<unsigned int i> Point barycenter(Dart_const_handle dh) const;

/*!
Inserts a point, copy of `p`, in the <I>i</I>-cell containing `dh`.
Returns a handle on one dart of this cell.
\pre <I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension"\f$ \leq\f$ 2 and `*dh`\f$ \in\f$\ref CombinatorialMap::darts "darts()".

If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`,
if <I>i</I>-attributes are non void,
 \ref CellAttribute::On_split "Attribute_type<i>::type::On_split"(<I>a</I>,<I>a'</I>) is called,
with <I>a</I> the original <I>i</I>-attribute associated
with <I>dh</I> and <I>a'</I> each new <I>i</I>-attribute created during the operation.

\cgalAdvancedBegin
If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are
not updated; thus the combinatorial map can be no more valid after this operation.
\cgalAdvancedEnd

*/
template <unsigned int i> Dart_handle insert_point_in_cell(Dart_handle dh, Point p);

/*!
Inserts a point in the barycenter of the <I>i</I>-cell containing `dh`.
Returns a handle on one dart of this cell.
\pre <I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "dimension"\f$ \leq\f$ 2 and `*dh`\f$ \in\f$\ref CombinatorialMap::darts "darts()".

If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`,
if <I>i</I>-attributes are non void,
\ref CellAttribute::On_split "Attribute_type<i>::type::On_split"(<I>a</I>,<I>a'</I>) is called,
with <I>a</I> the original <I>i</I>-attribute associated
with <I>dh</I> and <I>a'</I> each new <I>i</I>-attribute created during the operation.

\cgalAdvancedBegin
If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are
not updated; thus the combinatorial map can be no more valid after this operation.
\cgalAdvancedEnd

*/
template <unsigned int i> Dart_handle insert_barycenter_in_cell(Dart_handle dh);

/*!
Inserts a 1-cell in the 2-cell containing `dh`, the 1-cell
being attached only by one of its vertex to the 0-cell containing `dh`.
The second vertex is associated with a new 0-attribute containing a copy of
`p` as point. Returns a handle on one dart belonging to the new 0-cell.
\pre 2\f$ \leq\f$\ref CombinatorialMap::dimension "dimension" and `*dh`\f$ \in\f$\ref CombinatorialMap::darts "darts()".

\cgalAdvancedBegin
If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are
not updated; thus the combinatorial map can be no more valid after this operation.
\cgalAdvancedEnd

*/
Dart_handle insert_dangling_cell_1_in_cell_2(Dart_handle dh, Point p);

/// @}

/// \name Constructions
/// @{

/*!
Creates an isolated segment in this linear cell complex (two darts linked by \f$ \beta_2\f$)
having `p0`, `p1` as points.
Returns a handle on the dart associated with `p0`.
\pre \ref CombinatorialMap::dimension "dimension"\f$ \geq\f$ 2.

\image html make_segment.png "Example of r=lcc.make_segment(p0,p1)."
\image latex make_segment.png "Example of r=lcc.make_segment(p0,p1)."
*/
Dart_handle make_segment(const Point& p0, const Point& p1);

/*!
Creates an isolated triangle in this linear cell complex having `p0`, `p1`, `p2` as points.
Returns a handle on the dart associated with `p0`.
\pre \ref CombinatorialMap::dimension "dimension"\f$ \geq\f$ 1.

\image html make_triangle.png "Example of r=lcc.make_triangle(p0,p1,p2)."
\image latex make_triangle.png "Example of r=lcc.make_triangle(p0,p1,p2)."
*/
Dart_handle make_triangle(const Point& p0, const Point& p1, const Point& p2);

/*!
Creates an isolated quadrangle in this linear cell complex having `p0`, `p1`,
`p2`, `p3` as points.
Returns a handle on the dart associated with `p0`.
\pre \ref CombinatorialMap::dimension "dimension"\f$ \geq\f$ 1.

\image html make_quadrilateral.png "Example of r=lcc.make_quadrangle(p0,p1,p2,p3)."
\image latex make_quadrilateral.png "Example of r=lcc.make_quadrangle(p0,p1,p2,p3)."
*/
Dart_handle make_quadrangle(const Point& p0,const Point& p1,const Point& p2,const Point& p3);

/*!
Creates an isolated tetrahedron in this linear cell complex having `p0`,
`p1`,`p2`,`p3` as points. Returns a handle on the dart
associated with `p0` and belonging to the 2-cell having
`p0`, `p1`, `p2` as points.
\pre \ref CombinatorialMap::dimension "dimension"\f$ \geq\f$ 2.

\image html make_tetrahedron.png "Example of r=lcc.make_tetrahedron(p0,p1,p2,p3)."
\image latex make_tetrahedron.png "Example of r=lcc.make_tetrahedron(p0,p1,p2,p3)."
*/
Dart_handle make_tetrahedron(const Point& p0,const Point& p1,const Point& p2,const Point& p3);

/*!
Creates an isolated hexahedron in this linear cell complex having `p0`, `p1`,
`p2`, `p3`, `p4`, `p5`, `p6`, `p7` as points.
Returns a handle on the dart associated with `p0` and
belonging to the 2-cell having `p0`, `p5`, `p6`, `p1`
as points.
\pre \ref CombinatorialMap::dimension "dimension" \f$ \geq \f$ 2.

\image html make_hexahedron.png "Example of r=lcc.make_hexahedron(p0,p1,p2,p3,p4,p5,p6,p7)."
\image latex make_hexahedron.png "Example of r=lcc.make_hexahedron(p0,p1,p2,p3,p4,p5,p6,p7)."
*/
Dart_handle make_hexahedron(const Point& p0,const Point& p1,const Point& p2,
                            const Point& p3,const Point& p4,const Point& p5,
                            const Point& p6,const Point& p7);
/// @}

}; /* end Linear_cell_complex */
} /* end namespace CGAL */
