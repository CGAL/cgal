/*!
\ingroup PkgCombinatorialMapsConcepts
\cgalConcept

The concept `CombinatorialMap` defines a <I>d</I>-dimensional combinatorial map.

\cgalRefines{GenericMap}

\cgalHasModelsBegin
\cgalHasModelsBare{\link CGAL::Combinatorial_map `CGAL::Combinatorial_map<d,Items,Alloc>`\endlink}
\cgalHasModelsEnd

For a combinatorial map, the function \link GenericMap::next `next`\endlink is equal to \f$ \beta_1\f$, \link GenericMap::previous `previous`\endlink is equal to \f$ \beta_0\f$  and the function \link GenericMap::opposite `opposite<i>`\endlink is equal to \f$ \beta_i\f$.
*/

class CombinatorialMap {
public:

/*!
The null dart descriptor constant.
A dart `d` is <I>i</I>-free if `beta(d, i)==null_dart_descriptor`.
Note that `*null_dart_descriptor`\f$ \notin\f$`darts()`.
*/
Dart_descriptor null_dart_descriptor;

/// \name Access Member Functions
/// @{

/*!
Returns \f$ \beta_j\f$(\f$ \beta_i\f$(`d`)).
Overloads of this member function are defined that take from one to nine integer as arguments.
For each function, betas are applied in the same order as their indices are given as parameters.

For example `beta(d,1)`=\f$ \beta_1\f$(`d`),
and `beta(d,1,2,3,0)`=\f$ \beta_0\f$(\f$ \beta_3\f$(\f$ \beta_2\f$(\f$ \beta_1\f$(`d`)))).
\pre 0 \f$ \leq\f$ <I>i</I> \f$ \leq\f$ \link GenericMap::dimension `dimension`\endlink,
  0 \f$ \leq\f$ <I>j</I> \f$ \leq\f$ \link GenericMap::dimension `dimension`\endlink
  and `d`\f$ \in\f$ \link GenericMap::darts `darts()`\endlink.
*/
Dart_descriptor beta(Dart_descriptor d, int i, int j);

/*!
Returns \f$ \beta_j\f$(\f$ \beta_i\f$(`d`)).
Overloads of this member function are defined that take from one to nine integer as arguments.
\pre 0 \f$ \leq\f$ <I>i</I> \f$ \leq\f$ \link GenericMap::dimension `dimension`\endlink,
     0 \f$ \leq\f$ <I>j</I> \f$ \leq\f$ \link GenericMap::dimension `dimension`\endlink
     and `d`\f$ \in\f$ \link GenericMap::darts `darts()`\endlink.

*/
Dart_const_descriptor beta(Dart_const_descriptor d, int i, int j) const;

/*!
Returns \f$ \beta_j\f$(\f$ \beta_i\f$(`d`)).
Overloads of this member function are defined that take from one to nine integer as template arguments.
For each function, betas are applied in the same order as their indices are given as template arguments.

For example `beta<1>(d)`=\f$ \beta_1\f$(`d`),
and `beta<1,2,3,0>(d)`=\f$ \beta_0\f$(\f$ \beta_3\f$(\f$ \beta_2\f$(\f$ \beta_1\f$(`d`)))).
\pre 0 \f$ \leq\f$ <I>i</I> \f$ \leq\f$ \link GenericMap::dimension `dimension`\endlink,
     0 \f$ \leq\f$ <I>j</I> \f$ \leq\f$ \link GenericMap::dimension `dimension`\endlink
  and `d`\f$ \in\f$ \link GenericMap::darts `darts()`\endlink.
*/
template<int i, int j>
Dart_descriptor beta(Dart_descriptor d);

/*!
Returns \f$ \beta_j\f$(\f$ \beta_i\f$(`d`)).
Overloads of this member function are defined that take from one to nine integer as template arguments.
\pre 0 \f$ \leq\f$ <I>i</I> \f$ \leq\f$ \link GenericMap::dimension `dimension`\endlink,
     0 \f$ \leq\f$ <I>j</I> \f$ \leq\f$ \link GenericMap::dimension `dimension`\endlink
     and `d`\f$ \in\f$ \link GenericMap::darts `darts()`\endlink.

*/
template<int i, int j>
Dart_const_descriptor beta(Dart_const_descriptor d) const;

/*!
Returns a descriptor to a dart belonging to the same edge than dart `d`, and not to the same vertex. `null_descriptor` if such a dart does not exist.
*/
Dart_descriptor opposite(Dart_descriptor d);

/*!
Returns a const descriptor to a dart belonging to the same edge than dart `d`, and not to the same vertex, when the dart is const. `null_descriptor` if such a dart does not exist.
*/
Dart_const_descriptor opposite(Dart_const_descriptor d) const;

/*!
Returns true iff the combinatorial map is valid.

A combinatorial map is valid (see Sections \ref sec_definition and \ref sseccombimapvalidity) if for all its darts `d` \f$ \in\f$ `darts()`:

- `d` is 0-free, or \f$ \beta_1(\beta_0(d))=d\f$;
- `d` is 1-free, or \f$ \beta_0(\beta_1(d))=d\f$;
- \f$ \forall\f$<I>i</I>, 2 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink: `d` is <i>i</i>-free, or \f$ \beta_i(\beta_i(d))=d\f$;

- \f$ \forall\f$<I>i</I>, \f$ \forall\f$ <I>j</I>, such that 0 \f$ \leq\f$ <I>i</I> \f$ \leq\f$ \link GenericMap::dimension `dimension`\endlink-2
and <I>i</I>+2 \f$ \leq\f$ <I>j</I> \f$ \leq\f$ \link GenericMap::dimension `dimension`\endlink:
\f$ \beta_j(\beta_i(d))=\varnothing\f$ or ; \f$ \beta_j(\beta_i(\beta_j(\beta_i(d))))=d\f$;

- \f$ \forall\f$<I>i</I>, 0 \f$ \leq\f$ <I>i</I> \f$ \leq\f$ \link GenericMap::dimension `dimension`\endlink
  such that <I>i</I>-attributes are non void:
  + \f$ \forall\f$<I>d2</I> in the same <I>i</I>-cell than <I>d</I>: <I>d</I> and <I>d2</I> have the same <I>i</I>-attribute;
  + \f$ \forall\f$<I>d2</I>  in a different <I>i</I>-cell than <I>d</I>: <I>d</I> and <I>d2</I> have different <I>i</I>-attributes.
*/
bool is_valid() const;

/// @}

/// \name Modifiers
/// @{

/*!
Links `d1` and `d2` by \f$ \beta_i\f$.
The combinatorial map can be no longer valid after this operation. If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`, non void attributes of `d1` and `d2` are updated: if one dart has an attribute and the second dart not, the non null attribute is associated to the dart having a null attribute.
If both darts have an attribute, the attribute of `d1` is associated to `d2`.
\pre 0 \f$ \leq\f$ <I>i</I> \f$ \leq\f$ \link GenericMap::dimension `dimension`\endlink,
    `d1`\f$ \in\f$ \link GenericMap::darts `darts()`\endlink, `d2`\f$ \in\f$ \link GenericMap::darts `darts()`\endlink and (<I>i</I>\f$ <\f$ 2 or `d1`\f$ \neq\f$`d2`).
*/
template <unsigned int i> void link_beta(Dart_descriptor d1, Dart_descriptor d2);

/*!
Unlinks `d` and \f$ \beta_i\f$(`d`) by \f$ \beta_i\f$.
The combinatorial map can be no longer valid after this operation.
Attributes of `d` and \f$ \beta_i\f$(`d`) are not modified.
\pre 0 \f$ \leq\f$ <I>i</I> \f$ \leq\f$ \link GenericMap::dimension `dimension`\endlink,
     `d`\f$ \in\f$ \link GenericMap::darts `darts()`\endlink, and `d` is not <I>i</I>-free.
*/
template <unsigned int i> void unlink_beta(Dart_descriptor d);

/// @}

/// \name Operations
/// @{
/*!
Returns true iff `d1` can be <I>i</I>-sewn with `d2` by keeping the generic map valid.

This is true if there is a bijection <I>f</I> between all the darts of the orbit <I>D1</I>=\f$ \langle{}\f$\f$ \beta_1\f$,\f$ \ldots\f$,\f$ \beta_{i-2}\f$,\f$ \beta_{i+2}\f$,\f$ \ldots\f$,\f$ \beta_d\f$\f$ \rangle{}\f$(<I>d1</I>) and <I>D2</I>=\f$ \langle{}\f$\f$ \beta_1\f$,\f$ \ldots\f$,\f$ \beta_{i-2}\f$,\f$ \beta_{i+2}\f$,\f$ \ldots\f$,\f$ \beta_d\f$\f$ \rangle{}\f$(<I>d2</I>) satisfying:
- <I>f</I>(<I>d1</I>)=<I>d2</I>,
- \f$ \forall \f$ <I>e</I>\f$ \in\f$ <I>D1</I>, \f$ \forall \f$ <I>j</I>\f$ \in\f$ {1,\f$ \ldots\f$,<I>i</I>-2,<I>i</I>+2,\f$ \ldots\f$,<I>d</I>}, <I>f</I>(\f$ \beta_j\f$(<I>e</I>))=\f$ \beta_j^{-1}\f$(<I>f</I>(<I>e</I>)).

\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink,
  `d1`\f$ \in \f$ `darts()`, and `d2`\f$ \in \f$ `darts()`.
*/
template <unsigned int i> bool is_sewable(Dart_const_descriptor d1, Dart_const_descriptor d2) const;

/*!
  <I>i</I>-sew darts `d1` and `d2`, by keeping the generic map valid.

Links by \f$ \beta_i \f$ two by two all the darts of the orbit <I>D1</I>=\f$ \langle{}\f$\f$ \beta_1\f$,\f$ \ldots\f$,\f$ \beta_{i-2}\f$,\f$ \beta_{i+2}\f$,\f$ \ldots\f$,\f$ \beta_d\f$\f$ \rangle{}\f$(`d1`) and <I>D2</I>=\f$ \langle{}\f$\f$ \beta_0\f$,\f$ \beta_2\f$,\f$ \ldots\f$,\f$ \beta_{i-2}\f$,\f$ \beta_{i+2}\f$,\f$ \ldots\f$,\f$ \beta_d\f$\f$ \rangle{}\f$(`d2`) such that <I>d2</I>=<I>f</I>(<I>d1</I>), where <I>f</I> is the bijection between <I>D1</I> and <I>D2</I> satisfying: <I>f</I>(<I>d1</I>)=<I>d2</I>, and for all <I>e</I> \f$ \in \f$ <I>D1</I>, for all <I>j</I> \f$ \in \f$ {1,\f$ \ldots\f$,<I>i</I>-2,<I>i</I>+2,\f$ \ldots\f$,<I>d</I>}, <I>f</I>(\f$ \beta_j\f$(<I>e</I>))=\f$ \beta_j^{-1}\f$(<I>f</I>(<I>e</I>)).

If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`, when necessary, non void attributes are updated to ensure the validity of the generic map: for each <I>j</I>-cells <I>c1</I> and <I>c2</I> which are merged into one <I>j</I>-cell during the sew, the two associated attributes <I>attr1</I> and <I>attr2</I> are considered. If one attribute is `null_descriptor` and the other not, the non `null_descriptor` attribute is associated to all the darts of the resulting cell. When the two attributes are non `null_descriptor`, functor \link CellAttribute::On_merge `Attribute_type<i>::type::On_merge`\endlink is called on the two attributes <I>attr1</I> and <I>attr2</I>. If set, the dynamic onmerge function of <i>i</i>-attributes is also called on <I>attr1</I> and <I>attr2</I>. Then, the attribute <I>attr1</I> is associated to all darts of the resulting <I>j</I>-cell. Finally, attribute <I>attr2</I> is removed from the generic map.
\pre \link CombinatorialMap::is_sewable `is_sewable<i>(d1,d2)`\endlink.

\cgalAdvancedBegin
If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are not updated; thus the generic map can be no more valid after this operation.
\cgalAdvancedEnd

*/
template <unsigned int i> void sew(Dart_descriptor d1,Dart_descriptor d2);

/*!
  <I>i</I>-unsew darts `d` and `opposite<i>(d)`, by keeping the generic map valid.

Unlinks by \f$ \beta_i\f$ all the darts in the orbit \f$ \langle{}\f$\f$ \beta_1\f$,\f$ \ldots\f$,\f$ \beta_{i-2}\f$,\f$ \beta_{i+2}\f$,\f$ \ldots\f$,\f$ \beta_d\f$\f$ \rangle{}\f$(`d`).

If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`, when necessary, non void attributes are updated to ensure the validity of the generic map: for each <I>j</I>-cell <I>c</I> split in two <I>j</I>-cells <I>c1</I> and <I>c2</I> by the operation, if <I>c</I> is associated to a <I>j</I>-attribute <I>attr1</I>, then this attribute is duplicated into <I>attr2</I>, and all the darts belonging to <I>c2</I> are associated with this new attribute. Finally, the functor \link CellAttribute::On_split `Attribute_type<i>::type::On_split`\endlink is called on the two attributes <I>attr1</I> and <I>attr2</I>. If set, the dynamic onsplit function of <i>i</i>-attributes is also called on <I>attr1</I> and <I>attr2</I>.
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink, `d`\f$ \in \f$ `darts()` and `d` is not <I>i</I>-free.

\cgalAdvancedBegin
If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are not updated thus the generic map can be no more valid after this operation.
\cgalAdvancedEnd
*/
template <unsigned int i> void unsew(Dart_descriptor d);

/*!
  Reverse the orientation (swap \f$ \beta_0\f$ and \f$ \beta_1\f$ links) of the entire map.
*/
void reverse_orientation();

/*!
    Reverse the orientation (swap \f$ \beta_0\f$ and \f$ \beta_1\f$ links) of the connected component containing the given dart.
*/
void reverse_orientation_connected_component(Dart_descriptor d);

/// @}

}; /* end CombinatorialMap */
