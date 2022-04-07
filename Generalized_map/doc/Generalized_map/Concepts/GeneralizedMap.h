/*!
\ingroup PkgGeneralizedMapsConcepts
\cgalConcept

The concept `GeneralizedMap` defines a <I>d</I>-dimensional generalized map.

\cgalRefines `GenericMap`

\cgalHasModel \link CGAL::Generalized_map `CGAL::Generalized_map<d,Items,Alloc>`\endlink

For a generalized map, the function \link GenericMap::next `next`\endlink is equal to \f$ \alpha_1\circ\alpha_0\f$, \link GenericMap::previous `previous`\endlink is equal to \f$ \alpha_0\circ\alpha_1\f$  and the function \link GenericMap::opposite `opposite<i>`\endlink is equal to \f$ \alpha_i\circ\alpha_0\f$.
*/

//@{
class GeneralizedMap {
public:

/// \name Access Member Functions
/// @{

/*!
Returns \f$ \alpha_j\f$(\f$ \alpha_i\f$(`*dh`)).
Overloads of this member function are defined that take from one to nine integer as arguments. For each function, alphas are applied in the same order as their indices are given as parameters.

For example `alpha(dh,1)`=\f$ \alpha_1\f$(`*dh`),
and `alpha(dh,1,2,3,0)`=\f$ \alpha_0\f$(\f$ \alpha_3\f$(\f$ \alpha_2\f$(\f$ \alpha_1\f$(`*dh`)))).
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink, 0 \f$ \leq \f$ <I>j</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink and `*dh`\f$ \in\f$ \link GenericMap::darts `darts()`\endlink.
*/
Dart_handle alpha(Dart_handle dh, int i, int j);

/*!
Returns \f$ \alpha_j\f$(\f$ \alpha_i\f$(`*dh`)). Overloads of this member function are defined that take from one to nine integer as arguments.
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink, 0 \f$ \leq \f$ <I>j</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink and `*dh`\f$ \in\f$ \link GenericMap::darts `darts()`\endlink.

*/
Dart_const_handle alpha(Dart_const_handle dh, int i, int j) const;

/*!
Returns \f$ \alpha_j\f$(\f$ \alpha_i\f$(`*dh`)). Overloads of this member function are defined that take from one to nine integer as template arguments. For each function, alphas are applied in the same order as their indices are given as template arguments.

For example `alpha<1>(dh)`=\f$ \alpha_1\f$(`*dh`), and `alpha<1,2,3,0>(dh)`=\f$ \alpha_0\f$(\f$ \alpha_3\f$(\f$ \alpha_2\f$(\f$ \alpha_1\f$(`*dh`)))). \pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink, 0 \f$ \leq \f$ <I>j</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink and `*dh`\f$ \in\f$ \link GenericMap::darts `darts()`\endlink.
*/
template<int i, int j>
Dart_handle alpha(Dart_handle dh);

/*!
Returns \f$ \alpha_j\f$(\f$ \alpha_i\f$(`*dh`)). Overloads of this member function are defined that take from one to nine integer as template arguments.
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink, 0 \f$ \leq \f$ <I>j</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink and `*dh`\f$ \in\f$ \link GenericMap::darts `darts()`\endlink.

*/
template<int i, int j>
Dart_const_handle alpha(Dart_const_handle dh) const;

/*!
Returns true iff the current generalized map is orientable.
*/
bool is_orientable() const;

/*!
Returns true iff the generalized map is valid.

A generalized map is valid (see Sections \ref sec_definition_gmap and \ref ssecgenmapvalidity) if for all its darts `d` \f$\in\f$ \link GenericMap::darts `darts()`\endlink:

- \f$ \forall \f$ <I>i</I>, 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink: \f$ \alpha_i(\alpha_i(d))=d\f$;
- \f$ \forall \f$ <I>i</I>, \f$ \forall \f$ <I>j</I> such that
0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink-2
and <I>i</I>+2 \f$ \leq \f$ <I>j</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink:
\f$ \alpha_j(\alpha_i(\alpha_j(\alpha_i(d))))=d\f$;
- \f$ \forall \f$ <I>i</I>, 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink such that <I>i</I>-attributes are non void:
  + \f$ \forall \f$ <I>d2</I> in the same <I>i</I>-cell than <I>d</I>: <I>d</I> and <I>d2</I> have the same <I>i</I>-attribute;
  + \f$ \forall \f$ <I>d2</I>  in a different <I>i</I>-cell than <I>d</I>: <I>d</I> and <I>d2</I> have different <I>i</I>-attributes.
*/
bool is_valid() const;

/// @}

/// \name Modifiers
/// @{

/*!
Links `*dh1` and `*dh2` by \f$ \alpha_i\f$. The generalized map can be no longer valid after this operation. If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`, non void attributes of `*dh1` and `*dh2` are updated: if one dart has an attribute and the second dart not, the non null attribute is associated to the dart having a null attribute. If both darts have an attribute, the attribute of `*dh1` is associated to `*dh2`.
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink, `*dh1`\f$ \in\f$ \link GenericMap::darts `darts()`\endlink, `*dh2`\f$ \in\f$ \link GenericMap::darts `darts()`\endlink and (<I>i</I>\f$ <\f$ 2 or `dh1`\f$ \neq\f$`dh2`).
*/
template <unsigned int i> void link_alpha(Dart_handle dh1, Dart_handle dh2);

/*!
Unlinks `*dh` and \f$ \alpha_i\f$(`*dh`) by \f$ \alpha_i\f$. The generalized map can be no longer valid after this operation. Attributes of `*dh` and \f$ \alpha_i\f$(`*dh`) are not modified.
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink, `*dh`\f$ \in\f$ \link GenericMap::darts `darts()`\endlink, and `*dh` is not <I>i</I>-free.
*/
template <unsigned int i> void unlink_alpha(Dart_handle dh);

/// @}

/// \name Operations
/// @{
/*!
Returns true iff `*dh1` can be <I>i</I>-sewn with `*dh2` by keeping the generic map valid.

This is true if there is a bijection <I>f</I> between all the darts of the orbit <I>D1</I>=\f$ \langle{}\f$\f$ \alpha_0\f$,\f$ \ldots\f$,\f$ \alpha_{i-2}\f$,\f$ \alpha_{i+2}\f$,\f$ \ldots\f$,\f$ \alpha_d\f$\f$ \rangle{}\f$(<I>*dh1</I>) and <I>D2</I>=\f$ \langle{}\f$\f$ \alpha_0\f$,\f$ \ldots\f$,\f$ \alpha_{i-2}\f$,\f$ \alpha_{i+2}\f$,\f$ \ldots\f$,\f$ \alpha_d\f$\f$ \rangle{}\f$(<I>*dh2</I>) satisfying:
- <I>f</I>(<I>*dh1</I>)=<I>*dh2</I>,
- \f$ \forall \f$ <I>e</I>\f$ \in\f$ <I>D1</I>, \f$ \forall \f$ <I>j</I>\f$ \in\f$ {1,\f$ \ldots\f$,<I>i</I>-2,<I>i</I>+2,\f$ \ldots\f$,<I>d</I>}, <I>f</I>(\f$ \alpha_j\f$(<I>e</I>))=\f$ \alpha_j\f$(<I>f</I>(<I>e</I>)).

\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink,
  `*dh1`\f$ \in \f$ `darts()`, and `*dh2`\f$ \in \f$ `darts()`.
*/
template <unsigned int i> bool is_sewable(Dart_const_handle dh1, Dart_const_handle dh2) const;

/*!
  <I>i</I>-sew darts `*dh1` and `*dh2`, by keeping the generic map valid.

Links by \f$ \alpha_i\f$ two by two all the darts of the orbit <I>D1</I>=\f$ \langle{}\f$\f$ \alpha_0\f$,\f$ \ldots\f$,\f$ \alpha_{i-2}\f$,\f$ \alpha_{i+2}\f$,\f$ \ldots\f$,\f$ \alpha_d\f$\f$ \rangle{}\f$(`*dh1`) and <I>D2</I>=\f$ \langle{}\f$\f$ \alpha_0\f$,\f$ \alpha_2\f$,\f$ \ldots\f$,\f$ \alpha_{i-2}\f$,\f$ \alpha_{i+2}\f$,\f$ \ldots\f$,\f$ \alpha_d\f$\f$ \rangle{}\f$(`*dh2`) such that <I>d2</I>=<I>f</I>(<I>d1</I>), where <I>f</I> is the bijection between <I>D1</I> and <I>D2</I> satisfying: <I>f</I>(<I>*dh1</I>)=<I>*dh2</I>, and for all <I>e</I>\f$ \in\f$ <I>D1</I>, for all <I>j</I>\f$ \in\f$ {1,\f$ \ldots\f$,<I>i</I>-2,<I>i</I>+2,\f$ \ldots\f$,<I>d</I>}, <I>f</I>(\f$ \alpha_j\f$(<I>e</I>))=\f$ \alpha_j\f$(<I>f</I>(<I>e</I>)).

If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`, when necessary, non void attributes are updated to ensure the validity of the generic map: for each <I>j</I>-cells <I>c1</I> and <I>c2</I> which are merged into one <I>j</I>-cell during the sew, the two associated attributes <I>attr1</I> and <I>attr2</I> are considered. If one attribute is `nullptr` and the other not, the non `nullptr` attribute is associated to all the darts of the resulting cell. When the two attributes are non `nullptr`, functor \link CellAttribute::On_merge `Attribute_type<i>::type::On_merge`\endlink is called on the two attributes <I>attr1</I> and <I>attr2</I>. If set, the dynamic onmerge function of <i>i</i>-attributes is also called on <I>attr1</I> and <I>attr2</I>. Then, the attribute <I>attr1</I> is associated to all darts of the resulting <I>j</I>-cell. Finally, attribute <I>attr2</I> is removed from the generic map.
\pre \link GeneralizedMap::is_sewable `is_sewable<i>(dh1,dh2)`\endlink.

\cgalAdvancedBegin
If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are not updated; thus the generic map can be no more valid after this operation.
\cgalAdvancedEnd

*/
template <unsigned int i> void sew(Dart_handle dh1,Dart_handle dh2);

/*!
  <I>i</I>-unsew darts `*dh` and `*opposite<i>(*dh)`, by keeping the generic map valid.

Unlinks by \f$ \alpha_i\f$ all the darts in the orbit \f$ \langle{}\f$\f$ \alpha_0\f$,\f$ \ldots\f$,\f$ \alpha_{i-2}\f$,\f$ \alpha_{i+2}\f$,\f$ \ldots\f$,\f$ \alpha_d\f$\f$ \rangle{}\f$(`*dh`).

If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`, when necessary, non void attributes are updated to ensure the validity of the generic map: for each <I>j</I>-cell <I>c</I> split in two <I>j</I>-cells <I>c1</I> and <I>c2</I> by the operation, if <I>c</I> is associated to a <I>j</I>-attribute <I>attr1</I>, then this attribute is duplicated into <I>attr2</I>, and all the darts belonging to <I>c2</I> are associated with this new attribute. Finally, the functor \link CellAttribute::On_split `Attribute_type<i>::type::On_split`\endlink is called on the two attributes <I>attr1</I> and <I>attr2</I>. If set, the dynamic onsplit function of <i>i</i>-attributes is also called on <I>attr1</I> and <I>attr2</I>.
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ \link GenericMap::dimension `dimension`\endlink, `*dh`\f$ \in \f$ `darts()` and `*dh` is not <I>i</I>-free.

\cgalAdvancedBegin
If \link GenericMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==false`, non void attributes are not updated thus the generic map can be no more valid after this operation.
\cgalAdvancedEnd
*/
template <unsigned int i> void unsew(Dart_handle dh);

/// @}

}; /* end GeneralizedMap */
//@}

