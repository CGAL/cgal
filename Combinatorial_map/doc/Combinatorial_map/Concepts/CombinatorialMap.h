/*!
\ingroup PkgCombinatorialMapsConcepts
\cgalConcept

The concept `CombinatorialMap` defines a <I>d</I>-dimensional combinatorial map.

\cgalRefines `BasicMap`

\cgalHasModel \link CGAL::Combinatorial_map `CGAL::Combinatorial_map<d,Items,Alloc>`\endlink

\cgalHeading{Basic functions}

For a combinatorial map, the function \link BasicMap::next `next`\endlink is equal to \f$ \beta_1\f$, \link BasicMap::previous `previous`\endlink is equal to \f$ \beta_0\f$  and the function \link BasicMap::opposite `opposite<i>`\endlink is equal to \f$ \beta_i\f$.

\cgalHeading{Validity}

The function \link BasicMap::is_valid `is_valid`\endlink returns true iff the combinatorial map is valid.
A combinatorial map is valid (see Sections \ref sseccombimapanddarts and \ref sseccombimapvalidity) if for all its darts `d` \f$\in\f$`darts()`:

- `d` is 0-free, or \f$ \beta_1(\beta_0(d))=d\f$;
- `d` is 1-free, or \f$ \beta_0(\beta_1(d))=d\f$;
- \f$ \forall\f$<I>i</I>, 2\f$ \leq\f$<I>i</I>\f$ \leq\f$\link BasicMap::dimension `dimension`\endlink: `d` is <i>i</i>-free, or \f$ \beta_i(\beta_i(d))=d\f$;
- \f$ \forall\f$<I>i</I>, <I>j</I>, 0\f$ \leq\f$<I>i</I>\f$ <\f$<I>i</I>+2\f$ \leq\f$<I>j</I>\f$ \leq\f$\link BasicMap::dimension `dimension`\endlink such that <I>j</I>\f$ \geq\f$ 3: \f$ \beta_j(\beta_i(d))=\varnothing\f$ or ; \f$ \beta_j(\beta_i(\beta_j(\beta_i(d))))=d\f$;
- \f$ \forall\f$<I>i</I>, 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\link BasicMap::dimension `dimension`\endlink
  such that <I>i</I>-attributes are non void:
  + \f$ \forall\f$<I>d2</I> in the same <I>i</I>-cell than <I>d</I>: <I>d</I> and <I>d2</I> have the same <I>i</I>-attribute;
  + \f$ \forall\f$<I>d2</I>  in a different <I>i</I>-cell than <I>d</I>: <I>d</I> and <I>d2</I> have different <I>i</I>-attributes.


\cgalHeading{Sew and unsew}

The function \link BasicMap::is_sewable `is_sewable`\endlink returns true iff `*dh1` can be <I>i</I>-sewn with `*dh2` by keeping the combinatorial map valid.
This is true if there is a bijection <I>f</I> between all the darts of the orbit <I>D1</I>=\f$ \langle{}\f$\f$ \beta_1\f$,\f$ \ldots\f$,\f$ \beta_{i-2}\f$,\f$ \beta_{i+2}\f$,\f$ \ldots\f$,\f$ \beta_d\f$\f$ \rangle{}\f$(<I>*dh1</I>) and <I>D2</I>=\f$ \langle{}\f$\f$ \beta_1\f$,\f$ \ldots\f$,\f$ \beta_{i-2}\f$,\f$ \beta_{i+2}\f$,\f$ \ldots\f$,\f$ \beta_d\f$\f$ \rangle{}\f$(<I>*dh2</I>) satisfying: <I>f</I>(<I>*dh1</I>)=<I>*dh2</I>, and for all <I>e</I>\f$ \in\f$<I>D1</I>, for all <I>j</I>\f$ \in\f${1,\f$ \ldots\f$,<I>i</I>-2,<I>i</I>+2,\f$ \ldots\f$,<I>d</I>}, <I>f</I>(\f$ \beta_j\f$(<I>e</I>))=\f$ \beta_j^{-1}\f$(<I>f</I>(<I>e</I>)).

The function \link BasicMap::sew `sew`\endlink <I>i</I>-sew darts `*dh1` and `*dh2`, by keeping the combinatorial map valid.
Links by \f$ \beta_i\f$ two by two all the darts of the orbit <I>D1</I>=\f$ \langle{}\f$\f$ \beta_1\f$,\f$ \ldots\f$,\f$ \beta_{i-2}\f$,\f$ \beta_{i+2}\f$,\f$ \ldots\f$,\f$ \beta_d\f$\f$ \rangle{}\f$(`*dh1`) and <I>D2</I>=\f$ \langle{}\f$\f$ \beta_0\f$,\f$ \beta_2\f$,\f$ \ldots\f$,\f$ \beta_{i-2}\f$,\f$ \beta_{i+2}\f$,\f$ \ldots\f$,\f$ \beta_d\f$\f$ \rangle{}\f$(`*dh2`) such that <I>d2</I>=<I>f</I>(<I>d1</I>).
<I>f</I> is the bijection between <I>D1</I> and <I>D2</I> satisfying: <I>f</I>(<I>*dh1</I>)=<I>*dh2</I>, and for all <I>e</I>\f$ \in\f$<I>D1</I>, for all <I>j</I>\f$ \in\f${1,\f$ \ldots\f$,<I>i</I>-2,<I>i</I>+2,\f$ \ldots\f$,<I>d</I>}, <I>f</I>(\f$ \beta_j\f$(<I>e</I>))=\f$ \beta_j^{-1}\f$(<I>f</I>(<I>e</I>)).

The function \link BasicMap::unsew `unsew`\endlink <I>i</I>-unsew darts `*dh` and \f$ \beta_i\f$`(*dh)`, by keeping the combinatorial map valid.
Unlinks by \f$ \beta_i\f$ all the darts in the orbit \f$ \langle{}\f$\f$ \beta_1\f$,\f$ \ldots\f$,\f$ \beta_{i-2}\f$,\f$ \beta_{i+2}\f$,\f$ \ldots\f$,\f$ \beta_d\f$\f$ \rangle{}\f$(`*dh`).

*/

class CombinatorialMap {
public:

/// \name Access Member Functions
/// @{

/*!
Returns \f$ \beta_j\f$(\f$ \beta_i\f$(`*dh`)).
Overloads of this member function are defined that take from one to nine integer as arguments.
For each function, betas are applied in the same order as their indices are given as parameters.

For example `beta(dh,1)`=\f$ \beta_1\f$(`*dh`),
and `beta(dh,1,2,3,0)`=\f$ \beta_0\f$(\f$ \beta_3\f$(\f$ \beta_2\f$(\f$ \beta_1\f$(`*dh`)))).
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\link CombinatorialMap::dimension `dimension`\endlink,
  0\f$ \leq\f$<I>j</I>\f$ \leq\f$\link CombinatorialMap::dimension `dimension`\endlink
  and `*dh`\f$ \in\f$`darts()`.
*/
Dart_handle beta(Dart_handle dh, int i, int j);

/*!
Returns \f$ \beta_j\f$(\f$ \beta_i\f$(`*dh`)).
Overloads of this member function are defined that take from one to nine integer as arguments.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\link CombinatorialMap::dimension `dimension`\endlink,
     0\f$ \leq\f$<I>j</I>\f$ \leq\f$\link CombinatorialMap::dimension `dimension`\endlink
     and `*dh`\f$ \in\f$`darts()`.

*/
Dart_const_handle beta(Dart_const_handle dh, int i, int j) const;

/*!
Returns \f$ \beta_j\f$(\f$ \beta_i\f$(`*dh`)).
Overloads of this member function are defined that take from one to nine integer as template arguments.
For each function, betas are applied in the same order as their indices are given as template arguments.

For example `beta<1>(dh)`=\f$ \beta_1\f$(`*dh`),
and `beta<1,2,3,0>(dh)`=\f$ \beta_0\f$(\f$ \beta_3\f$(\f$ \beta_2\f$(\f$ \beta_1\f$(`*dh`)))).
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\link CombinatorialMap::dimension `dimension`\endlink,
  0\f$ \leq\f$<I>j</I>\f$ \leq\f$\link CombinatorialMap::dimension `dimension`\endlink
  and `*dh`\f$ \in\f$`darts()`.
*/
template<int i, int j>
Dart_handle beta(Dart_handle dh);

/*!
Returns \f$ \beta_j\f$(\f$ \beta_i\f$(`*dh`)).
Overloads of this member function are defined that take from one to nine integer as template arguments.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\link CombinatorialMap::dimension `dimension`\endlink,
     0\f$ \leq\f$<I>j</I>\f$ \leq\f$\link CombinatorialMap::dimension `dimension`\endlink
     and `*dh`\f$ \in\f$`darts()`.

*/
template<int i, int j>
Dart_const_handle beta(Dart_const_handle dh) const;

/*!
Returns a handle to a dart belonging to the same edge than dart `*dh`, and not to the same vertex. `NULL` if such a dart does not exist.
*/
Dart_handle opposite(Dart_handle dh);

/*!
Returns a const handle to a dart belonging to the same edge than dart `*dh`, and not to the same vertex, when the dart is const. `NULL` if such a dart does not exist.
*/
Dart_const_handle opposite(Dart_const_handle dh) const;

/*!
Links `*dh1` and `*dh2` by \f$ \beta_i\f$.
The combinatorial map can be no more valid after this operation. If \link CombinatorialMap::are_attributes_automatically_managed `are_attributes_automatically_managed()`\endlink`==true`, non void attributes of `*dh1` and `*dh2` are updated: if one dart has an attribute and the second dart not, the non null attribute is associated to the dart having a null attribute.
If both darts have an attribute, the attribute of `*dh1` is associated to `*dh2`.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\link CombinatorialMap::dimension `dimension`\endlink,
    `*dh1`\f$ \in\f$`darts()`, `*dh2`\f$ \in\f$`darts()` and (<I>i</I>\f$ <\f$ 2 or `dh1`\f$ \neq\f$`dh2`).
*/
template <unsigned int i> void link_beta(Dart_handle dh1, Dart_handle dh2);

/*!
Unlinks `*dh` and \f$ \beta_i\f$(`*dh`) by \f$ \beta_i\f$.
The combinatorial map can be no more valid after this operation.
Attributes of `*dh` and \f$ \beta_i\f$(`*dh`) are not modified.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$\link CombinatorialMap::dimension `dimension`\endlink,
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
  
}; /* end CombinatorialMap */


