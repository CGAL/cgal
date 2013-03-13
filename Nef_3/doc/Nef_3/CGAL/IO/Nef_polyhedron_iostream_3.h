namespace CGAL {

/*!
\ingroup PkgNef3IOFunctions

This operator reads a Nef polyhedron, which is given in the proprietary file 
format written by the input operator <I>in</I> and assigns it to <I>N</I>. It includes the 
complete incidence structure, the geometric data, and the marks of each item. 

It is recommended to use the \cgal kernels `Homogeneous`, 
`Simple_homogeneous`, 
or `Extended_homogeneous` parametrized with any exact number type that models 
\f$ \mathbb{Z}\f$ (e.g.`Gmpz` or `leda_integer`). The input and 
output iterators of Nef polyhedra parametrized with 
either of these kernels are compatible as long as the Nef polyhedron is bounded. 
An unbounded Nef polyhedron can only be read by a Nef polyhedron parametrized with 
an extended kernel. It is also recommended to use the \cgal stream modifier 
`set_ascii_mode()`. 

\sa `CGAL::Nef_polyhedron_3<Traits>` 

*/
template <class Traits>
istream& operator>>( std::istream& in, CGAL::Nef_polyhedron_3<Traits>& N);

/*!

\ingroup PkgNef3IOFunctions

This operator writes the Nef polyhedron `N` to the output stream `out`
using a proprietary file format. It includes the complete incidence
structure, the geometric data, and the marks of each item.

Using \cgal stream modifiers the following output formats can be
chosen: ASCII (`set_ascii_mode()`), binary (`set_binary_mode()`) or
pretty (`set_pretty_mode()`). The mandatory format is the ASCII
format. It is recommended to use this format for file input and
output.

As the output depends on the output operators of the geometric
primitives provided by the traits class, it might not be possible that
the input operator and output operators of different traits classes
are not compatible. We recommend to use the \cgal kernels
`Homogeneous`, `Simple_homogeneous`, or `Extended_homogeneous`
parametrized with any exact number type that models \f$\mathbb{Z}\f$
(e.g. `Gmpz` or `leda_integer`).

A bounded `Nef_polyhedron_3<Extended_homogeneous>` is automatically
written as though `Nef_polyhedron_3<Homogeneous>` or
`Nef_polyhedron_3<Simple_homogeneous>` is used. As a result, the
input operator of each of these types can read the output.

\sa `CGAL::Nef_polyhedron_3<Traits>`

*/
template <class Traits>
ostream& operator<<( std::ostream& out, CGAL::Nef_polyhedron_3<Traits>& N);

} /* namespace CGAL */

