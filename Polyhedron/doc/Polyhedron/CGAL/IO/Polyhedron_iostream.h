namespace CGAL {

/*!
\relates Polyhedron_3
\ingroup PkgPolyhedronIOFunc

reads a polyhedral surface in %Object File Format, OFF,
with file extension <TT>.off</TT>, from the input stream `in` and
appends it to the polyhedral surface \f$ P\f$. Only the point coordinates
and facets from the input stream are used to build the polyhedral
surface. Neither normal vectors nor color attributes are evaluated.

This function overloads the generic function \link PkgBGLIoFuncsOFF `read_OFF(std::istream&, FaceGraph)`\endlink
(if no named parameters are used).

\note This function requires `Traits` to be a model of `Kernel` (a stronger requirement
than `PolyhedronTraits_3`).

\cgalHeading{Implementation}

This operator is implemented using the modifier mechanism for
polyhedral surfaces and the `Polyhedron_incremental_builder_3`
class, which allows the construction in a single, efficient scan pass
of the input and handles also all the possible flexibility of the
polyhedral surface.

\note Starting with \cgal 5.2, this function will no longer set the `ios::badbit`
      of the input stream `in` if the file contains 0 vertices.

\sa `CGAL::Polyhedron_3<Traits>`
\sa `CGAL::Polyhedron_incremental_builder_3<HDS>`
\sa \link PkgPolyhedronIOFunc `operator<<(std::ostream&, Polyhedron_3<Traits>&)`\endlink
*/
template <class Traits>
bool read_OFF( std::istream& in, Polyhedron_3<Traits>& P);

/*!
\relates Polyhedron_3
\deprecated This function is deprecated since \cgal 5.3,
            \link PkgPolyhedronIOFunc `CGAL::IO::read_OFF(std::ostream&, Polyhedron_3<Traits>&)` \endlink should be used instead.
*/
template <class Traits>
bool read_off( std::ostream& out, Polyhedron_3<Traits>& P);

/*!
\relates Polyhedron_3
\ingroup PkgPolyhedronIOFunc
calls \link read_OFF() `read_OFF(in, P)` \endlink.
*/
template <class Traits>
std::istream& operator>>( std::istream& in, Polyhedron_3<Traits>& P);

/*!
\relates Polyhedron_3
\ingroup PkgPolyhedronIOFunc

writes the polyhedral surface \f$P\f$ to the output
stream `out` using the %Object File Format, OFF, with file extension
<TT>.off</TT>. The output is in ASCII format. From the polyhedral surface, only the point
coordinates and facets are written. Neither normal vectors nor color
attributes are used.

For OFF an ASCII and a binary format exist. The format can be selected
with the \cgal modifiers for streams, `set_ascii_mode()` and
`set_binary_mode()` respectively. The modifier `set_pretty_mode()` can be used
to allow for (a few) structuring comments in the output. Otherwise,
the output would be free of comments. The default for writing is ASCII
without comments.

This function overloads the generic function \link PkgBGLIoFuncsOFF `write_OFF(std::istream&,FaceGraph)` \endlink
(if no named parameters are used).

\sa `CGAL::Polyhedron_3<Traits>`
\sa `CGAL::Polyhedron_incremental_builder_3<HDS>`
\sa \link PkgPolyhedronIOFunc `operator>>(std::istream& in, Polyhedron_3<Traits>& P)` \endlink
*/
template <class Traits>
bool write_OFF( std::ostream& out, Polyhedron_3<Traits>& P);

/*!
\relates Polyhedron_3
\deprecated This function is deprecated since \cgal 5.3,
            \link PkgPolyhedronIOFunc `CGAL::IO::write_OFF(std::ostream&, Polyhedron_3<Traits>&)` \endlink should be used instead.
*/
template <class Traits>
bool write_off( std::ostream& out, Polyhedron_3<Traits>& P);

/*!
\relates Polyhedron_3
\ingroup PkgPolyhedronIOFunc
calls \link write_OFF() `write_OFF(out, P)` \endlink.
*/
template <class Traits>
std::ostream& operator<<( std::ostream& out, Polyhedron_3<Traits>& P);

} /* namespace CGAL */

