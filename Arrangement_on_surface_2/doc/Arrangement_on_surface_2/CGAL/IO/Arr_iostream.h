namespace CGAL {

namespace IO {

/*!
  \defgroup PkgArrangementOnSurface2Read CGAL::IO::read()
  \ingroup PkgArrangementOnSurface2IO

Reads a given arrangement from a given input stream
using a specific input format.

\cgalHeading{Requirements}

<UL>
<LI>The instantiated `Formatter` class must model the
  `ArrangementInputFormatter` concept.
<LI>The instantiated `WithHistoryFormatter` class must model the
  `ArrangementWithHistoryInputFormatter` concept.
</UL>


\sa `PkgArrangementOnSurface2Write`


\sa `PkgArrangementOnSurface2op_left_shift`
\sa `PkgArrangementOnSurface2op_right_shift`
*/
/// @{

/*!
Reads the arrangement object `arr` from the given input stream `is`
using a specific input format defined by \"formatter\".
*/
template<typename Traits, typename Dcel, typename Formatter>
std::istream& read (Arrangement_2<Traits,Dcel>& arr,
                    std::istream& is,
                    Formatter& formatter);


/// @}

/*!
  \defgroup PkgArrangementOnSurface2Write CGAL::IO::write()
  \ingroup PkgArrangementOnSurface2IO

Writes a given arrangement into a given output stream
using a specific output format.

\cgalHeading{Requirements}

<UL>
<LI>The instantiated `Formatter` class must model the
  `ArrangementOutputFormatter` concept.
<LI>The instantiated `WithHistoryFormatter` class must model the
  `ArrangementWithHistoryOutputFormatter` concept.
</UL>

\sa `PkgArrangementOnSurface2Read`
\sa `PkgArrangementOnSurface2op_left_shift`
\sa `PkgArrangementOnSurface2op_right_shift`
*/
/// @{

/*!
Writes the arrangement object `arr` into the given output stream
`os` using a specific output format defined by `formatter`.
*/
template<typename Traits, typename Dcel, typename Formatter>
std::ostream& write (const Arrangement_2<Traits,Dcel>& arr,
                     std::ostream& os,
                     Formatter& formatter);

/// @}

} // namespace IO

/*!
\ingroup PkgArrangementOnSurface2op_left_shift
Inserts the arrangement object `arr` into the output stream
`os` using the output format defined by the
`Arr_text_formatter` class. Only the basic geometric and
topological features of the arrangement are inserted. Auxiliary data
that may be attached to the \dcel features is ignored.
*/
template<typename Traits, typename Dcel>
std::ostream& operator<< (std::ostream& os,
                          const Arrangement_2<Traits,Dcel>& arr);


/*!
\ingroup PkgArrangementOnSurface2op_right_shift
Extracts an arrangement from a given input stream using the input
format defined by the `Arr_text_formatter` class - that is, only the
basic geometric and topological features of the arrangement are read
and no auxiliary data is attached to the Dcel features.
*/
template<class Traits, class Dcel>
std::istream& operator>>(std::istream& is, Arrangement_2<Traits,Dcel>& arr);


} /* end namespace CGAL::IO*/
