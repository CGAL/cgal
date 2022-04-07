
namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2IO

`Arr_with_history_text_formatter` defines the format of an arrangement in an input or output stream
(typically a file stream), thus enabling reading and writing an
arrangement-with-history instance using a simple text format.

The `ArrFormatter` parameter servers as a base class for
`Arr_with_history_text_formatter` and must be a model of the `ArrangementInputFormatter`
and the `ArrangementOutputFormatter` concepts. It is used to read or write
the base arrangement, while the derived class is responsible for reading and
writing the set of curves inducing the arrangement and maintaining the
relations between these curves and the edges they induce.

\cgalModels `ArrangementWithHistoryInputFormatter`
\cgalModels `ArrangementWithHistoryOutputFormatter`

\sa `PkgArrangementOnSurface2Read`
\sa `PkgArrangementOnSurface2Write`

*/
template< typename ArrFormatter >
class Arr_with_history_text_formatter {
public:

}; /* end Arr_with_history_text_formatter */
} /* end namespace CGAL */
