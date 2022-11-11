
namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2IO

`Arr_extended_dcel_text_formatter` defines the format of an arrangement in an input or output stream
(typically a file stream), thus enabling reading and writing an `Arrangement`
instance using a simple text format. The `Arrangement` class should be
instantiated with a \dcel class which in turn instantiates the
`Arr_extended_dcel` template with the `VertexData`, `HalfedgeData` and
`FaceData` types.
The formatter supports reading and writing the data objects attached to the
arrangement vertices, halfedges and faces.

The `Arr_extended_dcel_text_formatter` class assumes that the nested `Point_2` and the `Curve_2` types
defined by the `Arrangement` template-parameter, as well as the `VertexData`,
`HalfedgeData` and `FaceData` types, can all be written to an input stream using
the `<<` operator and read from an input stream using the `>>` operator.

\cgalModels `ArrangementInputFormatter`
\cgalModels `ArrangementOutputFormatter`

\sa `PkgArrangementOnSurface2Read`
\sa `PkgArrangementOnSurface2Write`
\sa `Arr_extended_dcel<Traits,VData,HData,FData,V,H,F>`

*/
template< typename Arrangement >
class Arr_extended_dcel_text_formatter {
public:

}; /* end Arr_extended_dcel_text_formatter */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2IO

`Arr_face_extended_text_formatter` defines the format of an arrangement in an input or output stream
(typically a file stream), thus enabling reading and writing an `Arrangement`
instance using a simple text format. The `Arrangement` class should be
instantiated with a \dcel class which in turn instantiates the
`Arr_face_extended_dcel` template with a `FaceData` type.
The formatter supports reading and writing the data objects attached to the
arrangement faces as well.

The `Arr_face_extended_text_formatter` class assumes that the nested `Point_2` and the `Curve_2` types
defined by the `Arrangement` template-parameter and that the `FaceData` type
can all be written to an input stream using the `<<` operator and read from an input stream using the `>>` operator.

\cgalModels `ArrangementInputFormatter`
\cgalModels `ArrangementOutputFormatter`

\sa `PkgArrangementOnSurface2Read`
\sa `PkgArrangementOnSurface2Write`
\sa `Arr_face_extended_dcel<Traits,FData,V,H,F>`

*/
template< typename Arrangement >
class Arr_face_extended_text_formatter {
public:

}; /* end Arr_face_extended_text_formatter */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2IO

`Arr_text_formatter` defines the format of an arrangement in an input or output stream
(typically a file stream), thus enabling reading and writing an `Arrangement`
instance using a simple text format. The arrangement is assumed to store no auxiliary
data with its \dcel records (and if there are such records they will not be written
or read by the formatter).

The `Arr_text_formatter` class assumes that the nested `Point_2` and the `Curve_2` types
defined by the `Arrangement` template-parameter can both be written to an input
stream using the `<<` operator and read from an input stream using the `>>`
operator.

\cgalModels `ArrangementInputFormatter`
\cgalModels `ArrangementOutputFormatter`

\sa `PkgArrangementOnSurface2Read`
\sa `PkgArrangementOnSurface2Write`

*/
template< typename Arrangement >
class Arr_text_formatter {
public:

}; /* end Arr_text_formatter */
} /* end namespace CGAL */
