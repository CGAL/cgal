
/*!
\ingroup PkgArrangementOnSurface2Concepts
\cgalConcept

A model for the `ArrangementWithHistoryOutputFormatter` concept supports a set of functions that enable
writing an arrangement-with-history instance to an output stream using a
specific format.

\cgalRefines{ArrangementOutputFormatter}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Arr_with_history_text_formatter<ArrFormatter>}
\cgalHasModelsEnd

*/

class ArrangementWithHistoryOutputFormatter {
public:

/// \name Types
/// @{

/*!
the type of arrangement to output.
*/
typedef unspecified_type Arr_with_history_2;

/*!
the inducing curve type.
*/
typedef typename Arrangement_2::Curve_2 Curve_2;

/// @}

/// \name Formatted Output Functions
/// @{

/*!
writes a message indicating the beginning of the inducing curves.
*/
void write_curves_begin();

/*!
writes a message indicating the end of the inducing curves.
*/
void write_curves_end();

/*!
writes a message indicating the beginning of a single curve record.
*/
void write_curve_begin();

/*!
writes a message indicating the end of a single curve record.
*/
void write_curve_end();

/*!
writes a curve.
*/
void write_curve (const Curve_2& c);

/*!
writes a message indicating the beginning of the set of edges
induced by the current curve.
*/
void write_induced_edges_begin();

/*!
writes a message indicating the end of the induced edges set.
*/
void write_induced_edges_end();

/// @}

}; /* end ArrangementWithHistoryOutputFormatter */

