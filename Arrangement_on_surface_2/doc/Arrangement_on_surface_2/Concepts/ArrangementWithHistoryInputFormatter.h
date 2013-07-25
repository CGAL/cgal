
/*!
\ingroup PkgArrangement2Concepts
\cgalConcept

A model for the `ArrangementWithHistoryInputFormatter` concept supports a set of functions that enable 
reading an arrangement-with-history instance from an input stream using a 
specific format. 

\cgalRefines `ArrangementInputFormatter` 

\cgalHasModel `CGAL::Arr_with_history_text_formatter<ArrFormatter>`

*/

class ArrangementWithHistoryInputFormatter {
public:

/// \name Types 
/// @{

/*!
the type of arrangement to input. 
*/ 
typedef unspecified_type Arr_with_history_2; 

/*!
the inducing curve type. 
*/ 
typedef typename Arrangement_2::Curve_2 Curve_2; 

/// @} 

/// \name Formatted Input Functions 
/// @{

/*!
reads a message indicating the beginning of the inducing curves. 
*/ 
void read_curves_begin(); 

/*!
reads a message indicating the end of the inducing curves. 
*/ 
void read_curves_end(); 

/*!
reads a message indicating the beginning of a single curve record. 
*/ 
void read_curve_begin(); 

/*!
reads a message indicating the end of a single curve record. 
*/ 
void read_curve_end(); 

/*!
reads a curve. 
*/ 
void read_curve (Curve_2& c); 

/*!
reads a message indicating the beginning of the set of edges 
induced by the current curve. 
*/ 
void read_induced_edges_begin(); 

/*!
reads a message indicating the end of the induced edges set. 
*/ 
void read_induced_edges_end(); 

/// @}

}; /* end ArrangementWithHistoryInputFormatter */

