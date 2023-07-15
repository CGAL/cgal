
/*!
\ingroup PkgAlgebraicFoundationsConcepts
\cgalConcept

A model of the concept `FromIntConstructible` is required
to be constructible from int.

\cgalHasModelsBegin
\cgalModels{int}
\cgalModels{long}
\cgalModels{double}
\cgalHasModelsEnd

*/

class FromIntConstructible {
public:

/// \name Creation
/// @{

/*!

*/
FromIntConstructible(int& i);

/// @}

}; /* end FromIntConstructible */

