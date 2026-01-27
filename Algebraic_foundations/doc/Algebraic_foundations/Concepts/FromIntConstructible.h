
/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

A model of the concept `FromIntConstructible` is required
to be constructible from int.

\cgalHasModelsBegin
\cgalHasModels{int}
\cgalHasModels{long}
\cgalHasModels{double}
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
