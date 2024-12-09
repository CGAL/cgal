
/*!
\ingroup PkgHandlesAndCirculatorsConcepts
\cgalConcept

A constant handle. Refer to the `Handle` concept for more details.

\cgalRefines{Descriptor}

\cgalHasModelsBegin
\cgalHasModels{const T* (const pointers)}
\cgalHasModelsEnd

\sa `Handle`

*/
class ConstHandle {
public:

/// \name Dereference
/// @{

/*!
returns the object pointed to.
*/
const value_type& operator*();

/*!
returns a pointer to the object pointed to.
*/
const value_type* operator->();

/// @}

}; /* end Handle */

