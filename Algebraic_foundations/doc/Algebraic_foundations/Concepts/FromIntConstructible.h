
/*!
\ingroup PkgAlgebraicFoundationsConcepts
\cgalConcept

A model of the concept `FromIntConstructible` is required 
to be constructible from int. 

\cgalHasModel int 
\cgalHasModel long 
\cgalHasModel double 

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

