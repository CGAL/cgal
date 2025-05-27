/*!
\ingroup PkgHDVFConcepts
\cgalConcept

The concept `Ring` describes the requirements for the ring of coefficients used to compute homology in the `CGAL::HomologicalDiscreteVectorField` concept. Besides ring operators, it also specifies the functions needed to test invertibility in the ring.
 
\cgalHasModelsBegin
\cgalHasModelsBare{`CGAL::Z`}
\cgalHasModelsBare{`CGAL::Zp`}
\cgalHasModelsEnd

*/

class Ring
{
public:
/// \name Types
/// @{

/*!
 Type of elements of the ring.
 */
typedef unspecified_type Value ;
/// @}
    
    
/// \name Operators
/// @{

    
/// @}
    

/// \name Invertibility
/// @{

    /*!
     \brief Test if `v` is invertible in the ring.
     */
    bool is_invertible(Value v);
    
/// @}
};
