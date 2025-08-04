/*!
\ingroup PkgHDVFConcepts
\cgalConcept

The concept `Ring` describes the requirements for the ring of coefficients used to compute homology in the `HomologicalDiscreteVectorField` concept. Besides ring operators, it also specifies the functions needed to test invertibility in the ring.

\cgalHasModelsBegin
\cgalHasModelsBare{`CGAL::HDVF::Z2`}
\cgalHasModelsBare{`CGAL::HDVF::Zp`}
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


    /// \name Input/output
    /// @{

    /*!
     * \brief Output operator
     */
    std::ostream& operator<< (std::ostream& out, const Ring& v);

    /*!
     * \brief Input operator
     */
    std::istream& operator>> (std::istream& out, const Ring& v);

    /// @}

/// \name Invertibility
/// @{

    /*!
     \brief Test if `v` is invertible in the ring.
     */
    bool is_invertible(Value v);

/// @}
};
