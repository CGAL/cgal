/*!
\ingroup PkgHDVFConcepts
\cgalConcept

The concept `Ring` describes the requirements for the ring of coefficients used to compute homology in the `HDVF` concept. Besides ring operators, it also specifies the functions needed to test invertibility in the ring.

\cgalHasModelsBegin
\cgalHasModelsBare{`CGAL::HDVF::Z2`}
\cgalHasModelsBare{`CGAL::HDVF::Zp`}
\cgalHasModelsEnd

*/

class Ring
{
public:

/// \name Operators
/// @{

    /*! TODO */

/// @}


    /// \name Input/output
    /// @{

    /*!
     * \brief Output operator
     */
    std::ostream& operator<< (std::ostream& out, const Ring& r);

    /*!
     * \brief Input operator
     */
    std::istream& operator>> (std::istream& out, const Ring& r);

    /// @}

/// \name Invertibility
/// @{

    /*!
     \brief Test if the current element is invertible in the ring.
     */
    bool is_invertible();

/// @}
};
