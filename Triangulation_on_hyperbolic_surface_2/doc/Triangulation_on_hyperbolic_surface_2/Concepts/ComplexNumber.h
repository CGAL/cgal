/*!
  \ingroup PkgHyperbolicSurfaceTriangulation2Concepts
  \cgalConcept

  Describes a complex number type over a `FieldNumberType` for its real and imaginary parts.

  \cgalRefines{Field}

  \cgalHasModelsBegin
  \cgalHasModels{CGAL::Complex_number}
  \cgalHasModelsEnd
*/
class ComplexNumber
{
public:
  /// \name Types
  /// @{

  /*!
    Number type for real and imaginary parts: must be a model of `FieldNumberType`.
  */
  typedef unspecified_type FT;

  /// @}

  /// \name Creation
  /// @{

  /*!
    Default constructor, sets the both the real part and the imaginary part to \f$ 0 \f$.
  */
  ComplexNumber();

  /*!
    Constructor, sets the real part to <code> real_part </code> and the imaginary part to \f$ 0 \f$.
  */
  ComplexNumber(const FT& real_part);

  /*!
    Constructor, sets the real part to <code> real_part </code> and the imaginary part to <code> imaginary_part </code>.
  */
  ComplexNumber(const FT& real_part, const FT& imaginary_part);

  /*!
    Constructor, sets the real part to <code> real_part </code> and the
    imaginary part to <code> imaginary_part </code>. FT must be
    constructible from U and V.
  */
  template<class U,class V>
  ComplexNumber(U&& real_part, V&& imaginary_part);

  /// @}

  /// \name Getter and Setter
  /// @{

  /*!
    sets the real part to <code> real_part </code>.
  */
  void real(const FT& real_part);

  /*!
    sets the imaginary part to <code> imaginary_part </code>.
  */
  void imag(const FT& imaginary_part);

  /*!
    returns the real part.
  */
  FT real() const;

  /*!
    returns the imaginary part.
  */
  FT imag() const;

  /// @}

  /// \name Operations
  /// @{

  /*!
    Copy operator.
  */
  ComplexNumber operator=(const ComplexNumber& other) const;

  /*!
    writes the complex in a stream.
  */
  std::ostream& operator<<(std::ostream& s, const ComplexNumber& z);

  /*!
    reads the complex from a stream.
  */
  void operator>>(std::istream& s, ComplexNumber& z);

  /// @}

  /// \relates ComplexNumber
  /// @{

  /*!
    returns the square of the modulus.
  */
  FT norm(ComplexNumber z) const;

  /*!
    returns the conjugate.
  */
  ComplexNumber conj(ComplexNumber z) const;

  /// @}
};
