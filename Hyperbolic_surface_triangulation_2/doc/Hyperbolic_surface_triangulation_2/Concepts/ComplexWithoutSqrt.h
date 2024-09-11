// Copyright (c) 2024 INRIA Nancy - Grand Est (France). LIGM Marne-la-Vallée (France)
// All rights reserved.

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2Concepts
\cgalConcept

Describes a complex number type that does not use square root.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Complex_without_sqrt}
\cgalHasModelsEnd
*/
class ComplexWithoutSqrt {
public:
  /// \name Types
  /// @{
  /*!
    Field number type: must be a model of `FieldNumberType`.
  */
  typedef unspecified_type FT;
  /// \name Creation
  /// @{
  /*!
    Default constructor, sets the both the real part and the imaginary part to \f$ 0 \f$.
  */
  ComplexWithoutSqrt();

  /*!
    Constructor, sets the real part to <code> real_part </code> and the imaginary part to \f$ 0 \f$.
  */
  ComplexWithoutSqrt(const FT& real_part);

  /*!
    Constructor, sets the real part to <code> real_part </code> and the imaginary part to <code> imaginary_part </code>.
  */
  ComplexWithoutSqrt(const FT& real_part, const FT& imaginary_part);
  /// @}

  /// \name Get and Set
  /// @{
  /*!
    sets the real part to <code> real_part </code>.
  */
  void set_real_part(const FT& real_part);

  /*!
    sets the imaginary part to <code> imaginary_part </code>.
  */
  void set_imaginary_part(const FT& imaginary_part);

  /*!
    returns the real part.
  */
  FT real_part() const;

  /*!
    returns the imaginary part.
  */
  FT imaginary_part() const;
    /// @}

  /// \name Operations
  /// @{
  /*!
    returns the square of the modulus.
  */
  FT squared_modulus() const;

  /*!
    returns the conjugate.
  */
  ComplexWithoutSqrt<FT> conjugate() const;

  /*!
    Complex addition.
  */
  ComplexWithoutSqrt<FT> operator+(const ComplexWithoutSqrt<FT>& other) const;

  /*!
    Complex substraction.
  */
  ComplexWithoutSqrt<FT> operator-(const ComplexWithoutSqrt<FT>& other) const;

  /*!
    returns the opposite.
  */
  ComplexWithoutSqrt<FT> operator-() const;

  /*!
    Complex multiplication.
  */
  ComplexWithoutSqrt<FT> operator*(const ComplexWithoutSqrt<FT>& other) const;

  /*!
    Complex division.
  */
  ComplexWithoutSqrt<FT> operator/(const ComplexWithoutSqrt<FT>& other) const;

  /*!
    Equality test.
  */
  bool operator==(const Complex_without_sqrt<FT>& z1, const Complex_without_sqrt<FT>& z2);
  /*!
    Inequality test.
  */
  bool operator!=(const Complex_without_sqrt<FT>& z1, const Complex_without_sqrt<FT>& z2);
  /*!
    writes the complex in a stream.
  */
  std::ostream& operator<<(std::ostream& s, const Complex_without_sqrt<FT>& z);
  /*!
    reads the complex from a stream.
  */
  void operator>>(std::istream& s, Complex_without_sqrt<FT>& z);
  /// @}
};
