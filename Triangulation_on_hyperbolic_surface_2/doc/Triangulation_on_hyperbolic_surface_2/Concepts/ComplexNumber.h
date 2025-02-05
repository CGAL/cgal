// Copyright (c) 2024 INRIA Nancy - Grand Est (France). LIGM Marne-la-Vallée (France)
// All rights reserved.

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2Concepts
\cgalConcept

Describes a complex number type over a `FieldNumberType` for its real and
imaginary parts.

\cgalRefines{Field}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Complex_number}
\cgalHasModelsEnd
*/
class ComplexNumber {
public:
  /// \name Types
  /// @{
  /*!
    Number type for real and imaginary parts: must be a model of `FieldNumberType`.
  */
  typedef unspecified_type FT;
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

  /// \name Get and Set
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
  /* /\*! */
  /*   returns +z. */
  /* *\/ */
  /* //ComplexNumber operator+(const ComplexNumber& z) const; */

  /* /\*! */
  /*    returns -z. */
  /* *\/ */
  /* //ComplexNumber operator-(const ComplexNumber& z) const; */

  /* /\*! */
  /*   Unary complex addition. */
  /* *\/ */
  /* //ComplexNumber operator+=(const ComplexNumber& other) const; */

  /* /\*! */
  /*   Unary complex substraction. */
  /* *\/ */
  /* //ComplexNumber operator-=(const ComplexNumber& other) const; */

  /* /\*! */
  /*   Unary complex multiplication. */
  /* *\/ */
  /* //ComplexNumber operator*=(const ComplexNumber& other) const; */

  /* /\*! */
  /*   Unary complex division. */
  /* *\/ */
  /* //ComplexNumber operator/=(const ComplexNumber& other) const; */

  /*!
    Copy operator.
  */
  ComplexNumber operator=(const ComplexNumber& other) const;

  /* /\*! */
 /*    Equality test. */
 /*  *\/ */
 /*  //bool operator==(const ComplexNumber& z1, const ComplexNumber& z2); */
 /*  /\*! */
 /*    Inequality test. */
 /*  *\/ */
 /*  // bool operator!=(const ComplexNumber& z1, const ComplexNumber& z2); */

 /*  /\*! */
 /*    Binary complex addition. */
 /*  *\/ */
 /*  //ComplexNumber operator+(const ComplexNumber& z1, const ComplexNumber& z2); */

 /*  /\*! */
 /*    Binary complex substraction. */
 /*  *\/ */
 /*  //ComplexNumber operator-(const ComplexNumber& z1, const ComplexNumber& z2); */

 /* /\*! */
 /*    Binary complex multiplication. */
 /*  *\/ */
 /*  //ComplexNumber operator*(const ComplexNumber& z1, const ComplexNumber& z2); */

 /*  /\*! */
 /*    Binary complex division. */
 /*  *\/ */
 /*  //ComplexNumber operator/(const ComplexNumber& z1, const ComplexNumber& z2); */

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
