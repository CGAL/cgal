// Copyright (c) 2024 INRIA Nancy - Grand Est (France). LIGM Marne-la-Vallée (France)
// All rights reserved.

namespace CGAL {

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

Represents a complex number over a field.

\cgalModels{CGAL::ComplexWithoutSqrt}

\tparam FT is the field type and must be a model of `FieldNumberType`.
*/
template <class FT>
class Complex_without_sqrt {
public:
  /// \name Creation
  /// @{
  /*!
    Default constructor, sets the complex number to 0 + 0 * i.
  */
  Complex_without_sqrt();

  /*!
    Constructor from the real part, sets the complex number to 'real_part' + 0 * i.
  */
  Complex_without_sqrt(const FT& real_part);

  /*!
    Constructor from the real and imaginary parts, sets the complex number to 'real_part' + 'imaginary_part' * i.
  */
  Complex_without_sqrt(const FT& real_part, const FT& imaginary_part);
  /// @}

  /// \name Access functions
  /// @{
  /*!
    Sets the real part to 'real_part'
  */
  void set_real_part(const FT& real_part);

  /*!
    Sets the imaginary part to 'imaginary_part'
  */
  void set_imaginary_part(const FT& imaginary_part);

  /*!
    Returns the real part
  */
  FT real_part() const;

  /*!
    Returns the imaginary part
  */
  FT imaginary_part() const;
  /// @}

  /// \name Operations
  /// @{
  /*!
    Returns the square of the modulus
  */
  FT squared_modulus() const;

  /*!
    Returns the conjugate
  */
  Complex_without_sqrt<FT> conjugate() const;

  /*!
    Returns the sum of itself and other
  */
  Complex_without_sqrt<FT> operator+(const Complex_without_sqrt<FT>& other) const;

  /*!
    Returns the difference of itself and other
  */
  Complex_without_sqrt<FT> operator-(const Complex_without_sqrt<FT>& other) const;

  /*!
    Returns the opposite of itself
  */
  Complex_without_sqrt<FT> operator-() const;

  /*!
    Returns the multiplication of itself and other
  */
  Complex_without_sqrt<FT> operator*(const Complex_without_sqrt<FT>& other) const;

  /*!
    Returns the division of itself by other
  */
  Complex_without_sqrt<FT> operator/(const Complex_without_sqrt<FT>& other) const;
  /// @}

  /// \name Equality test
  /// @{
  /*!
    Equality test operator.
  */
  template<class FT> bool operator==(const Complex_without_sqrt<FT>& z1, const Complex_without_sqrt<FT>& z2);
  /*!
    Inequality test operator.
  */
  template<class FT> bool operator!=(const Complex_without_sqrt<FT>& z1, const Complex_without_sqrt<FT>& z2);
  /// @}

  /// \name Input/output
  /// @{
  /*!
    Outputs the complex number in a stream.
  */
  template<class FT>std::ostream& operator<<(std::ostream& s, const Complex_without_sqrt<FT>& z);

  /*!
    Reads the complex number from a stream.
  */
  template<class FT>void operator>>(std::istream& s, Complex_without_sqrt<FT>& z);
  /// @}
};


} //namespace CGAL
