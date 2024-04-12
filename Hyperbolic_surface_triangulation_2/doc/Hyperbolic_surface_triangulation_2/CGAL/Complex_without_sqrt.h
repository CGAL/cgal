// Copyright (c) 2024 INRIA Nancy - Grand Est (France). LIGM Marne-la-Vallée (France)
// All rights reserved.

namespace CGAL {

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

Templated by a field type FT, represents a complex number over FT.

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
    Constructor from the real part, sets the complex number to 'real' + 0 * i.
  */
  Complex_without_sqrt(const FT& real);

  /*!
    Constructor from the real and imaginary parts, sets the complex number to 'real' + 'imaginary' * i.
  */
  Complex_without_sqrt(const FT& real, const FT& imag);
  /// @}

  /// \name Get and set
  /// @{
  /*!
    Sets the real part to 'real'
  */
  void set_real(const FT& real);

  /*!
    Sets the imaginary part to 'imag'
  */
  void set_imag(const FT& imag);

  /*!
    Returns the real part
  */
  FT real() const;

  /*!
    Returns the imaginary part
  */
  FT imag() const;
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
};

} //namespace CGAL
