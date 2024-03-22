// Copyright (c) 2024 INRIA Nancy - Grand Est (France). LIGM Marne-la-Vallée (France)
// All rights reserved.

namespace CGAL {

template<class ComplexType>
class Hyperbolic_point_2 {
  public:
    /// \name Types
    /// @{
    typedef ComplexType                           Complex;
    /// @}

    /// \name Creation
    /// @{
      /*!
        %Default constructor
      */
    Hyperbolic_point_2();

      /*!
        %Constructor from the complex number representing the point
      */
    Hyperbolic_point_2(const Complex& z);
    /// @}

    /// \name Get and set
    /// @{
      /*!
        %Get the complex number representing the point
      */
    Complex get_z() const;

      /*!
        %Set the complex number representing the point
      */
    void set_z(const Complex& z);
    /// @}
};

template<class ComplexType>
std::ostream& operator<<(std::ostream& s, const Hyperbolic_point_2<ComplexType>& point);

template<class ComplexType>
void operator>>(std::istream& s, Hyperbolic_point_2<ComplexType>& point);

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2TraitsClasses

The class `Hyperbolic_surfaces_traits_2` is designed as one of the
default models for the traits concept `HyperbolicSurfacesTraits_2`
offered by \cgal.

\tparam FieldType must be a model of `Field`.

This class provides exact constructions and predicates. The default value for `FieldType` is `CGAL::Gmpq`,
which guarantees exact constructions when the input points have rational coordinates.

\sa `Hyperbolic_surfaces_traits_2`

\cgalModels{HyperbolicSurfacesTraits_2}
*/
template<class FieldType>
class Hyperbolic_surfaces_traits_2 {
public:
  /// \name Types
  /// @{
  typedef FieldType                          FT;
  typedef Complex_without_sqrt<FieldType>    Complex;
  typedef Hyperbolic_point_2<Complex>        Point_2;
    /// @}
};

} //namespace CGAL
