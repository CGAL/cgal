// Copyright (c) 2024 INRIA Nancy - Grand Est (France). LIGM Marne-la-Vallée (France)
// All rights reserved.

namespace CGAL {

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

Represents an isometry in the Poincaré disk model.
The isometry f is represented by a list (c0, c1, c2, c3) of complex numbers,
so that f(z) = (c0 z + c1) / (c2 z + c3) holds on every complex z in the open unit disk.

Facilities are offered to compose isometries, and apply an isometry to a point.

\tparam Traits is the traits class and must be a model of `HyperbolicSurfacesTraits`.
*/
template<class Traits>
class Hyperbolic_isometry_2{
  public:
    /// \name Types
    /// @{
    /*!
    Complex number type.
    */
    typedef typename Traits::Complex                         ComplexNumber;
    /*!
    Point type.
    */
    typedef typename Traits::Hyperbolic_point_2              Point;
    /// @}
    /// \name Creation
    /// @{
    /*!
      Default constructor.
    */
    Hyperbolic_isometry_2();
    /// @}

    /// \name Access functions
    /// @{
    /*!
      Set the isometry to the identity.
    */
    void set_to_identity();

    /*!
      Can be used to set the coefficients of the isometry manually, be careful when doing so : the implementation does not check that the resulting moebius transform fixes the unit circle.
    */
    void set_coefficients(const ComplexNumber& c0, const ComplexNumber& c1, const ComplexNumber& c2, const ComplexNumber& c3);

    /*!
      Can be used to set one coefficient of the isometry manually, be careful when doing so : the implementation does not check that the resulting moebius transform fixes the unit circle.
    */
    void set_coefficient(int index, const ComplexNumber& coefficient);

    /*!
      Returns the index-th coefficient.
    */
    ComplexNumber get_coefficient(int index) const;
    /// @}

    /// \name Operations
    /// @{
    /*!
      Evaluates the isometry at the point.
    */
    Point evaluate(const Point& point) const;

    /*!
      Returns the composition of itself and other.
    */
    Hyperbolic_isometry_2<Traits> compose(const Hyperbolic_isometry_2<Traits>& other) const;
    /// @}

    /// \name Input/output
    /// @{
    /*!
    Writes the isometry in a stream.
    */
    template<class Traits> std::ostream& operator<<(std::ostream& s, const Hyperbolic_isometry_2<Traits>& isometry);
    /// @}


};

} //namespace CGAL
