// Copyright (c) 2024 INRIA Nancy - Grand Est (France). LIGM Marne-la-Vallée (France)
// All rights reserved.

namespace CGAL {

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

Represents an isometry in the Poincaré disk model.
The isometry f is represented by a list (c0, c1, c2, c3) of complex numbers,
so that f(z) = (c0 z + c1) / (c2 z + c3) holds on every complex z in the open unit disk.

Facilities are offered to compose isometries, and apply an isometry to a point.

\tparam Traits_2 is the traits class and must be a model of `HyperbolicSurfacesTraits_2`.
*/
template<class Traits_2>
class Hyperbolic_isometry_2{
  public:
    /// \name Creation
    /// @{
    /*!
      Default constructor.
    */
    Hyperbolic_isometry_2();
    /// @}

    /// \name Get and set
    /// @{
    /*!
      Set the isometry to the identity.
    */
    void set_to_identity();

    /*!
      Can be used to set the coefficients of the isometry manually, be careful when doing so : the implementation does not check that the resulting moebius transform fixes the unit circle.
    */
    void set_coefficients(const Traits_2::Complex& c0, const Traits_2::Complex& c1, const Traits_2::Complex& c2, const Traits_2::Complex& c3);

    /*!
      Can be used to set one coefficient of the isometry manually, be careful when doing so : the implementation does not check that the resulting moebius transform fixes the unit circle.
    */
    void set_coefficient(int index, const Traits_2::Complex& coefficient);

    /*!
      Returns the index-th coefficient.
    */
    Traits_2::Complex get_coefficient(int index) const;
    /// @}

    /// \name Operations
    /// @{
    /*!
      Evaluates the isometry at the point.
    */
    Traits_2::Hyperbolic_point_2 evaluate(const Traits_2::Hyperbolic_point_2& point) const;

    /*!
      Returns the composition of itself and other.
    */
    Hyperbolic_isometry_2<Traits_2> compose(const Hyperbolic_isometry_2<Traits_2>& other) const;
    /// @}

};

} //namespace CGAL
