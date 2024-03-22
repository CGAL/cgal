// Copyright (c) 2024 INRIA Nancy - Grand Est (France). LIGM Marne-la-Vallée (France)
// All rights reserved.

namespace CGAL {

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

Represents an isometry in the Poincare disk model.
The isometry f is stored as list (c0, c1, c2, c3) of four complex numbers,
so that f(z) = (c0 z + c1) / (c2 z + c3) holds on every complex z in the open unit disk.

\tparam GeometricTraits_2 is the geometric traits class and must be a model of `HyperbolicSurfacesTraits_2`.
*/
template<class GeometricTraits_2>
class Hyperbolic_isometry_2{
  private:
    typedef Hyperbolic_isometry_2<GeometricTraits_2>                    _Self;
    typedef Complex_without_sqrt<_FT>                                   _Cmplx;
    typedef typename GeometricTraits_2::Point_2                         _Point;

  public:
    /// \name Types
    /// @{
    typedef GeometricTraits_2 Geometric_traits_2;
    /// @}

    /// \name Creation
    /// @{
    /*!
      Default constructor
    */
    Hyperbolic_isometry_2();
    /// @}

    /// \name Get and set
    /// @{
    /*!
      Set the isometry to the identity
    */
    void set_to_identity();

    /*!
      Can be used to set the coefficients of the isometry manually.
      Be careful when doing so : the implementation does not check that the resulting moebius transform fixes the unit circle.
    */
    void set_coefficients(const _Cmplx& c0, const _Cmplx& c1, const _Cmplx& c2, const _Cmplx& c3);

    /*!
      Can be used to set one coefficient of the isometry smanually.
      Be careful when doing so : the implementation does not check that the resulting moebius transform fixes the unit circle.
    */
    void set_coefficient(int index, const _Cmplx& coefficient);

    /*!
      Returns the index-th coefficient
    */
    _Cmplx get_coefficient(int index) const;
    /// @}

    /// \name Operations
    /// @{
    /*!
      Evaluates the isometry at point
    */
    _Point evaluate(const _Point& point) const;

    /*!
      Returns the composition of itself and other
    */
    _Self compose(const _Self& other) const;
    /// @}

};

} //namespace CGAL
