namespace CGAL{

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

Represents an isometry in the Poincaré disk model.
The isometry \f$ f \f$ is represented by a list \f$ (c_0, c_1, c_2, c_3) \f$ of complex numbers,
so that \f$ f(z) = (c_0 z + c_1) / (c_2 z + c_3) \f$ holds on every complex \f$ z \f$ in the open unit disk.

Facilities are offered to compose isometries, and apply an isometry to a point.

\tparam Traits is the traits class and must be a model of `HyperbolicSurfaceTraits_2` (default model: `Hyperbolic_surface_traits_2`).
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

    /*!
      sets the isometry to the identity.
    */
    void set_to_identity();

    /*!
      sets the coefficients of the isometry manually. \note  Be
      careful when doing so: the implementation does not check that the
      resulting Möbius transform fixes the unit circle.

    */
    void set_coefficients(const ComplexNumber& c0, const ComplexNumber& c1, const ComplexNumber& c2, const ComplexNumber& c3);

    /*!
      sets a particular coefficient of the isometry manually. \note  Be
      careful when doing so: the implementation does not check that the
      resulting Möbius transform fixes the unit circle.

    */
    void set_coefficient(int index, const ComplexNumber& coefficient);

    /// \name Access Functions
    /// @{
    /*!
      returns the index-th coefficient.
    */
    ComplexNumber get_coefficient(int index) const;
    /// @}

    /// \name Operations
    /// @{
    /*!
      evaluates the isometry at the point.
    */
    Point evaluate(const Point& point) const;

    /// @{
    /*!
      returns the composition of two isometries.
    */
    template<class Traits>
      Hyperbolic_isometry_2<Traits>  operator*(const Hyperbolic_isometry_2<Traits>& iso1, const Hyperbolic_isometry_2<Traits>& iso2);
    /// @}

    /// \name Input/Output
    /// @{
    /*!
    writes the isometry in a stream.
    */
    std::ostream& operator<<(std::ostream& s, const Hyperbolic_isometry_2<Traits>& isometry);
    /// @}
};

}; // namespace CGAL
