 namespace CGAL {

/*!
  \ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

  represents an isometry in the Poincaré disk model.
  The isometry \f$ f \f$ is represented by a list \f$ (c_0, c_1, c_2, c_3) \f$ of complex numbers,
  so that \f$ f(z) = (c_0 z + c_1) / (c_2 z + c_3) \f$ holds on every complex \f$ z \f$ in the open unit disk.

  Functionalities are offered to compose isometries, and apply an isometry to a point.

  \tparam Traits must be a model of `HyperbolicSurfaceTraits_2`.
*/
template<class Traits>
class Hyperbolic_isometry_2
{
public:
  /// \name Types
  /// @{

  /*!
    Complex number type.
  */
  typedef typename Traits::Complex                         Complex_number;

  /*!
    Point type.
  */
  typedef typename Traits::Hyperbolic_point_2              Point;

  /// @}

  /// \name Creation
  /// @{

  /*!
    Default constructor to the identity.
  */
  Hyperbolic_isometry_2();

  /*!
    Constructor from coefficients.
  */
  Hyperbolic_isometry_2(const Complex_number& c0,
                        const Complex_number& c1,
                        const Complex_number& c2,
                        const Complex_number& c3);
  /// @}

  /*!
    sets the isometry to the identity.
  */
  void set_to_identity();

  /*!
    sets the coefficients of the isometry.

    \warning The implementation does not check that the resulting transformation is an isometry.
  */
  void set_coefficients(const Complex_number& c0,
                        const Complex_number& c1,
                        const Complex_number& c2,
                        const Complex_number& c3);

  /*!
    sets a particular coefficient of the isometry.

    \warning The implementation does not check that the resulting transformation is an isometry.
  */
  void set_coefficient(int index, const Complex_number& coefficient);

  /// \name Access Functions
  /// @{

  /*!
    returns the index-th coefficient.
  */
  const Complex_number& get_coefficient(int index) const;

  /// @}

  /// \name Operations
  /// @{

  /*!
    evaluates the isometry at point \f$ p \f$.
  */
  Point evaluate(const Point& p) const;

  /*!
    evaluates the isometry at point \f$ p \f$.
  */
  Point operator()(const Point& p) const;

  /*!
    returns the composition of two isometries.
  */
  template<class Traits>
  Hyperbolic_isometry_2<Traits> operator*(const Hyperbolic_isometry_2<Traits>& iso1,
                                          const Hyperbolic_isometry_2<Traits>& iso2);

  /// @}
};

} // namespace CGAL
