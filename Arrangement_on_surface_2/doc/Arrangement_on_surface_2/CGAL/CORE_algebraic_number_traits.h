namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2TraitsClasses
 *
 * `CORE_algebraic_number_traits` is a traits class for CORE's algebraic
 * number types.
 *
 * \sa `Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>`
 */
class CORE_algebraic_number_traits {
public:
  /// \name Types
  /// @{

  //! The integer number type.
  typedef CORE::BigInt                    Integer;

  //! The rational number type.
  typedef CORE::BigRat                    Rational;

  //! The polynomial type.
  typedef CORE::Polynomial<Integer>       Polynomial;

  //! The algebraic number type.
  typedef CORE::Expr                      Algebraic;

  /// @}

  /// \name Utility Functions
  /// @{

  /*! obtains the numerator of a rational number.
   *
   * \param q The rational number.
   * \return The numerator of `q`.
   */
  Integer numerator(const Rational& q) const;

  /*! obtains the denominator of a rational number.
   *
   * \param q The rational number.
   * \return The denominator of `q`.
   */
  Integer denominator(const Rational& q) const;

  /*! converts an integer to an algebraic number.
   *
   * \param z The integer.
   * \return The algebraic number equivalent to `z`.
   */
  Algebraic convert(const Integer& z) const;

  /*! converts a rational number to an algebraic number.
   *
   * \param q A rational number.
   * \return The algebraic number equivalent to `q`.
   */
  Algebraic convert(const Rational& q) const;

  /*! constructs a rational number that lies strictly between two algebraic
   * values.
   *
   * \param x1 The first algebraic value.
   * \param x2 The second algebraic value.
   * \pre The two values are not equal.
   * \return The rational number that lies in the open interval (`x1`, `x2`).
   */
  Rational rational_in_interval(const Algebraic& x1, const Algebraic& x2) const;

  /*! obtains a range of double-precision floats that contains the given
   * algebraic number.
   *
   * \param x The given number.
   * \return The range of double-precision floats that contain `x`.
   */
  std::pair<double, double> double_interval(const Algebraic& x) const;

  /*! converts a sequence of rational coefficients to an equivalent sequence
   * of integer coefficients. If the input coefficients are
   * \f$q(1),\ldots,q(k)\f$, where \f$q(i) = n(i)/d(i)\f$, then the output
   * coefficients will be of the form:
   * \f$a(i) = \frac{n(i) \cdot \mathrm{lcm}(d(1),\ldots,d(k))}{d(i) \cdot \mathrm{gcd}(n(1),\ldots, n(k))}\f$.
   * It inserts the output sequence into an output container given through an
   * output iterator.
   *
   * \param begin The begin iterator of the rational coefficients input
   *        container.
   * \param end The past-the-end iterator of the rational coefficients input
   *        container.
   * \param oi The output iterator of the integer coefficients output container.
   * \return The past-the-end iterator of the output container.
   *
   * \pre The value type of `InputIterator` is `Rational`.
   * \pre Dereferencing `oi` must yield an object convertible to `Integer`.
   */
  template <typename InputIterator, typename OutputIterator>
  OutputIterator convert_coefficients(InputIterator begin,
                                      InputIterator end,
                                      OutputIterator oi) const;

  /*! computes the square root of an algebraic number.
   *
   * \param x The number.
   * \return The square root of `x`.
   * \pre `x` is non-negative.
   */
  Algebraic sqrt(const Algebraic& x) const;

  /*! computes the roots of a quadratic equations \f$a*x^2+ b*x + c = 0\f$
   * with integer coefficients, and inserts them into an output container given
   * through an output iterator.
   *
   * \param a The coefficient of \f$x^2\f$
   * \param b The coefficient of \f$x\f$
   * \param c The free term.
   * \param oi The output iterator of the output container of real-valued
   *        solutions of the quadratic equation.
   * \return The past-the-end iterator of the output container.
   *
   * \pre Dereferencing `oi` must yield an object convertible to `Algebraic`.
   */
  template <typename NT, typename OutputIterator>
  OutputIterator solve_quadratic_equation(const NT& a, const NT& b, const NT& c,
                                          OutputIterator oi) const;

  /*! constructs a polynomial with integer coefficients.
   *
   * \param coeffs The coefficients of the input polynomial.
   * \param degree The degree of the input polynomial.
   * \return The polynomial.
   */
  Polynomial construct_polynomial(const Integer* coeffs,
                                  unsigned int degree) const;

  /*! constructs a polynomial with integer coefficients given rational
   * coefficients.
   *
   * \param coeffs The coefficients of the input polynomial.
   * \param degree The degree of the input polynomial.
   * \param poly Output: The resulting polynomial with integer coefficients.
   * \param poly_denom Output: The denominator for the polynomial.
   * \return Whether this polynomial is non-zero (false if the polynomial is
   *         zero).
   */
  bool construct_polynomial(const Rational *coeffs,
                            unsigned int degree,
                            Polynomial& poly, Integer& poly_denom) const;

  /*! constructs two polynomials with integer coefficients such that
   * \f$P(x)/Q(x)\f$ is a rational function equivalent to the one represented
   * by the two given vectors of rational coefficients. It is guaranteed that
   * the GCD of \f$P(x)\f$ and \f$Q(x)\f$ is trivial.
   *
   * \param p_coeffs The coefficients of the input numerator polynomial.
   * \param p_degree The degree of the input numerator polynomial.
   * \param q_coeffs The coefficients of the input denominator polynomial.
   * \param q_degree The degree of the input denominator polynomial.
   * \param p_poly Output: The resulting numerator polynomial with integer
   *                       coefficients.
   * \param q_poly Output: The resulting denominator polynomial with integer
   *                       coefficients.
   * \return `true` on success; `false` if the denominator is 0.
   */
  bool construct_polynomials(const Rational* p_coeffs,
                             unsigned int p_degree,
                             const Rational* q_coeffs,
                             unsigned int q_degree,
                             Polynomial& p_poly, Polynomial& q_poly) const;

  /*! Compute the degree of a polynomial.
   */
  int degree(const Polynomial& poly) const;

  /*! evaluates a polynomial at a given \f$x\f$-value.
   *
   * \param poly The polynomial.
   * \param x The value to evaluate at.
   * \return The value of the polynomial at `x`.
   */
  template <typename NT>
  NT evaluate_at(const Polynomial& poly, NT& x) const;

  /*! computes the derivative of the given polynomial.
   *
   * \param poly The polynomial \f$p(x)\f$.
   * \return The derivative \f$p'(x)\f$.
   */
  Polynomial derive(const Polynomial& poly) const;

  /*! multiplies a polynomial by some scalar coefficient.
   *
   * \param poly The polynomial \f$P(x)\f$.
   * \param a The scalar value.
   * \return The scalar multiplication \f$a \cdot P(x)\f$.
   */
  Polynomial scale(const Polynomial& poly, const Integer& a) const;

  /*! performs "long division" of two polynomials: Given \f$A(x)\f$ and
   * \f$B(x)\f$ compute two polynomials \f$Q(x)\f$ and \f$R(x)\f$ such that:
   * \f$A(x) = Q(x) \cdot B(x) + R(x)\f$ and \f$R(x)\f$ has minimal degree.
   *
   * \param poly_a The first polynomial \f$A(x)\f$.
   * \param poly_b The second polynomial \f$A(x)\f$.
   * \param rem Output: The remainder polynomial \f$R(x)\f$.
   * \return The quontient polynomial \f$Q(x)\f$.
   */
  Polynomial divide(const Polynomial& poly_a,
                    const Polynomial& poly_b,
                    Polynomial& rem) const;

  /*! computes the real-valued roots of a polynomial with integer coefficients,
   * and inserts them in ascending order into an output container given through
   * an output iterator.
   *
   * \param poly The input polynomial.
   * \param oi The output iterator of the output container of real-valued root
   *        of the polynomial.
   * \return The past-the-end iterator of the output container.
   *
   * \pre Dereferencing `oi` must yield an object convertible to `Algebraic`.
   */
  template <typename OutputIterator>
  OutputIterator compute_polynomial_roots(const Polynomial& poly,
                                          OutputIterator oi) const;

  /*! computes the real-valued roots of a polynomial with integer coefficients
   * within a given interval, and inserts them in ascending order into an output
   * container given through an output iterator.
   *
   * \param poly The input polynomial.
   * \param x_min The left bound of the interval.
   * \param x_max The right bound of the interval.
   * \param oi The output iterator of the output container of the real-valued
   *        root of the polynomial.
   * \return The past-the-end iterator of the output container.
   *
   * \pre Dereferencing `oi` must yield an object convertible to `Algebraic`.
   */
  template <typename OutputIterator>
  OutputIterator compute_polynomial_roots(const Polynomial& poly,
                                          double x_min, double x_max,
                                          OutputIterator oi) const;
  /// @}
};

} /* end namespace CGAL */
