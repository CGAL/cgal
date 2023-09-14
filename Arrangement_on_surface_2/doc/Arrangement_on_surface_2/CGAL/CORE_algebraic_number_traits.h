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

  /*! Obtain the numerator of a rational number.
   * \param q A rational number.
   * \return The numerator of q.
   */
  Integer numerator(const Rational& q) const;

  /*! Obtain the denominator of a rational number.
   * \param q A rational number.
   * \return The denominator of q.
   */
  Integer denominator(const Rational& q) const;

  /*! Convert an integer to an algebraic number.
   * \param z An integer.
   * \return The algebraic number equivalent to z.
   */
  Algebraic convert(const Integer& z) const;

  /*! Convert a rational number to an algebraic number.
   * \param q A rational number.
   * \return The algebraic number equivalent to q.
   */
  Algebraic convert(const Rational& q) const;

  /*! Construct a rational number that lies strictly between two algebraic
   * values.
   * \param x1 The first algebraic value.
   * \param x2 The second algebraic value.
   * \pre The two values are not equal.
   * \return A rational number that lies in the open interval (x1, x2).
   */
  Rational rational_in_interval(const Algebraic& x1, const Algebraic& x2) const;

  /*! Obtain a range of double-precision floats that contains the given
   * algebraic number.
   * \param x The given number.
   * \return A pair <x_lo, x_hi> that contain x.
   */
  std::pair<double, double> double_interval(const Algebraic& x) const;

  /*! Convert a sequence of rational coefficients to an equivalent sequence
   * of integer coefficients. If the input coefficients are q(1), ..., q(k),
   * where q(i) = n(i)/d(i) then the output coefficients will be of the
   * form:
   *               n(i) * lcm {d(1), ... , d(k)}
   *       a(i) = -------------------------------
   *               d(i) * gcd {n(1), ... , n(k)}
   *
   * \param q_begin The begin iterator of the rational sequence.
   * \param q_end The past-the-end iterator of the rational sequence.
   * \param zoi An output iterator for the integer coefficients.
   * \return A past-the-end iterator for the output container.
   * \pre The value type of q_begin and q_end is `Rational`, and
   *      the value type of zoi is `Integer`.
   */
  template <typename InputIterator, typename OutputIterator>
  OutputIterator convert_coefficients(InputIterator q_begin,
                                      InputIterator q_end,
                                      OutputIterator zoi) const;

  /*! Compute the square root of an algebraic number.
   * \param x The number.
   * \return The square root of x.
   * \pre x is non-negative.
   */
  Algebraic sqrt(const Algebraic& x) const;

  /*! Compute the roots of a quadratic equations \f$a*x^2+ b*x + c = 0\f$
   * with integer coefficients.
   * \param a The coefficient of \f$x^2\f$
   * \param b The coefficient of \f$x\f$
   * \param c The free term.
   * \param oi An output iterator for the real-valued solutions of the
   *           quadratic equation.
   * \return A past-the-end iterator for the output container.
   * \pre The value type of oi is Algebraic.
   */
  template <typename NT, typename OutputIterator>
  OutputIterator solve_quadratic_equation(const NT& a, const NT& b, const NT& c,
                                          OutputIterator oi) const;

  /*! Construct a polynomial with integer coefficients.
   * \param coeffs The coefficients of the input polynomial.
   * \param degree The degree of the input polynomial.
   * \return The polynomial.
   */
  Polynomial construct_polynomial(const Integer* coeffs,
                                  unsigned int degree) const;

  /*! Construct a polynomial with integer coefficients given rational
   * coefficients.
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

  /*! Construct two polynomials with integer coefficients such that P(x)/Q(x)
   * is a rational function equivalent to the one represented by the two
   * given vectors of rational coefficients. It is guaranteed that the GCD
   * of P(x) and Q(x) is trivial.
   * \param p_coeffs The coefficients of the input numerator polynomial.
   * \param p_degree The degree of the input numerator polynomial.
   * \param q_coeffs The coefficients of the input denominator polynomial.
   * \param q_degree The degree of the input denominator polynomial.
   * \param p_poly Output: The resulting numerator polynomial with integer
   *                       coefficients.
   * \param q_poly Output: The resulting denominator polynomial with integer
   *                       coefficients.
   * \return (true) on success; (false) if the denominator is 0.
   */
  bool construct_polynomials(const Rational* p_coeffs,
                             unsigned int p_degree,
                             const Rational* q_coeffs,
                             unsigned int q_degree,
                             Polynomial& p_poly, Polynomial& q_poly) const;

  /*! Compute the degree of a polynomial.
   */
  int degree(const Polynomial& poly) const;

  /*! Evaluate a polynomial at a given x-value.
   * \param poly A polynomial.
   * \param x The value to evaluate at.
   * \return The value of the polynomial at x.
   */
  template <typename NT>
  NT evaluate_at(const Polynomial& poly, NT& x) const;

  /*! Compute the derivative of the given polynomial.
   * \param poly The polynomial p(x).
   * \return The derivative p'(x).
   */
  Polynomial derive(const Polynomial& poly) const;

  /*! Multiply a polynomial by some scalar coefficient.
   * \param poly The polynomial P(x).
   * \param a The scalar value.
   * \return The scalar multiplication a*P(x).
   */
  Polynomial scale(const Polynomial& poly, const Integer& a) const;

  /*! Perform "long division" of two polynomials: Given A(x) and B(x) compute
   * two polynomials Q(x) and R(x) such that: A(x) = Q(x)*B(x) + R(x) and
   * R(x) has minimal degree.
   * \param polyA The first polynomial A(x).
   * \param polyB The second polynomial A(x).
   * \param rem Output: The remainder polynomial R(x).
   * \return The quontient polynomial Q(x).
   */
  Polynomial divide(const Polynomial& polyA,
                    const Polynomial& polyB,
                    Polynomial& rem) const;

  /*! Compute the real-valued roots of a polynomial with integer coefficients,
   * sorted in ascending order.
   * \param poly The input polynomial.
   * \param oi An output iterator for the real-valued root of the polynomial.
   * \return A past-the-end iterator for the output container.
   * \pre The value type of oi is Algebraic.
   */
  template <typename OutputIterator>
  OutputIterator compute_polynomial_roots(const Polynomial& poly,
                                          OutputIterator oi) const;

  /*! Compute the real-valued roots of a polynomial with integer coefficients,
   * within a given interval. The roots are sorted in ascending order.
   * \param poly The input polynomial.
   * \param x_min The left bound of the interval.
   * \param x_max The right bound of the interval.
   * \param oi An output iterator for the real-valued root of the polynomial.
   * \return A past-the-end iterator for the output container.
   * \pre The value type of oi is Algebraic.
   */
  template <typename OutputIterator>
  OutputIterator compute_polynomial_roots(const Polynomial& poly,
                                          double x_min, double x_max,
                                          OutputIterator oi) const;
  /// @}
};

} /* end namespace CGAL */
