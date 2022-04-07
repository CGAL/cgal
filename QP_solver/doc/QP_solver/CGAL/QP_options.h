
namespace CGAL {

/*!
\ingroup PkgQPSolverFunctions

This is a class used for passing options to the linear and
quadratic programming solvers. Currently, we support only
options referring to
<OL>
<LI>the verbosity,
<LI>the pricing strategy (see `Quadratic_program_pricing_strategy`),
<LI>the validation mode (see the Validity section of
`Quadratic_program_solution`)
</OL>
The idea is that this list grows in the future.

\cgalHeading{Operations}

Here we just have set/get pairs for any option type.

\sa `Quadratic_program_solution`
\sa `solve_quadratic_program`
\sa `solve_linear_program`
\sa `solve_nonnegative_quadratic_program`
\sa `solve_nonnegative_linear_program`

*/

class Quadratic_program_options {
public:

/// \name Creation
/// @{

/*!
constructs an instance of `Quadratic_program_options` where all available options
are at their defaults.
*/
Quadratic_program_options();

/// @}

/// \name Verbosity
/// @{

/*!
sets the verbosity of the solver to the value `verbosity` when
`options` is passed to any of the four solution functions. The provided
value must be a number between \f$ 0\f$ and \f$ 5\f$. Verbosity \f$ 0\f$ is the default and
results in the solver running silently. Verbosity \f$ 1\f$ prints a short
summary of every iteration. Higher verbosity values print more information
about the solution process, but these are mainly for debugging
purposes and have no effect if you compile with
<TT>CGAL_QP_NO_ASSERTIONS</TT> or <TT>NDEBUG</TT>.
*/
void set_verbosity (int verbosity);

/*!
returns the verbosity level of `options`.
*/
int get_verbosity () const;

/// @}

/// \name Pricing strategy
/// @{

/*!
sets the pricing strategy of the solver to the value `pricing_strategy`
when `options` is passed to any of the four solution functions. The pricing
strategy controls how the solver proceeds from any intermediate solution.
For the available strategies and their behavior, see the documentation of the
class `Quadratic_program_pricing_strategy`.
*/
void set_pricing_strategy
(Quadratic_program_pricing_strategy pricing_strategy);

/*!
returns the pricing strategy of `options`.
*/
Quadratic_program_pricing_strategy get_pricing_strategy() const;

/// @}

/// \name Validation mode
/// @{

/*!
sets the automatic validation mode of the solver to the value `validate`.
The default is `false`. By providing value `true` you can
tell the solver to automatically check whether the program has
correctly been solved, see the Validity section of the class
`Quadratic_program_solution`.
*/
void set_auto_validation(bool validate);

/*!
returns the validation mode of `options`.
*/
bool get_auto_validation() const;

/// @}

}; /* end Quadratic_program_options */

/*!
  \ingroup PkgQPSolverFunctions

  This is an enumeration type containing the values
  `QP_CHOOSE_DEFAULT`, `QP_DANTZIG`,
  `QP_PARTIAL_DANTZIG`, `QP_FILTERED_DANTZIG`,
  `QP_PARTIAL_FILTERED_DANTZIG`, and`QP_BLAND`.

  It indicates the pricing strategy to be used in
  solving a linear or quadratic program. This strategy determines
  how the solver gets from one intermediate solution to the next
  during any of its iterations.

  Here we briefly describe when to choose which strategy.

  \sa `Quadratic_program_options`
*/
enum Quadratic_program_pricing_strategy {
  /*!
    This is the default value of the pricing strategy in
    `Quadratic_program_options`, and it lets the solver choose the
    strategy that it thinks is most appropriate for the problem at hand.
    There are only few reasons to deviate from this default, but you are
    free to experiment, of course.
  */
  QP_CHOOSE_DEFAULT,

  /*!
    If the input type is <B>not</B> `double`, this is usually the
    best choice for linear and quadratic programs of medium size.
  */
  QP_PARTIAL_DANTZIG,

  /*!
    If the input type is <B>not</B> `double`, this can sometimes
    make a difference (be faster or slowe) than `QP_PARTIAL_DANTZIG`
    for problems with a high variable/constraint or constraint/variable ratio.
  */
  QP_DANTZIG,

  /*!
    If the input type <B>is</B> `double`, this is usually the best choice
    for linear and quadratic programs of medium size.
    If the input type is not `double`, this choice is equivalent
    to `QP_PARTIAL_DANTZIG`.

    <B>Note:</B> filtered strategies may in rare cases fail due to double
    exponent overflows, see
    Section \ref secQPcustomizationfiltering.
    In this case, the slower fallback option is
    the non-filtered variant `QP_PARTIAL_DANTZIG` of this strategy.
  */
  QP_PARTIAL_FILTERED_DANTZIG,

  /*!
    If the input type <B>is</B> `double`, this can sometimes
    make a difference (be faster or slowe) than `QP_PARTIAL_FILTERED_DANTZIG`
    for problems with a high variable/constraint or constraint/variable ratio.
    If the input type is not `double`, this choice is equivalent
    to `QP_DANTZIG`.

    \note Filtered strategies may in rare cases fail due to double
    exponent overflows, see Section \ref secQPcustomizationfiltering. In
    this case, the slower fallback option is the non-filtered variant
    `QP_DANTZIG` of this strategy.
  */
  QP_FILTERED_DANTZIG,

  /*!
    This is hardly ever the most efficient choice, but it is guaranteed
    to avoid internal cycling of the solution algorithm, see
    Section \ref secQPcustomizationcycling.
  */
  QP_BLAND
}; /* end Quadratic_program_pricing_strategy */
} /* end namespace CGAL */
