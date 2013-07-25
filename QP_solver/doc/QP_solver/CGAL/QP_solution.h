
namespace CGAL {

/*!
\ingroup PkgQPSolverClasses

An object of class `Quadratic_program_solution` represents the solution of a linear or 
convex quadratic program of the general form 

\f$
\newcommand{\qprel}{\gtreqless}
\newcommand{\qpx}{\mathbf{x}}
\newcommand{\qpl}{\mathbf{l}}
\newcommand{\qpu}{\mathbf{u}}
\newcommand{\qpc}{\mathbf{c}}
\newcommand{\qpb}{\mathbf{b}}
\newcommand{\qpy}{\mathbf{y}}
\newcommand{\qpw}{\mathbf{w}}
\newcommand{\qplambda}{\mathbf{\lambda}}
\f$


\f{eqnarray*}{
\mbox{(QP)}& \mbox{minimize} 
& \qpx^{T}D\qpx+\qpc^{T}\qpx+c_0 \\ 
&\mbox{subject to} & A\qpx\qprel \qpb, \\ 
& & \qpl \leq \qpx \leq \qpu 
\f}
in \f$ n\f$ real variables \f$ \qpx=(x_0,\ldots,x_{n-1})\f$. 

If \f$ D=0\f$, the program is 
a linear program; if the variable bounds are \f$ \qpx\geq 0\f$, we have a 
nonnegative program. Objects of type `Quadratic_program_solution` are returned by any of 
the four functions 
`solve_quadratic_program`, `solve_linear_program`, 
`solve_nonnegative_quadratic_program`, and 
`solve_nonnegative_linear_program`. 

Example 
-------------- 

\ref QP_solver/first_qp.cpp 

Terminology 
-------------- 

If there is no \f$ \qpx\f$ that satisfies all the (in)equalities, 
the program is called <I>infeasible</I>, otherwise, it is <I>feasible</I>, 
and any \f$ \qpx\f$ that satisfies all (in)equalities is called a <I>feasible 
solution</I>. 

If the objective function value becomes arbitrarily small over the 
<I>feasible region</I> (the set of feasible solutions), the program 
is called <I>unbounded</I>, and <I>bounded</I> otherwise. 

Any program that is both feasible and bounded has at least one 
feasible solution \f$ \qpx^*\f$ whose objective function value is not larger 
than that of any other feasible solution. This is called an 
<I>optimal solution</I>. 

Every convex quadratic program (even if it is infeasible or unbounded) 
has a 'solution' in form of an object of the class `Quadratic_program_solution`. 

The program concepts 
\sa `QuadraticProgram`
\sa `LinearProgram`
\sa `NonnegativeQuadraticProgram`
\sa `NonnegativeLinearProgram`

and the functions that compute objects of class
`Quadratic_program_solution` from models of these concepts:

\sa `solve_quadratic_program`
\sa `solve_linear_program`
\sa `solve_nonnegative_quadratic_program`
\sa `solve_nonnegative_linear_program`

*/
template< typename ET >
class Quadratic_program_solution {
public:

/// \name Types 
/// @{

/*!
The exact number type that was used to solve the 
program. 
*/ 
typedef unspecified_type ET; 

/*!
An iterator type with value type 
`Quotient<ET>` to go over the values of the variables in the solution. 
*/ 
typedef unspecified_type Variable_value_iterator; 

/*!
An iterator type with value type 
`ET` to go over the numerators of the variable values 
with respect to a common denominator. 
*/ 
typedef unspecified_type Variable_numerator_iterator; 

/*!
An iterator type with value type `int` 
to go over the indices of the basic variables and the basic constraints. 
*/ 
typedef unspecified_type Index_iterator; 

/*!
An iterator type with 
value type `Quotient<ET>` to go over an \f$ m\f$-vector \f$ \qplambda\f$ that proves 
optimality of the solution. 
*/ 
typedef unspecified_type Optimality_certificate_iterator; 

/*!
An iterator type 
with value type `ET` to go over the numerators of the vector \f$ \qplambda\f$ 
with respect to a common denominator. 
*/ 
typedef unspecified_type Optimality_certificate_numerator_iterator; 

/*!
An iterator type with 
value type `ET` to go over an \f$ m\f$-vector \f$ \qplambda\f$ that proves 
infeasibility of the solution. 
*/ 
typedef unspecified_type Infeasibility_certificate_iterator; 

/*!
An iterator type with 
value type `ET` to go over an \f$ n\f$-vector \f$ \qpw\f$ that proves unboundedness 
of the solution. 
*/ 
typedef unspecified_type Unboundedness_certificate_iterator; 

/// @} 

/// \name Creation 
/// Objects of type `Quadratic_program_solution` can be copied and
/// assigned. Objects of type `Quadratic_program_solution` that are
/// associated to an actual program are returned by any of the four
/// functions `solve_quadratic_program`, `solve_linear_program`,
/// `solve_nonnegative_quadratic_program`, and
/// `solve_nonnegative_linear_program`. 
/// <b>Example:</b> \ref QP_solver/first_qp.cpp
/// @{

/*!
constructs a void instance of `Quadratic_program_solution` that is associated to 
no program. 
*/ 
Quadratic_program_solution(); 

/// @} 

/// \name Operations 
/// @{

/*!
returns `true` iff `sol` is not 
associated to a program. The condition 
!`sol``.is_void()` is a precondition for all access methods below. 
*/ 
bool is_void() const; 

/// @} 

/// \name Solution status 
/// Here are the access methods for the status of the solution.
/// @{

/*!
returns `true` iff `sol` is an optimal solution of the associated program. 
*/ 
bool is_optimal() const; 

/*!
returns `true` iff the 
associated program is infeasible. 
*/ 
bool is_infeasible() const; 

/*!
returns `true` iff the 
associated program is unbounded. 
*/ 
bool is_unbounded() const; 

/*!
returns the status of the solution; 
this is one of the values `QP_OPTIMAL`, `QP_INFEASIBLE`, and 
`QP_UNBOUNDED`, depending on whether the program asociated 
to `sol` has an optimal solution, is infeasible, or is unbounded. 
*/ 
Quadratic_program_status status() const; 

/*!
returns the number of iterations that it took to solve the 
program asociated to `sol`. 
*/ 
int number_of_iterations() const; 

/// @} 

/// \name Solution values 
/// The actual solution (variable and objective function values) can be accessed as follows.
/// @{

/*!
returns the objective 
function value of `sol`. 
*/ 
Quotient<ET> objective_value() const; 

/*!
returns the numerator 
of the objective function value of `sol`. 
*/ 
ET objective_value_numerator() const; 

/*!
returns the denominator 
of the objective function value of `sol`. 
*/ 
ET objective_value_denominator() const; 

/*!
returns 
a random-access iterator over the values of the 
variables in `sol`. The value type is `Quotient<ET>`, and the valid 
iterator range has length \f$ n\f$. 
*/ 
Variable_value_iterator variable_values_begin() const; 

/*!
returns 
the corresponding past-the-end iterator. 
*/ 
Variable_value_iterator variable_values_end() const; 

/*!
returns a random-access iterator `it` over the values 
of the variables in `sol`, with respect to a common 
denominator of all variables. The value type is `ET`, and the valid 
iterator range has length \f$ n\f$. 
*/ 
Variable_numerator_iterator variable_numerators_begin() const; 

/*!
returns the corresponding past-the-end iterator. 
*/ 
Variable_numerator_iterator variable_numerators_end() const; 

/*!
returns the common denominator of the variable values as referred to 
by the previous two methods. 
*/ 
const ET& variables_common_denominator() const; 

/// @} 

/// \name Basic variables and constraints 
/// The solution of a linear or quadratic program distinguishes
/// 'important' variables (the ones not attaining one of their
/// bounds), and 'important' constraints (the ones being satisfied
/// with equality). The following methods grant access to them. \ref
/// QP_solver/important_variables.cpp \ref QP_solver/first_qp_basic_constraints.cpp
/// @{

/*!
returns a random access iterator over the indices of the basic 
variables. The value type is `int`. It is guaranteed that any 
variable that is not basic in `sol` attains one of its bounds. 
In particular, if the bounds are of type \f$ \qpx\geq0\f$, all non-basic 
variables have value \f$ 0\f$. 
*/ 
Index_iterator basic_variable_indices_begin() const; 

/*!
returns the corresponding past-the-end iterator. 
*/ 
Index_iterator basic_variable_indices_end() const; 

/*!
returns the number of basic variables, equivalently the length 
of the range determined by the previous two iterators. 
*/ 
int number_of_basic_variables() const; 

/*!
returns a random access iterator over the indices of the basic 
constraints in the system \f$ A\qpx\qprel\qpb\f$. The value type is `int`. 
It is guaranteed that any basic constraint is satisfied with equality. 
In particular, if the system is of type \f$ A\qpx=\qpb\f$, all constraints are 
basic. 
*/ 
Index_iterator basic_constraint_indices_begin() const; 

/*!
returns the corresponding past-the-end iterator. 
*/ 
Index_iterator basic_constraint_indices_end() const; 

/*!
returns the number of basic constraint, equivalently the length 
of the range determined by the previous two iterators. 
*/ 
int number_of_basic_constraints() const; 

/*!
writes the status of `sol` to the stream `out`. In case the 
status is `QP_OPTIMAL`, the optimal objective value and the values 
of the variables at the optimal solution are output as well. For more 
detailed information about the solution (like basic variables/constraints) 
please use the dedicated methods of `Quadratic_program_solution<ET>`. 
\relates Quadratic_program_solution 
*/ 
template <typename ET> 
std::ostream& operator<<(std::ostream& out, const Quadratic_program_solution<ET>& sol); 

/// @} 

/// \name Validity 
/// The following four methods allow you to check whether `sol` indeed
/// solves the program that you intended to solve. The methods use the
/// certificates described in the advanced section below and thus save
/// you from validating the certificates yourself (if you believe in
/// the correctness of these methods; otherwise, you can look at their
/// implementation to convince yourself). By passing a suitable option
/// to the solution function, you can make sure that this check is
/// done automatically after the solution of the program, see
/// `Quadratic_program_options`. If the check fails, a logfile is
/// generated that contains the details, and an error message is
/// written to `std::cerr` (see \ref QP_solver/cycling.cpp for an example that
/// uses this option). \ref QP_solver/first_qp.cpp \ref QP_solver/first_lp.cpp \ref
/// QP_solver/first_nonnegative_qp.cpp \ref QP_solver/first_nonnegative_lp.cpp
/// @{

/*!
returns `true` iff `sol` solves the quadratic program `qp`. 
If the result is `false`, you can get a message that describes the 
problem, through the method `get_error()`. 
*/ 
template <class QuadraticProgram> 
bool solves_quadratic_program 
(const QuadraticProgram& qp); 

/*!
returns `true` iff `sol` solves the linear program `lp`. 
If the result is `false`, you can get a message that describes the 
problem, through the method `get_error()`. 
*/ 
template <class LinearProgram> 
bool solves_linear_program 
(const LinearProgram& lp); 

/*!
returns `true` iff `sol` solves the nonnegative 
quadratic program `qp`. 
If the result is `false`, you can get a message that describes the 
problem, through the method `get_error()`. 
*/ 
template <class NonnegativeQuadraticProgram> 
bool solves_nonnegative_quadratic_program 
(const NonnegativeQuadraticProgram& qp); 

/*!
returns `true` iff `sol` solves the nonnegative 
linear program `lp`. If the result is `false`, you can get a message that describes the 
problem, through the method `get_error()`. 
*/ 
template <class NonnegativeLinearProgram> 
bool solves_nonnegative_linear_program 
(const NonnegativeLinearProgram& lp); 

/*!
returns `false` iff the validation through one of the 
previous four functions has failed. 
*/ 
bool is_valid() const; 

/*!
returns an error message in case any of the 
previous four validation functions has returned `false`. 
*/ 
const std::string& get_error() const; 

/// @} 

/*!
\name Certificates 
A certificate is a vector that admits a simple proof for the
correctness of the solution. Any non-void object of
`Quadratic_program_solution` comes with such a certificate. 
*/
/// @{

/*!
returns a random access iterator over the optimality certificate 
\f$ \qplambda\f$ as given in Lemma 1, with respect to the solution \f$ \qpx^*\f$ 
obtained from `sol``.variable_values_begin()`. The value type 
is `Quotient<ET>`, and the valid iterator range has length \f$ m\f$. 
\pre `sol``.is_optimal()`. 

<B>Lemma 1(optimality certificate):</B> A feasible \f$ n\f$-vector \f$\qpx^*\f$ is an
optimal solution of (QP) if an \f$ m\f$-vector \f$ \qplambda\f$ with
the following properties exist. 
<OL>
 <LI>if the \f$ i\f$-th constraint
 is of type \f$ \leq\f$ (\f$ \geq\f$, respectively), then 
 \f$\lambda_i\geq 0\f$ (\f$\lambda_i\leq 0\f$, respectively). 
<LI>\f$\qplambda^T(A\qpx^*-\qpb) = 0\f$
<LI>\f[
  \begin{array}{llll}
  &&\geq 0, & \mbox{if $x^*_j = l_j < u_j$} \\
  (\qpc^T + \qplambda^T A + 2{\qpx^*}^TD)_j& \quad  &= 0, & \mbox{if $l_j < x^*_j < u_j$} \\
  &&\leq 0, & \mbox{if $l_j < u_j = x^*_j$.}
  \end{array}
  \f]
</OL>

<B>Proof:</B> Let \f$\qpx\f$ be any feasible solution. We need to prove that
\f[\qpc^T\qpx + \qpx^TD\qpx \geq \qpc^T\qpx^* + {\qpx^*}^TD\qpx^*.\f]

For this, we argue as follows.
\f[
\begin{array}{lcll}
\qpc^T\qpx + 2{\qpx^*}^TD\qpx &\geq& \qpc^T\qpx + 2{\qpx^*}^TD\qpx + \qplambda^T(A\qpx-\qpb) &  
\mbox{(by $A\qpx\qprel \qpb$ and 1.)} \\
                  &=& (\qpc^T + \qplambda^T A + 2{\qpx^*}^TD)\qpx - \qplambda^Tb \\
                  &\geq& (\qpc^T + \qplambda^T A + 2{\qpx^*}^TD)\qpx^* - \qplambda^Tb &
\mbox{(by $\qpl\leq \qpx \leq \qpu$ and 3.)} \\
                  &=& \qpc^T\qpx^* + 2{\qpx^*}^TD\qpx^* &
\mbox{(by 2.)}
\end{array}
\f]

After adding \f$\qpx^TD\qpx - \qpx^TD\qpx - {\qpx^*}^TD\qpx^* = -{\qpx^*}^TD\qpx^*\f$ to both sides of
this inequality, we get
\f[
\qpc^T\qpx + \qpx^TD\qpx - (\qpx-\qpx^*)^TD(\qpx-\qpx^*) \geq \qpc^T\qpx^* + {\qpx^*}^TD\qpx^*,
\f]
and since \f$D\f$ is positive semidefinite, we have
\f$(\qpx-\qpx^*)^TD(\qpx-\qpx^*)\geq 0\f$ and the lemma follows.

\sa \ref QP_solver/optimality_certificate.cpp
*/ 
Optimality_certificate_iterator 
optimality_certifcate_begin() const; 

/*!
returns the corresponding past-the-end iterator. 

\sa `optimality_certifcate_begin()`
*/ 
Optimality_certificate_iterator 
optimality_certificate_end() const; 

/*!
returns a random access iterator over the numerator values 
of the optimality certificate \f$ \qplambda\f$, with respect to the 
common denominator returned by `sol`.`variables_common_denominator()`. 
The value type is `ET`, and the valid iterator range has length \f$ m\f$. 

\sa `optimality_certifcate_begin()`
*/ 
Optimality_certificate_numerator_iterator 
optimality_certifcate_numerators_begin() const; 

/*!
returns the corresponding past-the-end iterator. 

\sa `optimality_certifcate_begin()`
*/ 
Optimality_certificate_numerator_iterator 
optimality_certificate_numerators_end() const; 

/*!
returns a random access iterator over the infeasibility certificate 
\f$ \qplambda\f$ as given in Lemma 2. The value type 
is `ET`, and the valid iterator range has length \f$ m\f$. 
\pre `sol``.is_infeasible()`. 


<B>Lemma 2 (infeasibility certificate):</B> The program (QP) is
infeasible if an \f$m\f$-vector \f$\qplambda\f$ with the
following properties exist.

<OL>
<LI> if the \f$i\f$-th constraint is of type \f$\leq\f$ (\f$\geq\f$, respectively), 
then \f$\lambda_i\geq 0\f$ (\f$\lambda_i\leq 0\f$, respectively).
<LI>
\f[
\begin{array}{llll}
&&\geq 0 & \mbox{if \f$u_j=\infty\f$} \\
\qplambda^T A_j &\quad  \\
&&\leq 0 & \mbox{if \f$l_j=-\infty\f$.}
\end{array}
\f]
<LI> \f[\qplambda^T\qpb \quad<\quad \ccSum{j: \qplambda^TA_j <0}{}{ \qplambda^TA_j u_j }
\quad+\quad  \ccSum{j: \qplambda^TA_j >0}{}{ \qplambda^TA_j l_j}.\f]
</OL>

<B>Proof:</B> Let us assume for the purpose of obtaining a contradiction
that there is a feasible solution \f$\qpx\f$. Then we get
\f[
\begin{array}{lcll}
0 &\geq& \qplambda^T(A\qpx -\qpb) &  \mbox{(by \f$A\qpx\qprel \qpb\f$ and 1.)} \\
  &=& \ccSum{j: \qplambda^TA_j <0}{}{ \qplambda^TA_j x_j }
\quad+\quad  \ccSum{j: \qplambda^TA_j >0}{}{ \qplambda^TA_j x_j} - \qplambda^T \qpb \\
  &\geq& \ccSum{j: \qplambda^TA_j <0}{}{ \qplambda^TA_j u_j }
\quad+\quad  \ccSum{j: \qplambda^TA_j >0}{}{ \qplambda^TA_j l_j} - \qplambda^T \qpb &
\mbox{(by \f$\qpl\leq \qpx \leq \qpu\f$ and 2.)} \\
  &>& 0 & \mbox{(by 3.)}, 
\end{array}
\f]
and this is the desired contradiction \f$0>0\f$.

\sa \ref QP_solver/infeasibility_certificate.cpp
*/ 
Infeasibility_certificate_iterator 
infeasibility_certificate_begin() const; 

/*!
returns the corresponding past-the-end iterator. 

\sa `infeasibility_certificate_begin()`
*/ 
Infeasibility_certificate_iterator 
infeasibility_certificate_end() const; 

/*!
returns a random acess iterator over the unbounded direction \f$ \qpw\f$ 
as given in Lemma 3,with respect to the solution \f$ \qpx^*\f$ 
obtained from `sol``.variable_values_begin()`. The value type 
is `ET`, and the valid iterator range has length \f$ n\f$. 
\pre `sol``.is_unbounded()`. 

<B>Lemma 3 (unboundedness certificate:)</B> Let \f$\qpx^*\f$ be a feasible
solution of (QP). The program (QP) is unbounded if an \f$n\f$-vector 
\f$\qpw\f$ with the following properties exist.
<OL>
<LI> if the \f$i\f$-th constraint is of type \f$\leq\f$ (\f$\geq, =\f$, respectively),
then \f$(A\qpw)_i\leq 0\f$ (\f$(A\qpw)_i\geq 0, (A\qpw)_i=0\f$, respectively).
<LI> \f[
\begin{array}{llll}
&&\geq 0 & \mbox{if \f$l_j\f$ is finite} \\
w_j &\quad  \\
&&\leq 0 & \mbox{if \f$u_j\f$ is finite.}
\end{array}
\f]
<LI> \f$\qpw^TD\qpw=0\f$ and \f$(\qpc^T+2{\qpx^*}^TD)\qpw<0\f$.
</OL>

The vector \f$\qpw\f$ is called an <I>unbounded direction</I>.

<B>Proof:</B> For a real number \f$t\f$, consider the vector \f$\qpx(t):=\qpx^*+t\qpw\f$. By 1.
and 2., \f$\qpx(t)\f$ is feasible for all \f$t\geq 0\f$. The objective function value
of \f$\qpx(t)\f$ is
\f{eqnarray*}{
\qpc^T \qpx(t) + \qpx(t)^TD \qpx(t) &=& 
\qpc^T\qpx^* + t\qpc^T\qpw + {\qpx^*}^TD\qpx^* +  2t{\qpx^*}^TD\qpw + t^2 \qpw^TD\qpw  \\
&=& \qpc^T\qpx^* + {\qpx^*}^TD\qpx^* + t(\qpc^T + 2{\qpx^*}^TD)\qpw + t^2\qpw^TD\qpw.
\f}
By condition 3., this tends to \f$-\infty\f$ for \f$t\rightarrow\infty\f$, so
the problem is indeed unbounded. 
*/ 
Unboundedness_certificate_iterator 
unboundedness_certificate_begin() const; 

/*!
returns the corresponding past-the-end iterator. 

\sa `unboundedness_certificate_begin()`
*/ 
Unboundedness_certificate_iterator 
unboundedness_certificate_end(); 

/// @}

}; /* end Quadratic_program_solution */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgQPSolverClasses

This is an enumeration type containing the values 
`QP_OPTIMAL`, `QP_INFEASIBLE`, and `QP_UNBOUNDED`. 
It indicates the status of a linear or quadratic program solution 
as represented by an object of type 
`Quadratic_program_solution<ET>`. 

\sa `Quadratic_program_solution<ET>` 

*/
enum Quadratic_program_status {
  QP_OPTIMAL,
  QP_INFEASIBLE,
  QP_UNBOUNDED;
}; /* end Quadratic_program_status */
} /* end namespace CGAL */
