namespace CGAL {

/// \addtogroup PkgQPSolverFunctions
/// @{

/*!
\tparam LinearProgram a model of `LinearProgram`.

This function writes a linear program to an output stream (in
`MPSFormat`). The time
complexity is \f$ \Theta (mn)\f$, even if the program is very sparse.

It writes the linear program `lp` to `out` in
`MPSFormat`. The name of the program will be the one provided
by `problem_name`.

\cgalHeading{Requirements}
Output operators are defined for all entry types of `lp`.

*/
template <LinearProgram>
void print_linear_program
(std::ostream& out, const LinearProgram &lp,
const std::string& problem_name = std::string("MY_MPS"));


/*!
\tparam NonnegativeLinearProgram a model of `NonnegativeLinearProgram`

This function writes a nonnegative linear program to an output stream
(in `MPSFormat`). The time
complexity is \f$ \Theta (mn)\f$, even if the program is very sparse.

Writes the nonnegative linear program `lp` to `out` in
`MPSFormat`. The name of the program will be the one provided
by `problem_name`.

\cgalHeading{Requirements}
Output operators are defined for all entry types of `lp`.

*/
template <NonnegativeLinearProgram>
void print_nonnegative_linear_program
(std::ostream& out, const NonnegativeLinearProgram &lp,
const std::string& problem_name = std::string("MY_MPS"));


/*!
\tparam NonnegativeQuadraticProgram a model of `NonnegativeQuadraticProgram`

This function writes a nonnegative quadratic program
to an output stream (in `MPSFormat`). The time
complexity is \f$ \Theta (n^2 + mn)\f$, even if the program is very sparse.

Writes the nonnegative quadratic program `qp` to `out` in
`MPSFormat`. The name of the program will be the one provided
by `problem_name`.

\cgalHeading{Requirements}
Output operators are defined for all entry types of `qp`.

*/
template <NonnegativeQuadraticProgram>
void print_nonnegative_quadratic_program
(std::ostream& out, const NonnegativeQuadraticProgram &qp,
const std::string& problem_name = std::string("MY_MPS"));


/*!
\tparam QuadraticProgram a model of `QuadraticProgram`


This function writes a quadratic program to an output stream (in
`MPSFormat`). The time complexity is \f$ \Theta (n^2 + mn)\f$, even
if the program is very sparse.

Writes the quadratic program `qp` to `out` in `MPSFormat`.
The name of the program will be the one provided by `problem_name`.

\cgalHeading{Requirements}
Output operators are defined for all entry types of `qp`.

*/
template <QuadraticProgram>
void print_quadratic_program
(std::ostream& out, const QuadraticProgram &qp,
const std::string& problem_name = std::string("MY_MPS"));


/*!
This function solves a linear program, using some exact
Integral Domain `ET` for its computations. Various
options may be provided, see `Quadratic_program_options`.

\tparam LinearProgram a model of `LinearProgram` such as `Quadratic_program<NT>`,
`Quadratic_program_from_mps<NT>`, or `Linear_program_from_iterators<A_it, B_it, R_it, FL_it, L_it, FU_it, U_it, C_it>`

\tparam ET a model of the concepts `IntegralDomain` and
`RealEmbeddable`; it must
be an exact type, and all entries of `qp` are convertible to
`ET`.

Here are some recommended combinations of input type (the type of
the `qp` entries) and `ET`.

input type          |  ET
----------          | -------------
`double`            | `MP_Float`, `Gmpzf`, or `Gmpq`
`int`               | `MP_Float`, or `Gmpz`
any exact type `NT` |  `NT`

\note By default, this function performs a large number of
runtime-checks to ensure consistency during the solution process.
However, these checks slow down the computations by a considerable
factor. For maximum efficiency, it is advisable to define the macros
<TT>CGAL_QP_NO_ASSERTIONS</TT> or <TT>NDEBUG</TT>.

\returns the solution of the linear program `lp`, solved
with exact number type `ET`.

*/
template <LinearProgram, ET>
Quadratic_program_solution<ET> solve_linear_program
(const LinearProgram& lp, const ET&,
const Quadratic_program_options& options = Quadratic_program_options());


/*!
This function solves a nonnegative linear program, using some exact
Integral Domain `ET` for its computations. Various
options may be provided, see `Quadratic_program_options`.

\tparam NonnegativeQuadraticProgram a model of `NonnegativeQuadraticProgram` such as
`Quadratic_program<NT>`, `Quadratic_program_from_mps<NT>`, or `Nonnegative_quadratic_program_from_iterators<A_it, B_it, R_it, D_it, C_it>`.

\tparam ET is a model of the concepts `IntegralDomain` and
`RealEmbeddable`; it must
be an exact type, and all entries of `qp` are convertible to
`ET`.

Here are some recommended combinations of input type (the type of
the `qp` entries) and `ET`.

input type          |  ET
----------          | -------------
`double`            | `MP_Float`, `Gmpzf`, or `Gmpq`
`int`               | `MP_Float`, or `Gmpz`
any exact type `NT` |  `NT`

\note By default, this function performs a large number of
runtime-checks to ensure consistency during the solution process.
However, these checks slow down the computations by a considerable
factor. For maximum efficiency, it is advisable to define the macros
<TT>CGAL_QP_NO_ASSERTIONS</TT> or <TT>NDEBUG</TT>.

\returns the solution of the nonnegative linear program `lp`, solved
with exact number type `ET`.

*/
template <NonnegativeLinearProgram, ET>
Quadratic_program_solution<ET> solve_nonnegative_linear_program
(const NonnegativeLinearProgram& lp, const ET&,
const Quadratic_program_options& options = Quadratic_program_options());


/*!
This function solves a nonnegative quadratic program, using some exact
Integral Domain `ET` for its computations. Various
options may be provided, see `Quadratic_program_options`.

\tparam NonnegativeQuadraticProgram a model of `NonnegativeQuadraticProgram` such as
`Quadratic_program<NT>`, `Quadratic_program_from_mps<NT>`, or `Nonnegative_quadratic_program_from_iterators<A_it, B_it, R_it, D_it, C_it>`.

\tparam ET a model of the concepts `IntegralDomain` and
`RealEmbeddable`; it must
be an exact type, and all entries of `qp` are convertible to
`ET`.

Here are some recommended combinations of input type (the type of
the `qp` entries) and `ET`.

input type          |  ET
----------          | -------------
`double`            | `MP_Float`, `Gmpzf`, or `Gmpq`
`int`               | `MP_Float`, or `Gmpz`
any exact type `NT` |  `NT`

\note By default, this function performs a large number of
runtime-checks to ensure consistency during the solution process.
However, these checks slow down the computations by a considerable
factor. For maximum efficiency, it is advisable to define the macros
<TT>CGAL_QP_NO_ASSERTIONS</TT> or <TT>NDEBUG</TT>.

\returns the solution of the nonnegative quadratic program `qp`, solved
with exact number type `ET`.

*/
template <NonnegativeQuadraticProgram, ET>
Quadratic_program_solution<ET> solve_nonnegative_quadratic_program
(const NonnegativeQuadraticProgram& qp, const ET&,
const Quadratic_program_options& options = Quadratic_program_options());


/*!
This function solves a quadratic program, using some exact
Integral Domain `ET` for its computations. Various
options may be provided, see `Quadratic_program_options`.


\tparam NonnegativeQuadraticProgram a model of `NonnegativeQuadraticProgram` such as
`Quadratic_program<NT>`, `Quadratic_program_from_mps<NT>`, or `Nonnegative_quadratic_program_from_iterators<A_it, B_it, R_it, D_it, C_it>`.

\tparam ET is a model of the concepts `IntegralDomain` and
`RealEmbeddable`; it must
be an exact type, and all entries of `qp` are convertible to
`ET`.

Here are some recommended combinations of input type (the type of
the `qp` entries) and `ET`.

input type          |  ET
----------          | -------------
`double`            | `MP_Float`, `Gmpzf`, or `Gmpq`
`int`               | `MP_Float`, or `Gmpz`
any exact type `NT` |  `NT`

\note By default, this function performs a large number of
runtime-checks to ensure consistency during the solution process.
However, these checks slow down the computations by a considerable
factor. For maximum efficiency, it is advisable to define the macros
<TT>CGAL_QP_NO_ASSERTIONS</TT> or <TT>NDEBUG</TT>.

\returns the solution of the quadratic program `qp`, solved
with exact number type `ET`.

*/
template <QuadraticProgram, ET>
Quadratic_program_solution<ET> solve_quadratic_program
(const QuadraticProgram& qp, const ET&,
const Quadratic_program_options& options = Quadratic_program_options());

/// @}

} /* namespace CGAL */

