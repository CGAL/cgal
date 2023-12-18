
namespace CGAL {

/*!
  \ingroup PkgQPSolverClasses

  An object of class `Linear_program_from_iterators` describes a linear program of the form

  \f{eqnarray*}{
  \mbox{(QP)}& \mbox{minimize}
  & \qpc^{T}\qpx+c_0 \\
  &\mbox{subject to} & A\qpx\qprel \qpb, \\
  & & \qpl \leq \qpx \leq \qpu
  \f}
  in \f$ n\f$ real variables \f$ \qpx=(x_0,\ldots,x_{n-1})\f$.

  Here,
  <UL>
  <LI>\f$ A\f$ is an \f$ m\times n\f$ matrix (the constraint matrix),
  <LI>\f$ \qpb\f$ is an \f$ m\f$-dimensional vector (the right-hand side),
  <LI>\f$ \qprel\f$ is an \f$ m\f$-dimensional vector of relations
  from \f$ \{\leq, =, \geq\}\f$,

  <LI>\f$ \qpl\f$ is an \f$ n\f$-dimensional vector of lower
  bounds for \f$ \qpx\f$, where \f$ l_j\in\mathbb{R}\cup\{-\infty\}\f$ for all \f$ j\f$
  <LI>\f$ \qpu\f$ is an \f$ n\f$-dimensional vector of upper bounds for
  \f$ \qpx\f$, where \f$ u_j\in\mathbb{R}\cup\{\infty\}\f$ for all \f$ j\f$

  <LI>\f$ \qpc\f$ is an \f$ n\f$-dimensional vector (the linear objective
  function), and
  <LI>\f$ c_0\f$ is a constant.

  </UL>

  This class is simply a wrapper for existing iterators, and it does not
  copy the program data.

  It frequently happens that all values in one of the vectors from
  above are the same, for example if the system \f$ Ax\qprel b\f$ is
  actually a system of equations \f$ Ax=b\f$. To get an iterator over such a
  vector, it is not necessary to store multiple copies of the value in
  some container; an instance of the class `Const_oneset_iterator<T>`,
  constructed from the value in question, does the job more efficiently.

  \cgalModels{QuadraticProgram,LinearProgram}

  \cgalHeading{Example}

  \ref QP_solver/first_lp_from_iterators.cpp

  The following example for the simpler model
  `Nonnegative_linear_program_from_iterators<A_it, B_it, R_it, C_it>`
  should give you a flavor of the use of this
  model in practice.

  \ref QP_solver/solve_convex_hull_containment_lp.h

  \ref QP_solver/convex_hull_containment.cpp

  \sa `LinearProgram`
  \sa `Quadratic_program<NT>`
  \sa `Quadratic_program_from_mps<NT>`

*/
template< typename A_it, typename B_it, typename R_it, typename FL_it, typename L_it, typename FU_it, typename U_it, typename C_it >
class Linear_program_from_iterators {
public:

  /// \name Creation
  /// @{

  /*!
    constructs `lp` from given random-access iterators and the constant `c0`. The passed iterators are merely stored, no copying of the program data takes place. How these iterators are supposed to encode the linear program is
    described in `LinearProgram`.
  */
  Linear_program_from_iterators(int n, int m,
                                const A_it& a,
                                const B_it& b,
                                const R_it& r,
                                const FL_it& fl,
                                const L_it& l,
                                const FU_it& fu,
                                const U_it& u,
                                const C_it& c,
                                const std::iterator_traits<C_it>value_type& c0 = 0
    );

  /// @}

}; /* end Linear_program_from_iterators */

/*!
  \ingroup PkgQPSolverFunctions

  This template function creates an instance of
  `Linear_program_from_iterators<A_it, B_it, R_it, FL_it, L_it, FU_it, U_it, C_it>` from given iterators. This function can be useful if the types of these
  iterators are too complicated (or of too little interest for you)
  to write them down explicitly.

  \returns an instance of `Linear_program_from_iterators<A_it, B_it, R_it, FL_it, L_it, FU_it, U_it, C_it>`, constructed from the given iterators.

  \cgalHeading{Example}

  The following example demonstrates the typical usage of makers
  with the simpler function `make_nonnegative_linear_program_from_iterators()`.

  `QP_solver/solve_convex_hull_containment_lp2.h`

  \ref QP_solver/convex_hull_containment2.cpp

  \sa `Linear_program_from_iterators<A_it, B_it, R_it, FL_it, L_it, FU_it, U_it, C_it>`
*/
template <
  typename A_it,
  typename B_it,
  typename R_it,
  typename FL_it,
  typename L_it,
  typename FU_it,
  typename U_it,
  typename C_it >
Linear_program_from_iterators<A_it, B_it, R_it, FL_it, L_it, FU_it, U_it, C_it>
make_linear_program_from_iterators (
  int n, int m,
  const A_it& a,
  const B_it& b,
  const R_it& r,
  const FL_it& fl,
  const L_it& l,
  const FU_it& fu,
  const U_it& u,
  const C_it& c,
  std::iterator_traits<C_it>::value_type c0 = std::iterator_traits<C_it>::value_type(0));


/*!
  \ingroup PkgQPSolverFunctions

  This template function creates an instance of
  `Nonnegative_linear_program_from_iterators<A_it, B_it, R_it, C_it>`
  from given iterators. This function can be useful if the types of these
  iterators are too complicated (or of too little interest for you)
  to write them down explicitly.

  \returns an instance of `Nonnegative_linear_program_from_iterators<A_it, B_it, R_it, C_it>`, constructed from the given iterators.


  \cgalHeading{Example}

  `QP_solver/solve_convex_hull_containment_lp2.h`

  \ref QP_solver/convex_hull_containment2.cpp

  \sa `Nonnegative_linear_program_from_iterators<A_it, B_it, R_it, C_it>`
*/
template <
  A_it,
  B_it,
  R_it,
  C_it >
Nonnegative_linear_program_from_iterators
<A_it, B_it, R_it, C_it>
make_nonnegative_linear_program_from_iterators (
  int n, int m,
  const A_it& a,
  const B_it& b,
  const R_it& r,
  const C_it& c,
  std::iterator_traits<C_it>::value_type c0 =
  std::iterator_traits<C_it>::value_type(0));


/*!
  \ingroup PkgQPSolverFunctions

  This template function creates an instance of
  `Nonnegative_quadratic_program_from_iterators<A_it, B_it, R_it, D_it, C_it>` from given iterators. This function can be useful if the types of these
  iterators are too complicated (or of too little interest for you)
  to write them down explicitly.

  \returns an instance of
  `Nonnegative_quadratic_program_from_iterators<A_it, B_it, R_it, D_it,C_it>`,
  constructed from the given iterators.

  \cgalHeading{Example}

  The following example demonstrates the typical usage of makers
  with the simpler function `make_nonnegative_linear_program_from_iterators()`.

  `QP_solver/solve_convex_hull_containment_lp2.h`

  \ref QP_solver/convex_hull_containment2.cpp

  \sa `Nonnegative_quadratic_program_from_iterators<A_it, B_it, R_it, D_it, C_it>`


*/
template <
  A_it,
  B_it,
  R_it,
  D_it,
  C_it >
Nonnegative_quadratic_program_from_iterators
<A_it, B_it, R_it, D_it, C_it>
make_nonnegative_quadratic_program_from_iterators (
  int n, int m,
  const A_it& a,
  const B_it& b,
  const R_it& r,
  const D_it& d,
  const C_it& c,
  std::iterator_traits<C_it>::value_type c0 =
  std::iterator_traits<C_it>::value_type(0));

/*!
  \ingroup PkgQPSolverFunctions

  This template function creates an instance of
  `Quadratic_program_from_iterators<A_it, B_it, R_it, FL_it, L_it, FU_it, U_it, D_it, C_it>` from given iterators. This function can be useful if the types of these
  iterators are too complicated (or of too little interest for you)
  to write them down explicitly.

  \returns an instance of `Quadratic_program_from_iterators<A_it, B_it, R_it, FL_it, L_it, FU_it, U_it, D_it, C_it>`, constructed from the given iterators.

  \cgalHeading{Example}

  The following example demonstrates the typical usage of makers
  with the simpler function `make_nonnegative_linear_program_from_iterators()`.

  `QP_solver/solve_convex_hull_containment_lp2.h`

  \ref QP_solver/convex_hull_containment2.cpp

  \sa `Quadratic_program_from_iterators<A_it, B_it, R_it, FL_it, L_it, FU_it, U_it, D_it, C_it>`
*/
template <
  typename A_it,
  typename B_it,
  typename R_it,
  typename FL_it,
  typename L_it,
  typename FU_it,
  typename U_it,
  typename D_it,
  typename C_it >
Quadratic_program_from_iterators<A_it, B_it, R_it, FL_it, L_it, FU_it, U_it, D_it, C_it>
make_quadratic_program_from_iterators (
  int n, int m,
  const A_it& a,
  const B_it& b,
  const R_it& r,
  const FL_it& fl,
  const L_it& l,
  const FU_it& fu,
  const U_it& u,
  const D_it& d,
  const C_it& c,
  std::iterator_traits<C_it>::value_type c0 =
  std::iterator_traits<C_it>::value_type(0));

/*!
  \ingroup PkgQPSolverClasses

  An object of class `Nonnegative_linear_program_from_iterators` describes a linear program of the form

  \f{eqnarray*}{
  \mbox{(QP)}& \mbox{minimize}
  &\qpc^{T}\qpx+c_0 \\
  &\mbox{subject to} & A\qpx\qprel \qpb, \\
  & & \qpx \geq 0
  \f}
  in \f$ n\f$ real variables \f$ \qpx=(x_0,\ldots,x_{n-1})\f$.

  Here,
  <UL>
  <LI>\f$ A\f$ is an \f$ m\times n\f$ matrix (the constraint matrix),
  <LI>\f$ \qpb\f$ is an \f$ m\f$-dimensional vector (the right-hand side),
  <LI>\f$ \qprel\f$ is an \f$ m\f$-dimensional vector of relations
  from \f$ \{\leq, =, \geq\}\f$,

  <LI>\f$ \qpc\f$ is an \f$ n\f$-dimensional vector (the linear objective
  function), and
  <LI>\f$ c_0\f$ is a constant.

  </UL>

  This class is simply a wrapper for existing iterators, and it does not
  copy the program data.

  It frequently happens that all values in one of the vectors from
  above are the same, for example if the system \f$ Ax\qprel b\f$ is
  actually a system of equations \f$ Ax=b\f$. To get an iterator over such a
  vector, it is not necessary to store multiple copies of the value in
  some container; an instance of the class `Const_oneset_iterator<T>`,
  constructed from the value in question, does the job more efficiently.

  \cgalModels{QuadraticProgram,LinearProgram,NonnegativeQuadraticProgram,NonnegativeLinearProgram}

  \cgalHeading{Example}

  \ref QP_solver/first_nonnegative_lp_from_iterators.cpp

  \ref QP_solver/solve_convex_hull_containment_lp.h

  \ref QP_solver/convex_hull_containment.cpp

  \sa `NonnegativeLinearProgram`
  \sa `Quadratic_program<NT>`
  \sa `Quadratic_program_from_mps<NT>`

*/
template< typename A_it, typename B_it, typename R_it, typename C_it >
class Nonnegative_linear_program_from_iterators {
public:

  /// \name Creation
  /// @{

  /*!
    constructs `lp` from given random-access iterators and the constant
    `c0`. The passed iterators are merely stored, no copying of the program
    data takes place. How these iterators are supposed to encode the nonnegative
    linear program is described in `NonnegativeLinearProgram`.
  */
  Nonnegative_linear_program_from_iterators(int n, int m,
                                            const A_it& a,
                                            const B_it& b,
                                            const R_it& r,
                                            const C_it& c,
                                            const std::iterator_traits<C_it>value_type& c0 = 0
    );

  /// @}

}; /* end Nonnegative_linear_program_from_iterators */

/*!
  \ingroup PkgQPSolverClasses

  An object of class `Nonnegative_quadratic_program_from_iterators` describes a convex quadratic program of the form

  \f{eqnarray*}{
  \mbox{(QP)}& \mbox{minimize}
  & \qpx^{T}D\qpx+\qpc^{T}\qpx+c_0 \\
  &\mbox{subject to} & A\qpx\qprel \qpb, \\
  & & \qpx \geq 0
  \f}
  in \f$ n\f$ real variables \f$ \qpx=(x_0,\ldots,x_{n-1})\f$.

  Here,
  <UL>
  <LI>\f$ A\f$ is an \f$ m\times n\f$ matrix (the constraint matrix),
  <LI>\f$ \qpb\f$ is an \f$ m\f$-dimensional vector (the right-hand side),
  <LI>\f$ \qprel\f$ is an \f$ m\f$-dimensional vector of relations
  from \f$ \{\leq, =, \geq\}\f$,

  <LI>\f$ D\f$ is a symmetric positive-semidefinite \f$ n\times n\f$ matrix (the
  quadratic objective function),

  <LI>\f$ \qpc\f$ is an \f$ n\f$-dimensional vector (the linear objective
  function), and
  <LI>\f$ c_0\f$ is a constant.

  </UL>

  This class is simply a wrapper for existing iterators, and it does not
  copy the program data.

  It frequently happens that all values in one of the vectors from
  above are the same, for example if the system \f$ Ax\qprel b\f$ is
  actually a system of equations \f$ Ax=b\f$. To get an iterator over such a
  vector, it is not necessary to store multiple copies of the value in
  some container; an instance of the class `Const_oneset_iterator<T>`,
  constructed from the value in question, does the job more efficiently.

  \cgalModels{QuadraticProgram,NonnegativeQuadraticProgram}

  \cgalHeading{Example}

  \ref QP_solver/first_nonnegative_qp_from_iterators.cpp

  The following example for the simpler model
  `Nonnegative_linear_program_from_iterators<A_it, B_it, R_it, C_it>`
  should give you a flavor of the use of this
  model in practice.

  \ref QP_solver/solve_convex_hull_containment_lp.h

  \ref QP_solver/convex_hull_containment.cpp

  \sa `NonnegativeQuadraticProgram`
  \sa `Quadratic_program<NT>`
  \sa `Quadratic_program_from_mps<NT>`

*/
template< typename A_it, typename B_it, typename R_it, typename D_it, typename C_it >
class Nonnegative_quadratic_program_from_iterators {
public:

  /// \name Creation
  /// @{

  /*!
    constructs `qp` from given random-access iterators and the constant `c0`. The passed iterators are merely stored, no copying of the program data takes place. How these iterators are supposed to encode the nonnegative
    quadratic program is described in `NonnegativeQuadraticProgram`.
  */
  Nonnegative_quadratic_program_from_iterators(int n, int m,
                                               const A_it& a,
                                               const B_it& b,
                                               const R_it& r,
                                               const D_it& d,
                                               const C_it& c,
                                               const std::iterator_traits<C_it>value_type& c0 = 0
    );

  /// @}

}; /* end Nonnegative_quadratic_program_from_iterators */

/*!
  \ingroup PkgQPSolverClasses

  An object of class `Quadratic_program_from_iterators` describes a convex quadratic program of the form

  \f{eqnarray*}{
  \mbox{(QP)}& \mbox{minimize}
  & \qpx^{T}D\qpx+\qpc^{T}\qpx+c_0 \\
  &\mbox{subject to} & A\qpx\qprel \qpb, \\
  & & \qpl \leq \qpx \leq \qpu
  \f}
  in \f$ n\f$ real variables \f$ \qpx=(x_0,\ldots,x_{n-1})\f$.

  Here,
  <UL>
  <LI>\f$ A\f$ is an \f$ m\times n\f$ matrix (the constraint matrix),
  <LI>\f$ \qpb\f$ is an \f$ m\f$-dimensional vector (the right-hand side),
  <LI>\f$ \qprel\f$ is an \f$ m\f$-dimensional vector of relations
  from \f$ \{\leq, =, \geq\}\f$,

  <LI>\f$ \qpl\f$ is an \f$ n\f$-dimensional vector of lower
  bounds for \f$ \qpx\f$, where \f$ l_j\in\mathbb{R}\cup\{-\infty\}\f$ for all \f$ j\f$
  <LI>\f$ \qpu\f$ is an \f$ n\f$-dimensional vector of upper bounds for
  \f$ \qpx\f$, where \f$ u_j\in\mathbb{R}\cup\{\infty\}\f$ for all \f$ j\f$

  <LI>\f$ D\f$ is a symmetric positive-semidefinite \f$ n\times n\f$ matrix (the
  quadratic objective function),

  <LI>\f$ \qpc\f$ is an \f$ n\f$-dimensional vector (the linear objective
  function), and
  <LI>\f$ c_0\f$ is a constant.

  </UL>

  This class is simply a wrapper for existing iterators, and it does not
  copy the program data.

  It frequently happens that all values in one of the vectors from
  above are the same, for example if the system \f$ Ax\qprel b\f$ is
  actually a system of equations \f$ Ax=b\f$. To get an iterator over such a
  vector, it is not necessary to store multiple copies of the value in
  some container; an instance of the class `Const_oneset_iterator<T>`,
  constructed from the value in question, does the job more efficiently.

  \cgalModels{QuadraticProgram}

  \cgalHeading{Example}

  \ref QP_solver/first_qp_from_iterators.cpp

  The following example for the simpler model
  `Nonnegative_linear_program_from_iterators<A_it, B_it, R_it, C_it>`
  should give you a flavor of the use of this
  model in practice.

  \ref QP_solver/solve_convex_hull_containment_lp.h

  \ref QP_solver/convex_hull_containment.cpp

  \sa `QuadraticProgram`
  \sa `Quadratic_program<NT>`
  \sa `Quadratic_program_from_mps<NT>`

*/
template< typename A_it, typename B_it, typename R_it, typename FL_it, typename L_it, typename FU_it, typename U_it, typename D_it, typename C_it >
class Quadratic_program_from_iterators {
public:

  /// \name Creation
  /// @{

  /*!
    constructs `qp` from given random-access iterators and the constant `c0`. The passed iterators are merely stored, no copying of the program data takes place. How these iterators are supposed to encode the quadratic program is
    described in `QuadraticProgram`.
  */
  Quadratic_program_from_iterators(int n, int m,
                                   const A_it& a,
                                   const B_it& b,
                                   const R_it& r,
                                   const FL_it& fl,
                                   const L_it& l,
                                   const FU_it& fu,
                                   const U_it& u,
                                   const D_it& d,
                                   const C_it& c,
                                   const std::iterator_traits<C_it>value_type& c0 = 0
    );

  /// @}

}; /* end Quadratic_program_from_iterators */

/*!
  \ingroup PkgQPSolverClasses

  An object of class `Quadratic_program_from_mps` describes a convex quadratic program of the
  general form

  \f{eqnarray*}{
  \mbox{(QP)}& \mbox{minimize}
  & \qpx^{T}D\qpx+\qpc^{T}\qpx+c_0 \\
  &\mbox{subject to} & A\qpx\qprel \qpb, \\
  & & \qpl \leq \qpx \leq \qpu
  \f}
  in \f$ n\f$ real variables \f$ \qpx=(x_0,\ldots,x_{n-1})\f$.

  Here,
  <UL>
  <LI>\f$ A\f$ is an \f$ m\times n\f$ matrix (the constraint matrix),
  <LI>\f$ \qpb\f$ is an \f$ m\f$-dimensional vector (the right-hand side),
  <LI>\f$ \qprel\f$ is an \f$ m\f$-dimensional vector of relations
  from \f$ \{\leq, =, \geq\}\f$,

  <LI>\f$ \qpl\f$ is an \f$ n\f$-dimensional vector of lower
  bounds for \f$ \qpx\f$, where \f$ l_j\in\mathbb{R}\cup\{-\infty\}\f$ for all \f$ j\f$
  <LI>\f$ \qpu\f$ is an \f$ n\f$-dimensional vector of upper bounds for
  \f$ \qpx\f$, where \f$ u_j\in\mathbb{R}\cup\{\infty\}\f$ for all \f$ j\f$

  <LI>\f$ D\f$ is a symmetric positive-semidefinite \f$ n\times n\f$ matrix (the
  quadratic objective function),

  <LI>\f$ \qpc\f$ is an \f$ n\f$-dimensional vector (the linear objective
  function), and
  <LI>\f$ c_0\f$ is a constant.

  </UL>

  If \f$ D=0\f$, the program is
  a linear program; if the variable bounds are \f$ x\geq 0\f$, we have a
  nonnegative program.

  The program data are read from an input stream in `MPSFormat`. This is
  a commonly used format for encoding linear and quadratic programs that
  is understood by many solvers. All values are expected to be readable
  into type `NT`. The constructed program can be further manipulated
  by using the set-methods below.

  \cgalModels{QuadraticProgram,LinearProgram,NonnegativeQuadraticProgram,NonnegativeLinearProgram}

  \cgalHeading{Example}

  \ref QP_solver/first_qp_from_mps.cpp

  \ref QP_solver/first_lp_from_mps.cpp

  \ref QP_solver/first_nonnegative_qp_from_mps.cpp

  \ref QP_solver/first_nonnegative_lp_from_mps.cpp

  \sa `Quadratic_program<NT>`
  \sa `Quadratic_program_from_iterators<A_it, B_it, R_it, FL_it, L_it, FU_it, U_it, D_it, C_it>`
  \sa `Linear_program_from_iterators<A_it, B_it, R_it, FL_it, L_it, FU_it, U_it, C_it>`
  \sa `Nonnegative_quadratic_program_from_iterators<A_it, B_it, R_it, D_it, C_it>`
  \sa `Nonnegative_linear_program_from_iterators<A_it, B_it, R_it, C_it>`

*/
template< typename NT >
class Quadratic_program_from_mps {
public:

  /// \name Types
  /// @{

  /*!
    The number type of the program entries.
  */
  typedef unspecified_type NT;

  /// @}

  /// \name Creation
  /// @{

  /*!
    reads `qp` from the input stream `in`.
  */
  Quadratic_program_from_mps(std::istream& in);

  /// @}

  /// \name Operations
  /// @{

  /*!
    returns `true` if and only if an
    MPS-encoded quadratic program could be extracted from the input stream.
  */
  bool is_valid() const;

  /*!
    if !`qp``.is_valid()`,
    this method returns an error message explaining why the input does not
    conform to the `MPSFormat`.
  */
  const std::string& get_error() const;

  /*!
    returns the name of the \f$ j\f$-th variable. \pre \f$ j\f$ must not refer to a variable that has been added later, using one of the set methods below.
  */
  const std::string& variable_name_by_index (int j) const;

  /*!
    returns the index of the variable with name `name`. If there is
    no variable with this name, the result is \f$ -1\f$.
  */
  int variable_index_by_name (const std::string& name) const;

  /*!
    returns the name of the \f$ i\f$-th constraint. \pre \f$ i\f$ must not refer to a constraint that has been added later, using one of the set methods below.
  */
  const std::string& constraint_name_by_index (int i) const;

  /*!
    returns the index of the constraint with name `name`. If there is
    no constraint with this name, the result is \f$ -1\f$.
  */
  int constraint_index_by_name (const std::string& name) const;

  /*!
    returns `true` if and only if
    `qp` is a linear program.
  */
  bool is_linear() const;

  /*!
    returns `true` if and only if
    `qp` is a nonnegative program.
  */
  bool is_nonnegative() const;

  /*!
    sets the entry \f$ A_{ij}\f$
    in column \f$ j\f$ and row \f$ i\f$ of the constraint matrix \f$ A\f$ of `qp` to
    `val`. An existing entry is overwritten. `qp` is enlarged if
    necessary to accommodate this entry.
  */
  void set_a (int j, int i, const NT& val);

  /*!
    sets the entry \f$ b_i\f$
    of `qp` to `val`. An existing entry is overwritten.
    `qp` is enlarged if necessary to accommodate this entry.
  */
  void set_b (int i, const NT& val);

  /*!
    sets the entry \f$ \qprel_i\f$ of `qp` to `rel`. `CGAL::SMALLER`
    means that the \f$ i\f$-th constraint is of type "\f$ \leq\f$", `CGAL::EQUAL`
    means "\f$ =\f$", and `CGAL::LARGER` encodes "\f$ \geq\f$". An existing entry
    is overwritten. `qp` is enlarged if necessary to accommodate this entry.
  */
  void set_r (int i, CGAL::Comparison_result rel);

  /*!
    if `is_finite`, this sets the entry \f$ l_j\f$ of `qp` to `val`,
    otherwise it sets \f$ l_j\f$ to \f$ -\infty\f$. An existing entry is overwritten.
    `qp` is enlarged if necessary to accommodate this entry.
  */
  void set_l (int j, bool is_finite, const NT& val = NT(0));

  /*!
    if `is_finite`, this sets the entry \f$ u_j\f$ of `qp` to `val`,
    otherwise it sets \f$ u_j\f$ to \f$ \infty\f$. An existing entry is overwritten.
    `qp` is enlarged if necessary to accommodate this entry.
  */
  void set_u (int j, bool is_finite, const NT& val = NT(0));

  /*!
    sets the entry \f$ c_j\f$
    of `qp` to `val`. An existing entry is overwritten.
    `qp` is enlarged if necessary to accommodate this entry.
  */
  void set_c (int j, const NT& val);

  /*!
    sets the entry \f$ c_0\f$
    of `qp` to `val`. An existing entry is overwritten.
  */
  void set_c0 (const NT& val);

  /*!
    sets the entries
    \f$ 2D_{ij}\f$ and \f$ 2D_{ji}\f$ of `qp` to `val`. Existing entries are
    overwritten. `qp` is enlarged if necessary to accommodate these entries.
    \pre `j <= i`
  */
  void set_d (int i, int j, const NT& val);

  /// @}

}; /* end Quadratic_program_from_mps */
} /* end namespace CGAL */

namespace CGAL {

/*!
  \ingroup PkgQPSolverClasses

  An object of class `Quadratic_program` describes a convex quadratic program of the form
  \f{eqnarray*}{
  \mbox{(QP)}& \mbox{minimize}
  & \qpx^{T}D\qpx+\qpc^{T}\qpx+c_0 \\
  &\mbox{subject to} & A\qpx\qprel \qpb, \\
  & & \qpl \leq \qpx \leq \qpu
  \f}
  in \f$ n\f$ real variables \f$ \qpx=(x_0,\ldots,x_{n-1})\f$.

  Here,
  <UL>
  <LI>\f$ A\f$ is an \f$ m\times n\f$ matrix (the constraint matrix),
  <LI>\f$ \qpb\f$ is an \f$ m\f$-dimensional vector (the right-hand side),
  <LI>\f$ \qprel\f$ is an \f$ m\f$-dimensional vector of relations
  from \f$ \{\leq, =, \geq\}\f$,

  <LI>\f$ \qpl\f$ is an \f$ n\f$-dimensional vector of lower
  bounds for \f$ \qpx\f$, where \f$ l_j\in\mathbb{R}\cup\{-\infty\}\f$ for all \f$ j\f$
  <LI>\f$ \qpu\f$ is an \f$ n\f$-dimensional vector of upper bounds for
  \f$ \qpx\f$, where \f$ u_j\in\mathbb{R}\cup\{\infty\}\f$ for all \f$ j\f$

  <LI>\f$ D\f$ is a symmetric positive-semidefinite \f$ n\times n\f$ matrix (the
  quadratic objective function),

  <LI>\f$ \qpc\f$ is an \f$ n\f$-dimensional vector (the linear objective
  function), and
  <LI>\f$ c_0\f$ is a constant.

  </UL>

  If \f$ D=0\f$, the program is
  a linear program; if the variable bounds are \f$ x\geq 0\f$, we have a
  nonnegative program.

  This class allows you to build your program entry by entry, using
  the set-methods below.

  If you only need to wrap existing (random-access)
  iterators over your own data, then you may use any of the four models
  `Quadratic_program_from_iterators<A_it, B_it, R_it, FL_it, L_it, FU_it, U_it, D_it, C_it>`,
  `Linear_program_from_iterators<A_it, B_it, R_it, FL_it, L_it, FU_it, U_it, C_it>`,
  `Nonnegative_quadratic_program_from_iterators<A_it, B_it, R_it, D_it, C_it>`, and
  `Nonnegative_linear_program_from_iterators<A_it, B_it, R_it, C_it>`.

  If you want to read a quadratic program in `MPSFormat` from a stream,
  please use the model `Quadratic_program_from_mps<NT>`.

  \cgalModels{QuadraticProgram,LinearProgram,NonnegativeQuadraticProgram,NonnegativeLinearProgram}

  \cgalHeading{Example}

  \ref QP_solver/first_qp.cpp

  \ref QP_solver/first_lp.cpp

  \ref QP_solver/first_nonnegative_qp.cpp

  \ref QP_solver/first_nonnegative_lp.cpp

  \ref QP_solver/invert_matrix.cpp

  \sa `Quadratic_program_from_iterators<A_it, B_it, R_it, FL_it, L_it, FU_it, U_it, D_it, C_it>`
  \sa `Linear_program_from_iterators<A_it, B_it, R_it, FL_it, L_it, FU_it, U_it, C_it>`
  \sa `Nonnegative_quadratic_program_from_iterators<A_it, B_it, R_it, D_it, C_it>`
  \sa `Nonnegative_linear_program_from_iterators<A_it, B_it, R_it, C_it>`
  \sa `Quadratic_program_from_mps<NT>`

*/
template< typename NT >
class Quadratic_program {
public:

  /// \name Types
  /// @{

  /*!
    The number type of the program entries.
  */
  typedef unspecified_type NT;

  /// @}

  /// \name Creation
  /// @{

  /*!
    constructs a quadratic program with no variables and no constraints, ready
    for data to be added. Unless relations are explicitly set, they will
    be of type `default_r`. Unless bounds are explicitly set, they
    will be as specified by `default_fl` (finite lower bound?),
    `default_l` (lower bound value if lower bound is finite),
    `default_fu` (finite upper bound?), and
    `default_l` (upper bound value if upper bound is finite). If all
    parameters take their default values, we thus get equality constraints
    and bounds \f$ x\geq0\f$ by default. Numerical entries that are not
    explicitly set will default to \f$ 0\f$.\pre if `default_fl == default_fu == true`, then `default_l <= default_u`.
  */
  Quadratic_program
  (CGAL::Comparison_result default_r = CGAL::EQUAL,
   bool default_fl = true,
   const NT& default_l = 0,
   bool default_fu = false,
   const NT& default_u = 0);

  /// @}

  /// \name Operations
  /// @{

  /*!
    returns `true` if and only if
    `qp` is a linear program.
  */
  bool is_linear() const;

  /*!
    returns `true` if and only if
    `qp` is a nonnegative program.
  */
  bool is_nonnegative() const;

  /*!
    sets the entry \f$ A_{ij}\f$
    in column \f$ j\f$ and row \f$ i\f$ of the constraint matrix \f$ A\f$ of `qp` to
    `val`. An existing entry is overwritten. `qp` is enlarged if
    necessary to accommodate this entry.
  */
  void set_a (int j, int i, const NT& val);

  /*!
    sets the entry \f$ b_i\f$
    of `qp` to `val`. An existing entry is overwritten.
    `qp` is enlarged if necessary to accommodate this entry.
  */
  void set_b (int i, const NT& val);

  /*!
    sets the entry \f$ \qprel_i\f$ of `qp` to `rel`. `CGAL::SMALLER`
    means that the \f$ i\f$-th constraint is of type "\f$ \leq\f$", `CGAL::EQUAL`
    means "\f$ =\f$", and `CGAL::LARGER` encodes "\f$ \geq\f$". An existing entry
    is overwritten. `qp` is enlarged if necessary to accommodate this entry.
  */
  void set_r (int i, CGAL::Comparison_result rel);

  /*!
    if `is_finite`, this sets the entry \f$ l_j\f$ of `qp` to `val`,
    otherwise it sets \f$ l_j\f$ to \f$ -\infty\f$. An existing entry is overwritten.
    `qp` is enlarged if necessary to accommodate this entry.
  */
  void set_l (int j, bool is_finite, const NT& val = NT(0));

  /*!
    if `is_finite`, this sets the entry \f$ u_j\f$ of `qp` to `val`,
    otherwise it sets \f$ u_j\f$ to \f$ \infty\f$. An existing entry is overwritten.
    `qp` is enlarged if necessary to accommodate this entry.
  */
  void set_u (int j, bool is_finite, const NT& val = NT(0));

  /*!
    sets the entry \f$ c_j\f$
    of `qp` to `val`. An existing entry is overwritten.
    `qp` is enlarged if necessary to accommodate this entry.
  */
  void set_c (int j, const NT& val);

  /*!
    sets the entry \f$ c_0\f$
    of `qp` to `val`. An existing entry is overwritten.
  */
  void set_c0 (const NT& val);

  /*!
    sets the entries
    \f$ 2D_{ij}\f$ and \f$ 2D_{ji}\f$ of `qp` to `val`. Existing entries are
    overwritten. `qp` is enlarged if necessary to accommodate these entries.
    \pre `j <= i`
  */
  void set_d (int i, int j, const NT& val);

  /// @}

}; /* end Quadratic_program */
} /* end namespace CGAL */
