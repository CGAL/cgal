/*!
\ingroup PkgSolverInterfaceConcepts
\cgalConcept

A concept that describes the set of methods used to define and solve a
quadratic programming (`qp`) problem of the general form:
<center>
\f{eqnarray*}{
& \mbox{minimize} & \frac{1}{2}\mathbf{x}^{T}P\mathbf{x} + \mathbf{q}^{T}\mathbf{x} + r \\
& \mbox{subject to} & \mathbf{l} \leq A\mathbf{x} \leq \mathbf{u}
\f}
</center>
in \f$ n \f$ real variables \f$ \mathbf{x} = (x_0, \ldots, x_{n-1}) \f$ and \f$ m \f$ constraints.

Here,
<UL>
<LI>\f$ P \f$ is a symmetric positive-semidefinite \f$ n \times n\f$ matrix (the quadratic objective function),
<LI>\f$ \mathbf{q} \f$ is an \f$ n \f$-dimensional vector (the linear objective function),
<LI>\f$ r \f$ is a constant,
<LI>\f$ A \f$ is an \f$ m\times n\f$ matrix (the constraint matrix),
<LI>\f$ \mathbf{l} \f$ is an \f$ m \f$-dimensional vector of lower constraint bounds,
where \f$ l_i \in \mathbb{R} \cup \{-\infty\} \f$ for all \f$ i \f$,
<LI>\f$ \mathbf{u} \f$ is an \f$ m \f$-dimensional vector of upper constraint bounds,
where \f$ u_i \in \mathbb{R} \cup \{+\infty\} \f$ for all \f$ i \f$.
</UL>

\cgalHasModelsBegin
\cgalHasModels{CGAL::OSQP_quadratic_program_traits<T>}
\cgalHasModelsEnd
*/
class QuadraticProgramTraits {

public:

  /// \name Memory
  /// @{

  /*!
    Allocates memory for `n` variables and `m` constraints in `qp`.
  */
  void resize(const std::size_t n, const std::size_t m) { }

  /// @}

  /// \name Initialization
  /// @{

  /*!
    Sets the entries `Pij` and `Pji` of `qp` to `value`.
  */
  void set_P(const std::size_t i, const std::size_t j, const FT value) { }

  /*!
    Sets the entry `qi` of `qp` to `value`.
  */
  void set_q(const std::size_t i, const FT value) { }

  /*!
    Sets the entry `r` of `qp` to `value`.
  */
  void set_r(const FT value) { }

  /*!
    Sets the entry `Aij` in the row `i` and column `j` of
    the constraint matrix `A` of `qp` to `value`.
  */
  void set_A(const std::size_t i, const std::size_t j, const FT value) { }

  /*!
    Sets the entry `li` of `qp` to `value`.
  */
  void set_l(const std::size_t i, const FT value) { }

  /*!
    Sets the entry `ui` of `qp` to `value`.
  */
  void set_u(const std::size_t i, const FT value) { }

  /// @}

  /// \name Solution
  /// @{

  /*!
    \brief solves the quadratic program.

    Number of values in `solution` equals to the number `n` of values in the vector `x`.

    \tparam OutIterator
    a model of `OutputIterator` that accepts values of type `FieldNumberType`

    \param solution
    an output iterator with the solution

    \returns a status of the computation `success == true`
  */
  template<typename OutIterator>
  bool solve(OutIterator solution) { }

  /// @}
};
