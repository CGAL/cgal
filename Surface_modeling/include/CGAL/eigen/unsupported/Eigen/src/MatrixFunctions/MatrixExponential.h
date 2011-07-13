// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2009, 2010 Jitse Niesen <jitse@maths.leeds.ac.uk>
//
// Eigen is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// Alternatively, you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of
// the License, or (at your option) any later version.
//
// Eigen is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License and a copy of the GNU General Public License along with
// Eigen. If not, see <http://www.gnu.org/licenses/>.

#ifndef EIGEN_MATRIX_EXPONENTIAL
#define EIGEN_MATRIX_EXPONENTIAL

#ifdef _MSC_VER
  template <typename Scalar> Scalar log2(Scalar v) { return std::log(v)/std::log(Scalar(2)); }
#endif


/** \ingroup MatrixFunctions_Module
  * \brief Class for computing the matrix exponential.
  * \tparam MatrixType type of the argument of the exponential,
  * expected to be an instantiation of the Matrix class template.
  */
template <typename MatrixType>
class MatrixExponential {

  public:

    /** \brief Constructor.
      * 
      * The class stores a reference to \p M, so it should not be
      * changed (or destroyed) before compute() is called.
      *
      * \param[in] M  matrix whose exponential is to be computed.
      */
    MatrixExponential(const MatrixType &M);

    /** \brief Computes the matrix exponential.
      *
      * \param[out] result  the matrix exponential of \p M in the constructor.
      */
    template <typename ResultType> 
    void compute(ResultType &result);

  private:

    // Prevent copying
    MatrixExponential(const MatrixExponential&);
    MatrixExponential& operator=(const MatrixExponential&);

    /** \brief Compute the (3,3)-Pad&eacute; approximant to the exponential.
     *
     *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
     *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
     *
     *  \param[in] A   Argument of matrix exponential
     */
    void pade3(const MatrixType &A);

    /** \brief Compute the (5,5)-Pad&eacute; approximant to the exponential.
     *
     *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
     *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
     *
     *  \param[in] A   Argument of matrix exponential
     */
    void pade5(const MatrixType &A);

    /** \brief Compute the (7,7)-Pad&eacute; approximant to the exponential.
     *
     *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
     *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
     *
     *  \param[in] A   Argument of matrix exponential
     */
    void pade7(const MatrixType &A);

    /** \brief Compute the (9,9)-Pad&eacute; approximant to the exponential.
     *
     *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
     *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
     *
     *  \param[in] A   Argument of matrix exponential
     */
    void pade9(const MatrixType &A);

    /** \brief Compute the (13,13)-Pad&eacute; approximant to the exponential.
     *
     *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
     *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
     *
     *  \param[in] A   Argument of matrix exponential
     */
    void pade13(const MatrixType &A);

    /** \brief Compute Pad&eacute; approximant to the exponential.
     *
     * Computes \c m_U, \c m_V and \c m_squarings such that
     * \f$ (V+U)(V-U)^{-1} \f$ is a Pad&eacute; of
     * \f$ \exp(2^{-\mbox{squarings}}M) \f$ around \f$ M = 0 \f$. The
     * degree of the Pad&eacute; approximant and the value of
     * squarings are chosen such that the approximation error is no
     * more than the round-off error.
     *
     * The argument of this function should correspond with the (real
     * part of) the entries of \c m_M.  It is used to select the
     * correct implementation using overloading.
     */
    void computeUV(double);

    /** \brief Compute Pad&eacute; approximant to the exponential.
     *
     *  \sa computeUV(double);
     */
    void computeUV(float);

    typedef typename internal::traits<MatrixType>::Scalar Scalar;
    typedef typename NumTraits<Scalar>::Real RealScalar;

    /** \brief Reference to matrix whose exponential is to be computed. */
    typename internal::nested<MatrixType>::type m_M;

    /** \brief Even-degree terms in numerator of Pad&eacute; approximant. */
    MatrixType m_U;

    /** \brief Odd-degree terms in numerator of Pad&eacute; approximant. */
    MatrixType m_V;

    /** \brief Used for temporary storage. */
    MatrixType m_tmp1;

    /** \brief Used for temporary storage. */
    MatrixType m_tmp2;

    /** \brief Identity matrix of the same size as \c m_M. */
    MatrixType m_Id;

    /** \brief Number of squarings required in the last step. */
    int m_squarings;

    /** \brief L1 norm of m_M. */
    float m_l1norm;
};

template <typename MatrixType>
MatrixExponential<MatrixType>::MatrixExponential(const MatrixType &M) :
  m_M(M),
  m_U(M.rows(),M.cols()),
  m_V(M.rows(),M.cols()),
  m_tmp1(M.rows(),M.cols()),
  m_tmp2(M.rows(),M.cols()),
  m_Id(MatrixType::Identity(M.rows(), M.cols())),
  m_squarings(0),
  m_l1norm(static_cast<float>(M.cwiseAbs().colwise().sum().maxCoeff()))
{
  /* empty body */
}

template <typename MatrixType>
template <typename ResultType> 
void MatrixExponential<MatrixType>::compute(ResultType &result)
{
  computeUV(RealScalar());
  m_tmp1 = m_U + m_V;	// numerator of Pade approximant
  m_tmp2 = -m_U + m_V;	// denominator of Pade approximant
  result = m_tmp2.partialPivLu().solve(m_tmp1);
  for (int i=0; i<m_squarings; i++)
    result *= result;		// undo scaling by repeated squaring
}

template <typename MatrixType>
EIGEN_STRONG_INLINE void MatrixExponential<MatrixType>::pade3(const MatrixType &A)
{
  const Scalar b[] = {120., 60., 12., 1.};
  m_tmp1.noalias() = A * A;
  m_tmp2 = b[3]*m_tmp1 + b[1]*m_Id;
  m_U.noalias() = A * m_tmp2;
  m_V = b[2]*m_tmp1 + b[0]*m_Id;
}

template <typename MatrixType>
EIGEN_STRONG_INLINE void MatrixExponential<MatrixType>::pade5(const MatrixType &A)
{
  const Scalar b[] = {30240., 15120., 3360., 420., 30., 1.};
  MatrixType A2 = A * A;
  m_tmp1.noalias() = A2 * A2;
  m_tmp2 = b[5]*m_tmp1 + b[3]*A2 + b[1]*m_Id;
  m_U.noalias() = A * m_tmp2;
  m_V = b[4]*m_tmp1 + b[2]*A2 + b[0]*m_Id;
}

template <typename MatrixType>
EIGEN_STRONG_INLINE void MatrixExponential<MatrixType>::pade7(const MatrixType &A)
{
  const Scalar b[] = {17297280., 8648640., 1995840., 277200., 25200., 1512., 56., 1.};
  MatrixType A2 = A * A;
  MatrixType A4 = A2 * A2;
  m_tmp1.noalias() = A4 * A2;
  m_tmp2 = b[7]*m_tmp1 + b[5]*A4 + b[3]*A2 + b[1]*m_Id;
  m_U.noalias() = A * m_tmp2;
  m_V = b[6]*m_tmp1 + b[4]*A4 + b[2]*A2 + b[0]*m_Id;
}

template <typename MatrixType>
EIGEN_STRONG_INLINE void MatrixExponential<MatrixType>::pade9(const MatrixType &A)
{
  const Scalar b[] = {17643225600., 8821612800., 2075673600., 302702400., 30270240.,
  		      2162160., 110880., 3960., 90., 1.};
  MatrixType A2 = A * A;
  MatrixType A4 = A2 * A2;
  MatrixType A6 = A4 * A2;
  m_tmp1.noalias() = A6 * A2;
  m_tmp2 = b[9]*m_tmp1 + b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*m_Id;
  m_U.noalias() = A * m_tmp2;
  m_V = b[8]*m_tmp1 + b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*m_Id;
}

template <typename MatrixType>
EIGEN_STRONG_INLINE void MatrixExponential<MatrixType>::pade13(const MatrixType &A)
{
  const Scalar b[] = {64764752532480000., 32382376266240000., 7771770303897600.,
  		      1187353796428800., 129060195264000., 10559470521600., 670442572800.,
  		      33522128640., 1323241920., 40840800., 960960., 16380., 182., 1.};
  MatrixType A2 = A * A;
  MatrixType A4 = A2 * A2;
  m_tmp1.noalias() = A4 * A2;
  m_V = b[13]*m_tmp1 + b[11]*A4 + b[9]*A2; // used for temporary storage
  m_tmp2.noalias() = m_tmp1 * m_V;
  m_tmp2 += b[7]*m_tmp1 + b[5]*A4 + b[3]*A2 + b[1]*m_Id;
  m_U.noalias() = A * m_tmp2;
  m_tmp2 = b[12]*m_tmp1 + b[10]*A4 + b[8]*A2;
  m_V.noalias() = m_tmp1 * m_tmp2;
  m_V += b[6]*m_tmp1 + b[4]*A4 + b[2]*A2 + b[0]*m_Id;
}

template <typename MatrixType>
void MatrixExponential<MatrixType>::computeUV(float)
{
  if (m_l1norm < 4.258730016922831e-001) {
    pade3(m_M);
  } else if (m_l1norm < 1.880152677804762e+000) {
    pade5(m_M);
  } else {
    const float maxnorm = 3.925724783138660f;
    m_squarings = std::max(0, (int)ceil(log2(m_l1norm / maxnorm)));
    MatrixType A = m_M / std::pow(Scalar(2), Scalar(static_cast<RealScalar>(m_squarings)));
    pade7(A);
  }
}

template <typename MatrixType>
void MatrixExponential<MatrixType>::computeUV(double)
{
  if (m_l1norm < 1.495585217958292e-002) {
    pade3(m_M);
  } else if (m_l1norm < 2.539398330063230e-001) {
    pade5(m_M);
  } else if (m_l1norm < 9.504178996162932e-001) {
    pade7(m_M);
  } else if (m_l1norm < 2.097847961257068e+000) {
    pade9(m_M);
  } else {
    const double maxnorm = 5.371920351148152;
    m_squarings = std::max(0, (int)ceil(log2(m_l1norm / maxnorm)));
    MatrixType A = m_M / std::pow(Scalar(2), Scalar(m_squarings));
    pade13(A);
  }
}

/** \ingroup MatrixFunctions_Module
  *
  * \brief Proxy for the matrix exponential of some matrix (expression).
  *
  * \tparam Derived  Type of the argument to the matrix exponential.
  *
  * This class holds the argument to the matrix exponential until it
  * is assigned or evaluated for some other reason (so the argument
  * should not be changed in the meantime). It is the return type of
  * MatrixBase::exp() and most of the time this is the only way it is
  * used.
  */
template<typename Derived> struct MatrixExponentialReturnValue
: public ReturnByValue<MatrixExponentialReturnValue<Derived> >
{
    typedef typename Derived::Index Index;
  public:
    /** \brief Constructor.
      *
      * \param[in] src %Matrix (expression) forming the argument of the
      * matrix exponential.
      */
    MatrixExponentialReturnValue(const Derived& src) : m_src(src) { }

    /** \brief Compute the matrix exponential.
      *
      * \param[out] result the matrix exponential of \p src in the
      * constructor.
      */
    template <typename ResultType>
    inline void evalTo(ResultType& result) const
    {
      const typename Derived::PlainObject srcEvaluated = m_src.eval();
      MatrixExponential<typename Derived::PlainObject> me(srcEvaluated);
      me.compute(result);
    }

    Index rows() const { return m_src.rows(); }
    Index cols() const { return m_src.cols(); }

  protected:
    const Derived& m_src;
  private:
    MatrixExponentialReturnValue& operator=(const MatrixExponentialReturnValue&);
};

namespace internal {
template<typename Derived>
struct traits<MatrixExponentialReturnValue<Derived> >
{
  typedef typename Derived::PlainObject ReturnType;
};
}

template <typename Derived>
const MatrixExponentialReturnValue<Derived> MatrixBase<Derived>::exp() const
{
  eigen_assert(rows() == cols());
  return MatrixExponentialReturnValue<Derived>(derived());
}

#endif // EIGEN_MATRIX_EXPONENTIAL
