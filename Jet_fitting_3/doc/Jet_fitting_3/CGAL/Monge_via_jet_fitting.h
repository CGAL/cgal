namespace CGAL {

/*!
\ingroup PkgJetFitting3Ref

The class `Monge_via_jet_fitting` is designed to perform the estimation of the
local differential quantities at a given point. The point range is
given by a pair of input iterators, and it is assumed that the point
where the calculation is carried out is the point that the begin
iterator refers to.
The results are stored in an instance of the nested class `Monge_form`,
the particular information returned depending on the degrees specified
for the polynomial fitting and for the Monge form.

If `CGAL_EIGEN3_ENABLED` is defined, `LocalKernel` and `SvdTraits`
template parameters have defaults, `Simple_cartesian<double>` and `Eigen_svd` respectively.

\tparam DataKernel provides the geometric classes and tools
corresponding to the input points, and also members of the
`Monge_form` class.

\tparam LocalKernel provides
the geometric classes and tools required by local
computations.

\tparam SvdTraits features the linear
algebra algorithm required by the fitting method.  The scalar type, `SvdTraits::FT`, must be the same as that of   the `LocalKernel` concept : `LocalKernel::FT`.

\sa `Eigen_svd`
\sa `Monge_form`

\note This class requires the \ref thirdpartyEigen library.
*/
template< typename DataKernel, typename LocalKernel, typename SvdTraits >
class Monge_via_jet_fitting {
public:

/// \name Types
/// @{

/*!
\ingroup PkgJetFitting3Ref

The class `Monge_form` stores the Monge representation, i.e., the Monge
coordinate system and the coefficients of the Monge form in this
system.

The usual insert operator (`operator<<`) is overloaded for
`Monge_form`, it gives the Monge coordinate system (the origin and an
orthonormal basis) and the coefficients of the Monge form in this
system.


\sa `Monge_via_jet_fitting`

*/
class Monge_form {
public:

/// \name Types
/// @{

/*!

*/
typedef typename DataKernel::FT FT;

/*!

*/
typedef typename DataKernel::Point_3 Point_3;

/*!

*/
typedef typename DataKernel::Vector_3 Vector_3;

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
Monge_form();

/// @}

/// \name Access Functions
/// @{

/*!
Point on the fitted surface where
differential quantities are computed.
*/
Point_3 origin() const;

/// @}

/// \name Monge Basis
/// @{

/*!

*/
Vector_3 maximal_principal_direction() const;

/*!

*/
Vector_3 minimal_principal_direction() const;

/*!

*/
Vector_3 normal_direction() const;

/// @}

/// \name Monge Coefficients
/// @{

/*!
\f$ i=0\f$ for the maximum and \f$ i=1\f$ for the minimum.
*/
FT principal_curvatures(size_t i)const;

/*!
\f$ 0 \leq i \leq3\f$
*/
FT third_order_coefficients(size_t i) const;

/*!
\f$ 0 \leq i \leq4\f$
*/
FT fourth_order_coefficients(size_t i) const;

/// @}

/// \name Operations
/// @{

/*!
change principal basis and Monge coefficients so that the
`given_normal` and the Monge normal make an acute angle.
If
`given_normal * normal_direction()  < 0` then change the orientation: if
\f$ z=g(x,y)\f$ in the basis (d1,d2,n) then in the basis (d2,d1,-n)
\f$ z=h(x,y)=-g(y,x)\f$.
*/
void comply_wrt_given_normal(const Vector_3& given_normal);

/// @}

}; /* end Monge_form */

/*!

*/
typedef DataKernel Data_kernel;

/*!

*/
typedef LocalKernel Local_kernel;

/*!

*/
typedef typename Local_kernel::FT FT;

/*!

*/
typedef typename Local_kernel::Vector_3 Vector_3;

/*!
see the page `Monge_via_jet_fitting::Monge_form`.
*/
typedef unspecified_type Monge_form;

/// @}

/// \name Creation
/// @{

/*!
default constructor
*/
Monge_via_jet_fitting();

/// @}

/// \name Operations
/// @{

/*!
This operator performs all the computations. The \f$ N\f$ input points are
given by the `InputIterator` parameters which value-type are
`Data_kernel::Point_3`, `d` is the degree of the fitted
polynomial, \c d' is the degree of the expected Monge
coefficients. \pre \f$ N \geq N_{d}:=(d+1)(d+2)/2\f$, \f$ 1 \leq d' \leq\min(d,4) \f$.
*/
template <class InputIterator> Monge_form
operator()(InputIterator begin, InputIterator end, size_t d,
size_t d');

/*!
condition number of the linear fitting system.
*/
FT condition_number();

/*!
pca eigenvalues and eigenvectors, the pca_basis has always 3 such pairs.
\pre \f$ i\f$ ranges from 0 to 2.
*/
std::pair<FT, Vector_3> pca_basis(size_t i);

/// @}

}; /* end Monge_via_jet_fitting */
} /* end namespace CGAL */
