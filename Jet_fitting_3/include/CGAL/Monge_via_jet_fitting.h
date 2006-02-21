#ifndef _MONGE_VIA_JET_FITTING_H_
#define _MONGE_VIA_JET_FITTING_H_

#include <CGAL/basic.h>
#include <CGAL/circulator.h>
#include <CGAL/Linear_algebraCd.h>
#include "jet_fitting_3_assertions.h"

//#include <CGAL/eigen.h> //for ALTERNATIVE  with CGAL eigen code

#include <math.h>

CGAL_BEGIN_NAMESPACE

int fact(int n)
{
  if (n == 0)
    return(1);
  return(n * fact(n-1));
}
////////////////////// CLASS Monge_rep ////////////////////////
template <class DataKernel>
class Monge_rep {
public: 
  typedef typename DataKernel::FT        DFT;
  typedef typename DataKernel::Point_3   DPoint;
  typedef typename DataKernel::Vector_3  DVector;
protected:
  //point on the fitted surface where diff quantities are computed
  DPoint m_origin_pt;
  //the monge trihedron (d1,d2,n) is orthonormal direct
  DVector m_d1;//maximal ppal dir
  DVector m_d2;//minimal ppal dir
  DVector m_n;//normal direction
  //coeff = (k1, k2, //ppal curv
  //         b0, b1, b2, b3, //third order
  //         c0, c1, c2, c3, c4) //fourth order
  //     if (degree==1) no coeff needed
  std::vector<DFT> m_coefficients;
  
public:
  //constructor
  Monge_rep() {
    m_origin_pt  = DPoint(0.,0.,0.); 
    m_d1 = DVector(0.,0.,0.);
    m_d2 = DVector(0.,0.,0.);
    m_n = DVector(0.,0.,0.);
    m_coefficients = std::vector<DFT>();
  }
  ~Monge_rep() {}
  //access
  const DPoint origin_pt() const { return m_origin_pt; }
  DPoint& origin_pt() { return m_origin_pt; }
  const DVector d1() const { return m_d1; }
  DVector& d1() { return m_d1; }
  const DVector d2() const { return m_d2; }
  DVector& d2() { return m_d2; }
  const DVector n() const { return m_n; }
  DVector& n() { return m_n; }
  const std::vector<DFT> coefficients() const { return m_coefficients; }
  std::vector<DFT>& coefficients() { return m_coefficients; }

  void set_up(int degree);
};

//if d>=1, number of coeffs = (d+1)(d+2)/2 -4. 
//we remove cst, linear and the xy coeff which vanish
template <class DataKernel>
void Monge_rep<DataKernel>::
set_up(int degree) {
  if ( degree >= 1 ) std::fill_n(back_inserter(m_coefficients),
				 (degree+1)*(degree+2)/2-4, 0.);
}

////////////////////// CLASS Monge_info ////////////////////////
template <class LocalKernel>
class Monge_info {  
public:
  typedef typename LocalKernel::FT        LFT;
  typedef typename LocalKernel::Vector_3  LVector;
protected:  
  LFT m_pca_eigen_vals[3];
  LVector m_pca_eigen_vecs[3];
  LFT m_cond_nb;//of the least square system

public:
  //constructor
  Monge_info() {
    m_cond_nb = 0.;
    std::fill_n(m_pca_eigen_vals, 3, 0.);
    std::fill_n(m_pca_eigen_vecs, 3, LVector()); 
  }
  //access
  const LFT* pca_eigen_vals() const { return m_pca_eigen_vals; }
  LFT* pca_eigen_vals() { return m_pca_eigen_vals; }
  const LVector* pca_eigen_vecs() const { return m_pca_eigen_vecs; }
  LVector* pca_eigen_vecs() { return m_pca_eigen_vecs; }
  const LFT cond_nb() const { return m_cond_nb; }
  LFT& cond_nb() { return m_cond_nb; }
};


////////////////////// CLASS Monge_via_jet_fitting ////////////////////////
template < class DataKernel, class LocalKernel, class LinAlgTraits>  
class Monge_via_jet_fitting {
public:
  typedef DataKernel   Data_Kernel;
  typedef LocalKernel  Local_Kernel;
  typedef typename std::vector<typename Data_Kernel::Point_3>::iterator Range_Iterator;
  typedef Monge_rep<Data_Kernel>   Monge_rep;
  typedef Monge_info<Local_Kernel> Monge_info;

public:
  Monge_via_jet_fitting(Range_Iterator begin, Range_Iterator end, 
			int d, int dprime, 
			Monge_rep &monge_rep, Monge_info &monge_info);

protected:
  typedef typename Local_Kernel::FT       LFT;
  typedef typename Local_Kernel::Point_3  LPoint;
  typedef typename Local_Kernel::Vector_3 LVector;
  typedef CGAL::Aff_transformation_3<Local_Kernel> Aff_transformation;

  typedef typename Data_Kernel::FT       DFT;
  typedef typename Data_Kernel::Point_3  DPoint;

  typedef typename LinAlgTraits::Vector LAVector;
  typedef typename LinAlgTraits::Matrix LAMatrix;

protected:
  int deg;
  int deg_monge;
  int nb_d_jet_coeff;
  int nb_input_pts;
  LFT preconditionning;
  CGAL::Sqrt<LFT> Lsqrt;

  //translate_p0 changes the origin of the world to p0 the first point 
  //  of the input data points
  //change_world2fitting (coord of a vector in world) = coord of this 
  //  vector in fitting. The matrix tranform has as lines the coord of
  //  the basis vectors of fitting in the world coord. 
  //idem for change_fitting2monge
  Aff_transformation translate_p0, change_world2fitting,
    change_fitting2monge;

  //eigen val and vect stored in monge_info, 
  // change_world2fitting is computed 
  void compute_PCA(Range_Iterator begin, Range_Iterator end,
		   Monge_info &monge_info); 

  //Coordinates of input points are computed in the fitting basis with 
  //  p0 as origin.
  //Preconditionning is computed, M and Z are filled
  void fill_matrix(Range_Iterator begin, Range_Iterator end,
		   int d, LAMatrix& M, LAVector& Z);

  //A is computed, solving MA=Z in the ls sense
  //Preconditionning is needed
  //the condition number of the matrix M is stored in monge_info 
  void solve_linear_system(LAMatrix &M, LAVector &A, const LAVector &Z,
			   Monge_info& monge_info);
  
  //Classical differential geometric calculus
  //change_fitting2monge is computed
  //if deg_monge =1 only 1st order info
  //if deg_monge >= 2 2nd order info are computed
  void compute_Monge_basis(const LAVector &A, Monge_rep& monge_rep);

  //if deg_monge >=3 then 3rd (and 4th) order info are computed
  void compute_Monge_coefficients(const LAVector &A, int dprime, 
				  Monge_rep& monge_rep);

  //for a trihedron (v1,v2,v3) switches v1 to -v1 if det(v1,v2,v3) < 0
  void switch_to_direct_orientation(LVector& v1, const LVector& v2,
				   const LVector& v3);
};

//-------------------------------------------------------------------
// Implementation
//------------------------------------------------------------------

template < class DataKernel, class LocalKernel, class LinAlgTraits>  
Monge_via_jet_fitting<DataKernel, LocalKernel, LinAlgTraits>::
Monge_via_jet_fitting(Range_Iterator begin, Range_Iterator end, 
		      int d, int  dprime, 
		      Monge_rep<Data_Kernel> &monge_rep,  
		      Monge_info<Local_Kernel> &monge_info)
{
  // precondition: on the degrees, jet and monge
  CGAL_precondition( (d >=1) && (dprime >= 1) 
		             && (dprime <= 4) && (dprime <= d) );
  this->deg = d;
  this->deg_monge = dprime;
  this->nb_d_jet_coeff = (d+1)*(d+2)/2;
  this->nb_input_pts = end - begin;
  // precondition: solvable ls system
  CGAL_precondition( nb_input_pts >= nb_d_jet_coeff );

  //Initialize
  monge_rep.set_up(dprime);
  //for the system MA=Z
  LAMatrix M(nb_input_pts, nb_d_jet_coeff);
  LAVector A(nb_d_jet_coeff);
  LAVector Z(nb_input_pts);

  compute_PCA(begin, end, monge_info);
  fill_matrix(begin, end, d, M, Z);//with precond
  solve_linear_system(M, A, Z, monge_info);  //correct with precond
  compute_Monge_basis(A, monge_rep);
  if ( dprime >= 3) compute_Monge_coefficients(A, dprime, monge_rep);
} 

template < class DataKernel, class LocalKernel, class LinAlgTraits>  
void Monge_via_jet_fitting<DataKernel, LocalKernel, LinAlgTraits>::
compute_PCA(Range_Iterator begin, Range_Iterator end,
	    Monge_info &monge_info)
{
  LAMatrix Cov(3,3);
  LAVector eval(3);
  LAMatrix evec(3,3);
  
  int n = this->nb_input_pts;
  LFT x, y, z,
    sumX = 0., sumY = 0., sumZ = 0.,
    sumX2 = 0., sumY2 = 0., sumZ2 = 0.,
    sumXY = 0., sumXZ = 0., sumYZ = 0., 
    xx, yy, zz, xy, xz, yz;
  
  for (; begin != end; begin++)
    {
      x = (*begin).x();
      y = (*begin).y();
      z = (*begin).z();   
      sumX += x / n;
      sumY += y / n;
      sumZ += z / n;
      sumX2 += x * x / n;
      sumY2 += y * y / n;
      sumZ2 += z * z / n;
      sumXY += x * y / n;
      sumXZ += x * z / n;
      sumYZ += y * z / n;
    }
  xx = sumX2 - sumX * sumX;
  yy = sumY2 - sumY * sumY;
  zz = sumZ2 - sumZ * sumZ;
  xy = sumXY - sumX * sumY;
  xz = sumXZ - sumX * sumZ;
  yz = sumYZ - sumY * sumZ;
  Cov[0][0] = xx;
  Cov[0][1] = xy;
  Cov[0][2] = xz;
  Cov[1][0] = xy;
  Cov[1][1] = yy;
  Cov[1][2] = yz;
  Cov[2][0] = xz;
  Cov[2][1] = yz;
  Cov[2][2] = zz;
  // solve for eigenvalues and eigenvectors.
  // eigen values are sorted in descending order, 
  // eigen vectors are sorted in accordance.
  LinAlgTraits::eigen_symm_algo(Cov, eval, evec);
 
  //store in monge_info
  for (int i=0; i<3; i++)
    {
      monge_info.pca_eigen_vals()[i] = eval[i];//implicit cast LAFT->LFT
      LVector temp_vect(evec[0][i],evec[1][i],evec[2][i]);
      monge_info.pca_eigen_vecs()[i] = temp_vect;
    }

//   //ALTERNATIVE  with CGAL eigen code
//   // assemble covariance matrix as a
//   // semi-definite matrix. 
//   // Matrix numbering:
//   // 0
//   // 1 2
//   // 3 4 5
//   LFT covariance[6] = {xx,xy,yy,xz,yz,zz};

//   // solve for eigenvalues and eigenvectors.
//   // eigen values are sorted in descending order, 
//   // eigen vectors are sorted in accordance.
//   LFT eigen_values[3];
//   LFT eigen_vectors[9];
//     CGAL::CGALi::eigen_symmetric<LFT>(covariance,3,eigen_vectors,eigen_values);
//   //store in monge_info
//   for (int i=0; i<3; i++)
//     {
//       monge_info.pca_eigen_vals()[i] = eigen_values[i];//implicit cast LAFT->LFT
//     }
//   LVector v1(eigen_vectors[0],eigen_vectors[1],eigen_vectors[2]);
//       monge_info.pca_eigen_vecs()[0] = v1;
//   LVector v2(eigen_vectors[3],eigen_vectors[4],eigen_vectors[5]);
//       monge_info.pca_eigen_vecs()[1] = v2;
//   LVector v3(eigen_vectors[6],eigen_vectors[7],eigen_vectors[8]);
//       monge_info.pca_eigen_vecs()[2] = v3;
//     ///end ALTERNATIVE with CGAL eigen code

  switch_to_direct_orientation(monge_info.pca_eigen_vecs()[0],
			       monge_info.pca_eigen_vecs()[1],
			       monge_info.pca_eigen_vecs()[2]);
 
  //Store the change of basis W->F
  const LVector* pca_vecs = monge_info.pca_eigen_vecs();
  Aff_transformation 
    change_basis (pca_vecs[0][0], pca_vecs[0][1], pca_vecs[0][2], 
		  pca_vecs[1][0], pca_vecs[1][1], pca_vecs[1][2],
		  pca_vecs[2][0], pca_vecs[2][1], pca_vecs[2][2]);
  this->change_world2fitting = change_basis;
}

template < class DataKernel, class LocalKernel, class LinAlgTraits>  
void Monge_via_jet_fitting<DataKernel, LocalKernel, LinAlgTraits>::
fill_matrix(Range_Iterator begin, Range_Iterator end,
	    int d, LAMatrix &M, LAVector &Z)
{
  //origin of fitting coord system = first input data point
  LPoint point0 = *begin;
  //transform coordinates of sample points with a
  //translation ($-p$) and multiplication by $ P_{W\rightarrow F}$.
  LPoint orig(0.,0.,0.);
  LVector v_point0_orig(orig - point0);
  Aff_transformation transl(CGAL::TRANSLATION, v_point0_orig);
  this->translate_p0 = transl;
  Aff_transformation transf_points = this->change_world2fitting *
    this->translate_p0;
  
  //compute and store transformed points
  std::vector<LPoint> pts_in_fitting_basis;
  CGAL_For_all(begin,end){//implicit cast DPoint->LPoint
    LPoint cur_pt = transf_points(*begin);
    pts_in_fitting_basis.push_back(cur_pt);
  }
  
  //Compute preconditionning
  LFT precond = 0.;
  typename std::vector<LPoint>::iterator itb = pts_in_fitting_basis.begin(),
    ite = pts_in_fitting_basis.end();
  CGAL_For_all(itb,ite) precond += std::fabs(itb->x()) + std::fabs(itb->y());
  precond /= 2*this->nb_input_pts;
  this->preconditionning = precond;
  //fill matrices M and Z
  itb = pts_in_fitting_basis.begin();
  int line_count = 0;
  LFT x, y;
  CGAL_For_all(itb,ite) {
    x = itb->x();
    y = itb->y();
    Z[line_count] = itb->z();
    for (int k=0; k <= d; k++) for (int i=0; i<=k; i++)
      M[line_count][k*(k+1)/2+i] =
	std::pow(x,k-i)*std::pow(y,i)
	/(fact(i)*fact(k-i)*std::pow(this->preconditionning,k));
    line_count++;
  }
}

template < class DataKernel, class LocalKernel, class LinAlgTraits>  
void Monge_via_jet_fitting<DataKernel, LocalKernel, LinAlgTraits>::
solve_linear_system(LAMatrix &M, LAVector &A, const LAVector &Z,
		    Monge_info& monge_info)
{
 LinAlgTraits::solve_ls_svd_algo(M, A, Z, monge_info.cond_nb()); 
  for (int k=0; k <= this->deg; k++) for (int i=0; i<=k; i++)
    A[k*(k+1)/2+i] /= std::pow(this->preconditionning,k);
}

template < class DataKernel, class LocalKernel, class LinAlgTraits>  
void Monge_via_jet_fitting<DataKernel, LocalKernel, LinAlgTraits>::
compute_Monge_basis(const LAVector &A, Monge_rep& monge_rep)
{
  // only 1st order info.
  if ( this->deg_monge == 1 ) {  
    LPoint orig_monge(0., 0., A[0]);
    LVector  normal(-A[1], -A[2], 1.);
    LFT norm2 = normal * normal;
    normal = normal / Lsqrt(norm2);
    monge_rep.origin_pt() = 
      (this->translate_p0.inverse() * 
       this->change_world2fitting.inverse()) (orig_monge );
    monge_rep.n() = this->change_world2fitting.inverse()(normal);
  }
  // else (deg_monge >= 2) then 2nd order info are computed
  else {
  //bi-index to uni-index conversion : A(i,j)=A[(i+j)(i+j+1)/2+j]
  LPoint orig_monge(0., 0., A[0]);
  //normal = Xu crossprod Xv
  LVector Xu(1.,0.,A[1]), Xv(0.,1.,A[2]), normal(-A[1], -A[2], 1.);
  LFT norm2 = normal * normal;
  normal = normal / Lsqrt(norm2);

  //Surface in fitting_basis : X(u,v)=(u,v,J_A(u,v))
  //in the basis Xu=(1,0,A[1]), Xv=(0,1,A[2]), Weingarten=-I^{-1}II
  //first fond form I=(e,f,f,g)
  //                 =(Xu.Xu, Xu.Xv, Xu.Xv, Xv.Xv)
  //second fond form II=(l,m,m,n)/norm2^(1/2)
  //                   =(n.Xuu, n.Xuv, n.Xuv, n.Xvv)
  //ppal curv are the opposite of the eigenvalues of Weingarten or the
  //  eigenvalues of weingarten = -Weingarten = I^{-1}II
  typedef typename CGAL::Linear_algebraCd<LFT>::Matrix Matrix;

  LFT e = 1+A[1]*A[1], f = A[1]*A[2], g = 1+A[2]*A[2],
    l = A[3], m = A[4], n = A[5];
  Matrix  weingarten(2,2,0.);
  weingarten(0,0) = (g*l-f*m)/ (Lsqrt(norm2)*norm2);
  weingarten(0,1) = (g*m-f*n)/ (Lsqrt(norm2)*norm2);
  weingarten(1,0) = (e*m-f*l)/ (Lsqrt(norm2)*norm2);
  weingarten(1,1) = (e*n-f*m)/ (Lsqrt(norm2)*norm2);
  // Y, Z are normalized GramSchmidt of Xu, Xv
  // Xu->Y=Xu/||Xu||;
  // Xv->Z=Xv-(Xu.Xv)Xu/||Xu||^2;
  // Z-> Z/||Z||
  LVector Y, Z;
  LFT normXu = Lsqrt( Xu*Xu );
  Y = Xu / normXu;
  LFT XudotXv = Xu * Xv;
  Z = Xv - XudotXv * Xu / (normXu*normXu);
  LFT normZ = Lsqrt( Z*Z );
  Z = Z / normZ;
  Matrix change_XuXv2YZ(2,2,0.);
  change_XuXv2YZ(0,0) = 1 / normXu;
  change_XuXv2YZ(0,1) = -XudotXv / (normXu * normXu * normZ);
  change_XuXv2YZ(1,0) = 0;
  change_XuXv2YZ(1,1) = 1 / normZ;
  LFT det = 0.;
  Matrix inv = CGAL::Linear_algebraCd<LFT>::inverse ( change_XuXv2YZ, det );
  //in the new orthonormal basis (Y,Z) of the tangent plane :
  weingarten = inv *(1/det) * weingarten * change_XuXv2YZ;
  
  //switch to LinAlgTraits for diagonalization of weingarten
  LAMatrix W(2,2);
  for (int i=0; i<=1; i++) for (int j=0; j<=1; j++) 
    W[i][j] = weingarten(i,j);
  LAVector eval(2);
  LAMatrix evec(2,2);

  LinAlgTraits::eigen_symm_algo(W, eval, evec);
  LVector d_max = evec[0][0]*Y + evec[1][0]*Z,
    d_min = evec[0][1]*Y + evec[1][1]*Z;

  switch_to_direct_orientation(d_max, d_min, normal);
  Aff_transformation change_basis (d_max[0], d_max[1], d_max[2], 
				   d_min[0], d_min[1], d_min[2],
				   normal[0], normal[1], normal[2]);
  this->change_fitting2monge = change_basis;

  //store the monge basis origin and vectors with their world coord
  //store ppal curv
  monge_rep.origin_pt() = 
    (this->translate_p0.inverse() * 
     this->change_world2fitting.inverse()) (orig_monge );
  monge_rep.d1() = this->change_world2fitting.inverse()(d_max);
  monge_rep.d2() = this->change_world2fitting.inverse()(d_min);
  monge_rep.n()  = this->change_world2fitting.inverse()(normal);
  monge_rep.coefficients()[0] = eval[0];
  monge_rep.coefficients()[1] = eval[1];
  }
  //end else
}

template < class DataKernel, class LocalKernel, class LinAlgTraits>  
void Monge_via_jet_fitting<DataKernel, LocalKernel, LinAlgTraits>::
compute_Monge_coefficients(const LAVector &A, int dprime, 
			   Monge_rep& monge_rep)
{
  //One has the equation w=J_A(u,v) of the fitted surface S 
  // in the fitting_basis
  //Substituing (u,v,w)=change_fitting2monge^{-1}(x,y,z)
  //One has the equation f(x,y,z)=0 on this surface S in the monge
  //  basis
  //The monge form of the surface at the origin is the bivariate fct
  //   g(x,y) s.t. f(x,y,g(x,y))=0
  //voir les calculs Maple dans monge.mws
  //Notations are f123= d^3f/dxdydz
  //              g(x,y)=sum (gij x^i y^j/ i!j!) with 
  //              g00=g10=g01=g11=0, g20=kmax, g02=kmin
  //
  //g(x,y)= 1/2*(k1x^2 +k2y^2) 
  //       +1/6*(b0x^3 +3b1x^2y +3b2xy^2 +b3y^3) 
  //       +1/24*(c0x^4 +4c1x^3y +6c2x^2y^2 +4c3xy^3 +c4y^4)
  //       +...
  // p stores change_fitting2monge^{-1}=change_fitting2monge^{T}
  LFT p[3][3];
  p[0][0] = this->change_fitting2monge.m(0,0);
  p[1][0] = this->change_fitting2monge.m(0,1);
  p[2][0] = this->change_fitting2monge.m(0,2);
  p[0][1] = this->change_fitting2monge.m(1,0);
  p[1][1] = this->change_fitting2monge.m(1,1);
  p[2][1] = this->change_fitting2monge.m(1,2);
  p[0][2] = this->change_fitting2monge.m(2,0);
  p[1][2] = this->change_fitting2monge.m(2,1);
  p[2][2] = this->change_fitting2monge.m(2,2);

  //     note f1 = f2 = f12 = 0 
  //     LFT f1 = A[1] * p[0][0] + A[2] * p[1][0] - p[2][0];
  //     LFT f2 = A[2] * p[1][1] + A[1] * p[0][1] - p[2][1];
  //     LFT f12 = 
  //     2 * A[3] * p[0][0] * p[0][1]
  //     + 2 * A[5] * p[1][0] * p[1][1]
  //     + A[4] * p[0][1] * p[1][0] 
  //     + A[4] * p[0][0] * p[1][1];
  //         -f11 / f3 = kmax
  //         -f22 / f3 = kmin 
 
  LFT f3 = A[1] * p[0][2] + A[2] * p[1][2] - p[2][2];
  LFT f11 =
    2 * A[4] * p[0][0] * p[1][0]
    + 2 * A[5] * p[1][0] * p[1][0]
    + 2 * A[3] * p[0][0] * p[0][0];
  LFT f13 =
    A[4] * p[0][0] * p[1][2]
    + A[4] * p[0][2] * p[1][0]
    + 2 * A[5] * p[1][0] * p[1][2]
    + 2 * A[3] * p[0][0] * p[0][2];
  LFT f22 =
    2 * A[4] * p[0][1] * p[1][1]
    + 2 * A[5] * p[1][1] * p[1][1]
    + 2 * A[3] * p[0][1] * p[0][1];
  LFT f23 =
    A[4] * p[0][1] * p[1][2]
    + 2 * A[5] * p[1][1] * p[1][2]
    + A[4] * p[0][2] * p[1][1]
    + 2 * A[3] * p[0][1] * p[0][2];
  LFT f33 =
    2 * A[5] * p[1][2] * p[1][2]
    + 2 * A[3] * p[0][2] * p[0][2]
    + 2 * A[4] * p[0][2] * p[1][2];
  LFT f111 =
    6 * A[8] * p[0][0] * p[1][0] * p[1][0]
    + 6 * A[7] * p[0][0] * p[0][0] * p[1][0]
    + 6 * A[6] * p[0][0] * p[0][0] * p[0][0]
    + 6 * A[9] * p[1][0] * p[1][0] * p[1][0];
  LFT f222 =
    6 * A[7] * p[0][1] * p[0][1] * p[1][1]
    + 6 * A[8] * p[0][1] * p[1][1] * p[1][1]
    + 6 * A[9] * p[1][1] * p[1][1] * p[1][1]
    + 6 * A[6] * p[0][1] * p[0][1] * p[0][1];
  LFT f112 =
    2 * A[7] * p[0][0] * p[0][0] * p[1][1]
    + 6 * A[6] * p[0][0] * p[0][0] * p[0][1]
    + 2 * A[8] * p[0][1] * p[1][0] * p[1][0]
    + 4 * A[8] * p[0][0] * p[1][0] * p[1][1]
    + 6 * A[9] * p[1][0] * p[1][0] * p[1][1]
    + 4 * A[7] * p[0][0] * p[0][1] * p[1][0];
  LFT f122 =
    4 * A[8] * p[0][1] * p[1][0] * p[1][1]
    + 2 * A[8] * p[0][0] * p[1][1] * p[1][1]
    + 6 * A[6] * p[0][0] * p[0][1] * p[0][1]
    + 2 * A[7] * p[0][1] * p[0][1] * p[1][0]
    + 4 * A[7] * p[0][0] * p[0][1] * p[1][1]
    + 6 * A[9] * p[1][0] * p[1][1] * p[1][1];
  LFT f113 = 
    6*A[6]*p[0][0]*p[0][0]*p[0][2]
    +6*A[9]*p[1][0]*p[1][0]*p[1][2]
    +2*A[7]*p[0][0]*p[0][0]*p[1][2]
    +2*A[8]*p[0][2]*p[1][0]*p[1][0]
    +4*A[7]*p[0][0]*p[0][2]*p[1][0]
    +4*A[8]*p[0][0]*p[1][0]*p[1][2];
  LFT f223 = 
    2*A[8]*p[0][2]*p[1][1]*p[1][1]
    +6*A[6]*p[0][1]*p[0][1]*p[0][2]
    +6*A[9]*p[1][1]*p[1][1]*p[1][2]
    +2*A[7]*p[0][1]*p[0][1]*p[1][2]
    +4*A[7]*p[0][1]*p[0][2]*p[1][1]
    +4*A[8]*p[0][1]*p[1][1]*p[1][2];
  LFT f123 = 
    2*A[8]*p[0][2]*p[1][0]*p[1][1]
    +2*A[7]*p[0][0]*p[0][1]*p[1][2]
    +2*A[7]*p[0][0]*p[0][2]*p[1][1]
    +6*A[9]*p[1][0]*p[1][1]*p[1][2]
    +2*A[7]*p[0][1]*p[0][2]*p[1][0]
    +6*A[6]*p[0][0]*p[0][1]*p[0][2]
    +2*A[8]*p[0][0]*p[1][1]*p[1][2]
    +2*A[8]*p[0][1]*p[1][0]*p[1][2];

  LFT b0 = 1/(f3*f3)*(-f111*f3+3*f13*11);
  LFT b1 = 1/(f3*f3)*(-f112*f3+f23*f11);
  LFT b2 = 1/(f3*f3)*(-f122*f3+f13*f22);
  LFT b3 = -1/(f3*f3)*(f222*f3-3*f23*f22);
  
  monge_rep.coefficients()[2] = b0;
  monge_rep.coefficients()[3] = b1;
  monge_rep.coefficients()[4] = b2;
  monge_rep.coefficients()[5] = b3;

  if ( dprime == 4 )
    {
      LFT f1111 = 
	24*A[13]*p[0][0]*p[1][0]*p[1][0]*p[1][0]
	+24*A[12]*p[0][0]*p[0][0]*p[1][0]*p[1][0]
	+24*A[11]*p[0][0]*p[0][0]*p[0][0]*p[1][0]
	+24*A[14]*p[1][0]*p[1][0]*p[1][0]*p[1][0]
	+24*A[10]*p[0][0]*p[0][0]*p[0][0]*p[0][0];
      LFT f1112 = 
	6*A[13]*p[0][1]*p[1][0]*p[1][0]*p[1][0]
	+18*A[13]*p[0][0]*p[1][0]*p[1][0]*p[1][1]
	+24*A[10]*p[0][0]*p[0][0]*p[0][0]*p[0][1]
	+12*A[12]*p[0][0]*p[0][1]*p[1][0]*p[1][0]
	+18*A[11]*p[0][0]*p[0][0]*p[0][1]*p[1][0]
	+24*A[14]*p[1][0]*p[1][0]*p[1][0]*p[1][1]
	+6*A[11]*p[0][0]*p[0][0]*p[0][0]*p[1][1]
	+12*A[12]*p[0][0]*p[0][0]*p[1][0]*p[1][1];
      LFT f1122 = 
	12*A[11]*p[0][0]*p[0][0]*p[0][1]*p[1][1]
	+12*A[13]*p[0][0]*p[1][0]*p[1][1]*p[1][1]
	+12*A[13]*p[0][1]*p[1][0]*p[1][0]*p[1][1]
	+16*A[12]*p[0][0]*p[0][1]*p[1][0]*p[1][1]
	+12*A[11]*p[0][0]*p[0][1]*p[0][1]*p[1][0]
	+24*A[10]*p[0][0]*p[0][0]*p[0][1]*p[0][1]
	+4*A[12]*p[0][1]*p[0][1]*p[1][0]*p[1][0]
	+4*A[12]*p[0][0]*p[0][0]*p[1][1]*p[1][1]
	+24*A[14]*p[1][0]*p[1][0]*p[1][1]*p[1][1];
      LFT f1222 = 
	6*A[13]*p[0][0]*p[1][1]*p[1][1]*p[1][1]
	+24*A[10]*p[0][0]*p[0][1]*p[0][1]*p[0][1]
	+24*A[14]*p[1][0]*p[1][1]*p[1][1]*p[1][1]
	+6*A[11]*p[0][1]*p[0][1]*p[0][1]*p[1][0]
	+18*A[11]*p[0][0]*p[0][1]*p[0][1]*p[1][1]
	+12*A[12]*p[0][0]*p[0][1]*p[1][1]*p[1][1]
	+12*A[12]*p[0][1]*p[0][1]*p[1][0]*p[1][1]
	+18*A[13]*p[0][1]*p[1][0]*p[1][1]*p[1][1];
      LFT f2222 =
	24*A[13]*p[0][1]*p[1][1]*p[1][1]*p[1][1]
	+24*A[11]*p[0][1]*p[0][1]*p[0][1]*p[1][1]
	+24*A[12]*p[0][1]*p[0][1]*p[1][1]*p[1][1]
	+24*A[10]*p[0][1]*p[0][1]*p[0][1]*p[0][1]
	+24*A[14]*p[1][1]*p[1][1]*p[1][1]*p[1][1];

      LFT c0 =
	-1/(f3*f3*f3)*(f1111*(f3*f3)-4*f13*f3*f111+12*f13*f13*f11-6*f113*f3*f11+3*f33*f11*f11);
      LFT c1 =
	1/(f3*f3*f3)*(f23*f3*f111+3*f3*f123*f11+3*f13*f3*f112-f1112*(f3*f3)-6*f13*f23*f11); 
      LFT c2 =
	1/(f3*f3*f3)*(-f33*f22*f11+f113*f3*f22+2*f13*f3*f122-2*f13*f13*f22+f223*f3*f11+2*f23*f3*f112-2*f23*f23*f11-f1122*(f3*f3)); 
      LFT c3 =
	1/(f3*f3*f3)*(-f1222*(f3*f3)-6*f13*f23*f22+3*f123*f3*f22+f13*f3*f222+3*f23*f3*f122); 
      LFT c4 =
	-1/(f3*f3*f3)*(f2222*(f3*f3)+3*f33*f22*f22-6*f223*f3*f22-4*f23*f3*f222+12*f23*f23*f22) ; 
      
      monge_rep.coefficients()[6] = c0;
      monge_rep.coefficients()[7] = c1;
      monge_rep.coefficients()[8] = c2;
      monge_rep.coefficients()[9] = c3;
      monge_rep.coefficients()[10] = c4;
    }
}

template < class DataKernel, class LocalKernel, class LinAlgTraits>  
void Monge_via_jet_fitting<DataKernel, LocalKernel, LinAlgTraits>::
switch_to_direct_orientation(LVector& v1, const LVector& v2,
			    const LVector& v3) 
{
  CGAL::Sign orientation = CGAL::sign_of_determinant3x3(v1[0], v2[0], v3[0],
							v1[1], v2[1], v3[1],
							v1[2], v2[2], v3[2]);
  if (orientation == CGAL::NEGATIVE) v1 = -v1;
}

CGAL_END_NAMESPACE

#endif //_MONGE_VIA_JET_FITTING_H_



