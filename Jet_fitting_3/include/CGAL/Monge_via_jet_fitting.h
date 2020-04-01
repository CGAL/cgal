// Copyright (c) 2007  INRIA Sophia-Antipolis (France), INRIA Lorraine LORIA.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Marc Pouget and Frédéric Cazals
#ifndef CGAL_MONGE_VIA_JET_FITTING_H_
#define CGAL_MONGE_VIA_JET_FITTING_H_

#include <CGAL/license/Jet_fitting_3.h>


#include <CGAL/Simple_cartesian.h>
#include <CGAL/circulator.h>
#include <CGAL/Linear_algebraCd.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Default_diagonalize_traits.h>
#include <math.h>
#include <utility>
#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_svd.h>
#endif

namespace CGAL {

inline
unsigned int fact(unsigned int n){
  unsigned int i, p=1;
  for(i=2; i<=n; i++) p *= i;
  return p;
}

////////////////////// CLASS Monge_via_jet_fitting ////////////////////////
#ifdef CGAL_EIGEN3_ENABLED
template < class DataKernel, class LocalKernel = Simple_cartesian<double>, class SvdTraits = Eigen_svd >
#else
template < class DataKernel, class LocalKernel, class SvdTraits >
#endif
  class Monge_via_jet_fitting {
 public:
//////////////////////begin nested CLASS Monge_form ///////////////////
class Monge_form {
 public:
  typedef typename DataKernel::FT        FT;
  typedef typename DataKernel::Point_3   Point_3;
  typedef typename DataKernel::Vector_3  Vector_3;
 protected:
  //point on the fitted surface where diff quantities are computed
  Point_3 m_origin_pt;
  //the monge trihedron (d1,d2,n) is orthonormal direct
  Vector_3 m_d1;   //maximal ppal dir
  Vector_3 m_d2;   //minimal ppal dir
  Vector_3 m_n;    //normal direction
  //coeff = (k1, k2, //ppal curv
  //         b0, b1, b2, b3, //third order
  //         c0, c1, c2, c3, c4) //fourth order
  //     if (degree==1) no coeff needed
  std::vector<FT> m_coefficients;

public:
  //constructor
  Monge_form() {
    m_origin_pt  = Point_3(0.,0.,0.);
    m_d1 = Vector_3(0.,0.,0.);
    m_d2 = Vector_3(0.,0.,0.);
    m_n = Vector_3(0.,0.,0.);
    m_coefficients = std::vector<FT>();
  }
  ~Monge_form() {}
  //access
  const Point_3 origin() const { return m_origin_pt; }
  Point_3& origin() { return m_origin_pt; }
  const Vector_3 maximal_principal_direction() const { return m_d1; }
  Vector_3& maximal_principal_direction() { return m_d1; }
  const Vector_3 minimal_principal_direction() const { return m_d2; }
  Vector_3& minimal_principal_direction() { return m_d2; }
  const Vector_3 normal_direction() const { return m_n; }
  Vector_3& normal_direction() { return m_n; }
  const std::vector<FT> coefficients() const { return m_coefficients; }
  std::vector<FT>& coefficients() { return m_coefficients; }

  const FT principal_curvatures(size_t i) const {
    CGAL_precondition( (i == 0 || i == 1) && coefficients().size() >=2 );
    return coefficients()[i]; }
  const FT third_order_coefficients(size_t i) const {
    CGAL_precondition( i <= 3 && coefficients().size() >=6 );
    return coefficients()[i+2]; }
  const FT fourth_order_coefficients(size_t i) const {
    CGAL_precondition( i <= 4 && coefficients().size() >=11 );
    return coefficients()[i+6]; }

  //if d>=2, number of coeffs = (d+1)(d+2)/2 -4.
  //we remove cst, linear and the xy coeff which vanish
  void set_up(std::size_t degree);
  //switch min-max ppal curv/dir wrt a given normal orientation.
  // if given_normal.monge_normal < 0 then change the orientation
  // if z=g(x,y) in the basis (d1,d2,n) then in the basis (d2,d1,-n)
  // z=h(x,y)=-g(y,x)
  void comply_wrt_given_normal(const Vector_3& given_normal);
  void dump_verbose(std::ostream& out_stream) const;
  void dump_4ogl(std::ostream& out_stream, const FT scale);
};
//////////////////////end nested CLASS Monge_form /////////////////////

//continue main class Monge_via_jet_fitting ////////
 public:
 typedef  DataKernel   Data_kernel;
 typedef  LocalKernel  Local_kernel;

 //used to convert number types, points and vectors back and forth
 typedef NT_converter<typename Local_kernel::FT, typename  Data_kernel::FT> L2D_NTconverter;
 Cartesian_converter<Data_kernel, Local_kernel> D2L_converter;
 Cartesian_converter<Local_kernel, Data_kernel> L2D_converter;

 typedef typename Local_kernel::FT       FT;
 typedef typename Local_kernel::Point_3  Point_3;
 typedef typename Local_kernel::Vector_3 Vector_3;
 typedef CGAL::Aff_transformation_3<Local_kernel> Aff_transformation;

 typedef typename Data_kernel::FT       DFT;

 typedef typename SvdTraits::Vector LAVector;
 typedef typename SvdTraits::Matrix LAMatrix;


 public:
 Monge_via_jet_fitting();
 template <class InputIterator>
 Monge_form operator()(InputIterator begin, InputIterator end,
                       size_t d, size_t dprime);
 const FT condition_number() const {return condition_nb;}
 const std::pair<FT, Vector_3> pca_basis(size_t i) const {
   CGAL_precondition( i<3 );
   return m_pca_basis[i];}

 protected:
 int deg;
 int deg_monge;
  int nb_d_jet_coeff;
  int nb_input_pts;
  FT preconditionning;
  CGAL::Sqrt<FT> Lsqrt;
  FT condition_nb;

  std::vector< std::pair<FT, Vector_3> > m_pca_basis;

  //translate_p0 changes the origin of the world to p0 the first point
  //  of the input data points
  //change_world2fitting (coord of a vector in world) = coord of this
  //  vector in fitting. The matrix tranform has as lines the coord of
  //  the basis vectors of fitting in the world coord.
  //idem for change_fitting2monge
  Aff_transformation translate_p0, change_world2fitting,
    change_fitting2monge;

  //eigen val and vect stored in m_pca_basis
  // change_world2fitting is computed
 template <class InputIterator>
  void compute_PCA(InputIterator begin, InputIterator end);

  //Coordinates of input points are computed in the fitting basis with
  //  p0 as origin.
  //Preconditionning is computed, M and Z are filled
 template <class InputIterator>
  void fill_matrix(InputIterator begin, InputIterator end,
                   std::size_t d, LAMatrix& M, LAVector& Z);
  //A is computed, solving MA=Z in the ls sense, the solution A is stored in Z
  //Preconditionning is needed
  void solve_linear_system(LAMatrix &M, LAVector& Z);

  //Classical differential geometric calculus
  //change_fitting2monge is computed
  //if deg_monge =1 only 1st order info
  //if deg_monge >= 2 2nd order info are computed
  void compute_Monge_basis(const FT* A, Monge_form& monge_form);

  //if deg_monge >=3 then 3rd (and 4th) order info are computed
  void compute_Monge_coefficients(FT* A, std::size_t dprime,
                                  Monge_form& monge_form);

  //for a trihedron (v1,v2,v3) switches v1 to -v1 if det(v1,v2,v3) < 0
  void switch_to_direct_orientation(Vector_3& v1, const Vector_3& v2,
                                   const Vector_3& v3);

    friend
    std::ostream&
    operator<<(std::ostream& out_stream,
               const typename Monge_via_jet_fitting<DataKernel, LocalKernel, SvdTraits>::Monge_form& monge){
      monge.dump_verbose(out_stream);
      return out_stream;
    }
};

//-------------------------------------------------------------------
// Implementation
//------------------------------------------------------------------

// Implementation nested Monge_form //////////////////////////////
//template <class DataKernel>
template < class DataKernel, class LocalKernel, class SvdTraits>
  void Monge_via_jet_fitting<DataKernel, LocalKernel, SvdTraits>::
  Monge_form::
set_up(std::size_t degree) {
  if ( degree >= 2 ) std::fill_n(back_inserter(m_coefficients),
                                 (degree+1)*(degree+2)/2-4, 0.);
}


template < class DataKernel, class LocalKernel, class SvdTraits>
void Monge_via_jet_fitting<DataKernel, LocalKernel, SvdTraits>::Monge_form::
comply_wrt_given_normal(const Vector_3& given_normal)
{
  if ( given_normal*this->normal_direction() < 0 )
    {
      normal_direction() = -normal_direction();
      std::swap(maximal_principal_direction(), minimal_principal_direction());
      if ( coefficients().size() >= 2)
        std::swap(coefficients()[0],coefficients()[1]);
      if ( coefficients().size() >= 6) {
        std::swap(coefficients()[2],coefficients()[5]);
        std::swap(coefficients()[3],coefficients()[4]);}
      if ( coefficients().size() >= 11) {
        std::swap(coefficients()[6],coefficients()[10]);
        std::swap(coefficients()[7],coefficients()[9]);}
      typename std::vector<FT>::iterator itb = coefficients().begin(),
        ite = coefficients().end();
      for (;itb!=ite;itb++) { *itb = -(*itb); }
    }
}

template < class DataKernel, class LocalKernel, class SvdTraits>
void Monge_via_jet_fitting<DataKernel, LocalKernel, SvdTraits>::Monge_form::
dump_verbose(std::ostream& out_stream) const
{
  out_stream << "origin : " << origin() << std::endl
             << "n : " << normal_direction() << std::endl;
  if ( coefficients().size() >= 2)
    out_stream << "d1 : " << maximal_principal_direction() << std::endl
               << "d2 : " << minimal_principal_direction() << std::endl
               << "k1 : " << coefficients()[0] << std::endl
               << "k2 : " << coefficients()[1] << std::endl;
  if ( coefficients().size() >= 6)
    out_stream << "b0 : " << coefficients()[2] << std::endl
               << "b1 : " << coefficients()[3] << std::endl
                << "b2 : " << coefficients()[4] << std::endl
                << "b3 : " << coefficients()[5] << std::endl;
  if ( coefficients().size() >= 11)
    out_stream << "c0 : " << coefficients()[6] << std::endl
               << "c1 : " << coefficients()[7] << std::endl
                << "c2 : " << coefficients()[8] << std::endl
                << "c3 : " << coefficients()[9] << std::endl
                << "c4 : " << coefficients()[10] << std::endl
               << std::endl;
}

template < class DataKernel, class LocalKernel, class SvdTraits>
void Monge_via_jet_fitting<DataKernel, LocalKernel, SvdTraits>::Monge_form::
dump_4ogl(std::ostream& out_stream, const FT scale)
{
  CGAL_precondition( coefficients().size() >= 2 );
  out_stream << origin()  << " "
             << maximal_principal_direction() * scale << " "
             << minimal_principal_direction() * scale << " "
             << coefficients()[0] << " "
             << coefficients()[1] << " "
             << std::endl;
}
//////////////////////////////////////////////////////////////
// Implementation main Monge_via_jet_fiting

template < class DataKernel, class LocalKernel, class SvdTraits>
  Monge_via_jet_fitting<DataKernel, LocalKernel, SvdTraits>::
  Monge_via_jet_fitting()
{
  m_pca_basis = std::vector< std::pair<FT, Vector_3> >(3);
}

template < class DataKernel, class LocalKernel, class SvdTraits>
template <class InputIterator>
  typename  Monge_via_jet_fitting<DataKernel, LocalKernel, SvdTraits>::Monge_form
  Monge_via_jet_fitting<DataKernel, LocalKernel, SvdTraits>::
  operator()(InputIterator begin, InputIterator end,
             size_t d, size_t dprime)
{
  // precondition: on the degrees, jet and monge
  CGAL_precondition( (d >=1) && (dprime >= 1)
                     && (dprime <= 4) && (dprime <= d) );
  this->deg = static_cast<int>(d);
  this->deg_monge = static_cast<int>(dprime);
  this->nb_d_jet_coeff = static_cast<int>((d+1)*(d+2)/2);
  this->nb_input_pts = static_cast<int>(end - begin);
  // precondition: solvable ls system
  CGAL_precondition( nb_input_pts >= nb_d_jet_coeff );

  //Initialize
  Monge_form monge_form;
  monge_form.set_up(dprime);
  //for the system MA=Z
  LAMatrix M(nb_input_pts, nb_d_jet_coeff);
  LAVector Z(nb_input_pts);

  compute_PCA(begin, end);
  fill_matrix(begin, end, d, M, Z);//with precond
  solve_linear_system(M, Z);  //correct with precond
  compute_Monge_basis(Z.vector(), monge_form);
  if ( dprime >= 3) compute_Monge_coefficients(Z.vector(), dprime, monge_form);
  return monge_form;
}

template < class DataKernel, class LocalKernel, class SvdTraits>
template <class InputIterator>
void Monge_via_jet_fitting<DataKernel, LocalKernel, SvdTraits>::
compute_PCA(InputIterator begin, InputIterator end)
{
  int n = this->nb_input_pts;
  FT x, y, z,
    sumX = 0., sumY = 0., sumZ = 0.,
    sumX2 = 0., sumY2 = 0., sumZ2 = 0.,
    sumXY = 0., sumXZ = 0., sumYZ = 0.,
    xx, yy, zz, xy, xz, yz;

  for (; begin != end; begin++)
    {
      Point_3 lp = D2L_converter(*begin);
      x = lp.x();
      y = lp.y();
      z = lp.z();
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

  // assemble covariance matrix as a
  // semi-definite matrix.
  // Matrix numbering:
  // 0 1 2
  //   3 4
  //     5
  std::array<FT, 6> covariance = {{ xx,xy,xz,yy,yz,zz }};
  std::array<FT, 3> eigen_values = {{ 0., 0., 0. }};
  std::array<FT, 9> eigen_vectors = {{ 0., 0., 0. }};

  // solve for eigenvalues and eigenvectors.
  // eigen values are sorted in ascending order,
  // eigen vectors are sorted in accordance.
  CGAL::Default_diagonalize_traits<FT,3>::diagonalize_selfadjoint_covariance_matrix
    (covariance, eigen_values, eigen_vectors);

  //store in m_pca_basis
  for (int i=0; i<3; i++)
    {
      m_pca_basis[i].first =  eigen_values[2-i];
    }

  Vector_3 v1(eigen_vectors[6],eigen_vectors[7],eigen_vectors[8]);
  m_pca_basis[0].second = v1;
  Vector_3 v2(eigen_vectors[3],eigen_vectors[4],eigen_vectors[5]);
  m_pca_basis[1].second = v2;
  Vector_3 v3(eigen_vectors[0],eigen_vectors[1],eigen_vectors[2]);
  m_pca_basis[2].second = v3;
  switch_to_direct_orientation(m_pca_basis[0].second,
                               m_pca_basis[1].second,
                               m_pca_basis[2].second);

  //Store the change of basis W->F
  Aff_transformation
    change_basis (m_pca_basis[0].second[0], m_pca_basis[0].second[1], m_pca_basis[0].second[2],
                  m_pca_basis[1].second[0], m_pca_basis[1].second[1], m_pca_basis[1].second[2],
                  m_pca_basis[2].second[0], m_pca_basis[2].second[1], m_pca_basis[2].second[2]);

   this->change_world2fitting = change_basis;
}

template < class DataKernel, class LocalKernel, class SvdTraits>
template <class InputIterator>
void Monge_via_jet_fitting<DataKernel, LocalKernel, SvdTraits>::
fill_matrix(InputIterator begin, InputIterator end,
            std::size_t d, LAMatrix &M, LAVector& Z)
{
  //origin of fitting coord system = first input data point
  Point_3 point0 = D2L_converter(*begin);
  //transform coordinates of sample points with a
  //translation ($-p$) and multiplication by $ P_{W\rightarrow F}$.
  Point_3 orig(0.,0.,0.);
  Vector_3 v_point0_orig(orig - point0);
  Aff_transformation transl(CGAL::TRANSLATION, v_point0_orig);
  this->translate_p0 = transl;
  Aff_transformation transf_points = this->change_world2fitting *
    this->translate_p0;

  //compute and store transformed points
  std::vector<Point_3> pts_in_fitting_basis;
  pts_in_fitting_basis.reserve(this->nb_input_pts);
  CGAL_For_all(begin,end){
    Point_3 cur_pt = transf_points(D2L_converter(*begin));
    pts_in_fitting_basis.push_back(cur_pt);
  }

  //Compute preconditionning
  FT precond = 0.;
  typename std::vector<Point_3>::iterator itb = pts_in_fitting_basis.begin(),
    ite = pts_in_fitting_basis.end();
  CGAL_For_all(itb,ite) precond += CGAL::abs(itb->x()) + CGAL::abs(itb->y());
  precond /= 2*this->nb_input_pts;
  this->preconditionning = precond;
  //fill matrices M and Z
  itb = pts_in_fitting_basis.begin();
  int line_count = 0;
  FT x, y;
  CGAL_For_all(itb,ite) {
    x = itb->x();
    y = itb->y();
    //  Z[line_count] = itb->z();
    Z.set(line_count,itb->z());
    for (std::size_t k=0; k <= d; k++) {
      for (std::size_t i=0; i<=k; i++) {
        M.set(line_count, k*(k+1)/2+i,
              std::pow(x,static_cast<int>(k-i))
              * std::pow(y,static_cast<int>(i))
              /( fact(static_cast<unsigned int>(i)) *
                 fact(static_cast<unsigned int>(k-i))
                 *std::pow(this->preconditionning,static_cast<int>(k))));
      }
    }
    line_count++;
  }
}

template < class DataKernel, class LocalKernel, class SvdTraits>
void Monge_via_jet_fitting<DataKernel, LocalKernel, SvdTraits>::
solve_linear_system(LAMatrix &M, LAVector& Z)
{
  condition_nb = SvdTraits::solve(M, Z);
  for (int k=0; k <= this->deg; k++) for (int i=0; i<=k; i++)
    // Z[k*(k+1)/2+i] /= std::pow(this->preconditionning,k);
    Z.set( k*(k+1)/2+i, Z(k*(k+1)/2+i) / std::pow(this->preconditionning,k) );
}

template < class DataKernel, class LocalKernel, class SvdTraits>
void Monge_via_jet_fitting<DataKernel, LocalKernel, SvdTraits>::
compute_Monge_basis(const FT* A, Monge_form& monge_form)
{
  // only 1st order info.
  if ( this->deg_monge == 1 ) {
    Point_3 orig_monge(0., 0., A[0]);
    Vector_3  normal(-A[1], -A[2], 1.);
    FT norm2 = normal * normal;
    normal = normal / Lsqrt(norm2);
    monge_form.origin() = L2D_converter(
      (this->translate_p0.inverse() *
       this->change_world2fitting.inverse()) (orig_monge) );
    monge_form.normal_direction() = L2D_converter(this->change_world2fitting.inverse()(normal));
  }
  // else (deg_monge >= 2) then 2nd order info are computed
  else {
  //bi-index to uni-index conversion : A(i,j)=A[(i+j)(i+j+1)/2+j]
  Point_3 orig_monge(0., 0., A[0]);
  //normal = Xu crossprod Xv
  Vector_3 Xu(1.,0.,A[1]), Xv(0.,1.,A[2]), normal(-A[1], -A[2], 1.);
  FT norm2 = normal * normal;
  normal = normal / Lsqrt(norm2);

  //Surface in fitting_basis : X(u,v)=(u,v,J_A(u,v))
  //in the basis Xu=(1,0,A[1]), Xv=(0,1,A[2]), Weingarten=-I^{-1}II
  //first fond form I=(e,f,f,g)
  //                 =(Xu.Xu, Xu.Xv, Xu.Xv, Xv.Xv)
  //second fond form II=(l,m,m,n)/norm2^(1/2)
  //                   =(n.Xuu, n.Xuv, n.Xuv, n.Xvv)
  //ppal curv are the opposite of the eigenvalues of Weingarten or the
  //  eigenvalues of weingarten = -Weingarten = I^{-1}II
  typedef typename CGAL::Linear_algebraCd<FT>::Matrix Matrix;

  FT e = 1+A[1]*A[1], f = A[1]*A[2], g = 1+A[2]*A[2],
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
  Vector_3 Y, Z;
  FT normXu = Lsqrt( Xu*Xu );
  Y = Xu / normXu;
  FT XudotXv = Xu * Xv;
  Z = Xv - XudotXv * Xu / (normXu*normXu);
  FT normZ = Lsqrt( Z*Z );
  Z = Z / normZ;
  Matrix change_XuXv2YZ(2,2,0.);
  change_XuXv2YZ(0,0) = 1 / normXu;
  change_XuXv2YZ(0,1) = -XudotXv / (normXu * normXu * normZ);
  change_XuXv2YZ(1,0) = 0;
  change_XuXv2YZ(1,1) = 1 / normZ;
  FT det = 0.;
  Matrix inv = CGAL::Linear_algebraCd<FT>::inverse ( change_XuXv2YZ, det );
  //in the new orthonormal basis (Y,Z) of the tangent plane :
  weingarten = inv *(1/det) * weingarten * change_XuXv2YZ;

  // diagonalization of weingarten
  std::array<FT,3> W = {{ weingarten(0,0), weingarten(1,0), weingarten(1,1) }};
  std::array<FT,2> eval = {{ 0., 0. }};
  std::array<FT,4> evec = {{ 0., 0., 0., 0. }};

  //eval in increasing order
  CGAL::Default_diagonalize_traits<FT,2>::diagonalize_selfadjoint_covariance_matrix
    (W, eval, evec);

  Vector_3 d_max = evec[2]*Y + evec[3]*Z,
    d_min = evec[0]*Y + evec[1]*Z;

  switch_to_direct_orientation(d_max, d_min, normal);
  Aff_transformation change_basis (d_max[0], d_max[1], d_max[2],
                                   d_min[0], d_min[1], d_min[2],
                                   normal[0], normal[1], normal[2]);
  this->change_fitting2monge = change_basis;

  //store the monge basis origin and vectors with their world coord
  //store ppal curv
  monge_form.origin() = L2D_converter(
    (this->translate_p0.inverse() *
     this->change_world2fitting.inverse()) (orig_monge ));
  monge_form.maximal_principal_direction() = L2D_converter(this->change_world2fitting.inverse()(d_max));
  monge_form.minimal_principal_direction() = L2D_converter(this->change_world2fitting.inverse()(d_min));
  monge_form.normal_direction()  = L2D_converter(this->change_world2fitting.inverse()(normal));
  monge_form.coefficients()[0] = L2D_NTconverter()(eval[1]);
  monge_form.coefficients()[1] = L2D_NTconverter()(eval[0]);
  }
  //end else
}

template < class DataKernel, class LocalKernel, class SvdTraits>
void Monge_via_jet_fitting<DataKernel, LocalKernel, SvdTraits>::
compute_Monge_coefficients(FT* A, std::size_t dprime,
                           Monge_form& monge_form)
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
  FT p[3][3];
  p[0][0] = this->change_fitting2monge.m(0,0);
  p[1][0] = this->change_fitting2monge.m(0,1);
  p[2][0] = this->change_fitting2monge.m(0,2);
  p[0][1] = this->change_fitting2monge.m(1,0);
  p[1][1] = this->change_fitting2monge.m(1,1);
  p[2][1] = this->change_fitting2monge.m(1,2);
  p[0][2] = this->change_fitting2monge.m(2,0);
  p[1][2] = this->change_fitting2monge.m(2,1);
  p[2][2] = this->change_fitting2monge.m(2,2);

  // formula are designed for w=sum( Aij ui vj), but we have J_A = sum( Aij/i!j! ui vj)
  for (int k=0; k <= this->deg; k++) for (int i=0; i<=k; i++)
    A[k*(k+1)/2+i] /= fact(k-i)*fact(i);//this is A(k-i;i)

/*   //debug */
/*   std::cout << "coeff of A" << std::endl */
/*             << A[0] << " "<< A[1] << " "<< A[2] << std::endl */
/*             << A[3] << " "<< A[4] << " "<< A[5] << std::endl */
/*             << A[6] << " "<< A[7] << " "<< A[8] << " "<< A[9]<< std::endl */
/*             << A[10] << " "<< A[11] << " "<< A[12] << " "<< A[13]<< " " << A[14] << std::endl; */



  //     note f1 = f2 = f12 = 0
  //     FT f1 = A[1] * p[0][0] + A[2] * p[1][0] - p[2][0];
  //     FT f2 = A[2] * p[1][1] + A[1] * p[0][1] - p[2][1];
  //     FT f12 =
  //     2 * A[3] * p[0][0] * p[0][1]
  //     + 2 * A[5] * p[1][0] * p[1][1]
  //     + A[4] * p[0][1] * p[1][0]
  //     + A[4] * p[0][0] * p[1][1];
  //         -f11 / f3 = kmax
  //         -f22 / f3 = kmin

  FT f3 = A[1] * p[0][2] + A[2] * p[1][2] - p[2][2];
  FT f11 =
    2 * A[4] * p[0][0] * p[1][0]
    + 2 * A[5] * p[1][0] * p[1][0]
    + 2 * A[3] * p[0][0] * p[0][0];
  FT f13 =
    A[4] * p[0][0] * p[1][2]
    + A[4] * p[0][2] * p[1][0]
    + 2 * A[5] * p[1][0] * p[1][2]
    + 2 * A[3] * p[0][0] * p[0][2];
  FT f22 =
    2 * A[4] * p[0][1] * p[1][1]
    + 2 * A[5] * p[1][1] * p[1][1]
    + 2 * A[3] * p[0][1] * p[0][1];
  FT f23 =
    A[4] * p[0][1] * p[1][2]
    + 2 * A[5] * p[1][1] * p[1][2]
    + A[4] * p[0][2] * p[1][1]
    + 2 * A[3] * p[0][1] * p[0][2];
  FT f33 =
    2 * A[5] * p[1][2] * p[1][2]
    + 2 * A[3] * p[0][2] * p[0][2]
    + 2 * A[4] * p[0][2] * p[1][2];
  FT f111 =
    6 * A[8] * p[0][0] * p[1][0] * p[1][0]
    + 6 * A[7] * p[0][0] * p[0][0] * p[1][0]
    + 6 * A[6] * p[0][0] * p[0][0] * p[0][0]
    + 6 * A[9] * p[1][0] * p[1][0] * p[1][0];
  FT f222 =
    6 * A[7] * p[0][1] * p[0][1] * p[1][1]
    + 6 * A[8] * p[0][1] * p[1][1] * p[1][1]
    + 6 * A[9] * p[1][1] * p[1][1] * p[1][1]
    + 6 * A[6] * p[0][1] * p[0][1] * p[0][1];
  FT f112 =
    2 * A[7] * p[0][0] * p[0][0] * p[1][1]
    + 6 * A[6] * p[0][0] * p[0][0] * p[0][1]
    + 2 * A[8] * p[0][1] * p[1][0] * p[1][0]
    + 4 * A[8] * p[0][0] * p[1][0] * p[1][1]
    + 6 * A[9] * p[1][0] * p[1][0] * p[1][1]
    + 4 * A[7] * p[0][0] * p[0][1] * p[1][0];
  FT f122 =
    4 * A[8] * p[0][1] * p[1][0] * p[1][1]
    + 2 * A[8] * p[0][0] * p[1][1] * p[1][1]
    + 6 * A[6] * p[0][0] * p[0][1] * p[0][1]
    + 2 * A[7] * p[0][1] * p[0][1] * p[1][0]
    + 4 * A[7] * p[0][0] * p[0][1] * p[1][1]
    + 6 * A[9] * p[1][0] * p[1][1] * p[1][1];
  FT f113 =
    6*A[6]*p[0][0]*p[0][0]*p[0][2]
    +6*A[9]*p[1][0]*p[1][0]*p[1][2]
    +2*A[7]*p[0][0]*p[0][0]*p[1][2]
    +2*A[8]*p[0][2]*p[1][0]*p[1][0]
    +4*A[7]*p[0][0]*p[0][2]*p[1][0]
    +4*A[8]*p[0][0]*p[1][0]*p[1][2];
  FT f223 =
    2*A[8]*p[0][2]*p[1][1]*p[1][1]
    +6*A[6]*p[0][1]*p[0][1]*p[0][2]
    +6*A[9]*p[1][1]*p[1][1]*p[1][2]
    +2*A[7]*p[0][1]*p[0][1]*p[1][2]
    +4*A[7]*p[0][1]*p[0][2]*p[1][1]
    +4*A[8]*p[0][1]*p[1][1]*p[1][2];
  FT f123 =
    2*A[8]*p[0][2]*p[1][0]*p[1][1]
    +2*A[7]*p[0][0]*p[0][1]*p[1][2]
    +2*A[7]*p[0][0]*p[0][2]*p[1][1]
    +6*A[9]*p[1][0]*p[1][1]*p[1][2]
    +2*A[7]*p[0][1]*p[0][2]*p[1][0]
    +6*A[6]*p[0][0]*p[0][1]*p[0][2]
    +2*A[8]*p[0][0]*p[1][1]*p[1][2]
    +2*A[8]*p[0][1]*p[1][0]*p[1][2];

  FT b0 = 1/(f3*f3)*(-f111*f3+3*f13*f11);
  FT b1 = 1/(f3*f3)*(-f112*f3+f23*f11);
  FT b2 = 1/(f3*f3)*(-f122*f3+f13*f22);
  FT b3 = -1/(f3*f3)*(f222*f3-3*f23*f22);

  monge_form.coefficients()[2] = L2D_NTconverter()(b0);
  monge_form.coefficients()[3] = L2D_NTconverter()(b1);
  monge_form.coefficients()[4] = L2D_NTconverter()(b2);
  monge_form.coefficients()[5] = L2D_NTconverter()(b3);

  if ( dprime == 4 )
    {
      FT f1111 =
        24*A[13]*p[0][0]*p[1][0]*p[1][0]*p[1][0]
        +24*A[12]*p[0][0]*p[0][0]*p[1][0]*p[1][0]
        +24*A[11]*p[0][0]*p[0][0]*p[0][0]*p[1][0]
        +24*A[14]*p[1][0]*p[1][0]*p[1][0]*p[1][0]
        +24*A[10]*p[0][0]*p[0][0]*p[0][0]*p[0][0];
      FT f1112 =
        6*A[13]*p[0][1]*p[1][0]*p[1][0]*p[1][0]
        +18*A[13]*p[0][0]*p[1][0]*p[1][0]*p[1][1]
        +24*A[10]*p[0][0]*p[0][0]*p[0][0]*p[0][1]
        +12*A[12]*p[0][0]*p[0][1]*p[1][0]*p[1][0]
        +18*A[11]*p[0][0]*p[0][0]*p[0][1]*p[1][0]
        +24*A[14]*p[1][0]*p[1][0]*p[1][0]*p[1][1]
        +6*A[11]*p[0][0]*p[0][0]*p[0][0]*p[1][1]
        +12*A[12]*p[0][0]*p[0][0]*p[1][0]*p[1][1];
      FT f1122 =
        12*A[11]*p[0][0]*p[0][0]*p[0][1]*p[1][1]
        +12*A[13]*p[0][0]*p[1][0]*p[1][1]*p[1][1]
        +12*A[13]*p[0][1]*p[1][0]*p[1][0]*p[1][1]
        +16*A[12]*p[0][0]*p[0][1]*p[1][0]*p[1][1]
        +12*A[11]*p[0][0]*p[0][1]*p[0][1]*p[1][0]
        +24*A[10]*p[0][0]*p[0][0]*p[0][1]*p[0][1]
        +4*A[12]*p[0][1]*p[0][1]*p[1][0]*p[1][0]
        +4*A[12]*p[0][0]*p[0][0]*p[1][1]*p[1][1]
        +24*A[14]*p[1][0]*p[1][0]*p[1][1]*p[1][1];
      FT f1222 =
        6*A[13]*p[0][0]*p[1][1]*p[1][1]*p[1][1]
        +24*A[10]*p[0][0]*p[0][1]*p[0][1]*p[0][1]
        +24*A[14]*p[1][0]*p[1][1]*p[1][1]*p[1][1]
        +6*A[11]*p[0][1]*p[0][1]*p[0][1]*p[1][0]
        +18*A[11]*p[0][0]*p[0][1]*p[0][1]*p[1][1]
        +12*A[12]*p[0][0]*p[0][1]*p[1][1]*p[1][1]
        +12*A[12]*p[0][1]*p[0][1]*p[1][0]*p[1][1]
        +18*A[13]*p[0][1]*p[1][0]*p[1][1]*p[1][1];
      FT f2222 =
        24*A[13]*p[0][1]*p[1][1]*p[1][1]*p[1][1]
        +24*A[11]*p[0][1]*p[0][1]*p[0][1]*p[1][1]
        +24*A[12]*p[0][1]*p[0][1]*p[1][1]*p[1][1]
        +24*A[10]*p[0][1]*p[0][1]*p[0][1]*p[0][1]
        +24*A[14]*p[1][1]*p[1][1]*p[1][1]*p[1][1];

      FT c0 =
        -1/(f3*f3*f3)*(f1111*(f3*f3)-4*f13*f3*f111+12*f13*f13*f11-6*f113*f3*f11+3*f33*f11*f11);
      FT c1 =
        1/(f3*f3*f3)*(f23*f3*f111+3*f3*f123*f11+3*f13*f3*f112-f1112*(f3*f3)-6*f13*f23*f11);
      FT c2 =
        1/(f3*f3*f3)*(-f33*f22*f11+f113*f3*f22+2*f13*f3*f122-2*f13*f13*f22+f223*f3*f11+2*f23*f3*f112-2*f23*f23*f11-f1122*(f3*f3));
      FT c3 =
        1/(f3*f3*f3)*(-f1222*(f3*f3)-6*f13*f23*f22+3*f123*f3*f22+f13*f3*f222+3*f23*f3*f122);
      FT c4 =
        -1/(f3*f3*f3)*(f2222*(f3*f3)+3*f33*f22*f22-6*f223*f3*f22-4*f23*f3*f222+12*f23*f23*f22) ;

      monge_form.coefficients()[6] = L2D_NTconverter()(c0);
      monge_form.coefficients()[7] = L2D_NTconverter()(c1);
      monge_form.coefficients()[8] = L2D_NTconverter()(c2);
      monge_form.coefficients()[9] = L2D_NTconverter()(c3);
      monge_form.coefficients()[10] = L2D_NTconverter()(c4);
    }
}

template < class DataKernel, class LocalKernel, class SvdTraits>
void Monge_via_jet_fitting<DataKernel, LocalKernel, SvdTraits>::
switch_to_direct_orientation(Vector_3& v1, const Vector_3& v2,
                            const Vector_3& v3)
{
  if (CGAL::orientation (v1, v2, v3) == CGAL::NEGATIVE)
    v1 = -v1;
}


// template < class DataKernel, class LocalKernel, class SvdTraits>
// inline
// std::ostream&
// operator<<(std::ostream& out_stream,
//           const typename Monge_via_jet_fitting<DataKernel, LocalKernel, SvdTraits>::Monge_form& monge)
// {
//   monge.dump_verbose(out_stream);
//   return out_stream;
// }

} //namespace CGAL

#endif //CGAL_MONGE_VIA_JET_FITTING_H_
