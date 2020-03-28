// This file has been adopted from CGAL, to integrate the below mentioned paper
//
// Paper         : Learning to Reconstruct Symmetric Shapes using Planar Parameterization of 3D Surface
// Author(s)     : Hardik Jain, Manuel Wöllhaf, Olaf Hellwich
// Conference    : IEEE International Conference on Computer Vision Workshops (ICCVW) 2019
//

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_ITERATIVE_AUTHALIC_PARAMETERIZER_3_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_ITERATIVE_AUTHALIC_PARAMETERIZER_3_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Surface_mesh_parameterization/internal/angles.h>
#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>
#include <CGAL/Surface_mesh_parameterization/Error_code.h>

#include <CGAL/Surface_mesh_parameterization/Fixed_border_iterative_parameterizer_3.h>

#include <CGAL/Default.h>

#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>
#endif

#include <iostream>
#include <iomanip>
#include <math.h>
#include <boost/foreach.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/unordered_set.hpp>
#include <CGAL/squared_distance_2.h> //for 2D functions
#include <CGAL/squared_distance_3.h> //for 3D functions

/// \file Iterative_authalic_parameterizer_3.h
namespace CGAL {

namespace Surface_mesh_parameterization {

/// \ingroup  PkgSurfaceMeshParameterizationMethods
///
/// The class `Iterative_authalic_parameterizer_3`
/// implements the *Iterative Parameterization* algorithm.
///
/// A one-to-one mapping is guaranteed if the surface's border is mapped onto a convex polygon.
///
/// This class is a strategy called by the main
/// parameterization algorithm `Fixed_border_iterative_parameterizer_3::parameterize()` and it:
/// - provides the template parameters `BorderParameterizer_` and `SolverTraits_`.
/// - implements `compute_w_ij()` to compute w_ij = (i, j), coefficient of the matrix A
///   for j neighbor vertex of i, based on Iterative Parameterization algorithm.
///
/// \cgalModels `Parameterizer_3`
///
/// \tparam TriangleMesh_ must be a model of `FaceGraph`.
///
/// \tparam BorderParameterizer_ is a Strategy to parameterize the surface border
///         and must be a model of `Parameterizer_3`.<br>
///         <b>%Default:</b>
/// \code
///   Square_border_arc_length_parameterizer_3<TriangleMesh_>
/// \endcode
///
/// \tparam SolverTraits_ must be a model of `SparseLinearAlgebraTraits_d`.<br>
///         Note that the system is *not* symmetric because `Fixed_border_iterative_parameterizer_3`
///         does not remove border vertices from the system.<br>
///         <b>%Default:</b> If \ref thirdpartyEigen "Eigen" 3.1 (or greater) is available
///         and `CGAL_EIGEN3_ENABLED` is defined, then an overload of `Eigen_solver_traits`
///         is provided as default parameter:
/// \code
///   CGAL::Eigen_solver_traits<
///           Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType,
///                           Eigen::IncompleteLUT< double > > >
/// \endcode
///
/// \sa `CGAL::Surface_mesh_parameterization::Fixed_border_iterative_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
///
template < class TriangleMesh_,
class BorderParameterizer_ = Default,
class SolverTraits_ = Default>
class Iterative_authalic_parameterizer_3
    : public Fixed_border_iterative_parameterizer_3<
      TriangleMesh_,
      typename Default::Get<
      BorderParameterizer_,
      Square_border_arc_length_parameterizer_3<TriangleMesh_> >::type,
      typename Default::Get<
      SolverTraits_,
#if defined(CGAL_EIGEN3_ENABLED)
      CGAL::Eigen_solver_traits<
      Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType,
      Eigen::IncompleteLUT<double> > > >::type >
#else
#pragma message("Error: You must either provide 'SolverTraits_' or link CGAL with the Eigen library")
SolverTraits_>::type > // no parameter provided, and Eigen is not enabled: don't compile
#endif
{
public:
#ifndef DOXYGEN_RUNNING
  typedef typename Default::Get<
      BorderParameterizer_,
      Square_border_arc_length_parameterizer_3<TriangleMesh_> >::type  Border_parameterizer;

  typedef typename Default::Get<
      SolverTraits_,
#if defined(CGAL_EIGEN3_ENABLED)
      CGAL::Eigen_solver_traits<
      Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType,
      Eigen::IncompleteLUT<double> > >
#else
  SolverTraits_ // no parameter provided, and Eigen is not enabled: so don't compile!
#endif
  >::type                                                     Solver_traits;
#else
  typedef Border_parameterizer_                               Border_parameterizer;
  typedef SolverTraits_                                       Solver_traits;
#endif

  typedef TriangleMesh_                                       TriangleMesh;

  // Private types
private:
  // Superclass
  typedef Fixed_border_iterative_parameterizer_3<TriangleMesh,
      Border_parameterizer,
      Solver_traits>         Base;

  // Private types
private:
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  typedef CGAL::Vertex_around_target_circulator<TriangleMesh> vertex_around_target_circulator;
  typedef CGAL::Face_around_target_circulator<TriangleMesh> face_around_target_circulator;
  typedef CGAL::Vertex_around_face_circulator<TriangleMesh> vertex_around_face_circulator;

  // Traits subtypes:
  typedef typename Base::PPM                                   PPM;
  typedef typename Base::Kernel                                Kernel;
  typedef typename Base::NT                                    NT;
  typedef typename Base::Point_3                               Point_3;
  typedef typename Base::Point_2                               Point_2;
  typedef typename Base::Vector_3                              Vector_3;
  typedef boost::unordered_set<vertex_descriptor>              Vertex_set;

  // Solver traits subtypes:
  typedef typename Solver_traits::Vector                       Vector;
  typedef typename Solver_traits::Matrix                       Matrix;

  typedef CGAL::dynamic_face_property_t<double>                                     Face_double_tag;
  typedef typename boost::property_map<TriangleMesh, Face_double_tag>::type         Face_Double_map;
  typedef CGAL::dynamic_vertex_property_t<double>                                   Vertex_double_tag;
  typedef typename boost::property_map<TriangleMesh, Vertex_double_tag>::type       Vertex_Double_map;
  typedef CGAL::dynamic_vertex_property_t<int>                                      Vertex_int_tag;
  typedef typename boost::property_map<TriangleMesh, Vertex_int_tag>::type          Vertex_Int_map;
  typedef CGAL::dynamic_vertex_property_t<Point_2>                                  Vertex_point2_tag;
  typedef typename boost::property_map<TriangleMesh, Vertex_point2_tag>::type       Vertex_point2_map;

  // Public operations
public:
  /// Constructor
  Iterative_authalic_parameterizer_3(Border_parameterizer border_param = Border_parameterizer(),
      ///< %Object that maps the surface's border to 2D space.
      Solver_traits sparse_la = Solver_traits())
///< Traits object to access a sparse linear system.
: Fixed_border_iterative_parameterizer_3<TriangleMesh,
  Border_parameterizer,
  Solver_traits>(border_param, sparse_la)
  { }

  // Default copy constructor and operator =() are fine

  // Protected operations
protected:

  virtual NT compute_faceArea(TriangleMesh& mesh) {
    areaMap = get(Face_double_tag(), mesh);
    BOOST_FOREACH(face_descriptor fd, faces(mesh))
    put(areaMap, fd, Polygon_mesh_processing::face_area(fd,mesh));
    return 1.0;
  }

  virtual NT compute_borderLength_3D(TriangleMesh& mesh)  {
    return 1.0;
  }

  virtual NT compute_faceWise_L2(TriangleMesh& mesh, Vertex_point2_map &uvmap) {
    fL2Map = get(Face_double_tag(), mesh);
    BOOST_FOREACH(face_descriptor fd, faces(mesh))  {
      Point_2 uv_points[3];
      Point_3 mesh_points[3];
      int i = 0;
      BOOST_FOREACH(vertex_descriptor vd, CGAL::vertices_around_face(halfedge(fd,mesh),mesh)) {
        uv_points[i] = get(uvmap, vd);
        mesh_points[i] = mesh.point(vd);
        i++;
      }
      put(fL2Map, fd, getfL2(mesh_points, uv_points));
    }
    return 1.0;
  }

  virtual NT compute_vertexWise_L2(TriangleMesh& mesh, Vertex_set& vertices) {
    // update weights for vertices
    vL2Map = get(Vertex_double_tag(), mesh);
    BOOST_FOREACH(vertex_descriptor v, vertices)
    put(vL2Map, v, getvL2(mesh, v));
    return 1.0;
  }

  virtual Error_code setup_inner_vertex_relations(Matrix& A,
      Matrix& A_prev,
      Vector&,
      Vector&,
      const TriangleMesh& mesh,
      vertex_descriptor vertex,
      Vertex_Int_map& vimap)
  {

    // circulate over vertices around 'vertex' to compute w_ii and w_ijs
    vertex_around_target_circulator v_j(halfedge(vertex, mesh), mesh), end_v_j = v_j;
    std::vector<NeighborList> NeighborList_;
    int neighborsCounter = 0;
    double theta_sum = 0.0;
    // create NeighborList vector with vertex and vector
    CGAL_For_all(v_j, end_v_j)  {
      NeighborList NL;
      NL.vertex = *v_j;
      NL.vector = Vector_3(mesh.point(vertex),mesh.point(*v_j));
      NL.length = sqrt(NL.vector.squared_length());
      NeighborList_.push_back(NL);
      neighborsCounter++;
    }

    if (neighborsCounter < 2)
      return ERROR_NON_TRIANGULAR_MESH;

    if(neighborsCounter==2 && mesh.is_border(vertex)) {
      std::cout << "Encountered inner border with valency-2 vertex (" << vertex << "), initializing with Tutte weights, this can affect optimization" << std::endl;
      // Tutte weights
      for(int k=0; k<neighborsCounter; k++)
        NeighborList_[k].weight = 1.0;
    }
    else  {
      for(int n=0; n<neighborsCounter; n++) {
        int n_prev = (n==0 ? neighborsCounter-1 : n-1);
        double theta = angle(NeighborList_[n].vector,NeighborList_[n_prev].vector);
        NeighborList_[n].angle = theta;
        theta_sum += theta;
      }

      // Normalise the angle
      double factor = 2.0  / theta_sum;
      factor *= M_PI;
      for(int n=0; n<neighborsCounter; n++)
        NeighborList_[n].angle *= factor;

      NeighborList_[0].angle = 0.0;
      for(int n=1; n<neighborsCounter; n++)
        NeighborList_[n].angle += NeighborList_[n-1].angle;
      for(int n=0; n<neighborsCounter; n++)
        NeighborList_[n].uv = Point_2(NeighborList_[n].length*cos(NeighborList_[n].angle), NeighborList_[n].length*sin(NeighborList_[n].angle));

      for(int j=0; j<neighborsCounter; j++)
      {
        /* Given the j-th neighbour of node i,
            find the two neighbours by intersecting the
            line through nodes i and j with all segments of the polygon
            made by the neighbours. Take the two neighbours on
            either side. Only one segment intersects this line. */
        for(int k=0; k<neighborsCounter; k++)
        {
          int kk = (k == neighborsCounter-1 ? 0 : k+1);
          if(k == j || kk == j) continue;

          double cross1 = determinant(NeighborList_[j].uv,NeighborList_[k].uv);
          double cross2 = determinant(NeighborList_[j].uv,NeighborList_[kk].uv);

          if(cross1 * cross2 <= 0.0)
          {
            double tau0,tau1,tau2;
            baryCoords0(NeighborList_[j].uv, NeighborList_[k].uv, NeighborList_[kk].uv, tau0, tau1, tau2);
            NeighborList_[j].weight += tau0;
            NeighborList_[k].weight  += tau1;
            NeighborList_[kk].weight += tau2;
            break;
          }
        }
      }
    }
    // Scale the weights so that they sum to 1.
    //  double ratio = 1.0 / (double)n;
    double sum = 0;
    for (int j = 0; j < neighborsCounter; ++j)
      sum +=NeighborList_[j].weight;
    double ratio = 1.0 / sum;
    for(int j=0; j<neighborsCounter; j++)
      NeighborList_[j].weight *= ratio;


    // assign these weights to the edge pairs now
    NT w_ii = 0;
    const int i = get(vimap,vertex);
    for(int n=0; n<neighborsCounter; n++) {
      NT w_ij = -1.0 * NeighborList_[n].weight;
      if(w_ij > 0)
        w_ij *= -1.0;
      w_ii -= w_ij;

      // Get j index
      const int j = get(vimap, NeighborList_[n].vertex);

      // Set w_ij in matrix
      A.set_coef(i,j, w_ij, true /*new*/);
      A_prev.set_coef(i,j, w_ij, true);
    }
    // Set w_ii in matrix
    A.set_coef(i,i, w_ii, true /*new*/);
    return OK;

  }

  virtual Error_code setup_inner_vertex_relations_cotangent(Matrix& A,
      Vector&,
      Vector&,
      const TriangleMesh& mesh,
      vertex_descriptor vertex,
      Vertex_Int_map& vimap)
  {
    int i = get(vimap,vertex);

    // circulate over vertices around 'vertex' to compute w_ii and w_ijs
    NT w_ii = 0;
    int vertexIndex = 0;
    int neighborsCounter = 0;
    vertex_around_target_circulator v_j(halfedge(vertex, mesh), mesh), end = v_j;

    CGAL_For_all(v_j, end){
      neighborsCounter++;
    }

    bool bcompute_w_ij = true;
    if (neighborsCounter < 2)
      return ERROR_NON_TRIANGULAR_MESH;
    else if(neighborsCounter==2 && mesh.is_border(vertex))  {
      bcompute_w_ij = false;
      std::cout << "Encountered inner border with valency-2 vertex (" << vertex << "), initializing with Tutte weights, this can affect optimization" << std::endl;
    }

    CGAL_For_all(v_j, end){
      // Call to virtual method to do the actual coefficient computation
      NT w_ij = -1.0;
      if (bcompute_w_ij)
        w_ij *= compute_w_ij(mesh, vertex, v_j);
      // w_ii = - sum of w_ijs
      w_ii -= w_ij;

      // Get j index
      int j = get(vimap, *v_j);

      // Set w_ij in matrix
      A.set_coef(i,j, w_ij, true /*new*/);
      vertexIndex++;
    }


    // Set w_ii in matrix
    A.set_coef(i,i, w_ii, true /*new*/);
    return OK;
  }

  virtual double compute_sig_ij(TriangleMesh& mesh, Vertex_point2_map &uvmap, vertex_descriptor v_i, vertex_descriptor v_j, double gamma) {
    double out = (pow(get(vL2Map,v_i),gamma)+pow(get(vL2Map,v_j),gamma))/2.0;
    if(out <= 0.0)
      std::cout << "compute_sig_ij <= 0.0" << std::endl;
    return out;
  }

  virtual double distError(TriangleMesh& mesh)  {
    double varphi = 0.0;
    double localArea = 0.0;
    BOOST_FOREACH(face_descriptor fd, faces(mesh))  {
      varphi += get(fL2Map,fd)*get(areaMap,fd);
      localArea += get(areaMap,fd);
    }
    double err = sqrt(varphi/localArea);
    return err;
  }

private:
  TriangleMesh mesh;
  Vertex_Double_map vL2Map;
  Face_Double_map areaMap;
  Face_Double_map fL2Map;

  double getfL2(Point_3 mesh_points[], Point_2 uv_points[]) const {
    const double A = getA(uv_points);
    Point_3 Ss = getSs(mesh_points, uv_points, A);
    Point_3 St = getSt(mesh_points, uv_points, A);

    double a = innerProduct(Ss,Ss);
    double c = innerProduct(St,St);
    return sqrt((a+c)/2.0);
  }

  double getvL2(const TriangleMesh& mesh, vertex_descriptor &vertex)  const {
    bool debug = false;
    if(vertex == 711111)
      debug = true;
    if(debug)
      std::cout << "\t" << vertex << "\t" << std::flush;
    halfedge_descriptor hf = halfedge(vertex, mesh);
    if(debug)
      std::cout << hf << "\t" << std::flush;
    face_around_target_circulator f_j(hf, mesh), end_f_j = f_j;
    double varphi = 0.0;
    double localArea = 0.0;
    int i=0;
    CGAL_For_all(f_j, end_f_j)  {
      if(debug)
        std::cout << *f_j << "\t" << std::flush;
      if(*f_j > mesh.number_of_faces())
        continue;
      varphi += get(fL2Map,*f_j)*get(areaMap,*f_j);
      localArea += get(areaMap,*f_j);
      i++;
    }

    if(debug) {
      std::cout << varphi << "\t" << std::flush;
      std::cout << localArea << "\t" << std::flush;
    }

    if(mesh.is_border(vertex) && mesh.degree(vertex) != i+1)
      std::cerr << std::endl;
    else if(!mesh.is_border(vertex) && mesh.degree(vertex) != i)
      std::cerr << std::endl;
    return sqrt(varphi/localArea);
  }

  /// Compute w_ij = (i, j), coefficient of matrix A for j neighbor vertex of i.
  // These weights are shape preserving weights proposed in Floater1997
  NT compute_w_ij(const TriangleMesh& mesh,
      vertex_descriptor main_vertex_v_i,
      vertex_around_target_circulator neighbor_vertex_v_j) const
  {
    const PPM ppmap = get(vertex_point, mesh);

    const Point_3& position_v_i = get(ppmap, main_vertex_v_i);
    const Point_3& position_v_j = get(ppmap, *neighbor_vertex_v_j);

    // Compute the square norm of v_j -> v_i vector
    Vector_3 edge = position_v_i - position_v_j;
    double square_len = edge*edge;

    // Compute cotangent of (v_k,v_j,v_i) corner (i.e. cotan of v_j corner)
    // if v_k is the vertex before v_j when circulating around v_i
    vertex_around_target_circulator previous_vertex_v_k = neighbor_vertex_v_j;
    previous_vertex_v_k--;
    const Point_3& position_v_k = get(ppmap, *previous_vertex_v_k);
    NT cotg_psi_ij = internal::cotangent<Kernel>(position_v_k, position_v_j, position_v_i);
    NT cotg_beta_ij = internal::cotangent<Kernel>(position_v_i, position_v_k, position_v_j);

    // Compute cotangent of (v_i,v_j,v_l) corner (i.e. cotan of v_j corner)
    // if v_l is the vertex after v_j when circulating around v_i
    vertex_around_target_circulator next_vertex_v_l = neighbor_vertex_v_j;
    next_vertex_v_l++;
    const Point_3& position_v_l = get(ppmap,*next_vertex_v_l);
    NT cotg_theta_ij = internal::cotangent<Kernel>(position_v_i, position_v_j, position_v_l);
    NT cotg_alpha_ij = internal::cotangent<Kernel>(position_v_j, position_v_l, position_v_i);

    NT weight = 0.0;
    CGAL_assertion(square_len != 0.0); // two points are identical!
    if(square_len != 0.0) {
      weight = cotg_beta_ij + cotg_alpha_ij;
    }
    return weight;
  }

  struct NeighborList  {
    vertex_descriptor vertex;
    double angle;
    double length;
    Vector_3 vector;
    Point_2 uv;
    double weight;
  };

  void baryCoords0(Point_2& uv0, Point_2& uv1, Point_2& uv2, double& tau0, double& tau1, double& tau2)  {
    double det0 = determinant(uv1,uv2);
    double det1 = determinant(uv2,uv0);
    double det2 = determinant(uv0,uv1);
    double det3 = CGAL::determinant(Vector_2<Kernel>(uv1.x()-uv0.x(), uv1.y()-uv0.y()), Vector_2<Kernel>(uv2.x()-uv0.x(), uv2.y()-uv0.y()));
    if(det3 <= 0.0)
      det3 = 1.0;

    tau0 = det0 / det3;
    tau1 = det1 / det3;
    tau2 = det2 / det3;
  }

  double determinant(Point_2& v0, Point_2& v1)  {
    return (v0.x()*v1.y() - v1.x()*v0.y());
  }

  double angle(Vector_3& v0, Vector_3& v1)  {
    return std::acos(v0*v1/(CGAL::sqrt(v0*v0)*CGAL::sqrt(v1*v1)));
  }

  double getA(Point_2 uv_points[]) const {
    double A = (((uv_points[1].x()-uv_points[0].x())*(uv_points[2].y()-uv_points[0].y()))-((uv_points[2].x()-uv_points[0].x())*(uv_points[1].y()-uv_points[0].y())))/2;
    if (A == 0.0)
      return 1.0;
    return A;
  }

  Point_3 getSs(Point_3 mesh_points[3], Point_2 uv_points[3],const double &A) const {
    double dt0 = uv_points[1].y()-uv_points[2].y();
    double dt1 = uv_points[2].y()-uv_points[0].y();
    double dt2 = uv_points[0].y()-uv_points[1].y();
    Point_3 Ss (
        (mesh_points[0].x()*dt0 + mesh_points[1].x()*dt1 + mesh_points[2].x()*dt2 )/(2.0*A),
        (mesh_points[0].y()*dt0 + mesh_points[1].y()*dt1 + mesh_points[2].y()*dt2 )/(2.0*A),
        (mesh_points[0].z()*dt0 + mesh_points[1].z()*dt1 + mesh_points[2].z()*dt2 )/(2.0*A));
    return Ss;
  }

  Point_3 getSt(Point_3 mesh_points[3], Point_2 uv_points[3],const double &A) const  {
    double ds0 = uv_points[2].x()-uv_points[1].x();
    double ds1 = uv_points[0].x()-uv_points[2].x();
    double ds2 = uv_points[1].x()-uv_points[0].x();
    Point_3 St (
        (mesh_points[0].x()*ds0 + mesh_points[1].x()*ds1 +mesh_points[2].x()*ds2)/(2.0*A),
        (mesh_points[0].y()*ds0 + mesh_points[1].y()*ds1 +mesh_points[2].y()*ds2)/(2.0*A),
        (mesh_points[0].z()*ds0 + mesh_points[1].z()*ds1 +mesh_points[2].z()*ds2)/(2.0*A));
    return St;
  }

  double innerProduct(const Point_3& pointA, const Point_3& pointB) const {
    return ((pointA.x())*(pointB.x()) + (pointA.y())*(pointB.y()) + (pointA.z())*(pointB.z()));
  }

};

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_ITERATIVE_AUTHALIC_PARAMETERIZER_3_H
