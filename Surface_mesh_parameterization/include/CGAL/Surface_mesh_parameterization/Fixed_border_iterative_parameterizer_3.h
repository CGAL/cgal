// This file has been adopted from CGAL, to integrate the below mentioned paper
//
// Paper         : Learning to Reconstruct Symmetric Shapes using Planar Parameterization of 3D Surface
// Author(s)     : Hardik Jain, Manuel Wöllhaf, Olaf Hellwich
// Conference    : IEEE International Conference on Computer Vision Workshops (ICCVW) 2019
//

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_FIXED_BORDER_ITERATIVE_PARAMETERIZER_3_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_FIXED_BORDER_ITERATIVE_PARAMETERIZER_3_H

#define DEBUG_L0 1
#define DEBUG_L1 0

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Surface_mesh_parameterization/internal/angles.h>
#include <CGAL/Surface_mesh_parameterization/internal/Containers_filler.h>
#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>

#include <CGAL/Surface_mesh_parameterization/Error_code.h>
#include <CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <CGAL/boost/graph/properties.h>
#include <CGAL/circulator.h>
#include <CGAL/Default.h>

#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>
#endif

#include <iostream>
#include <iomanip>
#include <boost/foreach.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/unordered_set.hpp>
#include <CGAL/squared_distance_2.h> //for 2D functions
#include <CGAL/squared_distance_3.h> //for 3D functions

/// \file Fixed_border_iterative_parameterizer_3.h

namespace CGAL {

namespace Surface_mesh_parameterization {

// ------------------------------------------------------------------------------------
// Declaration
// ------------------------------------------------------------------------------------

/// \ingroup PkgSurfaceMeshParameterizationMethods
///
/// The class `Fixed_border_iterative_parameterizer_3`
/// is the base class of fixed border parameterization methods (Tutte, Floater, ...).
///
/// A one-to-one mapping is guaranteed if the border of the surface is mapped onto a convex polygon.
///
/// This class is a pure virtual class and thus cannot be instantiated.
/// Nevertheless, it implements most of the parameterization algorithm `parameterize()`.
/// Subclasses are *Strategies* \cgalCite{cgal:ghjv-dpero-95} that modify the behavior of this algorithm:
/// - They provide the template parameters `BorderParameterizer_` and `SolverTraits_`.
/// - They implement `compute_w_ij()` to compute w_ij = (i, j), coefficient of matrix A
///   for j neighbor vertex of i.
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
/// \sa `CGAL::Surface_mesh_parameterization::Iterative_authalic_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
///
template < typename TriangleMesh_,
class BorderParameterizer_ = Default,
class SolverTraits_ = Default>
class Fixed_border_iterative_parameterizer_3
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
#pragma message("Error: You must either provide 'SolverTraits_' or link CGAL with the Eigen library");
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
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  typedef CGAL::Vertex_around_target_circulator<TriangleMesh> vertex_around_target_circulator;
  typedef CGAL::Face_around_target_circulator<TriangleMesh> face_around_target_circulator;
  typedef CGAL::Vertex_around_face_circulator<TriangleMesh> vertex_around_face_circulator;

  // Protected types
protected:
  // Traits subtypes:
  typedef typename internal::Kernel_traits<TriangleMesh>::Kernel    Kernel;
  typedef typename internal::Kernel_traits<TriangleMesh>::PPM       PPM;
  typedef typename Kernel::FT                                       NT;
  typedef typename Kernel::Point_2                                  Point_2;
  typedef typename Kernel::Point_3                                  Point_3;
  typedef typename Kernel::Vector_3                                 Vector_3;

  // Solver traits subtypes:
  typedef typename Solver_traits::Vector                            Vector;
  typedef typename Solver_traits::Matrix                            Matrix;

  typedef boost::unordered_set<vertex_descriptor>                                   Vertex_set;
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
  Fixed_border_iterative_parameterizer_3(Border_parameterizer border_param = Border_parameterizer(),
      ///< %Object that maps the surface's border to 2D space
      Solver_traits sparse_la = Solver_traits())
///< Traits object to access a sparse linear system
: m_borderParameterizer(border_param), m_linearAlgebra(sparse_la), LScounter(0)
{ }

  /// Destructor of base class should be virtual.
  virtual ~Fixed_border_iterative_parameterizer_3() { }

  // Default copy constructor and operator =() are fine

  /// Compute a one-to-one mapping from a triangular 3D surface mesh
  /// to a piece of the 2D space.
  /// The mapping is piecewise linear (linear in each triangle).
  /// The result is the (u,v) pair image of each vertex of the 3D surface.
  ///
  /// \tparam VertexUVmap must be a model of `ReadWritePropertyMap` with
  ///         `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and
  ///         %Point_2 (type deduced from `TriangleMesh` using `Kernel_traits`)
  ///         as value type.
  /// \tparam VertexIndexMap must be a model of `ReadablePropertyMap` with
  ///         `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and
  ///         a unique integer as value type.
  /// \tparam VertexParameterizedMap must be a model of `ReadWritePropertyMap` with
  ///         `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and
  ///         a Boolean as value type.
  ///
  /// \param mesh a triangulated surface.
  /// \param bhd a halfedge descriptor on the boundary of `mesh`.
  /// \param uvmap an instanciation of the class `VertexUVmap`.
  /// \param vimap an instanciation of the class `VertexIndexMap`.
  /// \param vpmap an instanciation of the class `VertexParameterizedMap`.
  /// \param iterations an integer number of iterations to run the parameterization.
  /// \param error return error value of the iterative process.
  ///
  /// \pre `mesh` must be a triangular mesh.
  /// \pre The mesh border must be mapped onto a convex polygon.
  /// \pre The vertices must be indexed (`vimap` must be initialized)
  ///
  template <typename VertexUVmap,
  typename VertexIndexMap,
  typename VertexParameterizedMap>
  Error_code parameterize(TriangleMesh& mesh,
      halfedge_descriptor bhd,
      VertexUVmap uvmap1,
      VertexIndexMap vimap1,
      VertexParameterizedMap vpmap,
      int& iterations,
      double& error)
  {
    Error_code status = OK;

    Vertex_set vertices;

    internal::Containers_filler<TriangleMesh> fc(mesh, vertices);
    Polygon_mesh_processing::connected_component(
        face(opposite(bhd, mesh), mesh),
        mesh,
        boost::make_function_output_iterator(fc));

    int nbVertices = static_cast<int>(vertices.size());

    if (nbVertices == 0)
      return ERROR_EMPTY_MESH;

    // Compute (u,v) for border vertices and mark them as "parameterized"
    status = get_border_parameterizer().parameterize(mesh, bhd, uvmap1, vimap1, vpmap);
    if (status != OK)
      return status;

    // Create two sparse linear systems "A*Xu = Bu" and "A*Xv = Bv" (one line/column per vertex)
    Matrix A(nbVertices, nbVertices);
    Matrix A_prev(nbVertices, nbVertices);
    Vector Xu(nbVertices), Xv(nbVertices), Bu(nbVertices), Bv(nbVertices);
    double err[iterations];

    // Initialize A, Xu, Xv, Bu and Bv after border parameterization
    // Fill the border vertices' lines in both linear systems:
    // "u = constant" and "v = constant"
    //
    initialize_system_from_mesh_border(A, Bu, Bv, mesh, bhd, uvmap1, vimap1);

    if(DEBUG_L1) {
      printMatrix(vertices, vimap1, A, "A");
      printVector(vertices, vimap1, Bu, "Bu");
      printVector(vertices, vimap1, Bv, "Bv");
    }
    // Fill the matrix for the inner vertices v_i: compute A's coefficient
    // w_ij for each neighbor j; then w_ii = - sum of w_ijs
    boost::unordered_set<vertex_descriptor> main_border;
    BOOST_FOREACH(vertex_descriptor v, vertices_around_face(bhd,mesh))
    main_border.insert(v);

    // create new uvmap & vimap
    Vertex_point2_map uvmap = get(Vertex_point2_tag(),mesh);
    Vertex_Int_map vimap = get(Vertex_int_tag(),mesh);

    BOOST_FOREACH(vertex_descriptor v, vertices)  {
      put(uvmap, v, get(uvmap1, v));
      put(vimap, v, get(vimap1, v));
    }

    // update last best uvmap wrt the border
    lastBestuvmap = get(Vertex_point2_tag(), mesh);
    BOOST_FOREACH(vertex_descriptor v, vertices)  {
      // border vertices only
      if(main_border.find(v) != main_border.end())
        put(lastBestuvmap, v, get(uvmap, v));
    }

    int lastBesti = 0;
    compute_faceArea(mesh);
    compute_borderLength_3D(mesh);
    double gamma = 1;
    bool isChanged = false;
    // iterate it with the new weights
    if(DEBUG_L0)
      std::cout << std::endl;
    int i=0;
    while (i<=iterations) {
      if(DEBUG_L0)
        std::cout << "Iteration " << i << std::flush;
      // update weights for inner vertices
      BOOST_FOREACH(vertex_descriptor v, vertices)  {
        // inner vertices only
        if(main_border.find(v) == main_border.end())  {
          // Compute the line i of matrix A for i inner vertex
          if (i==0) {
            //status = setup_inner_vertex_relations(A, A_prev, Bu, Bv, mesh, v, vimap);
            status = setup_inner_vertex_relations_cotangent(A, Bu, Bv, mesh, v, vimap);
            if(status != OK)
              return status;
          }
          else  {
            status = setup_iter_inner_vertex_relations(A, A_prev, Bu, Bv, mesh, v, vimap, uvmap, gamma);
            if(status != OK)
              return status;
          }
        }
      }

      if(DEBUG_L1) {
        print(A, Xu, Bu);
        std::cout << std::endl;
        print(A, Xv, Bv);
      }

      // solve linear equations
      // Solve "A*Xu = Bu". On success, solution is (1/Du) * Xu.
      // Solve "A*Xv = Bv". On success, solution is (1/Dv) * Xv.
      NT Du = 0, Dv = 0;
      if(!get_linear_algebra_traits().linear_solver(A, Bu, Xu, Du) ||
          !get_linear_algebra_traits().linear_solver(A, Bv, Xv, Dv))
      {
        status = ERROR_CANNOT_SOLVE_LINEAR_SYSTEM;
      }
      else
        LScounter = 0;

      if(status != OK)  {
        if(LScounter <4)  {
          // modify the weights and re-try the linear solver
          LScounter++;
          gamma /= 2;
          continue;
        }
        else  {
          status = OK;
          break;
        }
      }

      // WARNING: this package does not support homogeneous coordinates!
      CGAL_assertion(Du == 1.0);
      CGAL_assertion(Dv == 1.0);

      // Copy A to A_prev, it is a computationally inefficient task but neccesary
      copySparseMatrix(mesh, vertices, vimap, A, A_prev);

      // Copy Xu and Xv coordinates into the (u,v) pair of each vertex
      BOOST_FOREACH(vertex_descriptor v, vertices)
      {
        // inner vertices only
        if(main_border.find(v) == main_border.end()){
          int index = get(vimap,v);
          put(uvmap,v,Point_2(Xu[index],Xv[index]));
          put(vpmap,v,true);
        }
      }

      compute_faceWise_L2(mesh, uvmap);
      compute_vertexWise_L2(mesh, vertices);

      err[i] = areaDist(mesh, vertices, main_border, uvmap);

      if(DEBUG_L0)
        std::cout << " err " << err[i] << std::flush;

      if(err[i] <= err[lastBesti]) {
        updateUVMAP(vertices, main_border, uvmap);
        lastBesti = i;
        isChanged = false;
        if(DEBUG_L0)
          std::cout << " *****" << std::flush;
      }
      else if (err[i]>100)  {
        break;
      }
      else  {
        if(!isChanged)  {
          gamma /= 2;
          isChanged = true;
        }
      }

      i++;
      std::cout << std::endl;
    }

    // Check postconditions
    if(status != OK)
      return status;
    if(i==0 && i!=iterations)
      // means that the computation terminated for the first iteration may be because system was unsolvable
      return ERROR_CANNOT_SOLVE_LINEAR_SYSTEM;

    BOOST_FOREACH(vertex_descriptor v, vertices)  {
      put(uvmap1, v, get(lastBestuvmap, v));
      put(vimap1, v, get(vimap, v));
    }
    iterations = lastBesti;
    error = err[lastBesti];
    return status;
  }

  void updateUVMAP(Vertex_set &vertices, boost::unordered_set<vertex_descriptor> &main_border, Vertex_point2_map uvmap) {
    BOOST_FOREACH(vertex_descriptor v, vertices)  {
      // inner vertices only
      if(main_border.find(v) == main_border.end())
        put(lastBestuvmap, v, get(uvmap, v));
    }
  }


  /// Initialize A, Bu and Bv after border parameterization.
  /// Fill the border vertices' lines in both linear systems:
  /// "u = constant" and "v = constant".
  ///
  /// \tparam VertexUVmap must be a model of `ReadWritePropertyMap` with
  ///         `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and
  ///         %Point_2 (type deduced from `TriangleMesh` using `Kernel_traits`)
  ///         as value type.
  /// \tparam VertexIndexMap must be a model of `ReadablePropertyMap` with
  ///         `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and
  ///         a unique integer as value type.
  ///
  /// \param A the matrix in both linear system
  /// \param Bu the right hand side vector in the linear system of x coordinates
  /// \param Bv the right hand side vector in the linear system of y coordinates
  /// \param mesh a triangulated surface.
  /// \param bhd a halfedge descriptor on the boundary of `mesh`.
  /// \param uvmap an instanciation of the class `VertexUVmap`.
  /// \param vimap an instanciation of the class `VertexIndexMap`.
  ///
  /// \pre Vertices must be indexed (`vimap` must be initialized).
  /// \pre A, Bu and Bv must be allocated.
  /// \pre Border vertices must be parameterized.
  template <typename VertexUVmap, typename VertexIndexMap>
  void initialize_system_from_mesh_border(Matrix& A, Vector& Bu, Vector& Bv,
      const TriangleMesh& mesh,
      halfedge_descriptor bhd,
      VertexUVmap uvmap,
      VertexIndexMap vimap) const
  {
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(bhd, mesh)){
      // Get vertex index in sparse linear system
      int index = get(vimap, target(hd, mesh));
      // Write a diagonal coefficient of A
      A.set_coef(index, index, 1, true /*new*/);
      // get the halfedge uv
      // Write constant in Bu and Bv
      const Point_2& uv = get(uvmap, target(hd, mesh));
      Bu[index] = uv.x();
      Bv[index] = uv.y();
    }
  }

  /// Compute w_ij, coefficient of matrix A for j neighbor vertex of i.
  /// Implementation note: Subclasses must at least implement compute_w_ij().
  ///
  /// \param mesh a triangulated surface.
  /// \param main_vertex_v_i the vertex of `mesh` with index `i`
  /// \param neighbor_vertex_v_j the vertex of `mesh` with index `j`
  virtual NT compute_w_ij(const TriangleMesh& mesh,
      vertex_descriptor main_vertex_v_i,
      vertex_around_target_circulator neighbor_vertex_v_j) const = 0;

  virtual Error_code setup_inner_vertex_relations(Matrix& A,
      Matrix& A_prev,
      Vector&,
      Vector&,
      const TriangleMesh& mesh,
      vertex_descriptor vertex,
      Vertex_Int_map& vimap)  {
    return OK; // TODO: Need to check
  }

  virtual Error_code setup_inner_vertex_relations_cotangent(Matrix& A,
      Vector&,
      Vector&,
      const TriangleMesh& mesh,
      vertex_descriptor vertex,
      Vertex_Int_map& vimap)  {
    return OK; // TODO: Need to check
  }

  template <typename VertexIndexMap>
  Error_code setup_iter_inner_vertex_relations(Matrix& A,
      Matrix& A_prev,
      Vector&,
      Vector&,
      TriangleMesh& mesh,
      vertex_descriptor vertex,
      VertexIndexMap vimap,
      Vertex_point2_map &uvmap,
      double gamma)
  {

    int i = get(vimap,vertex);
    // circulate over vertices around 'vertex' to compute w_ii and w_ijs
    NT w_ii = 0;
    int vertexIndex = 0;

    //const double sigma = getvL2(mesh, vertex, fL2Map, areaMap);
    vertex_around_target_circulator v_j(halfedge(vertex, mesh), mesh), end = v_j;
    CGAL_For_all(v_j, end){
      // Get j index
      int j = get(vimap, *v_j);
      // Call to virtual method to do the actual coefficient computation
      //NT w_ij = A_prev.get_coef(i,j) / compute_sig_ij(vertex, *v_j) / gamma;
      NT w_ij = A_prev.get_coef(i,j) / pow(compute_sig_ij(mesh, uvmap, vertex, *v_j, 1.0),gamma);

      // w_ii = - sum of w_ijs
      w_ii -= w_ij;

      // Set w_ij in matrix
      A.set_coef(i,j, w_ij, false);
      //A_prev.set_coef(i,j, w_ij, false);
      vertexIndex++;
    }

    if (vertexIndex < 2)
      return ERROR_NON_TRIANGULAR_MESH;

    // Set w_ii in matrix
    A.set_coef(i,i, w_ii, true /*new*/);
    return OK;
  }

  virtual NT compute_faceArea(TriangleMesh& mesh) = 0;

  virtual NT compute_faceWise_L2(TriangleMesh& mesh, Vertex_point2_map &uvmap) = 0;

  virtual NT compute_vertexWise_L2(TriangleMesh& mesh, Vertex_set& vertices) = 0;

  virtual double compute_sig_ij(TriangleMesh& mesh, Vertex_point2_map &uvmap, vertex_descriptor v_i, vertex_descriptor v_j, double gamma) = 0;

  virtual NT compute_borderLength_3D(TriangleMesh& mesh) = 0;

  virtual double distError(TriangleMesh& mesh) = 0;

  // Measure L2 stretch
  template <typename VertexUVmap>
  double distError1(const TriangleMesh& mesh, Vertex_set &vertices,
      boost::unordered_set<vertex_descriptor> &main_border,
      VertexUVmap &uvmap)  {

    //    double car
    BOOST_FOREACH(face_descriptor fd, faces(mesh))  {

    }

    // iterate fpr all inner vertices and for each vertex
    std::vector<double> area_3D;
    std::vector<double> area_2D;
    std::vector<double> area_dist;

    double A_3D = Polygon_mesh_processing::area(mesh);
    double A_2D = 1.;

    BOOST_FOREACH(vertex_descriptor v, vertices){
      // inner vertices only
      if(main_border.find(v) == main_border.end()){
        double a_2D = 0;
        double a_3D = 0;
        // find the area of all the adjacent faces to this vertex
        face_around_target_circulator f_j(halfedge(v, mesh), mesh), end = f_j;
        CGAL_For_all(f_j, end)  {
          // get area in original mesh
          a_3D += (Polygon_mesh_processing::face_area(*f_j, mesh)/A_3D);

          // get area in parameterised mesh
          // iterate for all the vertices of this face and compute area
          std::vector<Point_2> uv_points;
          BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(halfedge(v, mesh), mesh))  {
            uv_points.push_back(get(uvmap,vd));
          }
          a_2D += (abs(CGAL::area(uv_points[0],uv_points[1],uv_points[2]))/A_2D);
        }
        area_3D.push_back(a_3D);
        area_2D.push_back(a_2D);
        area_dist.push_back(square((a_3D/a_2D) - 1));
      }
    }
    return sqrt(sum_vector(area_dist));
  }

  template <typename VertexUVmap>
  double areaDist(TriangleMesh& mesh, Vertex_set &vertices,
      boost::unordered_set<vertex_descriptor> &main_border,
      VertexUVmap &uvmap)  {
    // iterate for all inner vertices and for each vertex

    Face_Double_map area_3DMap = get(Face_double_tag(), mesh);
    Face_Double_map area_2DMap = get(Face_double_tag(), mesh);;
    std::vector<double> area_dist;
    std::vector<face_descriptor> innerFaces;
    double A_3D = 0.0;
    double A_2D = 0.0;
    /*
    BOOST_FOREACH(face_descriptor fd, faces(mesh))  {
      if (!mesh.is_border(target(mesh,fd)))
        innerFaces.push_back(fd);
    }
     */  BOOST_FOREACH(face_descriptor fd, faces(mesh))  {
       put(area_3DMap, fd, Polygon_mesh_processing::face_area(fd, mesh));
       // get area in parameterised mesh
       std::vector<Point_2> uv_points;
       BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(halfedge(fd, mesh), mesh))  {
         uv_points.push_back(get(uvmap,vd));
       }
       put(area_2DMap, fd, abs(CGAL::area(uv_points[0],uv_points[1],uv_points[2])));
     }

     BOOST_FOREACH(face_descriptor fd, faces(mesh))  {
       A_3D += get(area_3DMap, fd);
       A_2D += get(area_2DMap, fd);
     }

     BOOST_FOREACH(face_descriptor fd, faces(mesh))  {
       double a_3D = get(area_3DMap, fd);
       double a_2D = get(area_2DMap, fd);
       area_dist.push_back(abs(a_3D/A_3D - a_2D/A_2D));
     }

     return sum_vector(area_dist);
  }

  template <typename T>
  T sum_vector(typename std::vector<T> vec)  {
    T sum = 0;
    for (typename std::vector<T>::iterator it=vec.begin(); it!= vec.end(); it++) {
      sum += *it;
    }
    return sum;
  }

  template <typename VertexIndexMap>
  void copySparseMatrix(TriangleMesh& mesh, Vertex_set &vertices, VertexIndexMap& vimap, Matrix &src, Matrix &dest) {
    assert(src.row_dimension()==dest.row_dimension());
    assert(src.column_dimension()==dest.column_dimension());
    BOOST_FOREACH(vertex_descriptor vertex, vertices)  {
      int i = get(vimap,vertex);
      vertex_around_target_circulator v_j(halfedge(vertex, mesh), mesh), end = v_j;
      CGAL_For_all(v_j, end){
        int j = get(vimap, *v_j);
        dest.set_coef(i,j, src.get_coef(i,j), false);
      }
    }
  }

  template <typename VertexIndexMap>
  void printMatrix(Vertex_set &vertices, VertexIndexMap& vimap, Matrix &A, std::string name)  {
    std::cout << "Matrix " << name << "(" << A.row_dimension() << "x" << A.column_dimension() << ")" << std::endl;
    Matrix A1(A.row_dimension(), A.column_dimension());
    int r=0,c=0;
    BOOST_FOREACH(vertex_descriptor v1,vertices) {
      int i = get(vimap, v1);
      BOOST_FOREACH(vertex_descriptor v2,vertices){
        int j = get(vimap, v2);
        A1.set_coef(r,c,A.get_coef(i,j));
        c++;
      }
      r++; c = 0;
    }

    for(int r = 0; r < A.row_dimension(); r++)  {
      for(int c = 0; c < A.column_dimension(); c++)  {
        std::cout << std::setw(10) << A1.get_coef(r,c) << "\t" << std::flush;
      }
      std::cout << std::endl;
    }

  }

  template <typename VertexIndexMap>
  void printVector(Vertex_set &vertices, VertexIndexMap& vimap, Vector &A, std::string name)  {
    std::cout << "Vector " << name << "(" << A.size() << ")" << std::endl;
    Vector A1(A.size());
    int r = 0;
    BOOST_FOREACH(vertex_descriptor v1,vertices) {
      int i = get(vimap,v1);
      A1.set(r, A(i));
      r++;
    }
    for(int r = 0; r < A.size(); r++)  {
      std::cout << A1(r) << std::endl;
    }
  }

  void print(Matrix& A, Vector &Xu, Vector&Bu)  {
    std::cout << "Matrix " << "(" << A.row_dimension() << "x" << A.column_dimension() << ")" << std::endl;
    for(int r = 0; r < A.row_dimension(); r++)  {
      for(int c = 0; c < A.column_dimension(); c++)  {
        std::cout << std::setw(10) << A.get_coef(r,c) << "\t" << std::flush;
      }
      std::cout << "\t\t"  << Xu(r) << "\t\t" << Bu(r) << std::endl;
    }
  }

  // Protected accessors
protected:
  /// Get the object that maps the surface's border onto a 2D space.
  Border_parameterizer& get_border_parameterizer() { return m_borderParameterizer; }

  /// Get the sparse linear algebra (traits object to access the linear system).
  Solver_traits& get_linear_algebra_traits() { return m_linearAlgebra; }

  // Fields
private:
  // Object that maps the surface's border onto a 2D space.
  Border_parameterizer m_borderParameterizer;

  // Traits object to solve a sparse linear system
  Solver_traits m_linearAlgebra;

  // Object that keeps the last best UV Map
  Vertex_point2_map lastBestuvmap;

  // Counter to keep track of failure of Linear solver
  int LScounter;

};

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_FIXED_BORDER_ITERATIVE_PARAMETERIZER_3_H
