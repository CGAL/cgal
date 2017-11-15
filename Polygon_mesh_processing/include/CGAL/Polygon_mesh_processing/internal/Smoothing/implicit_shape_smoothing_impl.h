#ifndef SMOOTHING_H
#define SMOOTHING_H

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/Weights.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Eigen_matrix.h>
#include <CGAL/barycenter.h>
#include <Eigen/Sparse>

#include <fstream>
#include <unordered_set>
#include <unordered_map>

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {
//Cotangent_value_Meyer<PolygonMesh, VertexPointMap>>

std::ostream& nl(std::ostream& out)
{
  return out << "\n";
}

template<typename PolygonMesh, typename VertexPointMap,
         typename CotangentValue = CGAL::internal::Cotangent_value_Meyer<PolygonMesh, VertexPointMap>>
struct Edge_cotangent_weight : CotangentValue
{
    Edge_cotangent_weight(PolygonMesh& pmesh_, VertexPointMap vpmap_)
      : CotangentValue(pmesh_, vpmap_)
    {}

    PolygonMesh& pmesh()
    {
      return CotangentValue::pmesh();
    }

    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor   halfedge_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor     vertex_descriptor;

    double operator()(halfedge_descriptor he)
    {
      if(is_border_edge(he, pmesh()))
      {
        halfedge_descriptor h1 = next(he, pmesh());
        vertex_descriptor vs = source(he, pmesh());
        vertex_descriptor vt = target(he, pmesh());
        vertex_descriptor v1 = target(h1, pmesh());

        return (CotangentValue::operator ()(vs, v1, vt));
      }
      else
      {
        halfedge_descriptor h1 = next(he, pmesh());
        halfedge_descriptor h2 = prev(opposite(he, pmesh()), pmesh());
        vertex_descriptor vs = source(he, pmesh());
        vertex_descriptor vt = target(he, pmesh());
        vertex_descriptor v1 = target(h1, pmesh());
        vertex_descriptor v2 = source(h2, pmesh());

        return ( CotangentValue::operator()(vs, v1, vt) + CotangentValue::operator()(vs, v2, vt) ) / 2.0;
      }
    }
};






template<typename PolygonMesh>
struct Incident_area
{

  Incident_area(PolygonMesh& mesh) : pmesh(mesh){}

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;

  double operator()(halfedge_descriptor he)
  {

    halfedge_descriptor hopp = opposite(he, pmesh);
    face_descriptor f1 = face(he, pmesh);
    face_descriptor f2 = face(hopp, pmesh);

    double A1 = f1 == boost::graph_traits<PolygonMesh>::null_face() ? 0 : face_area(f1, pmesh);
    double A2 = f2 == boost::graph_traits<PolygonMesh>::null_face() ? 0 : face_area(f2, pmesh);

    // todo: check degenerecies
    //CGAL_assertion(A1>0 && A2>0);

    return A1 + A2;
  }

  PolygonMesh& pmesh;

};





template<typename PolygonMesh, typename VertexPointMap>
class Shape_smoother{

// types
private:

  typedef typename GetGeomTraits<PolygonMesh>::type GeomTraits;
  typedef typename GeomTraits::FT NT;
  typedef typename GeomTraits::Point_3 Point;
  typedef typename GeomTraits::Triangle_3 Triangle;

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;

  // vertex index map
  typedef typename boost::property_map<PolygonMesh, boost::vertex_index_t>::type IndexMap;

  //CGAL Sparse solver
  typedef typename CGAL::Eigen_sparse_matrix<double>::EigenType CEigenMatrix;

  typedef CGAL::Eigen_solver_traits< Eigen::SimplicialLDLT< CEigenMatrix > > Solver_traits;
  //typedef CGAL::Eigen_solver_traits<> Solver_traits; //BicGSTAB

  typedef typename Solver_traits::Matrix Matrix;
  typedef typename Solver_traits::Vector Vector;

  // Eigen sparse matrix & vector
  typedef typename Eigen::SparseMatrix<double> Eigen_matrix;
  typedef typename Eigen::SparseVector<double> Eigen_vector;



// data
private:

    IndexMap vimap_ = get(boost::vertex_index, mesh_);

    // geometry data
    PolygonMesh& mesh_;
    VertexPointMap& vpmap_;
    Edge_cotangent_weight<PolygonMesh, VertexPointMap> weight_calculator_;
    Incident_area<PolygonMesh> inc_areas_calculator_;

    std::size_t nb_vert_;

    // linear solver
    Solver_traits solver_;

    // constrained vertices
    std::unordered_set<vertex_descriptor> constrained_vertices_;


public:

  Shape_smoother(PolygonMesh& mesh, VertexPointMap& vpmap) : mesh_(mesh), vpmap_(vpmap),
      weight_calculator_(mesh, vpmap),
      inc_areas_calculator_(mesh),
      nb_vert_(static_cast<int>(vertices(mesh).size()))
  {  }


  Eigen_matrix calc_stiff_matrix()
  {
    Eigen_matrix L = stiff_matrix();
    return L;
  }


  void solve_system(Eigen::SparseMatrix<double>& L, const double time)
  {

    //gather_constrained_vertices();

    Matrix A(nb_vert_, nb_vert_);
    Vector Bx(nb_vert_);
    Vector By(nb_vert_);
    Vector Bz(nb_vert_);
    Vector Xx(nb_vert_);
    Vector Xy(nb_vert_);
    Vector Xz(nb_vert_);


    compute_coeff_matrix(A, L, time);

    //apply_constraints(A);


    //extract_matrix(A);


    compute_rhs(Bx, By, Bz);

    //extract_vectors(Bx,By,Bz);


    NT dx, dy, dz;
    if(!solver_.linear_solver(A, Bx, Xx, dx) ||
       !solver_.linear_solver(A, By, Xy, dy) ||
       !solver_.linear_solver(A, Bz, Xz, dz) )
    {
        std::cerr<<"Could not solve linear system!"<<nl;
        bool solved = false;
        CGAL_assertion(solved);
    }

    //extract_vectors(Xx,Xy,Xz);

    //translate_centroid(Xx, Xy, Xz);
    //extract_vectors(Xx,Xy,Xz);

    NT surface_area = area(faces(mesh_), mesh_);
    std::cout << "area= " << surface_area << nl;
    std::cout << "area sqrt= " << CGAL::sqrt(surface_area) << nl;


    normalize_area(Xx, Xy, Xz);
    //extract_vectors(Xx, Xy, Xz);

    update_map(Xx, Xy, Xz);

    surface_area = area(faces(mesh_), mesh_);
    std::cout << "surface_area normalized= " << surface_area << nl;


  }


private:

     // ----------- LINEAR SYSTEM -------- //
    Eigen_matrix get_stiffness_matrix()
    {
      Eigen_matrix mat(nb_vert_, nb_vert_);
      for(vertex_descriptor vi : vertices(mesh_))
      {
        NT sum_Lik = 0;
        for(halfedge_descriptor h : halfedges_around_source(vi, mesh_))
        {
          NT Lij = weight_calculator_(h);
          sum_Lik -= Lij;
          vertex_descriptor vj = target(h, mesh_);
          mat.coeffRef(vimap_[vi], vimap_[vj]) = Lij;
        }

        mat.coeffRef(vimap_[vi], vimap_[vi]) = sum_Lik;
      }

      return mat;
    }

    Eigen_matrix get_mass_matrix()
    {
      Eigen_matrix mat(nb_vert_, nb_vert_);
      for(vertex_descriptor vi : vertices(mesh_))
      {
        NT sum_Dik = 0;
        for(halfedge_descriptor h : halfedges_around_source(vi, mesh_))
        {
          NT Dij = inc_areas_calculator_(h) / 12.0;
          sum_Dik += Dij;
          vertex_descriptor vj = target(h, mesh_);
          mat.coeffRef(vimap_[vi], vimap_[vj]) = Dij;
        }

        mat.coeffRef(vimap_[vi], vimap_[vi]) = sum_Dik;
      }

      return mat;
    }


    Eigen_matrix cot_entries()
    {
      Eigen_matrix C(nb_vert_, nb_vert_);

      for(vertex_descriptor v : vertices(mesh_))
      {
        for(halfedge_descriptor h : halfedges_around_source(v, mesh_))
        {
          NT Lij = weight_calculator_(h);

          vertex_descriptor v_source = source(h, mesh_);
          vertex_descriptor v_target = target(h, mesh_);
          auto i_source = vimap_[v_source];
          auto i_target = vimap_[v_target];

          C.coeffRef(i_source, i_target) = Lij;
        }
      }

      return C;
    }

    Eigen_matrix stiff_matrix()
    {
      Eigen_matrix mat(nb_vert_, nb_vert_);

      // cot values
      Eigen_matrix C = cot_entries();

      for(face_descriptor f : faces(mesh_))
      {
        for(halfedge_descriptor hi : halfedges_around_face(halfedge(f, mesh_), mesh_))
        {
          vertex_descriptor v_source = source(hi, mesh_);
          vertex_descriptor v_target = target(hi, mesh_);
          auto i_source = vimap_[v_source];
          auto i_target = vimap_[v_target];

          auto cot_val = C.coeff(i_source, i_target);
          // i,j
          mat.coeffRef(i_source, i_target) += cot_val; // why not (i,j)
          // j,i
          mat.coeffRef(i_target, i_source) += cot_val;
          //diagonal
          mat.coeffRef(i_source, i_source) -= cot_val;
          mat.coeffRef(i_target, i_target) -= cot_val;
        }
      }

      return mat;
    }


    Eigen_matrix mass_matrix()
    {
      Eigen_matrix mat(nb_vert_, nb_vert_);
      for(face_descriptor f : faces(mesh_))
      {

        double area = face_area(f, mesh_);

        for(vertex_descriptor v : vertices_around_face(halfedge(f, mesh_), mesh_))
        {
          auto indx = vimap_[v];
          mat.coeffRef(indx, indx) += 2.0 * area;
        }

      }

      //export_eigen_matrix(mat, "data/D");

      mat /= 12.0;

      //export_eigen_matrix(mat, "data/D12");

      return mat;
    }




    void compute_coeff_matrix(Matrix& A, Eigen_matrix& L, const double time)
    {
      //Eigen_matrix L = get_stiffness_matrix();

      //Eigen_matrix L = stiff_matrix();

      export_eigen_matrix(L, "data/L");

      //Eigen_matrix D = get_mass_matrix();
      Eigen_matrix D = mass_matrix();

      export_eigen_matrix(D, "data/D");


      Eigen_matrix Ae = D - time * L;


      fill_sparse_matrix(A, Ae);

      A.assemble_matrix();

    }


    void compute_rhs(Vector& Bx, Vector& By, Vector& Bz)
    {

      Eigen_vector Xx(nb_vert_);
      Eigen_vector Xy(nb_vert_);
      Eigen_vector Xz(nb_vert_);

      for(vertex_descriptor vi : vertices(mesh_))
      {
        int index = vimap_[vi];
        Point p = get(vpmap_, vi);
        Xx.coeffRef(index) = p.x();
        Xy.coeffRef(index) = p.y();
        Xz.coeffRef(index) = p.z();
      }


      //Eigen_matrix D = get_mass_matrix();
      Eigen_matrix D = mass_matrix();



      //Eigen_matrix D(Bx.rows(), Bx.rows());
      //D.setIdentity();


      Eigen_vector Bxe = D * Xx;
      Eigen_vector Bye = D * Xy;
      Eigen_vector Bze = D * Xz;


      fill_sparse_vector(Bx, Bxe);
      fill_sparse_vector(By, Bye);
      fill_sparse_vector(Bz, Bze);

    }


    // -------------- COPY TO SPARSE SYSTEM ---------- //
    void fill_sparse_matrix(Matrix&A, Eigen_matrix& Ae)
    {
      CGAL_assertion(A.row_dimension() == Ae.rows() &&
                     A.column_dimension() == Ae.cols());

      // this is not great
      for(std::size_t i=0; i<Ae.rows(); ++i)
      {
        for(std::size_t j=0; j<Ae.cols(); ++j)
        {
          A.set_coef(i, j, Ae.coeff(i, j), true);
        }
      }
    }

    void fill_sparse_vector(Vector& V, Eigen_vector& Ve)
    {
      CGAL_assertion(V.rows() == Ve.rows());
      for(std::size_t i= 0; i<Ve.rows(); ++i)
      {
        V[i] = Ve.coeff(i);
      }
    }


    // ---------------- CONSTRAINS ---------------- //
    void gather_constrained_vertices()
    {
      for(vertex_descriptor v : vertices(mesh_))
      {
        if(is_border(v, mesh_))
        {
          constrained_vertices_.insert(v);
        }
      }
    }

    void apply_constraints(Matrix& A)
    { // to be called after gather_constrained_vertices()

      typename std::unordered_set<vertex_descriptor>::iterator it;
      for(it = constrained_vertices_.begin(); it != constrained_vertices_.end(); ++it)
      {
        int i = get(vimap_, *it);
        A.set_coef(i, i, 1.0);

        // also set all cols of the same row = 0 - not needed
        //for(std::size_t j = 0; j<A.column_dimension(); ++j)
        //  A.set_coef(i, j, 0);
      }
    }


    // -------------- EXTRACT : to move to debug file----------------- //
    void extract_matrix(Matrix& A)
    {
        std::ofstream out("data/mat.dat");
        for(int i=0; i < A.row_dimension(); ++i)
        {
            for(int j=0; j < A.column_dimension(); ++j)
            {
                NT val = A.get_coef(i, j);
                out<<val<<" ";
            }
            out<<nl;
        }
        out.close();
    }

    void extract_vectors(Vector& Vx, Vector& Vy, Vector& Vz)
    {
      CGAL_assertion(Vx.rows() == Vy.rows());
      CGAL_assertion(Vx.rows() == Vz.rows());

      std::ofstream out("data/vecs.dat");
      for(int j=0; j < Vx.dimension(); ++j)
      {
          out<<Vx[j]<<"\t"<<Vy[j]<<"\t"<<Vz[j]<<"\t"<<nl;
      }
      out.close();
    }

    void export_eigen_matrix(Eigen_matrix& A, const char* filename)
    {

      std::ofstream out(filename);
      for(auto i=0; i < A.rows(); ++i)
      {
          for(auto j=0; j < A.cols(); ++j)
          {
              NT val = A.coeff(i, j);
              out<<val<<" ";
          }
          out<<nl;
      }
      out.close();
    }


    // -------------- UPDATE MESH ----------------- //

    Triangle triangle(face_descriptor f) const
    {
      halfedge_descriptor h = halfedge(f, mesh_);
      vertex_descriptor v1  = target(h, mesh_);
      vertex_descriptor v2  = target(next(h, mesh_), mesh_);
      vertex_descriptor v3  = target(next(next(h, mesh_), mesh_), mesh_);
      return Triangle(get(vpmap_, v1), get(vpmap_, v2), get(vpmap_, v3));
    }

    void translate_centroid(Vector& Xx, Vector& Xy, Vector& Xz)
    {

      std::vector<std::pair<Point, NT>> barycenters;

      for(face_descriptor f : faces(mesh_))
      {
        Point tr_centroid = CGAL::centroid(triangle(f));
        barycenters.push_back(std::make_pair(tr_centroid, face_area(f, mesh_)));
      }

      Point centroid = CGAL::barycenter(barycenters.begin(), barycenters.end());
      std::cout << "centroid= " << centroid << std::endl;

      for(std::size_t i=0; i < Xx.rows(); ++i)
      {
        Xx[i] -= centroid.x();
        Xy[i] -= centroid.y();
        Xz[i] -= centroid.z();
      }

    }

    void normalize_area(Vector& Xx, Vector& Xy, Vector& Xz)
    {

      NT surface_area = area(faces(mesh_), mesh_);

      for(std::size_t i=0; i < Xx.rows(); ++i)
      {
        Xx[i] /= CGAL::sqrt(surface_area);
        Xy[i] /= CGAL::sqrt(surface_area);
        Xz[i] /= CGAL::sqrt(surface_area);
      }


      /*
      Xx /= surface_area;
      Xy /= surface_area;
      Xz /= surface_area;
      */

    }


    void update_map(Vector& Xx, Vector& Xy, Vector& Xz)
    {
        for (vertex_descriptor v : vertices(mesh_))
        {
            int index = get(vimap_, v);
            NT x_new = Xx[index];
            NT y_new = Xy[index];
            NT z_new = Xz[index];
            put(vpmap_, v, Point(x_new, y_new, z_new));
        }
    }




};


}















} // PMP
} // CGAL


#endif // SMOOTHING_H
