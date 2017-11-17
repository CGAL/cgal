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

  // Eigen sparse system
  typedef typename Eigen::SparseMatrix<double> Eigen_matrix;
  typedef typename Eigen::VectorXd Eigen_vector;

  typedef typename Eigen::BiCGSTAB<Eigen_matrix, Eigen::IncompleteLUT<double> > Eigen_solver;
  // typedef typename Eigen::SimplicialLDLT<Eigen_matrix> Eigen_solver;



// data
private:

    IndexMap vimap_ = get(boost::vertex_index, mesh_);

    // geometry data
    PolygonMesh& mesh_;
    VertexPointMap& vpmap_;
    Edge_cotangent_weight<PolygonMesh, VertexPointMap> weight_calculator_;
    Incident_area<PolygonMesh> inc_areas_calculator_;

    std::size_t nb_vert_;

    // constrained vertices
    std::unordered_set<vertex_descriptor> constrained_vertices_;


public:

  Shape_smoother(PolygonMesh& mesh, VertexPointMap& vpmap) : mesh_(mesh), vpmap_(vpmap),
      weight_calculator_(mesh, vpmap),
      inc_areas_calculator_(mesh),
      nb_vert_(static_cast<int>(vertices(mesh).size()))
  {  }

  void setup_system(Eigen_matrix& A, Eigen_matrix& L, Eigen_matrix& D,
                    Eigen_vector& bx, Eigen_vector& by, Eigen_vector& bz,
                    const double& time)
  {

    calc_mass_matrix(D);

    // todo : assert L, D sizes

    //std::cout<<"mass matrix= "<<D<<std::endl;

    //std::cerr << "compute coefficient matrix...";
    compute_coeff_matrix(A, L, D, time);

    //std::cout<<"mass matrix= "<<A<<std::endl;

    //std::cerr << "done" << std::endl;

    std::cerr << "rhs...";
    compute_rhs(bx, by, bz, D);


    std::cerr << "done" << std::endl;
  }



  //void solve_system(Eigen_matrix& stiffness_matrix)
  void solve_system(Eigen_matrix& A,
                    Eigen_vector& Xx, Eigen_vector& Xy, Eigen_vector& Xz,
                    Eigen_vector& bx, Eigen_vector& by, Eigen_vector& bz)
  {
    // Cholesky factorization
    std::cerr << "compute ...";
    Eigen_solver solver;
    solver.compute(A);
    std::cerr << "done" << std::endl;
    // solver.analyzePattern(A);
      // solver.factorize(A);

      // back-substitution
    std::cerr << "solve...";
    Xx = solver.solve(bx);
    Xy = solver.solve(by);
    Xz = solver.solve(bz);
    std::cerr << "done" << std::endl;

    if(solver.info() != Eigen::Success)
    {
      std::cerr << "Not Solved!" << std::endl;
      return;
    }

  }

  void calc_stiff_matrix(Eigen_matrix& mat)
  {

    //mat.resize(nb_vert_, nb_vert_);

  // TOFIX: very slow way to fill the matrx
  // avoid matrix C!! and use triplets
  // read http://eigen.tuxfamily.org/dox/group__TutorialSparse.html#TutorialSparseFilling

    // cot values
    Eigen_matrix C = cot_entries(); //TODO

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

  }

  void update_mesh(Eigen_vector& Xx, Eigen_vector& Xy, Eigen_vector& Xz)
  {
    update_map(Xx, Xy, Xz);

    // normalize area
    NT surface_area = area(faces(mesh_), mesh_);
    std::cout << "area= " << surface_area << nl;
    std::cout << "area sqrt= " << CGAL::sqrt(surface_area) << nl;


    normalize_area(Xx, Xy, Xz);

    update_map(Xx, Xy, Xz);

    surface_area = area(faces(mesh_), mesh_);
    std::cout << "surface_area normalized= " << surface_area << nl;

  }


private:

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




    void calc_mass_matrix(Eigen_matrix& D)
    {
      //D.resize(nb_vert_, nb_vert_);

      for(face_descriptor f : faces(mesh_))
      {
        double area = face_area(f, mesh_);
        for(vertex_descriptor v : vertices_around_face(halfedge(f, mesh_), mesh_))
        {
          auto indx = vimap_[v];
          D.coeffRef(indx, indx) += 2.0 * area;
        }
      }

      D /= 12.0;
    }



    void compute_coeff_matrix(Eigen_matrix& A, Eigen_matrix& L, Eigen_matrix& D, const double time)
    {

      double delta = time;

      //Eigen_matrix D(L.rows(), L.cols());
      //D.setIdentity();

      //A.resize(nb_vert_, nb_vert_);
      A = D - delta * L;


    }


    void compute_rhs(Eigen_vector& bx, Eigen_vector& by, Eigen_vector& bz,
                     Eigen_matrix& D)
    {
      //bx.resize(nb_vert_);
      //by.resize(nb_vert_);
      //bz.resize(nb_vert_);

      for(vertex_descriptor vi : vertices(mesh_))
      {
        int index = vimap_[vi];
        Point p = get(vpmap_, vi);
        bx.coeffRef(index) = p.x();
        by.coeffRef(index) = p.y();
        bz.coeffRef(index) = p.z();
      }

      //calc_mass_matrix(D);

      //Eigen_matrix D(Bx.rows(), Bx.rows());
      //D.setIdentity();

      bx = D * bx;
      by = D * by;
      bz = D * bz;


      //std::cout<<"bx= "<<bx<<std::endl;
      //std::cout<<"by= "<<by<<std::endl;
      //std::cout<<"bz= "<<bz<<std::endl;
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

    /*
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
    */

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

    /*
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
    */

    void normalize_area(Eigen_vector& Xx, Eigen_vector& Xy, Eigen_vector& Xz)
    {

      NT surface_area = area(faces(mesh_), mesh_);

      Xx /= CGAL::sqrt(surface_area);
      Xy /= CGAL::sqrt(surface_area);
      Xz /= CGAL::sqrt(surface_area);

    }


    void update_map(Eigen_vector& Xx, Eigen_vector& Xy, Eigen_vector& Xz)
    {
        for (vertex_descriptor v : vertices(mesh_))
        {
            int index = get(vimap_, v);
            NT x_new = Xx.coeffRef(index);
            NT y_new = Xy.coeffRef(index);
            NT z_new = Xz.coeffRef(index);
            put(vpmap_, v, Point(x_new, y_new, z_new));
        }
    }




};


}















} // PMP
} // CGAL


#endif // SMOOTHING_H
