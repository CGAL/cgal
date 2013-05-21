// adapted from fill.cpp -IOY
#ifndef CGAL_FILL_HOLE_POLYHEDRON_H
#define CGAL_FILL_HOLE_POLYHEDRON_H

#include <CGAL/Mesh_3/dihedral_angle_3.h>
#include <CGAL/internal/Weights.h>
#include <vector>
#include <limits>
#include <set>
#include <map>

// for default parameters
#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>  // for sparse linear system solver
#if defined(CGAL_SUPERLU_ENABLED)
#include <Eigen/SuperLUSupport>
#else
#include <Eigen/SparseLU>
#endif
#endif

namespace CGAL {

template<class Polyhedron, class SparseLinearSolver, class WeightCalculator>
class Fill_hole_Polyhedron_3 {

  typedef typename Polyhedron::Traits::Point_3 Point_3;
  typedef typename Polyhedron::Vertex_handle Vertex_handle;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef typename Polyhedron::Facet_handle Facet_handle;
  typedef typename Polyhedron::Vertex_iterator Vertex_iterator;
  typedef typename Polyhedron::Halfedge_iterator Halfedge_iterator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator  Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator  Halfedge_around_vertex_circulator;
  typedef std::vector<Halfedge_handle> Polyline_3;

  typedef SparseLinearSolver Sparse_linear_solver;

//members
  std::map<Vertex_handle, double> scale_attribute;
  bool is_refining;
  bool is_fairing;
  double alpha;
  WeightCalculator weight_calculator;

public:
  Fill_hole_Polyhedron_3(bool is_refining, double alpha, bool is_fairing, WeightCalculator weight_calculator = WeightCalculator()) 
    : is_refining(is_refining), is_fairing(is_fairing), alpha(alpha), weight_calculator(weight_calculator)
  { }

private:
  struct Weight {

    std::pair<double,double> w;

    Weight()
      : w(std::make_pair<double>(0, 0))
    {}

    Weight(double angle, double area)
      : w(angle, area)
    {}

    friend Weight operator+(const Weight& w1, const Weight& w2)
    {
      return Weight((std::max)(w1.w.first, w2.w.first), w1.w.second + w2.w.second);
    }

    friend bool operator<(const Weight& w1, const Weight& w2)
    {
      if(w1.w.first < w2.w.first){
        return true;
      }else if(w1.w.first == w2.w.first){
        return w1.w.second < w2.w.second;
      }
      return false;
    }
  };

  Weight
  sigma(const Polyline_3& P, int i, int j, int k)
  {
    assert(i+1 ==j);
    assert(j+1 ==k);
    const Point_3& p = P[i]->opposite()->vertex()->point();
    const Point_3& q = P[j]->opposite()->vertex()->point();
    const Point_3& r = P[k]->opposite()->vertex()->point();
    // The CGAL::dihedral angle is measured between the oriented triangles, that is it goes from [-pi, pi]
    // What we need is the angle between the normals of the triangles between [0, pi]
    double ang1 = 180 - CGAL::abs( CGAL::Mesh_3::dihedral_angle(p,q,r,P[i]->opposite()->next()->vertex()->point()));
    double ang2 = 180 - CGAL::abs( CGAL::Mesh_3::dihedral_angle(q,r,p, P[j]->opposite()->next()->vertex()->point()));
    return Weight((std::max)(ang1, ang2), sqrt(CGAL::squared_area(p,q,r)));
  }

  Weight
  sigma(const Polyline_3& P, int i, int m, int k, const std::vector<int>& lambda)
  {
    assert(i < m);
    assert(m < k);
    int n = P.size() -1; // because the first and last point are equal
    const Point_3& pi = P[i]->opposite()->vertex()->point();
    const Point_3& pm = P[m]->opposite()->vertex()->point();
    const Point_3& pk = P[k]->opposite()->vertex()->point();
    double ang1=0, ang2=0;
    if(lambda[i*n+m] != -1){
      const Point_3& pim = P[lambda[i*n+m]]->opposite()->vertex()->point();
      ang1 = 180 - CGAL::abs(CGAL::Mesh_3::dihedral_angle(pi,pm,pk, pim));
    }
    if(lambda[m*n+k] != -1){
      const Point_3& pmk = P[lambda[m*n+k]]->opposite()->vertex()->point();
      ang2 = 180 - CGAL::abs( CGAL::Mesh_3::dihedral_angle(pm,pk,pi, pmk));
    }
    return Weight((std::max)(ang1, ang2), sqrt(CGAL::squared_area(pi,pm,pk)));
  }

  Halfedge_handle add_facets(const Polyline_3& P, const std::vector<int>& lambda, 
    int i, int k, Polyhedron& poly, std::set<Facet_handle>& facets, bool last = true)
  {
    Halfedge_handle h, g;
    int n = P.size() -1; // because the first and last point are equal
    if(i+2 == k){
      if(last){
        h = poly.fill_hole(P[i+1]);
        assert(h->facet() != Facet_handle());
        facets.insert(h->facet());
        return h;
      } else {
        h = poly.add_facet_to_border(P[i]->prev(), P[i+1]);
        assert(h->facet() != Facet_handle());
        facets.insert(h->facet());
        return h->opposite();
      }
    } else {
      int la = lambda[i*n + k];
      if(la != i+1){
        h = add_facets(P, lambda, i, la, poly, facets, false);
        assert(h->opposite()->facet() != Facet_handle());
        facets.insert(h->opposite()->facet());
      } else {
        h = P[i];
      }
      if(la != k-1){
        g = add_facets(P, lambda, la, k, poly, facets, false);
        assert(g->opposite()->facet() != Facet_handle());
        facets.insert(g->opposite()->facet());
      } else {
        g = P[la];
      }
      if(last){
        h = poly.fill_hole(g);
        assert(h->facet() != Facet_handle());
        facets.insert(h->facet());
        return h;
      } else {
        h = poly.add_facet_to_border(h->prev(), g);
        assert(h->facet() != Facet_handle());
        facets.insert(h->facet());
        return h->opposite();
      }
    }
  }

  void average_length(Polyhedron& poly, Vertex_handle vh)
  {
    const Point_3& vp = vh->point(); 
    Halfedge_around_vertex_circulator circ(vh->vertex_begin()), done(circ);
    int deg = 0;
    double sum = 0;
    do {
      const Point_3& vq = circ->opposite()->vertex()->point();
      sum += sqrt(CGAL::squared_distance(vp, vq));
      ++deg;
      ++circ;
    } while(circ != done);
    scale_attribute[vh] = sum/deg;
  }

  bool relax(Polyhedron& poly, Halfedge_handle h)
  {
    const Point_3& p = h->vertex()->point();
    const Point_3& q = h->opposite()->vertex()->point();
    const Point_3& r = h->next()->vertex()->point();
    const Point_3& s = h->opposite()->next()->vertex()->point();
    if( (CGAL::ON_UNBOUNDED_SIDE  != CGAL::side_of_bounded_sphere(p,q,r,s)) ||
      (CGAL::ON_UNBOUNDED_SIDE  != CGAL::side_of_bounded_sphere(p,q,s,r)) ){
        poly.flip_edge(h);
        return true;
    }
    return false;
  }

  bool subdivide(Polyhedron& poly, std::set<Facet_handle>& facets)
  {
    std::list<Facet_handle> new_facets;
    for(typename std::set<Facet_handle>::iterator it = facets.begin(); it!= facets.end(); ++it){
      if(*it  == Facet_handle()){
      } else {
        Halfedge_handle hh =  (*it)->halfedge();
        Vertex_handle vi = (*it)->halfedge()->vertex();
        Vertex_handle vj = (*it)->halfedge()->next()->vertex();
        Vertex_handle vk = (*it)->halfedge()->prev()->vertex();
        Point_3 c = CGAL::centroid(vi->point(), vj->point(), vk->point());
        double sac  = (scale_attribute[vi] + scale_attribute[vj] + scale_attribute[vk])/3.0;
        double dist_c_vi = sqrt(CGAL::squared_distance(c,vi->point()));
        double dist_c_vj = sqrt(CGAL::squared_distance(c,vj->point()));
        double dist_c_vk = sqrt(CGAL::squared_distance(c,vk->point()));
        if((alpha * dist_c_vi > sac) &&
          (alpha * dist_c_vj > sac) &&
          (alpha * dist_c_vk > sac) &&
          (alpha * dist_c_vi > scale_attribute[vi]) &&
          (alpha * dist_c_vj > scale_attribute[vj]) &&
          (alpha * dist_c_vk > scale_attribute[vk])){
            Halfedge_handle h = poly.create_center_vertex((*it)->halfedge());
            h->vertex()->point() = c;
            scale_attribute[h->vertex()] = sac;

            // collect 2 new facets for next round 
            Facet_handle h1 = h->next()->opposite()->face();
            Facet_handle h2 = h->opposite()->face();
            new_facets.push_back(h1);
            new_facets.push_back(h2);
            // relax edges of the  patching mesh 
            Halfedge_handle e_ij = h->prev();
            Halfedge_handle e_ik = h->opposite()->next();
            Halfedge_handle e_jk = h->next()->opposite()->prev();

            if(facets.find(e_ij->opposite()->face()) != facets.end()){
              relax(poly, e_ij);
            }
            if(facets.find(e_ik->opposite()->face()) != facets.end()){
              relax(poly, e_ik);
            }
            if(facets.find(e_jk->opposite()->face()) != facets.end()){
              relax(poly, e_jk);
            }
        }
      }
    }
    for(typename std::list<Facet_handle>::iterator it = new_facets.begin();
      it!= new_facets.end();
      ++it){
        facets.insert(*it);
    }
    return ! new_facets.empty();
  }

  bool relax(Polyhedron& poly, std::set<Facet_handle>& facets)
  {
    int flips = 0;
    std::list<Halfedge_handle> interior_edges;
    for(typename std::set<Facet_handle>::iterator it = facets.begin(); it!= facets.end(); ++it){
      Halfedge_around_facet_circulator  circ = (*it)->facet_begin(), done(circ);
      do {
        Halfedge_handle h = circ;
        Halfedge_handle oh = h->opposite();
        if(facets.find(oh->face()) != facets.end()){
          // it's an interior edge
          interior_edges.push_back((h < oh) ? h : oh );
        }
        ++circ;
      } while(circ != done);
    }
    std::cerr << "Test " << interior_edges.size() << " edges " << std::endl;
    for(typename std::list<Halfedge_handle>::iterator it = interior_edges.begin();
      it != interior_edges.end();
      ++it){
        if(relax(poly,*it)){
          ++flips;
        }
    }
    std::cerr << "|flips| = " << flips << std::endl;
    return flips > 0;
  }

  void refine(Polyhedron& poly, std::set<Facet_handle>& facets)
  {
    int i = 0;
    do {
      if(i == 10){
        break;
      }
      i++;
    } while( subdivide(poly, facets) && relax(poly,facets) );

    // according to paper it should be like below (?) IOY
    //while(true) {
    //  bool subdiv = subdivide(poly, facets);
    //  if(!subdiv) { break; }
    //  while(relax(poly,facets)) {}
    //}
  }

  void compute_lambda(const Polyline_3& P, std::vector<int>& lambda) {
    int n = P.size() - 1; // because the first and last point are equal
    lambda = std::vector<int>(n*n,-1);
    std::vector<Weight> W(n*n,Weight(0,0));

    for(int i=0; i < n-2; ++i){
      W[i*n + (i+2)] = sigma(P, i, i+1, i+2);
      lambda[i*n + (i+2)] = i+1;
    }

    for(int j = 3; j< n; ++j){
      for(int i=0; i<n-j; ++i){
        int k = i+j;
        int m_min = 0;
        Weight w_min((std::numeric_limits<double>::max)(), (std::numeric_limits<double>::max)());
        for(int m = i+1; m<k; ++m){
          Weight w = W[i*n + m] + W[m*n + k] + sigma(P,i,m,k, lambda);
          if(w < w_min){
            w_min = w;
            m_min = m;
          }
        }
        W[i*n+k] = w_min;
        lambda[i*n+k] = m_min;
      }
    }
  }

  double weight_sum(Vertex_handle v, Polyhedron& polyhedron) {
    double weight = 0;
    Halfedge_around_vertex_circulator circ = v->vertex_begin();
    do {
      weight += weight_calculator(circ, polyhedron);
    } while(++circ != v->vertex_begin());
    return weight;
  }

  //template<class InputIterator>
  //void fill_matrix(InputIterator vb, InputIterator ve, 
  //  typename Sparse_linear_solver::Matrix& A, 
  //  typename Sparse_linear_solver::Vector& Bx,
  //  typename Sparse_linear_solver::Vector& By,
  //  typename Sparse_linear_solver::Vector& Bz,
  //  const std::set<Vertex_handle>& interior_vertices,
  //  std::map<Vertex_handle, std::size_t>& vertex_id_map
  //  )
  //{
  //  
  //}

  void fair(std::set<Vertex_handle>& interior_vertices, Polyhedron& polyhedron)
  {
    const std::size_t nb_vertices = interior_vertices.size();
    typename Sparse_linear_solver::Vector X(nb_vertices), Bx(nb_vertices);
    typename Sparse_linear_solver::Vector Y(nb_vertices), By(nb_vertices);
    typename Sparse_linear_solver::Vector Z(nb_vertices), Bz(nb_vertices);

    std::map<Vertex_handle, std::size_t> vertex_id_map;
    std::size_t id = 0;
    for(typename std::set<Vertex_handle>::iterator it = interior_vertices.begin(); it != interior_vertices.end(); ++it, ++id) {
      vertex_id_map[*it] = id;      
    }
    
    typename Sparse_linear_solver::Matrix A(nb_vertices);
    // this one-to-one corresponds to "second-order weighted umbrella operator" equation in [Filling Holes in Meshes] paper
    // fixed vertex positions are transfered to right hand side
    for(typename std::set<Vertex_handle>::iterator vb = interior_vertices.begin(); vb != interior_vertices.end(); ++vb) {
      std::size_t v_id = vertex_id_map[*vb];
      double v_weight_sum = weight_sum(*vb, polyhedron);
      double x, y, z;
      x = y = z = 0.0;
      Halfedge_around_vertex_circulator circ = (*vb)->vertex_begin();
      do {
        Vertex_handle nv = circ->opposite()->vertex();        
        double nv_weight = weight_calculator(circ, polyhedron);
        double nv_normalized_weight = nv_weight / v_weight_sum;

        if(interior_vertices.find(nv) != interior_vertices.end()) {
          std::size_t nv_id = vertex_id_map[nv];
          A.add_coef(v_id, nv_id, -nv_normalized_weight);
        }
        else {
          x += nv_normalized_weight * nv->point().x();
          y += nv_normalized_weight * nv->point().y();
          z += nv_normalized_weight * nv->point().z();
        }

        double nv_weight_sum = weight_sum(nv, polyhedron);
        Halfedge_around_vertex_circulator circ2 = nv->vertex_begin();
        do {
          Vertex_handle nnv = circ2->opposite()->vertex();          
          double nnv_weight = weight_calculator(circ2, polyhedron);
          double nnv_normalized_weight = nnv_weight / nv_weight_sum;

          if(interior_vertices.find(nnv) != interior_vertices.end()) {
            std::size_t nnv_id = vertex_id_map[nnv];
            A.add_coef(v_id, nnv_id, nv_normalized_weight * nnv_normalized_weight);
          }
          else {
            x += - nv_normalized_weight * nnv_normalized_weight * nnv->point().x();
            y += - nv_normalized_weight * nnv_normalized_weight * nnv->point().y();
            z += - nv_normalized_weight * nnv_normalized_weight * nnv->point().z();
          }
        } while(++circ2 != nv->vertex_begin());

        if(interior_vertices.find(nv) != interior_vertices.end()) {
          std::size_t nv_id = vertex_id_map[nv];
          A.add_coef(v_id, nv_id, -nv_normalized_weight);
        }
        else {
          x += nv_normalized_weight * nv->point().x();
          y += nv_normalized_weight * nv->point().y();
          z += nv_normalized_weight * nv->point().z();
        }
      } while(++circ != (*vb)->vertex_begin());

      A.add_coef(v_id, v_id, 1);
      Bx[v_id] = x; By[v_id] = y; Bz[v_id] = z;
    }

    // factorize
    double D;
    Sparse_linear_solver m_solver;
    bool prefactor_ok = m_solver.pre_factor(A, D);
    if(!prefactor_ok) {
      CGAL_warning(false && "pre_factor failed!");
      return;
    }
    // solve
    bool is_all_solved = m_solver.linear_solver(Bx, X) && m_solver.linear_solver(By, Y) && m_solver.linear_solver(Bz, Z);
    if(!is_all_solved) {
      CGAL_warning(false && "linear_solver failed!"); 
      return; 
    }
    // update 
    id = 0;
    for(typename std::set<Vertex_handle>::iterator it = interior_vertices.begin(); it != interior_vertices.end(); ++it, ++id) {
      (*it)->point() = Point_3(X[id], Y[id], Z[id]);
    }

    //typedef Sparse_linear_solver::Matrix::EigenType EigenMatrix;
    //EigenMatrix A2 = A.eigen_object()* A.eigen_object();
    //Solver solver;
    //solver.compute(A.eigen_object());
    //if(solver.info() != Eigen::Success) { CGAL_warning(false); return; }
    //X = solver.solve(Bx);
    //Y = solver.solve(By);
    //Z = solver.solve(Bz);
  }

public:
  void operator()(Polyhedron& poly, Halfedge_handle it)
  {
    std::cerr << "fill" << std::endl;
    if(! it->is_border()){
      return;
    }
    
    std::set<Vertex_handle> boundary_vertices;
    Polyline_3 P;
    Halfedge_around_facet_circulator circ(it), done(circ);
    do{
      boundary_vertices.insert(circ->vertex());
      P.push_back(circ);
      average_length(poly, circ->vertex());
    } while (++circ != done);
    P.push_back(circ);

    std::vector<int> lambda;
    compute_lambda(P, lambda);

    int n = P.size() - 1; // because the first and last point are equal
    std::set<Facet_handle> facets;
    add_facets(P, lambda, 0, n-1, poly, facets);

    if(is_refining) {
      refine(poly, facets);
      std::cerr << "|facets| = " << facets.size() << std::endl;

      if(is_fairing) {
        std::set<Vertex_handle> interior_vertices;
        for(typename std::set<Facet_handle>::iterator it = facets.begin(); it != facets.end(); ++it) {
          Halfedge_around_facet_circulator circ = (*it)->facet_begin();
          do {
            if(boundary_vertices.find(circ->vertex()) == boundary_vertices.end()) {
              interior_vertices.insert(circ->vertex());
            }
          } while(++circ != (*it)->facet_begin());
        }

        std::cerr << "|boundary vertices| = " << boundary_vertices.size() << std::endl;
        std::cerr << "|interior vertices| = " << interior_vertices.size() << std::endl;

        fair(interior_vertices, poly);
      } // if(is_fairing)
    } // if(is_refining)
  }
};

enum Fairing_weight_type_tag {
  UNIFORM_WEIGHTING,
  SCALE_DEPENDENT_WEIGHTING,
  COTANGENT_WEIGHTING
};

template<class Sparse_linear_system, class Polyhedron>
void fill(Polyhedron& poly, 
  typename Polyhedron::Halfedge_handle it, 
  bool refine = true,
  double density_control_factor = 1.41 /* ~sqrt(2) */,
  bool fair = true,
  Fairing_weight_type_tag weight_tag = SCALE_DEPENDENT_WEIGHTING
  )
{
  typedef CGAL::internal::Uniform_weight<Polyhedron> Uniform_weight;
  typedef CGAL::internal::Cotangent_weight<Polyhedron, CGAL::internal::Cotangent_value_Meyer<Polyhedron> > Cotangent_weight;
  typedef CGAL::internal::Scale_dependent_weight<Polyhedron> Scale_dependent_weight;

  if(weight_tag == COTANGENT_WEIGHTING) {
    Fill_hole_Polyhedron_3<Polyhedron, Sparse_linear_system, Cotangent_weight>
      (refine, density_control_factor, fair)(poly, it);
  }
  else if(weight_tag == UNIFORM_WEIGHTING) {
    Fill_hole_Polyhedron_3<Polyhedron, Sparse_linear_system, Uniform_weight>
      (refine, density_control_factor, fair)(poly, it);
  }
  else if(weight_tag == SCALE_DEPENDENT_WEIGHTING) {
    Fill_hole_Polyhedron_3<Polyhedron, Sparse_linear_system, Scale_dependent_weight>
      (refine, density_control_factor, fair)(poly, it);
  }
}

template<class Polyhedron>
void fill(Polyhedron& poly, 
  typename Polyhedron::Halfedge_handle it, 
  bool refine = true,
  double density_control_factor = 1.41 /* ~sqrt(2) */,
  bool fair = true,
  Fairing_weight_type_tag weight_tag = SCALE_DEPENDENT_WEIGHTING
  )
{
  typedef
#if defined(CGAL_EIGEN3_ENABLED)
#if defined(CGAL_SUPERLU_ENABLED)
  CGAL::Eigen_solver_traits<Eigen::SuperLU<CGAL::Eigen_sparse_matrix<double>::EigenType> >
#else
  CGAL::Eigen_solver_traits<
    Eigen::SparseLU<
      CGAL::Eigen_sparse_matrix<double, Eigen::ColMajor>::EigenType,
      Eigen::COLAMDOrdering<int> >  >
#endif
#endif
  Sparse_linear_system;

  fill<Sparse_linear_system, Polyhedron>(poly, it, refine, density_control_factor, fair, weight_tag);
}

} // namespace CGAL

#endif















