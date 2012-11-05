#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Mesh_3/dihedral_angle_3.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <set>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef Polyhedron::Vertex_handle Vertex_handle;
typedef Polyhedron::Halfedge_handle Halfedge_handle;
typedef Polyhedron::Facet_handle Facet_handle;
typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator  Halfedge_around_facet_circulator;
typedef Polyhedron::Halfedge_around_vertex_circulator  Halfedge_around_vertex_circulator;

typedef std::vector<Halfedge_handle> Polyline_3;

std::map<Vertex_handle, double> scale_attribute;

struct Weight {

  std::pair<double,double> w;

  Weight()
    : w(std::make_pair<double>(0, 0))
  {}

  Weight(double angle, double area)
    : w(angle, area)
  {}
};

Weight operator+(const Weight& w1, const Weight& w2)
{
  return Weight((std::max)(w1.w.first, w2.w.first), w1.w.second + w2.w.second);
}

bool operator<(const Weight& w1, const Weight& w2)
{
  if(w1.w.first < w2.w.first){
    return true;
  }else if(w1.w.first == w2.w.first){
    return w1.w.second < w2.w.second;
  }
  return false;
}



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




Halfedge_handle add_facets(const Polyline_3& P, const std::vector<int>& lambda, int i, int k, Polyhedron& poly, std::set<Facet_handle>& facets, bool last = true)
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
  double alpha = sqrt(2.0);
  std::list<Facet_handle> new_facets;
  for(std::set<Facet_handle>::iterator it = facets.begin(); it!= facets.end(); ++it){
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
  for(std::list<Facet_handle>::iterator it = new_facets.begin();
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
  for(std::set<Facet_handle>::iterator it = facets.begin(); it!= facets.end(); ++it){
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
  for(std::list<Halfedge_handle>::iterator it = interior_edges.begin();
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
}


void fill(Polyhedron& poly, Halfedge_handle it)
{
  std::cerr << "fill" << std::endl;
  if(! it->is_border()){
    return;
  }
  std::set<Facet_handle> facets;
  Polyline_3 P;
  Halfedge_around_facet_circulator circ(it), done(circ);
  do{
    P.push_back(circ);
    average_length(poly, circ->vertex());
  } while (++circ != done);
  P.push_back(circ);

  int n = P.size() - 1; // because the first and last point are equal
  std::vector<Weight> W(n*n,Weight(0,0));
  std::vector<int> lambda(n*n,-1);

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
  add_facets(P, lambda, 0, n-1, poly, facets);

  refine(poly, facets);
  std::cerr << "|facets| = " << facets.size() << std::endl;
}



int main()
{
  Polyhedron poly;
  std::cin >> poly;

  for(Halfedge_iterator it = poly.halfedges_begin(); it != poly.halfedges_end(); ++it){
    if(it->is_border()){
      fill(poly, it);
    }
  }
  std::cout.precision(20);
  std::cout << poly << std::endl;

  std::cerr << "done" << std::endl;
  return 0;
}
