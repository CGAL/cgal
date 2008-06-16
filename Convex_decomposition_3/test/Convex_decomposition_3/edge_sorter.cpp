#include<CGAL/Gmpz.h>
#include<CGAL/Homogeneous.h>
#include<CGAL/Convex_decomposition_3/Edge_sorter.h>

typedef CGAL::Gmpz RT;
typedef CGAL::Homogeneous<RT> Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::FT FT;

class LVertex {
  Point_3 pnt;
public:
  LVertex() {}
  LVertex(Point_3 p) : pnt(p) {}
  Point_3 point() { return pnt; }
};

class LEdge {
  LVertex* v;
  LEdge* t;
public:
  LEdge(LVertex* v_) : v(v_) {}
  ~LEdge() { if(v != NULL) delete v; v = NULL; }
  LVertex* source() { return v; }
  LEdge* twin() { return t; }
  void set_twin(LEdge* t_) { t = t_; }
  void reset_source() { v = NULL; }
};

template<typename Edget, typename Vertext>
class Report_new_vertex {
  
public:
  Report_new_vertex() {}
  Edget* operator()(Edget* e, Point_3 ip) {
    Edget* x = new Edget(new Vertext(ip));
    Edget* y = new Edget(new Vertext(ip));
    x->set_twin(e->twin());
    y->set_twin(e);
    e->twin()->set_twin(x);
    e->set_twin(y);
    return x;
  }
};

template<typename Iterator, typename Compare>
bool check_sorting(Iterator begin, Iterator end, Compare& compare)
{
  if(begin == end) return true;
  Iterator next = begin;
  for(++next; next != end; ++begin, ++next) {
    if((*next)->source()->point() < (*begin)->source()->point())
      return false;
  }
  return true;
}

typedef CGAL::Need_to_split<Kernel, std::less<Point_3> > SplitTest;
typedef CGAL::Generic_edge_sorter<Point_3, std::less<Point_3>,  std::less<FT>,
				  SplitTest, Report_new_vertex<LEdge, LVertex>, std::deque<LEdge*> > GesSM;

int main(int argc, char* argv[]) {


  std::deque<LEdge*> edges;
  LEdge* e0 = new LEdge(new LVertex(Point_3(-1,0,-1)));
  LEdge* t0 = new LEdge(new LVertex(Point_3(-1,0,1)));
  LEdge* e0b = new LEdge(new LVertex(Point_3(-1,0,-1)));
  LEdge* t0b = new LEdge(new LVertex(Point_3(-1,0,1)));
  LEdge* e1 = new LEdge(new LVertex(Point_3(-2,-1,0)));
  LEdge* t1 = new LEdge(new LVertex(Point_3(2,1,0)));
  LEdge* e2 = new LEdge(new LVertex(Point_3(2,-1,0)));
  LEdge* t2 = new LEdge(new LVertex(Point_3(2,1,0)));
  LEdge* e3 = new LEdge(new LVertex(Point_3(-2,-1,-1)));
  LEdge* t3 = new LEdge(new LVertex(Point_3(2,1,-1)));
  LEdge* e4 = new LEdge(new LVertex(Point_3(-2,-1,1)));
  LEdge* t4 = new LEdge(new LVertex(Point_3(2,1,1)));
  LEdge* e5 = new LEdge(new LVertex(Point_3(-1,-1,0)));
  LEdge* t5 = new LEdge(new LVertex(Point_3(3,1,0)));
  LEdge* e6 = new LEdge(e0->source());
  LEdge* t6 = new LEdge(new LVertex(Point_3(0,0,-10)));
  LEdge* e7 = new LEdge(new LVertex(Point_3(-1,0,-1)));
  LEdge* t7 = new LEdge(new LVertex(Point_3(3,0,1)));

  e0->set_twin(t0);
  t0->set_twin(e0);
  e0b->set_twin(t0b);
  t0b->set_twin(e0b);
  e1->set_twin(t1);
  t1->set_twin(e1);
  e2->set_twin(t2);
  t2->set_twin(e2);
  e3->set_twin(t3);
  t3->set_twin(e3);
  e4->set_twin(t4);
  t4->set_twin(e4);
  e5->set_twin(t5);
  t5->set_twin(e5);
  e6->set_twin(t6);
  t6->set_twin(e6);
  e7->set_twin(t7);
  t7->set_twin(e7);

  std::less<Point_3> Less;
  SplitTest st;
  Report_new_vertex<LEdge, LVertex> rnv;
  GesSM gesSM;

  edges.clear();
  edges.push_back(e2);
  edges.push_back(e0);
  gesSM(edges, st, rnv);
  CGAL_assertion(check_sorting(edges.begin(), edges.end(), Less));
  CGAL_assertion(edges.size() == 2);
  CGAL_assertion(edges[0]->source()->point() == Point_3(-1,0,-1));

  edges.clear();
  edges.push_back(e3);
  edges.push_back(e0);
  gesSM(edges, st, rnv);
  CGAL_assertion(check_sorting(edges.begin(), edges.end(), Less));
  CGAL_assertion(edges.size() == 2);
  CGAL_assertion(edges[0]->source()->point() == Point_3(-2,-1,-1));

  edges.clear();
  edges.push_back(e0);
  edges.push_back(e4);
  gesSM(edges, st, rnv);
  CGAL_assertion(check_sorting(edges.begin(), edges.end(), Less));
  CGAL_assertion(edges.size() == 2);
  CGAL_assertion(edges[0]->source()->point() == Point_3(-2,-1,1)); 

  edges.clear();
  edges.push_back(e5);
  edges.push_back(e1);
  gesSM(edges, st, rnv);
  CGAL_assertion(check_sorting(edges.begin(), edges.end(), Less));
  CGAL_assertion(edges.size() == 2);
  CGAL_assertion(edges[0]->source()->point() == Point_3(-2,-1,0)); 

  edges.clear();
  edges.push_back(e6);
  edges.push_back(e0);
  edges.push_back(e1);
  gesSM(edges, st, rnv);
  CGAL_assertion(check_sorting(edges.begin(), edges.end(), Less));
  CGAL_assertion(edges.size() == 4);
  CGAL_assertion(edges[2]->twin()->source()->point() == Point_3(-1,0,0));
  CGAL_assertion(edges[3]->source()->point() == Point_3(-1,0,0));
  CGAL_assertion(edges[2] == e0);
  CGAL_assertion(edges[3]->twin() == t0);
  LEdge* en0 = edges[2]->twin();
  LEdge* en1 = edges[3];

  edges.clear();
  edges.push_back(e0b);
  edges.push_back(e1);
  edges.push_back(e6);
  gesSM(edges, st, rnv);
  CGAL_assertion(check_sorting(edges.begin(), edges.end(), Less));
  CGAL_assertion(edges.size() == 4);
  CGAL_assertion(edges[1]->twin()->source()->point() == Point_3(-1,0,0));
  CGAL_assertion(edges[3]->source()->point() == Point_3(-1,0,0));
  LEdge* en2 = edges[1]->twin();
  LEdge* en3 = edges[3];

  edges.clear();
  edges.push_back(e1);
  edges.push_back(e7);
  gesSM(edges, st, rnv);
  CGAL_assertion(check_sorting(edges.begin(), edges.end(), Less));
  CGAL_assertion(edges.size() == 2);
  CGAL_assertion(edges[0]->source()->point() == Point_3(-2,-1,0)); 

  delete e0;
  delete e0b;
  delete e1;
  delete e2;
  delete e3;
  delete e4;
  delete e5;
  e6->reset_source();
  delete e6;
  delete e7;
  delete t0;
  delete t0b;
  delete t1;
  delete t2;
  delete t3;
  delete t4;
  delete t5;
  delete t6;
  delete t7;
  delete en0;
  delete en1;
  delete en2;
  delete en3;
}
