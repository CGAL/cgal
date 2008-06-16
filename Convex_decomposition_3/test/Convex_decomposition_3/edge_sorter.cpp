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
  ~LEdge() { delete v; }
  LVertex* source() { return v; }
  LEdge* twin() { return t; }
  void set_twin(LEdge* t_) { t = t_; }
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
  LEdge* e1 = new LEdge(new LVertex(Point_3(-2,-1,0)));
  LEdge* t1 = new LEdge(new LVertex(Point_3(2,1,0)));
  e0->set_twin(t0);
  t0->set_twin(e0);
  e1->set_twin(t1);
  t1->set_twin(e1);

  SplitTest st;
  Report_new_vertex<LEdge, LVertex> rnv;
  GesSM gesSM;

  edges.push_back(e0);
  edges.push_back(e1);
  gesSM(edges, st, rnv);
  std::less<Point_3> Less;
  CGAL_assertion(check_sorting(edges.begin(), edges.end(), Less));
  CGAL_assertion(edges.size() == 3);
  CGAL_assertion(edges[1]->twin()->source()->point() == Point_3(-1,0,0));
  CGAL_assertion(edges[2]->source()->point() == Point_3(-1,0,0));

  
}
