#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/Erase_event.h>
#include <CGAL/Kinetic/Inexact_simulation_traits.h>
#include <CGAL/Kinetic/Exact_simulation_traits.h>
#include <CGAL/Kinetic/Insert_event.h>
#include <CGAL/Kinetic/Ref_counted.h>
#include <CGAL/Kinetic/listeners.h>
#include <CGAL/Kinetic/Sort.h>
#include <CGAL/Kinetic/Sort_visitor_base.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <list>
#include <map>

#ifdef _MSC_VER
#pragma warning(disable:4355) // complaint about using 'this' to
                                // initialize a member
#endif

template <class Arrangement>
struct Arrangement_visitor: public CGAL::Kinetic::Sort_visitor_base
{
  Arrangement_visitor(Arrangement *a):p_(a){}
  template <class Vertex_handle>
  void remove_vertex(Vertex_handle a) {
    p_->erase(a);
  }
  template <class Vertex_handle>
  void create_vertex(Vertex_handle a) {
    p_->insert(a);
  }
  template <class Vertex_handle>
  void after_swap(Vertex_handle a, Vertex_handle b) {
    p_->swap(a, b);
  }
  Arrangement *p_;
};



template <class TraitsT>
class Planar_arrangement:
  public CGAL::Kinetic::Sort<TraitsT,
			     Arrangement_visitor<Planar_arrangement<TraitsT> > > {
  typedef TraitsT Traits;
  typedef Planar_arrangement<TraitsT> This;
  typedef typename CGAL::Kinetic::Sort<TraitsT,
				       Arrangement_visitor<This> > Sort;
  typedef Arrangement_visitor<This> Visitor;
  typedef typename Traits::Active_points_1_table::Key Key;

public:
  typedef CGAL::Exact_predicates_inexact_constructions_kernel::Point_2 Approximate_point;
  typedef std::pair<int,int> Edge;
  typedef typename Sort::Vertex_handle Vertex_handle;

  // Register this KDS with the MovingObjectTable and the Simulator
  Planar_arrangement(Traits tr): Sort(tr, Visitor(this)), tr_(tr) {}

  Approximate_point vertex(int i) const
  {
    return approx_coords_[i];
  }

  size_t vertices_size() const
  {
    return approx_coords_.size();
  }

  typedef std::vector<Edge >::const_iterator Edges_iterator;
  Edges_iterator edges_begin() const
  {
    return edges_.begin();
  }
  Edges_iterator edges_end() const
  {
    return edges_.end();
  }

  void insert(Vertex_handle k) {
    last_points_[*k]=new_point(*k);
  }

  void swap(Vertex_handle a, Vertex_handle b) {
    int swap_point= new_point(*a);
    edges_.push_back(Edge(swap_point, last_points_[*a]));
    edges_.push_back(Edge(swap_point, last_points_[*b]));
    last_points_[*a]= swap_point;
    last_points_[*b]= swap_point;
  }

  void erase(Vertex_handle a) {
    edges_.push_back(Edge(last_points_[*a], new_point(*a)));
  }

  int new_point(typename Traits::Active_points_1_table::Key k) {
    double tv= CGAL::to_double(tr_.simulator_handle()->current_time());
    double dv= CGAL::to_double(tr_.active_points_1_table_handle()->at(k).x()(tv));
    approx_coords_.push_back(Approximate_point(tv, dv));
    return approx_coords_.size()-1;
  }

  std::vector<Approximate_point > approx_coords_;
  std::map<Key, int> last_points_;
  std::vector<Edge> edges_;
  Traits tr_;
};

template <class NT>
double snap(NT v) {
  double ival= std::floor(CGAL::to_double(v)*1024.0)/1024.0;
  if (ival > 10) return 20.0;
  else if (ival <-10) return 0.0;
  else return ival+10;
}


int main(int, char *[])
{
  typedef CGAL::Kinetic::Exact_simulation_traits Traits;
  typedef Traits::Kinetic_kernel::Point_1 Point;
  typedef Traits::Simulator::Time Time;
  typedef CGAL::Kinetic::Insert_event<Traits::Active_points_1_table> Insert_event;
  typedef CGAL::Kinetic::Erase_event<Traits::Active_points_1_table> Erase_event;
  typedef Planar_arrangement<Traits> Arrangement;

  Traits tr(0,1000000.0);
  Arrangement sort(tr);

  typedef Traits::Kinetic_kernel::Function_kernel::FT NT;

  Traits::Simulator::Handle sp= tr.simulator_handle();

  std::ifstream in("data/sweepline.input");

  int num=0;
  std::vector<std::pair<NT, NT> > extents;
  std::vector<Point> points;
  do {
    char buf[1000];
    in.getline(buf, 1000);
    if (!in) break;
    std::istringstream iss(buf);
    NT begin, end;
    Point pt;
    iss >> begin;
    iss >> end;
    iss >> pt;
    CGAL_assertion(begin < end);
    extents.push_back(std::make_pair(begin, end));
    points.push_back(pt);
    tr.simulator_handle()->new_event(Time(begin),
				      Insert_event(pt, tr.active_points_1_table_handle()));
    tr.simulator_handle()->new_event(Time(end),
				      Erase_event(Traits::Active_points_1_table::Key(num),
						  tr.active_points_1_table_handle()));
    ++num;
  } while (true);

  while (sp->next_event_time() < sp->end_time()) {
    sp->set_current_event_number(sp->current_event_number()+1);
    //std::cout << *sp << std::endl;
    //std::cout << sort << std::endl;
  }

  // output to metapost
  std::ofstream mpostfile("sweepline.mp");

  mpostfile << "beginfig(0)\n";
  mpostfile << "u:=20.0pt;\n";
  mpostfile << "pickup pencircle scaled 0.05u;\n";
  mpostfile << "%" << std::ios::fixed<< std::endl;

  for (Arrangement::Edges_iterator it = sort.edges_begin(); it != sort.edges_end(); ++it) {
    mpostfile << "draw(" << snap(sort.vertex(it->first).x()) << "u, "
	      << snap(sort.vertex(it->first).y()) << "u)";
    mpostfile << "--(" << snap(sort.vertex(it->second).x()) << "u, "
	      << snap(sort.vertex(it->second).y()) << "u) withcolor black;\n";
  }

  for (unsigned int i=0; i< points.size(); ++i) {
    mpostfile << "draw(";
    //std::cout << extents[i].first << " " << extents[i].second << std::endl;
    for (NT t= extents[i].first; t < extents[i].second; t += .01) {
      NT ntv= points[i].x()(t);
      double val= snap(ntv);

      mpostfile << snap(t) << "u, " << val << "u)\n--(";

    }
    mpostfile << snap(extents[i].second) << "u, "
	      << snap(points[i].x()(extents[i].second))
	      << "u) withcolor red;\n";
  }

  //for (double t =0; t<

  mpostfile << "endfig\n";
  return EXIT_SUCCESS;
}
