#include <CGAL/Arrangement_of_spheres_3/Slice_data_structure.h>
#include <map>
#include <set>
#include <vector>


Slice_data_structure::Slice_data_structure() {
  initialize();
  //std::cout << hds_.size_of_bfaces() << " " << hds_.size_of_halfedges() << " " << hds_.size_of_vertices() << std::endl;
  //hds_.clear();
}

Slice_data_structure::HDS& Slice_data_structure::hds() {
  return hds_;
}


bool Slice_data_structure::has_vertex(Face_const_handle fh, Vertex_const_handle vh) const {
  Halfedge_const_handle h= vh->halfedge();
  do {
    if (h->face() == fh) return true;
    h= h->opposite()->prev();
    CGAL_assertion(h->vertex() == vh);
  } while (h != vh->halfedge());
  return false;
}


void Slice_data_structure::audit_vertex(Vertex_const_handle v) const {
  Point pt = v->point();
  // check for tangency points
  std::vector<Curve> curves;
  HDS::Halfedge_const_handle c= v->halfedge();
  do {
    curves.push_back(c->curve());
    c= c->next()->opposite();
  } while (c != v->halfedge());
    
  CGAL_assertion(!curves.empty());
  CGAL_assertion(curves.size() != 1);
  if (curves.size() ==2) {
    CGAL_assertion(!curves[0].is_finite());
    CGAL_assertion(!curves[1].is_finite());
    CGAL_assertion(curves[0].is_rule());
    CGAL_assertion(curves[1].is_rule());
    CGAL_assertion(curves[0].is_vertical()
		   != curves[1].is_vertical());
  } else {
      
    //CGAL_assertion(curves.size() ==4);
    std::vector<bool> dirs(4, false);
    std::map<int,int> arcs;
    std::map<int,std::pair<int,bool> > rules;
    std::vector<Curve> ordered_arcs;
    for (unsigned int i=0; i< curves.size(); ++i){
      Curve c= curves[i];
      Curve nc= curves[(i+1)%curves.size()];
      int ind= c.index();
      int nind= nc.index();
      CGAL_assertion(c.is_rule() || nc.is_rule() || ind != nind);
      if (c.is_rule()){
	int bin=0;
	if (c.is_negative()) bin +=2;
	if (c.is_vertical()) bin +=1;
	CGAL_assertion(!dirs[bin]);
	dirs[bin]=true;
	
	if (rules.find(ind)== rules.end()){
	  rules[ind]= std::make_pair(0, c.is_vertical());
	} else {
	  CGAL_assertion(rules[ind].first==0);
	  ++rules[ind].first;
	  CGAL_assertion(rules[ind].second== c.is_vertical());
	}
      } else {
	if (arcs.find(ind) == arcs.end()){
	  arcs[ind]=1;
	} else {
	  ++arcs[ind];
	}
	ordered_arcs.push_back(c);
      }
    }
      
    for (std::map<int,int>::const_iterator it= arcs.begin();
	 it != arcs.end(); ++it){
      CGAL_assertion(it->second ==2);
    }
      
      
    // check that all arc pairs cross
    CGAL_assertion(ordered_arcs.size()%2==0);
    int half= ordered_arcs.size()/2;
    for (unsigned int i=0; i< ordered_arcs.size(); ++i){
      CGAL_assertion(ordered_arcs[i].index() 
		     == ordered_arcs[(i+half)%ordered_arcs.size()].index());
    }
  }
}



void Slice_data_structure::audit() const {
  std::set<HE_key> reachable;
  for (HDS::Vertex_const_iterator it= hds_.vertices_begin(); 
       it != hds_.vertices_end(); ++it){

    audit_vertex(it);

    HDS::Halfedge_const_handle c= it->halfedge();
    do {
      reachable.insert(c);
      c= c->next()->opposite();
    } while (c != it->halfedge());
  }

  for (HDS::Halfedge_iterator it= hds_.halfedges_begin(); it != hds_.halfedges_end(); ++it){
    HDS::Halfedge_const_handle h= it;
    HDS::Halfedge_const_handle ho= it->opposite();

    if (it->next() == HDS::Halfedge_handle()) {
      std::cerr<< "Invalid next for ";
      write(it, std::cerr) << std::endl;
      errors_.push_back(it->curve().index());
    }
    if (it->prev() == HDS::Halfedge_handle()) {
      std::cerr<< "Invalid prev for ";
      write(it, std::cerr) << std::endl;
      errors_.push_back(it->curve().index());
    } else if (it != it->prev()->next()){
      std::cerr<< "Invalid prev/next for ";
      write(it, std::cerr) << std::endl;
      errors_.push_back(it->curve().index());
    }
    if (reachable.find(h) == reachable.end()){
      std::cerr << "Non-vertex reachable halfedge ";
      write(it, std::cerr) << std::endl;
      errors_.push_back(it->curve().index());
    }
      
    CGAL_assertion(it->curve().is_valid());
    CGAL_assertion(it->curve().is_inside() != it->opposite()->curve().is_inside());
  }

  CGAL::HalfedgeDS_const_decorator<HDS> chds(hds_);
  if (!chds.is_valid(true, 3)) {
    std::cerr << "Not valid." << std::endl;
  }

    
}

void Slice_data_structure::set_is_building(bool tf) {
  if (tf==false) {
    for (std::map<Edge, Halfedge_handle>::iterator it= unmatched_hedges_.begin(); it != unmatched_hedges_.end(); ++it){
      //std::cout << "Searching for next for ";
      //write(it->second, std::cout) << std::endl;
      it->second->set_face(inf_);
      inf_->set_halfedge(it->second);
      Halfedge_handle c=it->second->opposite()->prev()->opposite();
      Vertex_handle v= it->second->vertex();
      while (c->prev() != Halfedge_handle()){
	//write( c, std::cout) << std::endl;
	Vertex_handle vo= c->opposite()->vertex();
	CGAL_assertion(v==vo);
	c= c->prev()->opposite();
      }
      //std::cout << "Found ";
      //write(c, std::cout) << std::endl;
      c->set_prev(it->second);
      it->second->set_next(c);
    }
    unmatched_hedges_.clear();
  }
}
  
std::ostream &Slice_data_structure::write(Halfedge_const_handle h, std::ostream &out) const {
  out << h->opposite()->vertex()->point() << " -- " << h->curve()
      << " -- " << h->vertex()->point();
  return out;
}

void Slice_data_structure::reserve(int nv, int ne, int nf) {
  hds_.reserve(nv, ne, nf);
}

void Slice_data_structure::initialize() {
  inf_= hds_.faces_push_back(HDS::Face());
}

void Slice_data_structure::clear() {
  hds_.vertices_clear();
  hds_.edges_clear();
  hds_.faces_clear();
  inf_=Face_handle();
  errors_.clear();
    
  unmatched_hedges_.clear();
  points_.clear();
  initialize();
}

//typedef std::pair<Point,Point>  ED;

Slice_data_structure::Vertex_handle Slice_data_structure::new_point(Point p) {
  CGAL_precondition(p.is_valid());
  //std::cout << "Creating point " << p << std::endl;
  HDS::Vertex v;
  v.point()=p;
  points_[p]=hds_.vertices_push_back(v);
  return points_[p];
}

Slice_data_structure::Halfedge_handle Slice_data_structure::new_hedge(Point s, Curve ff, Point f) {
  //std::cout << "Creating edge " << s << " -- " << ff << " -- " << f << std::endl;
  CGAL_precondition(points_.find(s) != points_.end());
  if (points_.find(f) == points_.end()) {
    new_point(f);
  }

  Halfedge_handle h;
  Edge ep(s, ff, f);
  if (unmatched_hedges_.find(ep) != unmatched_hedges_.end()){
    h= unmatched_hedges_[ep];
    unmatched_hedges_.erase(ep);
    CGAL_assertion(h->is_border());
    //std::cout << "matched" << std::endl;
  } else {
    //std::cout << "unmatched" << std::endl;
    h=hds_.edges_push_back(HDS::Halfedge(), HDS::Halfedge());
    h->set_curve(ff);
    h->opposite()->set_curve(ff.other_side());
    //h->set_inside(inside);
    //h->opposite()->set_inside(!inside);
    unmatched_hedges_[Edge(f, ff.other_side(), s)]= h->opposite();
  }
  points_[f]->set_halfedge(h);
  h->set_vertex(points_[f]);
  points_[s]->set_halfedge(h->opposite());
  h->opposite()->set_vertex(points_[s]);
  return h;
}


