#ifndef KINETIC_REGULAR_TRIANGULATION_3_H
#define KINETIC_REGULAR_TRIANGULATION_3_H
#include <CGAL/basic.h>


#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/KDS/internal/Delaunay_triangulation_base_3.h>
#include <CGAL/KDS/Regular_triangulation_vertex_base_3.h>
#include <CGAL/KDS/Regular_triangulation_cell_base_3.h>
#include <CGAL/KDS/Regular_triangulation_visitor_base_3.h>
#include <CGAL/KDS/Listener.h>
#include <CGAL/KDS/Ref_counted.h>

#include <CGAL/KDS/Simulator_kds_listener.h>
#include <CGAL/KDS/Notifying_table_listener_helper.h>

CGAL_KDS_BEGIN_INTERNAL_NAMESPACE

template <class KD, class Root_stack, class VH>
class Regular_3_pop_event: public Delaunay_event_base_3<KD, Root_stack> {
  typedef Delaunay_event_base_3<KD, Root_stack>  P;
public:
  Regular_3_pop_event(const Root_stack &s,
		      const VH &vh,
		      KD *kd): P(s, kd), vh_(vh){
  }

  void process(const typename KD::Simulator::Time&){
    P::kdel_->pop(vh_);
  }

  VH vertex() const {
    return vh_;
  }

  void write(std::ostream &out) const {
    out << "Pop " << vh_->point();
  }

  /*virtual bool is_of_type(int tp) const {
    return (tp &type())!=0 || P::is_of_type(tp);
    }
    static int type() {
    return 2;
    }*/

  virtual ~Regular_3_pop_event(){};

protected:
  const VH vh_;
};

template <class K, class S, class VH>
std::ostream& operator<<(std::ostream &out, const Regular_3_pop_event<K,S,VH> &e){
  e.write(out);
  return out;
}











template <class KD, class Root_stack, class K, class Cell>
class Regular_3_non_vertex_event: public Delaunay_event_base_3<KD, Root_stack>  {
  typedef Delaunay_event_base_3<KD, Root_stack>  P;
public:

  Regular_3_non_vertex_event(const Root_stack &s,
			     const K &k,
			     const Cell &c,
			     KD *kd): P(s,kd), k_(k), cell_(c){
  }

  Regular_3_non_vertex_event(){}
 
  void write(std::ostream &out) const {
    out << "Nothing " << P::vh_->point();
  }

  /* virtual bool is_of_type(int tp) const {
     return (tp &type())!=0 || P::is_of_type(tp);
     }
     static int type() {
     return 4;
     }*/

  K point() const {return k_;}

  Cell cell_handle() const {return cell_;}

  virtual ~Regular_3_non_vertex_event(){};

protected:
  const K k_;
  const Cell cell_;
};









template <class KD, class Root_stack, class K, class Cell>
class Regular_3_move_event: public Regular_3_non_vertex_event<KD, Root_stack, K, Cell>  {
  typedef Regular_3_non_vertex_event<KD, Root_stack, K, Cell>  P;
public:
  Regular_3_move_event(const Root_stack &s,
		       const K &k,
		       const Cell &c,
		       int dir,
		       KD *kd): P(s,k, c, kd), dir_(dir){
  }

  void process(const typename KD::Simulator::Time&){
    P::kdel_->move(P::k_, P::cell_, dir_);
  }
 
  void write(std::ostream &out) const {
    out << "Move " << P::point();
  }

  /*virtual bool is_of_type(int tp) const {
    return (tp &type())!=0 || P::is_of_type(tp);
    }
    static int type() {
    return 8;
    }*/

  virtual ~Regular_3_move_event(){};

protected:
  int dir_;
};

template <class K,  class S, class KK, class C>
std::ostream& operator<<(std::ostream &out, const Regular_3_move_event<K,S,KK,C> &e){
  e.write(out);
  return out;
}












template <class KD, class Root_stack, class K, class Cell>
class Regular_3_push_event: public Regular_3_non_vertex_event<KD, Root_stack, K, Cell>  {
  typedef Regular_3_non_vertex_event<KD, Root_stack, K, Cell>  P;
public:

  Regular_3_push_event(const Root_stack &s,
		       const K &k,
		       const Cell &c,
		       KD *kd): P(s,k, c, kd){
  }

  void process(const typename KD::Simulator::Time&){
    P::kdel_->push(P::k_, P::cell_);
  }
 
 
  void write(std::ostream &out) const {
    out << "Push " << P::point();
  }
  /*virtual bool is_of_type(int tp) const {
    return (tp &type())!=0 || P::is_of_type(tp);
    }
    static int type() {
    return 16;
    }*/

  virtual ~Regular_3_push_event(){};
};

template <class K, class S, class KK, class C>
std::ostream& operator<<(std::ostream &out, const Regular_3_push_event<K,S,KK,C> &e){
  e.write(out);
  return out;
}


template <class Traits> 
struct Regular_triangulation_3_types {
  typedef typename Traits::Moving_point_table MPT;
  typedef typename Traits::Kinetic_kernel KK;
  typedef CGAL::Triangulation_cell_base_3<typename Traits::Instantaneous_kernel> CFB;

  /*typedef CGAL::Triangulation_cell_base_with_info_3<Delaunay_cache_3<MPT, KK>, 
    typename Traits::Instantaneous_kernel, CFB> CFBI;*/

  //typedef Triangulation_labeled_edge_cell_base_3<CFBI, typename Traits::Simulator::Event_key> TFB;
  //typedef Triangulation_labeled_facet_cell_base_3<TFB, typename Traits::Simulator::Event_key> FB;
  /*typedef CGAL::Triangulation_vertex_base_3<typename Traits::Instantaneous_kernel> CVB;
  typedef CGAL::Triangulation_vertex_base_with_info_3<typename Traits::Simulator::Event_key, 
						      typename Traits::Instantaneous_kernel, 
						      CVB> LVB;*/
  typedef CGAL::KDS::Regular_triangulation_cell_base_3<Traits> FB;
  typedef CGAL::KDS::Regular_triangulation_vertex_base_3<Traits> LVB;
  
  typedef CGAL::Triangulation_data_structure_3<LVB, FB> TDS;
 
  typedef CGAL::Regular_triangulation_3<typename Traits::Instantaneous_kernel, TDS> Default_triangulation;
  
  //friend class CGAL::Delaunay_triangulation_3<typename P::Instantaneous_kernel, TDS>;

};


CGAL_KDS_END_INTERNAL_NAMESPACE


CGAL_KDS_BEGIN_NAMESPACE

/*!
  redundant_cells_ maps each cell with redundant points to the ids of the points in that cell

  redundant_points_ maps each redundant point to a certificate
*/
template <class TraitsT, 
	  class VisitorT= Regular_triangulation_visitor_base_3,
	  class TriangulationT= typename internal::Regular_triangulation_3_types<TraitsT>::Default_triangulation>
class Regular_triangulation_3: 
  public Ref_counted<Regular_triangulation_3<TraitsT, VisitorT, TriangulationT> > {
private:
  typedef Regular_triangulation_3<TraitsT, VisitorT, TriangulationT> This;
  
public:
  typedef TraitsT Traits;
  typedef typename Traits::Moving_point_table::Key Point_key;

protected:
  typedef typename Traits::Moving_point_table MPT;
  typedef typename Traits::Simulator Simulator;
  typedef typename Traits::Simulator::Event_key Event_key;
  typedef typename Traits::Simulator::Time Time;
  
  typedef typename Traits::Simulator::Root_stack Root_stack;
  typedef TriangulationT Delaunay;
public:
  typedef internal::Delaunay_3_edge_flip_event<This, Root_stack, 
					       typename Delaunay::Edge> Edge_flip;
  friend class internal::Delaunay_3_edge_flip_event<This, Root_stack, 
						    typename Delaunay::Edge>;
  typedef internal::Delaunay_3_facet_flip_event<This, Root_stack, 
						typename Delaunay::Facet> Facet_flip;
  friend class internal::Delaunay_3_facet_flip_event<This, Root_stack, 
						typename Delaunay::Facet>;
  

  typedef internal::Regular_3_pop_event<This, Root_stack,
					typename Delaunay::Vertex_handle> Pop_event;
  friend class internal::Regular_3_pop_event<This, Root_stack,
					typename Delaunay::Vertex_handle>;
  typedef internal::Regular_3_non_vertex_event<This, Root_stack,
					       Point_key,
					       typename Delaunay::Cell_handle> Non_vertex_event;
  typedef internal::Regular_3_move_event<This, Root_stack,
					 Point_key,
					 typename Delaunay::Cell_handle> Move_event;
  
  friend class internal::Regular_3_move_event<This, Root_stack,
					 Point_key,
					 typename Delaunay::Cell_handle> ;
  typedef internal::Regular_3_push_event<This, Root_stack,
					 Point_key,
					 typename Delaunay::Cell_handle> Push_event;
  friend class internal::Regular_3_push_event<This, Root_stack,
					 Point_key,
					 typename Delaunay::Cell_handle> ;
  
protected:
  typedef std::multimap<typename Delaunay::Cell_handle,
			Point_key> RCMap;
  typedef std::map<Point_key, Event_key> RPMap;
  
  struct Base_traits;
  friend struct Base_traits;

  struct Base_traits: public TraitsT {
    typedef This Wrapper;
    typedef TriangulationT Triangulation;
    typedef typename This::Edge_flip Edge_flip;
    typedef typename This::Facet_flip Facet_flip;
    typedef typename TraitsT::Kinetic_kernel::Power_test_3 Side_of_oriented_sphere_3;
    typedef typename TraitsT::Kinetic_kernel::Weighted_orientation_3 Orientation_3;
    
    Side_of_oriented_sphere_3 side_of_oriented_sphere_3_object() const {
      return TraitsT::kinetic_kernel_object().power_test_3_object();
    }

    Orientation_3 orientation_3_object() const {
      return TraitsT::kinetic_kernel_object().weighted_orientation_3_object();
    }

    Base_traits(This *t, const TraitsT &tr): TraitsT(tr), wr_(t) {}

    Wrapper* wrapper_pointer() {
      return wr_;
    }
    const Wrapper* wrapper_pointer() const {
      return wr_;
    }

    Wrapper *wr_;
  };

  typedef internal::Delaunay_triangulation_base_3<Base_traits, VisitorT> KDel;

  
  typedef typename CGAL::KDS::Simulator_kds_listener<typename TraitsT::Simulator::Listener, This> Simulator_listener;
  friend  class CGAL::KDS::Simulator_kds_listener<typename TraitsT::Simulator::Listener, This>;
  typedef typename CGAL::KDS::Notifying_table_listener_helper<typename TraitsT::Moving_point_table::Listener, This> Moving_point_table_listener;
  friend class CGAL::KDS::Notifying_table_listener_helper<typename TraitsT::Moving_point_table::Listener, This>;

public:
  typedef VisitorT Visitor;
  
  Regular_triangulation_3(Traits tr, Visitor v= Visitor()): kdel_(Base_traits(this, tr), v),
							    siml_(tr.simulator_pointer(), this),
							    motl_(tr.moving_point_table_pointer(), this),
							    listener_(NULL){
  }

   const Visitor &visitor() const {
    return kdel_.visitor();
  }
  
  typedef TriangulationT Triangulation;
  const Triangulation& triangulation() const {
    return kdel_.triangulation();
  }

 

  struct Listener_core{
    typedef typename This::Pointer Notifier_pointer;
    typedef enum {TRIANGULATION} Notification_type;
  };
  friend class Listener<Listener_core>;
  typedef Listener<Listener_core> Listener;
  
protected:

  
  void set_listener(Listener *l){
    listener_= l;
  }
  Listener* listener() const {
    return listener_;
  }
  void audit_structure() const {
    if (!has_certificates()) return;
    for (typename Triangulation::Finite_vertices_iterator vit= triangulation().finite_vertices_begin();
	 vit != triangulation().finite_vertices_end(); ++vit){
      if (triangulation().degree(vit) == 4){
	Point_key k= vit->point(); 
	CGAL_assertion(vit->info() || !vit->point()); // it could be infinite
      } else {
	CGAL_assertion(!vit->info());
      }
      CGAL_assertion(redundant_points_.find(vit->point())== redundant_points_.end());
    }
   
  }
  
  void audit() const {
    CGAL_KDS_LOG(LOG_LOTS, "Verifying regular.\n");
    if (!has_certificates()) return;
    CGAL_KDS_LOG(LOG_LOTS, *this << std::endl);
    //P::instantaneous_kernel().set_time(P::simulator()->audit_time());
    kdel_.audit();
    audit_structure();
    Delaunay dt(kdel_.triangulation().geom_traits());
    std::vector<typename Delaunay::Vertex_handle> nvhs(kdel_.moving_object_table()->size());
    for (typename MPT::Keys_iterator it= kdel_.moving_object_table()->keys_begin();
	 it != kdel_.moving_object_table()->keys_end(); ++it){
      nvhs[it->index()] = dt.insert(*it);
    }
    CGAL_KDS_LOG(LOG_LOTS, "Done building." << std::endl);
    for (typename Delaunay::Finite_vertices_iterator fit= dt.finite_vertices_begin();
	 fit != dt.finite_vertices_end(); ++fit){
      std::vector<typename Delaunay::Vertex_handle> neighbors;
      dt.incident_vertices(fit, back_inserter(neighbors));
      std::vector<Point_key> neighborsp;
      for (unsigned int i=0; i< neighbors.size(); ++i){
	neighborsp.push_back(neighbors[i]->point());
      }
      
      std::sort(neighborsp.begin(), neighborsp.end());
      typename Delaunay::Vertex_handle kvh= kdel_.vertex_handle(fit->point());
      if (kvh == NULL){
	std::cerr << "Missing vertex from kinetic: " << fit->point() << std::endl;
      }
      std::vector<typename Delaunay::Vertex_handle> kneighbors;
      triangulation().incident_vertices(kvh, back_inserter(kneighbors));
      std::vector<Point_key> kneighborsp;
      for (unsigned int i=0; i< kneighbors.size(); ++i){
	kneighborsp.push_back(kneighbors[i]->point());
      }
      
      std::sort(kneighborsp.begin(), kneighborsp.end());
      for (unsigned int i=0; i< std::max(neighborsp.size(), kneighborsp.size()); ++i){
	Point_key nr;
	if (i < neighborsp.size()) nr= neighborsp[i];
	Point_key nk;
	if (i < kneighborsp.size()) nk= kneighborsp[i];
	if (nr != nk){
	  std::cerr << "Mismatched neighbors " << nr << nk << std::endl;
	}
      }
    }
    for (typename Delaunay::Finite_cells_iterator cit= dt.finite_cells_begin(); 
	 cit != dt.finite_cells_end(); ++cit){
      std::vector<Point_key> ks;
      for (unsigned int i=0; i<4; ++i){
	ks.push_back(cit->vertex(i)->point());
      }
      std::vector<typename Delaunay::Vertex_handle> kvhs;
      for (unsigned int i=0; i<4; ++i){
	kvhs.push_back(kdel_.vertex_handle(ks[i]));
      }
      typename Triangulation::Cell_handle h;
      int i,j,k,l;
      if (!triangulation().is_cell(kvhs[0], kvhs[1], kvhs[2], kvhs[3], h,i,j,k,l)){
	std::cerr << "Missing cell " << ks[0] << ks[1] << ks[2] << ks[3] << std::endl;
      }
    }

    for (typename Triangulation::Finite_cells_iterator cit= triangulation().finite_cells_begin(); 
	 cit != triangulation().finite_cells_end(); ++cit){
      std::vector<Point_key> ks;
      for (unsigned int i=0; i<4; ++i){
	ks.push_back(cit->vertex(i)->point());
      }
      std::vector<typename Delaunay::Vertex_handle> kvhs;
      for (unsigned int i=0; i<4; ++i){
	kvhs.push_back(nvhs[ks[i].index()]);
      }
      typename Triangulation::Cell_handle h;
      int i,j,k,l;
      if (!dt.is_cell(kvhs[0], kvhs[1], kvhs[2], kvhs[3], h,i,j,k,l)){
	std::cerr << "Missing cell " << ks[0] << ks[1] << ks[2] << ks[3] << std::endl;
      }
    }
  }
public:
  void write(std::ostream &out) const {
    if (triangulation().dimension() != 3) return;
    kdel_.write(out);
    out << "Redundant points: ";
    for (typename RPMap::const_iterator it= redundant_points_.begin(); it != redundant_points_.end();
	 ++it){
      out << it->first << " ";
    }
    out << std::endl;
    typename Delaunay::Cell_handle last;
    out << "Redundant cells: ";
    for (typename RCMap::const_iterator it= redundant_cells_.begin(); it != redundant_cells_.end();
	 ++it){
      if (it->first != last){
	last= it->first;
	internal::write_cell(last, out);
	out << ": ";
      } 
      out << it->second << " ";
    }
    out << std::endl;
  }

  void push(Point_key k, typename Triangulation::Cell_handle h) {
    kdel_.visitor().pre_push(k,h);
    CGAL_KDS_LOG(LOG_LOTS, "Pushing " << k << " into cell ");
    CGAL_KDS_LOG_WRITE(LOG_LOTS, internal::write_cell(h, LOG_STREAM));
    CGAL_KDS_LOG(LOG_LOTS, std::endl);
    redundant_points_[k]= Event_key();
    std::vector<Point_key> redundant;
    add_cell(h, redundant);
    redundant_points_.erase(k);

    for (unsigned int i=0; i< 4; ++i){
      typename Triangulation::Vertex_handle ov= h->vertex(i);
      if (ov->info()){
	kdel_.simulator()->delete_event(ov->info());
	ov->info()=  Event_key();
      }
    }

    typename Triangulation::Vertex_handle vh= kdel_.push_vertex(k, h);

    std::vector<typename Triangulation::Cell_handle> ic;
    triangulation().incident_cells(vh, back_inserter(ic));
    for (unsigned int i=0; i< redundant.size(); ++i){
      if (redundant[i] != k) {
	handle_redundant(redundant[i], ic.begin(), ic.end());
      }
    }
    vh->info()= make_certificate(vh);
    
    on_geometry_changed();
    kdel_.visitor().post_push(vh);
  }

  void pop(typename Triangulation::Vertex_handle vh){
    kdel_.visitor().pre_pop(vh);
    CGAL_KDS_LOG(LOG_LOTS, "Popping " << vh->point() << std::endl);
    std::vector<Point_key> redundant;
    Point_key k= vh->point();
    //add_cell(h, redundant);
    std::vector<typename Triangulation::Cell_handle> ic;
    triangulation().incident_cells(vh, back_inserter(ic));
    for (unsigned int i=0; i< ic.size(); ++i){
      add_cell(ic[i], redundant);
    }

    typename Triangulation::Cell_handle h= kdel_.pop_vertex(vh);

    for (unsigned int i=0; i< redundant.size(); ++i){
      bool success= handle_redundant(redundant[i],h);
      CGAL_postcondition(success);
      if (!success){
	std::cerr << "dropped a vertex in pop.\n";
	redundant_points_[k]=kdel_.simulator()->null_event();
      }
    }
    bool success=handle_redundant(k, h, true);
    if (!success){
      std::cerr << "dropped a vertex when popped.\n";
      redundant_points_[k]=kdel_.simulator()->null_event();
    }
    CGAL_postcondition(success);
    on_geometry_changed();
    kdel_.visitor().post_pop(k,h);
  }

  void move(Point_key k, typename Triangulation::Cell_handle h, int dir){
    kdel_.visitor().pre_move(k,h);
    CGAL_KDS_LOG(LOG_LOTS, "Moving " << k << " from ");
    CGAL_KDS_LOG_WRITE(LOG_LOTS, internal::write_cell(h, LOG_STREAM)); 
    CGAL_KDS_LOG(LOG_LOTS, " to ");
    CGAL_KDS_LOG_WRITE(LOG_LOTS, internal::write_cell(h->neighbor(dir), LOG_STREAM ));
    CGAL_KDS_LOG(LOG_LOTS, std::endl);
    typename Triangulation::Cell_handle neighbor = h->neighbor(dir);
    typename RCMap::iterator it= redundant_cells_.equal_range(h).first;
    while (it->second != k){
      ++it;
      CGAL_assertion(it != redundant_cells_.equal_range(h).second);
    }
    redundant_cells_.erase(it);
    bool hinf=false;
    for (unsigned int i=0; i<4; ++i){
      if (neighbor->vertex(i)== triangulation().infinite_vertex()){
	hinf=true;
	break;
      }
    }
    if (hinf){ 
      push(k, neighbor);
    } else {
      bool success= handle_redundant(k, neighbor, true);
      CGAL_postcondition(success);
      if (!success){
	// fix this later
	std::cerr << "dropped a vertex in move?.\n";
	redundant_points_[k]=kdel_.simulator()->null_event();
      }
    }
    kdel_.visitor().post_move(k,neighbor);
  }

  //! remove an object
  /*!
    See if it is redundant. If so, remove it from the list and delete its certificate.
    Otherwise, pass it along.
  */
  void erase(Point_key ){
    CGAL_assertion(0);
    on_geometry_changed();
  }

  void set(Point_key k){
    if (!kdel_.has_certificates()) return;
    if (kdel_.vertex_handle(k) != NULL){
      typename Triangulation::Vertex_handle vh= kdel_.change_vertex(k);
      std::vector<typename Triangulation::Cell_handle> ic;
      triangulation().incident_cells(vh, std::back_insert_iterator<std::vector<typename Triangulation::Cell_handle> >(ic));
      for (unsigned int i=0; i< ic.size(); ++i){
	typename Triangulation::Cell_handle h=ic[i];
	std::vector<Point_key> objs;
	for (typename RCMap::iterator it= redundant_cells_.equal_range(h).first;
	     it != redundant_cells_.equal_range(h).second; ++it){
	  objs.push_back(it->second);
	}
	redundant_cells_.erase(redundant_cells_.equal_range(h).first, redundant_cells_.equal_range(h).second);
	for (unsigned int j=0; j< objs.size(); ++j){
	  bool success= handle_redundant(objs[j], h);
	  if (!success){
	    // fix this later
	    std::cerr << "dropped a vertex where there is no way I should.\n";
	    redundant_points_[k]=kdel_.simulator()->null_event();
	  }
	}
      }
    } else {
      // Hack
      if (redundant_points_[k] == kdel_.simulator()->null_event()
	  || redundant_points_.find(k) == redundant_points_.end()
	  || !redundant_points_[k] ){
	// error handling code, this should not happen
	std::cerr << "Hack in handling lost point.\n";
	triangulation().geom_traits().set_time(kdel_.simulator()->rational_current_time());
	typename Delaunay::Cell_handle h= triangulation().locate(k);
	redundant_points_[k]= make_certificate(k, h);
      } else {
	//kdel_.simulator()->template event<Non_vertex_event>(redundant_points_[k]);
	typename Triangulation::Cell_handle h= kdel_.simulator()->template event<Non_vertex_event>(redundant_points_[k]).cell_handle();
	kdel_.simulator()->delete_event(redundant_points_[k]);
	redundant_points_.erase(k);
	redundant_points_[k]= make_certificate(k, h);
      }
    }
  }

 
  void insert(Point_key k){
    bool had_certificates= has_certificates();
    set_has_certificates(false);
    typename Triangulation::Cell_handle h= kdel_.triangulation().locate(k);
    typename KDel::Triangulation::Vertex_handle vh= kdel_.new_vertex_regular(k, h);
    if (vh==NULL){
      //if (redundant_points_.size() <= k.index()) redudant_points_.resize(k.index()+1);
      //redundant_points[k]=Event_key::null();
      //redundant_cells_.insert(typename RCMap::value_type(h, create_non_vertex_event(k, h)));
    }/* else {
	int deg = triangulation().degree(vh);
	if (deg==4){
	
	}
	}*/
    if (had_certificates) set_has_certificates(true);
    on_geometry_changed();
  }
public:
  void set_has_certificates(bool tf){
    if (tf == has_certificates()) return;
    if (tf==false){
      for (typename Triangulation::Finite_vertices_iterator vit= triangulation().finite_vertices_begin();
	   vit != triangulation().finite_vertices_end(); ++vit){
	if (vit->info()){
	  kdel_.simulator()->delete_event(vit->info());
	  vit->info()=  Event_key();
	}
      }
      for (typename RPMap::iterator it = redundant_points_.begin(); it != redundant_points_.end(); ++it){
	kdel_.simulator()->delete_event(it->second);
	it->second= Event_key();
      }
      redundant_points_.clear();
      redundant_cells_.clear();
      kdel_.set_has_certificates(false);
    } else {
      kdel_.set_has_certificates(true);
      if (kdel_.triangulation().dimension()==3) {
	// must be first so the vertex handles are set
	CGAL_KDS_LOG(LOG_LOTS, "Setting up certificates.\n");
	for (typename Triangulation::Finite_vertices_iterator vit= triangulation().finite_vertices_begin();
	     vit != triangulation().finite_vertices_end(); ++vit){
	  if (internal::has_degree_4(triangulation(), vit)){
	    vit->info()= make_certificate(vit);
	  }
	}
	for (typename Base_traits::Moving_point_table::Keys_iterator kit= kdel_.moving_object_table()->keys_begin(); 
	     kit != kdel_.moving_object_table()->keys_end(); ++kit){
	  typename Triangulation::Vertex_handle vh= kdel_.vertex_handle(*kit);
	  if (vh == NULL){
	    CGAL_KDS_LOG(LOG_LOTS, "On init " << *kit << " is redundant" << std::endl);
	    typename Triangulation::Cell_handle h= kdel_.triangulation().locate(*kit);
	    redundant_points_[*kit]= make_certificate(*kit, h);
	    redundant_cells_.insert(make_pair(h, *kit));	  
	  }
	}
      } else {
	CGAL_KDS_LOG(LOG_LOTS, "Triangulation does not have dimension 3.\n");
      }
    }
    
  }

  bool has_certificates() const {
    return kdel_.has_certificates();
  }
protected:
  //! also much check for vertex_events
  void flip(const typename Triangulation::Edge &edge){
    for (int i=0; i<2; ++i) {
      if (internal::has_degree_4(triangulation(), internal::vertex_of_edge(edge, i))){
	typename Delaunay::Vertex_handle vh= internal::vertex_of_edge(edge, i);
	Event_key k= vh->info();
	if (k){
	  //typename P::Time t= kdel_.extract_time(k);
	  //CGAL_assertion(t==simulator()->current_time());
	  CGAL_KDS_LOG(LOG_SOME, "diverting edge flip to pop.\n");
	  internal::set_edge_label(triangulation(), edge, kdel_.simulator()->null_event());
	  return;
	}
      }
    }
    std::vector<Point_key> redundant;

    typename Triangulation::Cell_circulator cc= triangulation().incident_cells(edge), ce=cc;
    do {
      add_cell(cc, redundant);
      ++cc;
    } while (cc != ce);

    typename Triangulation::Facet_circulator fc= triangulation().incident_facets(edge), fe=fc;
    do {
      typename Triangulation::Vertex_handle vh= internal::other_vertex(*fc, edge);
      if (vh->info()){
	kdel_.simulator()->delete_event(vh->info());
	vh->info()=Event_key();
      }
      ++fc;
    } while (fc != fe);

    typename Triangulation::Facet f= kdel_.flip(edge);

    std::vector<typename Triangulation::Cell_handle> ncells;
    ncells.push_back(f.first);
    ncells.push_back(f.first->neighbor(f.second));
    for (unsigned int i=0; i< redundant.size(); ++i){
      handle_redundant(redundant[i], ncells.begin(), ncells.end());
    }

    for (unsigned int i=0; i< 3; ++i){
      CGAL_assertion(!internal::has_degree_4(triangulation(), internal::vertex_of_facet(f, i)));
    }
    typename Triangulation::Vertex_handle v=f.first->vertex(f.second);
    if (internal::has_degree_4(triangulation(), v)){
      v->info()= make_certificate(v);
    }
    typename Triangulation::Vertex_handle vm= triangulation().mirror_vertex(f.first, f.second);
    if (internal::has_degree_4(triangulation(), vm)){
      vm->info()=make_certificate(vm);
    }

    /*typename Triangulation::Cell_handle cells[2];
      cells[0]= f.first;
      cells[1]= f.first->neighbor(f.second);
      for (typename std::vector<Point_key>::iterator it= redundant.begin(); 
      it != redundant.end(); ++it){
      typename P::Simulator::Root_stack s[2];
      s[0]= root_stack(*it, cells[0]);
      typename P::Time f[2];
      f[0]= s.next_time_negative();
      s[1]= root_stack(*it, cells[1]);
      f[1]= s1.next_time_negative();
      if (f[0] < f[1]){
      create_non_vertex_event(*it, cells[0], f[0], s[0]);
      } else {
      create_non_vertex_event(*it, cells[1], f[1], s[1]);
      }
      }*/
    

    on_geometry_changed();
  }

  void flip(const typename KDel::Facet &flip_facet){
    std::vector<Point_key> redundant;

    typename Triangulation::Cell_handle ch= flip_facet.first;
    add_cell(ch, redundant);
    
    typename Triangulation::Cell_handle och= flip_facet.first->neighbor(flip_facet.second);
    add_cell(och, redundant);

    if (flip_facet.first->vertex(flip_facet.second)->info()){
      kdel_.simulator()->delete_event(flip_facet.first->vertex(flip_facet.second)->info());
      flip_facet.first->vertex(flip_facet.second)->info()=Event_key();
    }
    if (triangulation().mirror_vertex(flip_facet.first, flip_facet.second)->info()){
      kdel_.simulator()->delete_event(triangulation().mirror_vertex(flip_facet.first, flip_facet.second)->info());
      triangulation().mirror_vertex(flip_facet.first, flip_facet.second)->info()= Event_key();
    }
  
    typename Triangulation::Edge edge=  kdel_.flip(flip_facet);

    std::vector<typename Triangulation::Cell_handle> ncells;
    typename Triangulation::Cell_circulator cc= triangulation().incident_cells(edge), ce= cc;
    do {
      ncells.push_back(cc);
      ++cc;
    } while (cc != ce);

    for (unsigned int i=0; i< redundant.size(); ++i){
      handle_redundant(redundant[i], ncells.begin(), ncells.end());
    }

    typename Triangulation::Facet_circulator fc= triangulation().incident_facets(edge), fe= fc;
    do {
      typename Triangulation::Vertex_handle vh= internal::other_vertex(*fc, edge);
      if (internal::has_degree_4(triangulation(), vh)){
	vh->info()=make_certificate(vh);
      }
      ++fc;
    } while (fc != fe);
    on_geometry_changed();
  }
  
  bool handle_redundant(Point_key k, typename Triangulation::Cell_handle h, bool must_handle=false){
    CGAL_KDS_LOG(LOG_LOTS, "Handle redundant " << k << " ") ;
    CGAL_KDS_LOG_WRITE(LOG_LOTS, internal::write_cell( h, LOG_STREAM));
    CGAL_KDS_LOG(LOG_LOTS, std::endl);
    unsigned int i=0;
    if (!must_handle) {
      for (i=0; i< 4; ++i){
	typename Triangulation::Facet f(h, i);
	typename Base_traits::Simulator::Function_kernel::Function cf= kdel_.orientation_object()(point(vertex_of_facet(f,0)->point()),
												  point(vertex_of_facet(f,1)->point()),
												  point(vertex_of_facet(f,2)->point()),
												  point(k));
	typename Base_traits::Simulator::Function_kernel::Sign_at sar
	  = kdel_.simulator()->function_kernel_object().sign_at_object(cf);
	CGAL::Sign sn = CGAL::sign(sar(kdel_.simulator()->current_time()));
	
#ifndef NDEBUG
	{
	  if ( sn != CGAL::sign(sar(to_double(kdel_.simulator()->current_time())))) {
	    std::cerr <<"Difference of opinion on sign at" << std::endl;
	    std::cerr <<"Real root " <<kdel_.simulator()->current_time() << std::endl;
	    std::cerr <<"Approximation " << to_double(kdel_.simulator()->current_time()) << std::endl;
	    std::cerr <<"Polynomial " << cf << std::endl;
	  }
	  
	  typename Base_traits::Simulator::Function_kernel::Sign_at csar
	    = kdel_.simulator()->function_kernel_object().sign_at_object(kdel_.orientation_object()(point(internal::vertex_of_facet(f,0)->point()),
												    point(internal::vertex_of_facet(f,1)->point()),
												    point(internal::vertex_of_facet(f,2)->point()),
												    point(f.first->vertex(f.second)->point())));
	  CGAL::Sign csn = CGAL::sign(csar(kdel_.simulator()->current_time()));
	  CGAL_assertion(csn==CGAL::POSITIVE);
	}
#endif
	
	if (sn == CGAL::NEGATIVE) {
	  CGAL_KDS_LOG(LOG_LOTS, "rejected because of side ");
	  CGAL_KDS_LOG_WRITE(LOG_LOTS, internal::write_facet(f, LOG_STREAM));
	  CGAL_KDS_LOG(LOG_LOTS, std::endl << internal::vertex_of_facet(f,0)->point() << " : " 
		       << point(internal::vertex_of_facet(f,0)->point()) << std::endl);
	  CGAL_KDS_LOG(LOG_LOTS, internal::vertex_of_facet(f,1)->point() << " : " 
		       << point(internal::vertex_of_facet(f,1)->point()) << std::endl);
	  CGAL_KDS_LOG(LOG_LOTS, internal::vertex_of_facet(f,2)->point() << " : " 
		       << point(internal::vertex_of_facet(f,2)->point()) << std::endl);
	  CGAL_KDS_LOG(LOG_LOTS, k << " : " << point(k) << std::endl);
	  CGAL_KDS_LOG(LOG_LOTS, cf << " : " << kdel_.simulator()->current_time() << std::endl << std::endl);

	  break;
	}
      }
    }
    if (i == 4 || must_handle){
      redundant_points_[k]=make_certificate(k, h);
      redundant_cells_.insert(typename RCMap::value_type(h, k));
      return true;
    } else {
      if (must_handle){
	std::cerr << "Time is " << kdel_.simulator()->current_time() << std::endl;
	for (int i=0; i< 4; ++i){
	  typename Triangulation::Facet f(h, i);
	  std::cerr << "Facet " << i << std::endl;
	  std::cerr << internal::vertex_of_facet(f,0)->point() << " : " << point(internal::vertex_of_facet(f,0)->point()) << std::endl;
	  std::cerr << internal::vertex_of_facet(f,1)->point() << " : " << point(internal::vertex_of_facet(f,1)->point()) << std::endl;
	  std::cerr << internal::vertex_of_facet(f,2)->point() << " : " << point(internal::vertex_of_facet(f,2)->point()) << std::endl;
	  std::cerr << k << " : " << point(k)<<std::endl;
	}
	
      }
      return false;
    } 
  }

  template <class It>
  void handle_redundant(Point_key k, It beg, It end){
    for (It cur=beg; cur != end; ++cur){
      //bool is_ok=true;
      bool hinf=false;
      for (unsigned int i=0; i< 4; ++i){
	if (!(*cur)->vertex(i)->point()) {
	  hinf=true;
	  break;
	}
      }
      if (hinf) continue;

      bool handled= handle_redundant(k, *cur);
      if (handled) return;
    }
    std::cerr << "None of the cells could handle vertex " << k << std::endl;
    std::cerr << "Supplied cells are: " << std::endl;
    for (It cur=beg; cur != end; ++cur){
      std::cerr << (*cur)->vertex(0)->point() 
		<< " " << (*cur)->vertex(1)->point() 
		<< " " << (*cur)->vertex(2)->point() 
		<< " " << (*cur)->vertex(3)->point() << std::endl;
    }

    for (typename Triangulation::Finite_cells_iterator it = kdel_.triangulation().finite_cells_begin(); 
	 it != kdel_.triangulation().finite_cells_end(); ++it){
      bool hinf=false;
      for (unsigned int i=0; i< 4; ++i){
	if (!(it->vertex(i)->point())) {
	  hinf=true;
	  break;
	}
      }
      if (hinf) continue;

      bool handled= handle_redundant(k, *beg);
      if (handled) {
	std::cerr << "A full search found : ";
	std::cerr << (it)->vertex(0)->point() 
		  << " " << (it)->vertex(1)->point() 
		  << " " << (it)->vertex(2)->point() 
		  << " " << (it)->vertex(3)->point() << std::endl;
	return;
      }
    }
    CGAL_postcondition(0);
    std::cerr << "dropped a vertex in handle.\n";
    redundant_points_[k]=kdel_.simulator()->null_event();
  }

  template <class V>
  std::back_insert_iterator<V> make_inserter(V &v){
    
    return std::back_insert_iterator<V>(v);
  }

  template <class A, class B>
  std::pair<A, B> make_pair(const A &a, const B&b) const {
    return std::pair<A, B>(a, b);
  }

  Event_key make_certificate(typename Triangulation::Vertex_handle vh){
    CGAL_precondition( internal::has_degree_4(triangulation(), vh));
    
    typename Triangulation::Cell_handle ch= vh->cell();
    typename Triangulation::Facet f(ch, ch->index(vh));
    std::vector<typename Triangulation::Vertex_handle> n(4);
    if (vh== triangulation().infinite_vertex()) return kdel_.simulator()->null_event();
    for (int i=0; i<3; ++i){
      n[i]= internal::vertex_of_facet(f,i);
      if (n[i]== triangulation().infinite_vertex()) return kdel_.simulator()->null_event();
    }
    int ind= (f.second+1)%4;
    n[3] = triangulation().mirror_vertex(ch, ind); // some vertex on facet
    if (n[3]== triangulation().infinite_vertex()) return kdel_.simulator()->null_event();
   
    CGAL_KDS_LOG(LOG_LOTS, "Making D4 certificate for " << n[0]->point() << n[1]->point() 
		 << n[2]->point() << n[3]->point() << " around " << vh->point() << std::endl);

    // hack, kill incident edge events, should co-opt them
    std::vector<typename Triangulation::Cell_handle> cells;
    triangulation().incident_cells(vh, make_inserter(cells));
    /*for (unsigned int i=0; i< cells.size(); ++i){
      typename Triangulation::Cell_handle ch= cells[i];
      int vi= ch->index(vh);
      for (int j=0; j<4; ++j){
      if (j==vi) continue;
      typename Triangulation::Edge e(ch, vi, j);
      kdel_.suppress_event(e);
      }
      }*/

    //! The order is switched to invert the predicate since we want it to fail when it goes outside
    typename Simulator::Root_stack s
      = kdel_.simulator()->root_stack_object(kdel_.power_test_object()(point(n[1]->point()),
									    point(n[0]->point()),
									    point(n[2]->point()),
									    point(n[3]->point()),
									    point(vh->point())));
    if (!s.empty()){
      Time t= s.top();
      s.pop();
      return kdel_.simulator()->new_event(t, Pop_event(s, vh, this));
    } else {
      return kdel_.simulator()->null_event();
    }
  }


  Event_key make_certificate(Point_key k, typename Triangulation::Cell_handle h){
    typename Simulator::Root_stack ps
      = kdel_.simulator()->root_stack_object(kdel_.power_test_object()(point(h->vertex(0)->point()),
									    point(h->vertex(1)->point()),
									    point(h->vertex(2)->point()),
									    point(h->vertex(3)->point()),
									    point(k)));
    Time pst;
    if (!ps.empty()) pst = ps.top();
    else pst= std::numeric_limits<Time>::infinity();

    int first=0;
    for (unsigned int i=0; i< 4; ++i){
      typename Triangulation::Facet f(h, i);
      // order matters
      typename Simulator::Root_stack cs
	= kdel_.simulator()->root_stack_object(kdel_.orientation_object()(point(internal::vertex_of_facet(f,0)->point()),
									       point(internal::vertex_of_facet(f,1)->point()),
									       point(internal::vertex_of_facet(f,2)->point()),
									       point(k)));
      if (!cs.empty() && cs.top() < pst){
	pst= cs.top();
	ps=cs;
	first= i+1;
      }
    }
    if (pst != std::numeric_limits<Time>::infinity()){
      if (first==0 ){
	CGAL_KDS_LOG(LOG_LOTS, "Making push certificate for " << k << std::endl);
	ps.pop();
	return kdel_.simulator()->new_event(pst, Push_event(ps, k, h, this));
      } else {
	CGAL_KDS_LOG(LOG_LOTS, "Making move certificate for " << k << std::endl);
	ps.pop();
	return kdel_.simulator()->new_event(pst, Move_event(ps, k, h, first-1, this));
      }
    } else {
      return kdel_.simulator()->null_event();
    }
  }

  void on_geometry_changed(){
    if (listener_!= NULL){
      listener_->new_notification(Listener::TRIANGULATION);
    }
    CGAL_KDS_LOG(LOG_LOTS, *this);
    audit_structure();
  }

  void add_cell(typename Triangulation::Cell_handle cc, std::vector<Point_key> &redundant){
    std::pair<typename RCMap::iterator, typename RCMap::iterator> ip= redundant_cells_.equal_range(cc);
    for (typename RCMap::iterator it = ip.first;
	 it != ip.second; ++it){
      CGAL_KDS_LOG(LOG_LOTS, it->second << " is redundant " << std::endl);
      redundant.push_back(it->second);
      if (redundant_points_[it->second]) {
	kdel_.simulator()->delete_event(redundant_points_[it->second]);
	redundant_points_[it->second]= Event_key();
      }
    }
    redundant_cells_.erase(ip.first, ip.second);
  }

  typename MPT::Data point(Point_key k) {
    return kdel_.moving_object_table()->at(k);
  }

  KDel kdel_;
  Simulator_listener siml_; 
  Moving_point_table_listener motl_; 
  Listener *listener_;
  RPMap redundant_points_;
  RCMap redundant_cells_;
  //typename P::Instantaneous_kernel::Orientation_3 po_;
  // typename P::Kinetic_kernel::Weighted_orientation_3 por_;
};

template <class Traits, class Triang, class Visit>
std::ostream &operator<<(std::ostream &out, const Regular_triangulation_3<Traits, Triang, Visit> &rt){
  rt.write(out);
  return out;
}

CGAL_KDS_END_NAMESPACE

#endif
