#ifndef CGAL_KDS_KINETIC_DELAUNAY_3_H
#define CGAL_KDS_KINETIC_DELAUNAY_3_H

#include <CGAL/KDS/basic.h>

#include <CGAL/KDS/internal/Delaunay_triangulation_base_3.h>

#include <CGAL/KDS/Listener.h>
#include <CGAL/KDS/Ref_counted.h>

// Triangulations
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>


// Local helpers
#include <CGAL/KDS/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/KDS/Simulator_kds_listener.h>
#include <CGAL/KDS/Notifying_table_listener_helper.h>
#include <CGAL/KDS/Delaunay_triangulation_visitor_base_3.h>


CGAL_KDS_BEGIN_INTERNAL_NAMESPACE

template <class Traits> 
struct Delaunay_triangulation_3_types {
  typedef typename Traits::Moving_point_table MPT;
  typedef typename Traits::Kinetic_kernel KK;
  typedef CGAL::KDS::Delaunay_triangulation_cell_base_3<Traits> CFBI;
  /*typedef CGAL::Triangulation_cell_base_with_info_3<Delaunay_cache_3<MPT, KK>, 
    typename Traits::Instantaneous_kernel, CFB> CFBI;*/
  typedef CGAL::Triangulation_vertex_base_3<typename Traits::Instantaneous_kernel> CVB;
  typedef CGAL::Triangulation_data_structure_3<CVB, CFBI> TDS;
 
  typedef CGAL::Delaunay_triangulation_3<typename Traits::Instantaneous_kernel, TDS> Default_triangulation;
  
  //friend class CGAL::Delaunay_triangulation_3<typename P::Instantaneous_kernel, TDS>;

};

CGAL_KDS_END_INTERNAL_NAMESPACE

CGAL_KDS_BEGIN_NAMESPACE


//! A 3D kinetic Delaunay triangulation.
template <class TraitsT, 
	  class Visitor= Delaunay_triangulation_visitor_base_3,
	  class TriangulationT=typename internal::Delaunay_triangulation_3_types<TraitsT>::Default_triangulation>
class Delaunay_triangulation_3: public Ref_counted<Delaunay_triangulation_3<TraitsT, Visitor, TriangulationT> > {
private:
  typedef Delaunay_triangulation_3<TraitsT, Visitor, TriangulationT> This;

  typedef internal::Delaunay_3_edge_flip_event<This,
					       typename TraitsT::Simulator::Root_stack, 
					       typename TriangulationT::Edge> Edge_flip;
  friend class internal::Delaunay_3_edge_flip_event<This,
						    typename TraitsT::Simulator::Root_stack, 
						    typename TriangulationT::Edge>;
  
  typedef internal::Delaunay_3_facet_flip_event<This, typename TraitsT::Simulator::Root_stack, 
						typename TriangulationT::Facet> Facet_flip;
  friend class  internal::Delaunay_3_facet_flip_event<This,
						      typename TraitsT::Simulator::Root_stack, 
						      typename TriangulationT::Facet>;

  struct Base_traits: public TraitsT {
    typedef This Wrapper;
    typedef TriangulationT Triangulation;
    typedef typename This::Edge_flip Edge_flip;
    typedef typename This::Facet_flip Facet_flip;
    typedef typename TraitsT::Kinetic_kernel::Side_of_oriented_sphere_3 Side_of_oriented_sphere_3;
    typedef typename TraitsT::Kinetic_kernel::Orientation_3 Orientation_3;
    
    Side_of_oriented_sphere_3 side_of_oriented_sphere_3_object() const {
      return TraitsT::kinetic_kernel_object().side_of_oriented_sphere_3_object();
    }

    Orientation_3 orientation_3_object() const {
      return TraitsT::kinetic_kernel_object().orientation_3_object();
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

  typedef internal::Delaunay_triangulation_base_3<Base_traits, Visitor> KDel;
  typedef typename TraitsT::Moving_point_table::Key Point_key;
    
  struct Listener_core{
    typedef typename This::Pointer Notifier_pointer;
    typedef enum {TRIANGULATION} Notification_type;
  };

  
  typedef typename CGAL::KDS::Simulator_kds_listener<typename TraitsT::Simulator::Listener, This> Simulator_listener;
  friend  class CGAL::KDS::Simulator_kds_listener<typename TraitsT::Simulator::Listener, This>;
  typedef typename CGAL::KDS::Notifying_table_listener_helper<typename TraitsT::Moving_point_table::Listener, This> Moving_point_table_listener;
  friend class CGAL::KDS::Notifying_table_listener_helper<typename TraitsT::Moving_point_table::Listener, This>;

public:
  //! Initialize it. 
  Delaunay_triangulation_3(TraitsT tr, Visitor v= Visitor()): kdel_(Base_traits(this, tr), v),
							      siml_(tr.simulator_pointer(), this),
							      motl_(tr.moving_point_table_pointer(), this),
							      listener_(NULL){
  }


  //! The type of the underlying triangulation
  typedef TriangulationT Triangulation;
  //! access the underlying triangulation
  const Triangulation& triangulation() const {
    return kdel_.triangulation();
  }

  Visitor& visitor() {
    return kdel_.visitor();
  }

  const Visitor& visitor() const {
    return kdel_.visitor();
  }

  friend class Listener<Listener_core>;
  //! This listener allows a class to know when the triangulation changes.
  /*!
    There is one notifaction type, TRIANGULATION.
  */
  typedef Listener<Listener_core> Listener;

  

  void write(std::ostream &out) const {
    kdel_.write(out);
  }

  void set_listener(Listener *l){
    listener_= l;
  }
  //! make the structure have or not have certificates
  void set_has_certificates(bool tf) {
    kdel_.set_has_certificates(tf);
  }
  Listener* listener() const {
    return listener_;
  }
  void audit() const {
    kdel_.audit();
  }
  
  //! true if the structure has certificates
  bool has_certificates() const {
    return kdel_.has_certificates();
  }

  void erase(Point_key k){
    kdel_.delete_vertex(k);
    on_geometry_changed();
  }

  void set(Point_key k){
    kdel_.change_vertex(k);
  }

 
  void insert(Point_key k){
    kdel_.new_vertex(k);
    /*if (kdel_.triangulation()->dimension() ==3){
      kdel_.set_has_certificates(true);
      }*/
    on_geometry_changed();
  }

  void flip(const typename KDel::Edge &edge){
    kdel_.flip(edge);
    on_geometry_changed();
  }

  void flip(const typename KDel::Facet &flip_facet){
    kdel_.flip(flip_facet);
    on_geometry_changed();
  }

  void on_geometry_changed(){
    if (listener_!= NULL){
      listener_->new_notification(Listener::TRIANGULATION);
    }
  }


  KDel kdel_;
  Simulator_listener siml_; 
  Moving_point_table_listener motl_; 
  Listener *listener_;
};

CGAL_KDS_END_NAMESPACE

#endif

