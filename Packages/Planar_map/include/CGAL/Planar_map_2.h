// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.4-I-65 $
// release_date  : $CGAL_Date: 2002/03/19 $
//
// file          : include/CGAL/Planar_map_2.h
// package       : Planar_map (5.87)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Iddo Hanniel      <hanniel@math.tau.ac.il>
//                 Eyal Flato        <flato@post.tau.ac.il>
//                 Oren Nechushtan   <theoren@math.tau.ac.il>
//                 Eti Ezra          <estere@post.tau.ac.il>
//                 Shai Hirsch       <shaihi@post.tau.ac.il>
//                 Eugene Lipovetsky <eug@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_PLANAR_MAP_2_H
#define CGAL_PLANAR_MAP_2_H

/*! \file
 * The implementation of the Planar_map_2<Dcel,Traits> class.
 */

#include <CGAL/Planar_map_2/Planar_map_misc.h>
#include <CGAL/Planar_map_2/Pm_change_notification.h>
#include <CGAL/Topological_map.h>

#ifndef CGAL_NO_PM_DEFAULT_POINT_LOCATION
#include <CGAL/Pm_default_point_location.h>
#include <CGAL/Pm_walk_along_line_point_location.h>
#include <CGAL/Pm_naive_point_location.h>
#include <CGAL/Pm_point_location_base.h>
#endif // CGAL_NO_PM_DEFAULT_POINT_LOCATION

// for solving the dynamic cast in the copy constructor, 
// these lines will be removed after writing 
// copy construtor for point location.
#include <CGAL/Pm_walk_along_line_point_location.h>
#include <CGAL/Pm_naive_point_location.h>

// default bounding box for finite curves
#include <CGAL/Pm_unbounding_box.h>

// default bounding box for infinite curves
#include <CGAL/Pm_dynamic_open_bounding_box.h>

#include <CGAL/IO/Pm_file_scanner.h>

#include <list>

CGAL_BEGIN_NAMESPACE

////////////////////////////////////////////////////////////////////////////
//      PLANAR_MAP_2   

/*! An object $pm$ of the class Planar_map_2<Dcel,Traits> is the planar
 * subdivision induced by a set of $x$-monotone curves such that no curve
 * intersects the interior of any other curve.
 */
template <class PlanarMapDcel_2, class PlanarMapTraits_2> 
class Planar_map_2 : public Topological_map<PlanarMapDcel_2>
{
public:
  typedef PlanarMapDcel_2                       Dcel;
  typedef PlanarMapTraits_2                     Traits;
  typedef Planar_map_2<Dcel,Traits>             Self;
  typedef Planar_map_traits_wrap<Traits>        Traits_wrap;
  typedef typename Traits::X_curve              X_curve_2;
  typedef typename Traits::Point                Point_2;
  
  typedef Topological_map<Dcel> TPM;
  typedef typename TPM::Vertex_iterator         Vertex_iterator;
  typedef typename TPM::Halfedge_iterator       Halfedge_iterator;
  typedef typename TPM::Edge_iterator           Edge_iterator;

  typedef typename TPM::Face_iterator           Face_iterator;
  typedef typename TPM::Vertex_const_iterator   Vertex_const_iterator;
  typedef typename TPM::Halfedge_const_iterator Halfedge_const_iterator;
  typedef typename TPM::Edge_const_iterator     Edge_const_iterator;

  typedef typename TPM::Face_const_iterator     Face_const_iterator;
  typedef typename TPM::Vertex_handle           Vertex_handle;
  typedef typename TPM::Vertex_const_handle     Vertex_const_handle;
  typedef typename TPM::Halfedge_handle         Halfedge_handle;
  typedef typename TPM::Face_handle             Face_handle;
  typedef typename TPM::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename TPM::Face_const_handle       Face_const_handle;
  typedef typename TPM::Halfedge_around_vertex_circulator
                                             Halfedge_around_vertex_circulator;
  typedef typename TPM::Halfedge_around_vertex_const_circulator 
                                       Halfedge_around_vertex_const_circulator;
  typedef typename TPM::Holes_iterator          Holes_iterator;
  typedef typename TPM::Holes_const_iterator    Holes_const_iterator;
  typedef typename TPM::Ccb_halfedge_const_circulator 
                                                Ccb_halfedge_const_circulator;
  typedef typename TPM::Ccb_halfedge_circulator Ccb_halfedge_circulator;
  
  typedef typename TPM::Size                    Size;
  
  typedef Pm_point_location_base<Self>          Point_location_base;
  typedef Pm_bounding_box_base<Self>            Bounding_box_base;
  
  typedef Pm_change_notification<Self>          Change_notification;

  // Obsolete, for backward compatability
  typedef Point_2                               Point;
  typedef X_curve_2                             X_curve;

  // Implementation Types
  // --------------------
protected:
  typedef std::list<X_curve_2>                  X_curve_2_container;
  
public:
  typedef enum { 
    VERTEX = 1, 
    EDGE, 
    FACE , 
    UNBOUNDED_VERTEX, 
    UNBOUNDED_EDGE, 
    UNBOUNDED_FACE
  } Locate_type;

  //! A parameter-less constructor.
  Planar_map_2();

  //! A copy constructor.
  Planar_map_2(const Self & pm);

  //! A destructor
  virtual ~Planar_map_2();
    
  //! A constructor given a point-location parameter.
  Planar_map_2(Point_location_base * pl_ptr);
  
  //! A constructor given a copy traits, point-location, and
  //! Bounding-box parameters.
  Planar_map_2(const Traits & tr_, Point_location_base * pl_ptr, 
               Bounding_box_base * bb_ptr);
  
  //! A constructor given a traits, point-location, and Bounding-box
  //! parameters. set NULLs for defaults
  Planar_map_2(Traits_wrap * tr_ptr, Point_location_base * pl_ptr, 
               Bounding_box_base * bb_ptr);

  //! reads a planar map.
  bool read(std::istream &in);

  //! reads a planar map.
  template <class Scanner>
  bool read(std::istream & in, Scanner & scanner)
  {
    clear(); 
    return scan_planar_map(scanner);
  }

  //! inserts a new edge into the map. The new edge is represented by
  //! a given curve.

  /*! insert() inserts a new edge into the map. The new edge is represented by
   * a given curve.
   * \param cv the given curve
   * \param en the notification class. It's default value is NULL, which
   * implies no notification on insertion.
   * \return a handle to a new halfedge represented by the given curve.
   */
  Halfedge_handle insert(const X_curve_2     & cv, 
                         Change_notification * en = NULL);

  //! inserts a set of new edges into the map. The new edges are
  //! reprsented by a given range of cirves.

  /*! insert() iterates through a given range of curves, inserting a new couple
   * of halfedges into the planar map for each curve.
   * \param begin the input iterator that points to the first curve in the
   * range.
   * \param end the input past-the-end iterator of the range.
   * \param en the notification class.
   * \return a handle to the last halfedge inserted into the map. It is
   * represented by the last curve in the given range.
   */
  template <class X_curve_2_iterator>
  Halfedge_iterator insert(const X_curve_2_iterator & begin,
                           const X_curve_2_iterator & end,
                           Change_notification * en = NULL)
  {
    X_curve_2_iterator it = begin;
    Halfedge_iterator out;
    if (it!=end) {
      out=insert(*it, en);
      it++;
    }
    while (it != end) {
      insert(*it, en);
      it++;
    }
    return out;
  }

  //! inserts a new edge into the map. The new
  //! edge is represented by a given curve, and placed in the interior of a
  //! given face.

  /*! insert_in_face_interior() inserts a new edge into the map. The new edge
   * is represented by a given curve, and placed in the interior of a given
   * face.
   * \param cv the given curve.
   * \param f the given face.
   * \param en the notification class. It's default value is NULL, which
   * implies no notification on insertion.
   * \return a handle to a new halfedge directed in the same way as the curve
   * cv. That is, The curve source and target points coinside with the points
   * of the source and target vertices of the returned halfedge respectively.)
   */
  Halfedge_handle insert_in_face_interior(const X_curve_2 & cv, 
                                          Face_handle f, 
                                          Change_notification * en = NULL);

  //! inserts a new edge into the map. The new edge
  //! is represented by a given curve, for which one endpoint, held by the
  //! target vertex of a given halfedge, is already in the map.

  /*! insert_from_vertex() inserts a new edge into the map. The new edge
   * is represented by a given curve, for which one endpoint, held by the
   * target vertex v of a given halfedge, is already in the map.
   * The returened twin halfedge is inserted immediately after the given
   * halfedge in the circular list of halfedges that share the same target
   * vertex v.
   * This method is the quick version of insert_from_vertex(), as the search
   * for the previous halfedge in the circular list of halfedges that share the
   * same target vertex v is unnecessary, saving its computation time.
   * \param cv the given curve.
   * \param prev the reference halfedge. Its target vertex v, is already
   * in the map.
   * \param en the notification class. It's default value is NULL, which
   * implies no notification on insertion.
   * \return a handle to a new halfedge that has the vertex v as its source.
   */
  Halfedge_handle insert_from_vertex(const X_curve_2 & cv,
                                     Halfedge_handle prev, 
                                     Change_notification * en = NULL
#ifdef _MSC_VER
                                     ,int dummy = 0
#endif
                                     );

  //! inserts a new edge into the map. The new edge
  //! is represented by a given curve, for which one endpoint, a given vertex,
  //! is already in the map.

  /*! insert_from_vertex() inserts a new edge into the map. The new edge is
   * represented by a given curve, for which one endpoint, a given vertex, is
   * already in the map.
   * \param cv the given curve.
   * \param v1 the given vertex.
   * \param source indicates whether the target of the returned halfedge holds
   * the given curve (cv) target point.
   * \param en the notification class. It's default value is NULL, which
   * implies no notification on insertion.
   * \return a handle to a new halfedge that has v1 as its source vertex.
   */
  Halfedge_handle insert_from_vertex(const X_curve_2 & cv, 
                                     Vertex_handle v1, 
                                     Change_notification * en = NULL);

  // Obsolete
  Halfedge_handle insert_from_vertex(const X_curve_2 & cv, 
                                     Vertex_handle v1, bool source, 
                                     Change_notification * en = NULL);

  //! inserts a new edge into the map. The new edge
  //! is represented by a given curve, for which its two endpoints, held by
  //! the target vertices of two given halfedges, are already in the map.

  /*! insert_at_vertices() inserts a new edge into the map. The new edge
   * is represented by a given curve, for which its two endpoints, held by
   * the target vertices of two given halfedges, are already in the map.
   * Call the two target vertices of the two given halfedges v1 and v2
   * respectively. They will be the source and target vertices of the returned
   * halfedge respectively.
   * The returened twin halfedge is inserted immediately after the first given
   * halfedge in the circular list of halfedges that share the same target
   * vertex v1.
   * The returened halfedge is inserted immediately after the second given
   * halfedge in the circular list of halfedges that share the same target
   * vertex v2.
   * This method is the quick version of insert_at_vertices(), as the searches
   * for the previous halfedges in the circular lists of halfedges that share
   * the same target vertices v1 and v2 is unnecessary, saving its computation
   * time.
   * \param cv the given curve.
   * \param prev1 the given halfedge that its target vertex is already in
   * the map, and will be the source vertex of the returned halfedge.
   * \param prev2 the given halfedge that its target vertex is already in
   * the map, and will be the target vertex of the returned halfedge.
   * \param en the notification class. It's default value is NULL, which
   * implies no notification on insertion.
   * \return a handle to a new halfedge that has the target vertices of the
   * prev1 and prev2 as its source and target vertices respectively.
   */
  Halfedge_handle insert_at_vertices(const X_curve_2 & cv,
                                     Halfedge_handle prev1, 
                                     Halfedge_handle prev2,
                                     Change_notification * en = NULL
#ifdef _MSC_VER
                                     ,int dummy = 0
#endif
                                     );

  //! inserts a new edge into the map. The new edge
  //! is represented by a given curve, for which its two endpoints, held by two
  //! given vertices, are already in the map.

  /*! insert_at_vertices() inserts a new edge into the map. The new edge is
   * represented by a given curve, for which its two endpoints, held by two
   * given vertices, are already in the map.
   * \param cv the given curve.
   * \param v1 the vertex in the map that holds the curve source point.
   * \param v2 the vertex in the map that holds the curve target point.
   * \param en the notification class. It's default value is NULL, which
   * implies no notification on insertion.
   * \return a handle to a new halfedge that has v1 and v2 as its source and
   * target vertices respectively.
   */
  Halfedge_handle insert_at_vertices(const X_curve_2 & cv, 
                                     Vertex_handle v1, 
                                     Vertex_handle v2, 
                                     Change_notification * en = NULL);

  //! is the planar map empty?
  bool is_empty() const { return halfedges_begin() == halfedges_end(); }

  // Note that the return types are abstract base classes
  const Point_location_base * get_point_location() const { return pl; }
  const Bounding_box_base * get_bounding_box() const { return bb; }
  const Traits_wrap & get_traits() const { return *traits; }

private:

  //a private implementation which defines if prev1 is on an outer ccb of 
  //the new face (returns true) or on an inner ccb (returns false)
  bool prev1_inside_hole(Halfedge_const_handle prev1,
                         Halfedge_const_handle prev2,
                         const X_curve_2& cv);  
  
public:  
  Halfedge_handle split_edge(Halfedge_handle       e, 
                             const X_curve_2       & c1, 
                             const X_curve_2       & c2,
                             Change_notification * en = NULL);

  Halfedge_handle merge_edge(Halfedge_handle e1, 
                             Halfedge_handle e2, 
                             const X_curve_2 & cv, 
                             Change_notification * en = NULL);              

  Face_handle remove_edge(Halfedge_handle e);

  Halfedge_handle vertical_ray_shoot(const Point_2 & p, 
                                     Locate_type   & lt, 
                                     bool            up)
  {
    CGAL_precondition(pl);
    return pl->vertical_ray_shoot(p,lt,up);
  }
  
  Halfedge_const_handle vertical_ray_shoot(const Point_2& p,
                                           Locate_type &lt, bool up) const
  {
    CGAL_precondition(pl);
    return ((const Point_location_base*)pl)->vertical_ray_shoot(p,lt,up);
    // The type cast to const is there to ensure that the planar map 
    // is not changed.
  }

  Halfedge_handle locate(const Point_2& p, Locate_type &lt) {
    CGAL_precondition(pl);
    return pl->locate(p,lt);
  }

  Halfedge_const_handle locate(const Point_2& p, Locate_type &lt) const {
    CGAL_precondition(pl);
    return ((const Point_location_base*)pl)->locate(p,lt);
    // The type cast to const is there to ensure that the planar map 
    // is not changed.
  }

  //! determines whether a given point lies within the interior of a given face

  /*! is_point_in_face() is a predicate that determines whether a given point 
   * lies within the interior of a given face.
   * A point lies within a face interior, iff the number of intersections 
   * between the face boundary and a ray emanating from the point is even.
   * Note that if the given face is the unbounded face, and it has no holes,
   * the point must lie within the face interior.
   * \param p the given point.
   * \param f a handle to the given face.
   * \return true if the given point lies within the interior of the given 
   * face, and false otherwise.
   */

  bool is_point_in_face(const Point_2 & p, Face_const_handle f)
  {
    if (!f->is_unbounded()) {
      // f is bounded:
      Halfedge_const_handle h = f->halfedge_on_outer_ccb();
      return point_is_in(p, h, h->curve());
    return false;
    }
    // f is the unbounded face:
    if (f->holes_begin() == f->holes_end()) return true;
    // f has at least one hole:
    Halfedge_const_handle h = *(f->holes_begin());
    return point_is_in(p, h, h->curve());
  }

protected:

  //! determines whether a given point lies within the interior of a face
  //! incident to a given halfedge.

  /*! point_is_in() is a predicate that determines whether a given point lies
   * within the interior of a face incident to a given halfedge. The halfedge
   * curve is provided explicitly, in case the halfedge hasn't been fully
   * constructed yet.
   * A point lies within a face interior, iff the number of intersections 
   * between the face boundary and a ray emanating from the point is even.
   * This function counts the number of intersections with a vertical ray, by
   * counting the number of boundary halfedges that are above the input point,
   * and the input point is in their x-range. The functions carefuly handles
   * degenerate cases. For example, the vertical ray coinsides with a boundary
   * halfedge.
   * \param p the given point.
   * \param ne a handle to a halfedge incident to the face in question.
   * \param ncv the curve of the given halfedge (for cases where the
   * halfedge is premature curve-less)
   * \return true if the given point lies within the interior of the face
   * incident to the given halfedge, and false otherwise.
   */
  bool point_is_in(const Point_2       & p, 
                   Halfedge_const_handle ne,
                   const X_curve_2     & ncv) const;
  
  /////////////////////////////////////////////////////////
  // Assignment functions 
  // Those are temporary functions and it should not be used
public:
  void assign(const Self &pm)
  {
    TPM::assign(pm);
    // traits->assign(pm->traits);
    // bb->assign(pm->bb);
    // pl->assign(pm->pl);
  }

  Self& operator=(const Self& pm);

  // used in implementation of operator=(
  void clear();
  
protected:
  // used in implementation of operator=(
  void x_curve_container(X_curve_2_container &l) const;
  // default initializer for the bounding box.
#include <CGAL/Planar_map_2/Bounding_box_special_initializer.h>

#ifdef CGAL_PM_DEBUG // for private debugging use

public:
  void debug()
  {
    if (pl) (pl->debug());
  }

#endif
private:
  /////////////////////////////////////////////////////////////////////////
  //                 Scanning Arrangement.
  ///////////////////////////////////////////////////////////////////////// 

  template <class Scanner>
  bool  scan_planar_map (Scanner& scanner)
  {
    if (!build_dcel(scanner))
      return false;

    return true;
  }

  template <class Scanner>
  bool  build_dcel (Scanner& scanner)
  {
      typedef typename Dcel::Vertex                 D_vertex;
    typedef typename Dcel::Halfedge                 D_halfedge;
    typedef typename Dcel::Face                     D_face;
  
    typedef typename  Dcel::Vertex_iterator         D_vetrex_iterator;
    typedef typename  Dcel::Vertex_const_iterator   D_vetrex_const_iterator;
    typedef typename  Dcel::Halfedge_iterator       D_halfedge_iterator;
    typedef typename  Dcel::Halfedge_const_iterator D_halfedge_const_iterator;
    typedef typename  Dcel::Face_iterator           D_face_iterator;
    typedef typename  Dcel::Face_const_iterator     D_face_const_iterator;

    // keeping a vector of halfedges (to access them easily by their indices).
    std::vector<D_halfedge* >  halfedges_vec;  

    std::vector<D_vertex* >    vertices_vec; 
 
    if ( ! scanner.in())
      return 0;

    scanner.scan_pm_vhf_sizes();
    if ( ! scanner.in()){
      std::cerr << "can't read vhf values"<<std::endl;
      clear();
      return false;
    }

    // read in all vertices
    unsigned int  i;
    for (i = 0; i < scanner.number_of_vertices(); i++) {
      D_vertex* nv = d.new_vertex();
      Point_2 p;

      // scanner.scan_vertex_attributes (nv);
      // if ( ! scanner.in()){
      // std::cerr << "can't read vertex attributes"<<std::endl;
      // clear();
      // return false;
      // }

      scanner.scan_vertex (nv);
      if ( ! scanner.in()){
        std::cerr << "can't read vertex"<<std::endl;
        clear();
        return false;
      }
      //nv->set_point(p);

      // for debug.
      //std::cout<<"Reading vertex no " <<i<<" point is ";
      //std::cout<<nv->point()<<std::endl;

      vertices_vec.push_back(nv);

      bb->insert(nv->point());
    
      //scanner.skip_to_next_vertex();
      //if ( ! scanner.in()){
      //  std::cerr << "can't skip to next vertex"<<std::endl;
      //  scanner.in().clear( std::ios::badbit);
      //  clear();
      //  return false;
      // }
    }
  
    for (i = 0; i < scanner.number_of_halfedges(); i++, i++){
      D_halfedge *nh = NULL;
      void  *nv1, *nv2;
      std::size_t index1, index2;
      X_curve_2 cv;

      // std::cout<<"Reading Edge no " <<i<<std::endl;

      nh = d.new_edge();
    
      // scanner.scan_halfedge_attributes (nh);
      // if ( ! scanner.in()){
      // std::cerr << "can't read halfedge attributes"<<std::endl;
      // clear();
      // return false;
      // }

      scanner.scan_index(index1);
      if ( ! scanner.in()){
        std::cerr << "can't read source of halfedge"<<std::endl;
        clear();
        return false;
      }
      cv = scanner.scan_halfedge(nh);
      if ( ! scanner.in()){
        std::cerr << "can't read halfedge"<<std::endl;
        clear();
        return false;
      }
      //nh->set_curve(cv);

      // scanner.scan_halfedge_attributes (nh->opposite());
      // if ( ! scanner.in()){
      // std::cerr << "can't read halfedge attributes"<<std::endl;
      // clear();
      // return false;
      // }

      scanner.scan_index (index2);
      if ( ! scanner.in()){
        std::cerr << "can't read source of halfedge"<<std::endl;
        clear();
        return false;
      }
      scanner.scan_halfedge(nh->opposite());
      if ( ! scanner.in()){
        std::cerr << "can't read halfedge"<<std::endl;
        clear();
        return false;
      }
      //nh->opposite()->set_curve(cv);

      nv1 = vertices_vec[index1];
      ((D_vertex*) nv1)->set_halfedge(nh); 
      nh->set_vertex((D_vertex*) nv1);
      //for debug
      //std::cout<<((D_vertex*) nv1)->point()<<std::endl;
    
      nv2 = vertices_vec[index2];
      ((D_vertex*) nv2)->set_halfedge(nh->opposite()); 
      nh->opposite()->set_vertex((D_vertex*) nv2);
      //for debug
      //std::cout<<((D_vertex*) nv2)->point()<<std::endl;

      pl->insert(D_halfedge_iterator(nh->opposite()), cv);
      bb->insert(cv);
    
      halfedges_vec.push_back(nh);
      halfedges_vec.push_back(nh->opposite());

      //scanner.skip_to_next_halfedge();
      //if ( ! scanner.in()){
      //  std::cerr << "can't skip to next halfedge"<<std::endl;
      //  scanner.in().clear( std::ios::badbit);
      //  clear();
      //  return false;
      //} 
    }
  
    // read in all facets
    for (i = 0; i < scanner.number_of_faces(); i++) {
      //std::size_t  num_of_holes, num_halfedges_on_outer_ccb;
    
      //std::cout<<"Reading Face no " <<i<<std::endl;
   
      D_face* nf = u_face; //this is the unbounded face.
      if (i > 0)  // else - allocate the bounded face.
        nf = d.new_face();

      scanner.scan_face(nf);
      if ( ! scanner.in()){
        std::cerr << "can't read face"<<std::endl;
        clear();
        return false;
      }
    }
  
    if ( ! scanner.in() ) {
      scanner.in().clear( std::ios::badbit);
      return false;
    } 
  
    return true;
  }

  // Data Members
  // ------------

protected:
  Point_location_base * pl;
  Bounding_box_base   * bb;
  Traits_wrap         * traits;

private:
  bool use_delete_pl;
  bool use_delete_bb;
  bool use_delete_traits;
};

//-----------------------------------------------------------------------------
//  Member Function Definitions
//-----------------------------------------------------------------------------
template < class Dcel, class Traits >
Planar_map_2< Dcel, Traits >::Planar_map_2()
{
  traits = new Traits_wrap();
  use_delete_traits = true;
    
#ifndef CGAL_NO_PM_DEFAULT_POINT_LOCATION
  pl = new Pm_default_point_location<Self>;
  use_delete_pl = true;
  pl->init(*this,*traits);
#else
  CGAL_assertion_msg(false,
    "No default point location is defined; you must supply one.");
#endif
    
  bb=init_default_bounding_box((Traits*)traits);
  use_delete_bb=true;
  bb->init(*this,*traits);
}

//-----------------------------------------------------------------------------
template < class Dcel, class Traits >
Planar_map_2< Dcel, Traits >::
Planar_map_2(typename Planar_map_2<Dcel, Traits>::Point_location_base * pl_ptr)
{ 
  traits = new Traits_wrap();
  use_delete_traits = true;
    
  pl = pl_ptr;
  use_delete_pl = false;
  pl->init(*this,*traits);
    
  bb = init_default_bounding_box((Traits*)traits);
  use_delete_bb = true;
  bb->init(*this,*traits);
}
//-----------------------------------------------------------------------------
template < class Dcel, class Traits >
Planar_map_2< Dcel, Traits >::
Planar_map_2(
  const typename Planar_map_2< Dcel, Traits >::Traits & tr_, 
  typename Planar_map_2< Dcel, Traits >::Point_location_base   * pl_ptr, 
  typename Planar_map_2< Dcel, Traits >::Bounding_box_base     * bb_ptr )
{
  traits = new Traits_wrap(tr_);
  use_delete_traits = true;
    
  if (pl_ptr == NULL)
  {
#ifndef CGAL_NO_PM_DEFAULT_POINT_LOCATION
    pl = new Pm_default_point_location<Self>;
    use_delete_pl = true;
    pl->init(*this,*traits);
#else
    CGAL_assertion_msg( false,
   "No default point location is defined; you must supply one.");
#endif
  }
  else
  {
    pl = pl_ptr;
    use_delete_pl = false;
    pl->init(*this,*traits);
  }
    
  if (bb_ptr == NULL)
  {
    bb=init_default_bounding_box((Traits*)traits);
    use_delete_bb=true;
    bb->init(*this,*traits);
  }
  else
  {
    bb = bb_ptr;
    use_delete_bb = false;
    bb->init(*this,*traits);
  }
}
//-----------------------------------------------------------------------------
template < class Dcel, class Traits >
Planar_map_2< Dcel, Traits >::
Planar_map_2(
  typename Planar_map_2< Dcel, Traits >::Traits_wrap          *tr_ptr, 
  typename Planar_map_2< Dcel, Traits >::Point_location_base  *pl_ptr, 
  typename Planar_map_2< Dcel, Traits >::Bounding_box_base    *bb_ptr )
{
  if (tr_ptr == NULL)
  {
    traits = new Traits_wrap();
    use_delete_traits = true;
  }
  else
  {
    traits = tr_ptr;
    use_delete_traits = false;
  }
    
  if (pl_ptr == NULL)
  {
#ifndef CGAL_NO_PM_DEFAULT_POINT_LOCATION
    pl = new Pm_default_point_location<Self>;
    use_delete_pl = true;
    pl->init(*this,*traits);
#else
    CGAL_assertion_msg( false,
    "No default point location is defined; you must supply one.");
#endif
  }
  else
  {
    pl = pl_ptr;
    use_delete_pl = false;
    pl->init(*this,*traits);
  }
    
  if (bb_ptr == NULL)
  {
    bb=init_default_bounding_box((Traits*)traits);
    use_delete_bb=true;
    bb->init(*this,*traits);
  }
  else
  {
    bb = bb_ptr;
    use_delete_bb = false;
    bb->init(*this,*traits);
  }  
}
//-----------------------------------------------------------------------------
template < class Dcel, class Traits >
Planar_map_2< Dcel, Traits >::
Planar_map_2(const Planar_map_2< Dcel, Traits >& pm )
{
  // doing the same as Planar_map_2(pm.get_traits(),pm.get_point_location(),
  //                                pm.get_point_bbox());
    
  typedef Pm_naive_point_location<Planar_map_2<Dcel,Traits> >  Pm_naive;
  typedef Pm_naive*                                    Pm_naive_pointer;

  traits = new Traits_wrap();
  use_delete_traits = true;

  if (Pm_naive_pointer tmp_pl = dynamic_cast<Pm_naive_pointer>(pm.pl) ){
    //cout<<"Naive"<<std::endl;
    pl = new Pm_naive_point_location<Self>;
  }
  else if (Pm_walk_along_line_point_location<Self>* tmp_pl = 
           dynamic_cast<Pm_walk_along_line_point_location<Self>*>(pm.pl) ){
    pl = new Pm_walk_along_line_point_location<Self>;
    //cout<<"Walk"<<std::endl;
  }
  else{
    //cout<<"Default"<<std::endl;
#ifndef CGAL_NO_PM_DEFAULT_POINT_LOCATION
    pl = new Pm_default_point_location<Self>;
#else
    CGAL_assertion_msg( false,
    "No default point location is defined; you must supply one.");
#endif
  }
  use_delete_pl = true;
  pl->init(*this,*traits);
    
  bb=init_default_bounding_box((Traits*)traits);
  use_delete_bb=true;
  bb->init(*this,*traits);
    
  assign(pm);
    
  Halfedge_iterator h_iter;
  for (h_iter = halfedges_begin(); 
       h_iter != halfedges_end(); 
       h_iter++, h_iter++)
    pl->insert(h_iter, h_iter->curve());
    
  for (Vertex_iterator v_iter = vertices_begin(); 
       v_iter != vertices_end(); 
       v_iter++)
    bb->insert(v_iter->point());
    
  for (h_iter = halfedges_begin(); 
       h_iter !=  halfedges_end(); 
       h_iter++, h_iter++)
    bb->insert(h_iter->curve());
}
//-----------------------------------------------------------------------------
template < class Dcel, class Traits > 
Planar_map_2< Dcel, Traits >::~Planar_map_2()
{
  if (use_delete_pl) delete pl;
  if (use_delete_bb) delete bb;
  if (use_delete_traits) delete traits;
}
//-----------------------------------------------------------------------------
template < class Dcel, class Traits > 
bool 
Planar_map_2< Dcel, Traits >::
read(std::istream & in)
{
  clear();
  Pm_file_scanner<Self> scanner(in); 
  return scan_planar_map(scanner);
}

/*!
 */
template < class Dcel, class Traits >
typename Planar_map_2< Dcel, Traits >::Halfedge_handle 
Planar_map_2< Dcel, Traits >::
insert_in_face_interior(
  const typename Planar_map_2< Dcel, Traits >::X_curve_2 & cv, 
  typename Planar_map_2< Dcel, Traits >::Face_handle f, 
  Change_notification * en)
{
  Halfedge_handle h = Topological_map<Dcel>::insert_in_face_interior(f);
  h->set_curve(cv);  //should set the curve of the twin as well but for now
  h->twin()->set_curve(cv);
  
  //pl->insert(h);  //maybe should be above
  //iddo - for arrangement
  pl->insert(h,cv);

  h->source()->set_point(traits->curve_source(cv));
  h->target()->set_point(traits->curve_target(cv));

  if (en != NULL)
  {
    en->add_edge(cv, h, true, false);
    en->add_hole(f, h);
  }

  return h;

}

/*!
 */
template < class Dcel, class Traits >
typename Planar_map_2< Dcel, Traits >::Halfedge_handle 
Planar_map_2< Dcel, Traits >::
insert_from_vertex(const typename Planar_map_2<Dcel, Traits>::X_curve_2 & cv,
                   typename Planar_map_2<Dcel,Traits>::Halfedge_handle prev,
                   Change_notification * en
#ifdef _MSC_VER
  ,int dummy
#endif
  )
{
#ifdef _MSC_VER
  (void) dummy;
#endif
  CGAL_precondition_msg(traits->point_is_same(prev->target()->point(), 
                                              traits->curve_source(cv)) ||
                        traits->point_is_same(prev->target()->point(), 
                                              traits->curve_target(cv)),
  "Point of target vertex of input halfedge should be a curve endpoint.");

  Halfedge_handle h = Topological_map<Dcel>::insert_from_vertex(prev);  
  h->set_curve(cv);  
  h->twin()->set_curve(cv);

  pl->insert(h, cv);            // for arrangement

  bool source = traits->point_is_same(prev->target()->point(), 
                                      traits->curve_source(cv));
  (source) ?
    h->target()->set_point(traits->curve_target(cv)) :
    h->target()->set_point(traits->curve_source(cv));

  if (en != NULL) en->add_edge(cv, h, true, false);

  return h;
}

/*!
 */
template < class Dcel, class Traits >
typename Planar_map_2< Dcel, Traits >::Halfedge_handle 
Planar_map_2< Dcel, Traits >::
insert_from_vertex(const typename Planar_map_2< Dcel, Traits >::X_curve_2 & cv,
                   typename Planar_map_2< Dcel, Traits >::Vertex_handle v1, 
                   Change_notification * en)
{
  CGAL_precondition_msg(traits->point_is_same(v1->point(), 
                                              traits->curve_source(cv)) ||
                        traits->point_is_same(v1->point(), 
                                              traits->curve_target(cv)),
  "Point of input vertex should be a curve endpoint.");

  // Find the previous of cv:
  Halfedge_around_vertex_circulator prev = v1->incident_halfedges(),
      after = prev,
      infinite_loop = prev;
  ++after;

  if (after != prev) {
    while (!(traits->curve_is_between_cw(cv,prev->curve(),
                                         after->curve(),v1->point()))) {
      prev = after;
      ++after;
      if (prev == infinite_loop)  // infinite loop indication
      {
        std::cerr << std::endl << "Planar_map_2::insert_from_vertex("
                  << "const X_curve_2& cv, Vertex_handle v1, "
                  << "bool source) called with previously "
                  << "inserted curve " << std::endl;
        return Halfedge_handle();
      }
    }
  }

  return insert_from_vertex(cv, prev, en);
}

// Obsolete
template < class Dcel, class Traits >
typename Planar_map_2< Dcel, Traits >::Halfedge_handle 
Planar_map_2< Dcel, Traits >::
insert_from_vertex(const typename Planar_map_2< Dcel, Traits >::X_curve_2 & cv,
                   typename Planar_map_2< Dcel, Traits >::Vertex_handle v1,
                   bool source,
                   Change_notification * en)
{
  (void) source;
  // For some reason MSVC cannot handle the following call, even if the
  // definition is inlined in the class. Too many nested calls. Go figure...
#if 0
  return insert_from_vertex(cv, v1, en);
#else
  // Find the previous of cv:
  Halfedge_around_vertex_circulator prev = v1->incident_halfedges(),
      after = prev,
      infinite_loop = prev;
  ++after;

  if (after != prev) {
    while (!(traits->curve_is_between_cw(cv,prev->curve(),
                                         after->curve(),v1->point()))) {
      prev = after;
      ++after;
      if (prev == infinite_loop)  // infinite loop indication
      {
        std::cerr << std::endl << "Planar_map_2::insert_from_vertex("
                  << "const X_curve_2& cv, Vertex_handle v1, "
                  << "bool source) called with previously "
                  << "inserted curve " << std::endl;
        return Halfedge_handle();
      }
    }
  }

  return insert_from_vertex(cv, prev, en);
#endif
}

/*!
 */
template < class Dcel, class Traits >
typename Planar_map_2< Dcel, Traits >::Halfedge_handle 
Planar_map_2< Dcel, Traits >::
insert_at_vertices(const typename Planar_map_2<Dcel, Traits>::X_curve_2 & cv,
                   typename Planar_map_2<Dcel, Traits>::Halfedge_handle prev1, 
                   typename Planar_map_2<Dcel, Traits>::Halfedge_handle prev2,
                   Change_notification * en
#ifdef _MSC_VER
  ,int dummy
#endif
  )
{
#ifdef _MSC_VER
  (void) dummy;
#endif
  CGAL_precondition_msg(traits->point_is_same(prev1->target()->point(), 
                                              traits->curve_source(cv)) &&
                        traits->point_is_same(prev2->target()->point(), 
                                              traits->curve_target(cv)) ||
                        traits->point_is_same(prev2->target()->point(), 
                                              traits->curve_source(cv)) &&
                        traits->point_is_same(prev1->target()->point(), 
                                              traits->curve_target(cv)),
  "Points of target vertices of input halfedges should be curve endpoints.");

  Size num_before = number_of_faces();

  bool prev1_before_prev2 = prev1_inside_hole(prev1, prev2, cv);
  Halfedge_handle h = (prev1_before_prev2) ?
    Topological_map<Dcel>::insert_at_vertices(prev1, prev2) : 
    Topological_map<Dcel>::insert_at_vertices(prev2, prev1); 

  h->set_curve(cv);
  h->twin()->set_curve(cv);

  Size num_after = number_of_faces();
  if (num_after - num_before) {         // a face was added => move holes
    Face_handle nf = h->face();         // the new face is pointed at by h
    Face_handle of = h->twin()->face(); // old face

    Holes_iterator it = of->holes_begin();
    while (it != of->holes_end()) {
      // check if the hole is inside new face
      // new for arrangement
      if (point_is_in((*it)->target()->point(), h, cv)) {
        Holes_iterator tmp = it;  // deletion invalidates iterators so... 
        ++it;   // assumes only the erased iterator is invalidated (like stl
        // list) 

        move_hole(tmp, of, nf); 
      }
      else
        ++it;
    }
  }

  // v1 should be the source of h.
  if (!prev1_before_prev2) h = h->twin();

  //pl->insert(h);
  //iddo - for arrangement
  pl->insert(h, cv);
    
  // Notifying change.
  if (en != NULL) {
    Face_handle orig_face =
        (!prev1_before_prev2) ? h->face() : h->twin()->face();

    en->add_edge(cv, h, true, false);
      
    // After fixing the notifier we won't have to check that since
    // h->face() will be surely the new face.
    (h->face() == orig_face) ?
      en->split_face(h->face(), h->twin()->face()) :
      en->split_face(h->twin()->face(), h->face());
      
    // we surely know h->face() is the new face.
    // en->split_face(h->twin()->face(), h->face());
  }
    
  return h;
}

/*!
 */
template < class Dcel, class Traits >
typename Planar_map_2< Dcel, Traits >::Halfedge_handle 
Planar_map_2< Dcel, Traits >::
insert_at_vertices(const typename Planar_map_2<Dcel, Traits>::X_curve_2 & cv, 
                   typename Planar_map_2< Dcel, Traits >::Vertex_handle   v1, 
                   typename Planar_map_2< Dcel, Traits >::Vertex_handle   v2, 
                   Change_notification                                  * en)
{
  CGAL_precondition_msg(traits->point_is_same(v1->point(), 
                                              traits->curve_source(cv)) &&
                        traits->point_is_same(v2->point(), 
                                              traits->curve_target(cv)) ||
                        traits->point_is_same(v2->point(), 
                                              traits->curve_source(cv)) &&
                        traits->point_is_same(v1->point(), 
                                              traits->curve_target(cv)),
                        "Points of input vertices should be curve endpoints.");

  Halfedge_around_vertex_circulator prev1 = v1->incident_halfedges(),
    after = prev1,
    infinite_loop = prev1;
  ++after;

  if (after != prev1) {
    while (!(traits->curve_is_between_cw(cv, prev1->curve(),
                                         after->curve(), v1->point())))
    {
      prev1 = after;
      ++after;
      if (prev1 == infinite_loop)  // infinite loop indication
      {
        std::cerr << std::endl << "Planar_map_2::insert_at_vertices("
                  << "const X_curve_2 & cv, Vertex_const_handle v1, "
                  << "Vertex_const_handle v2) called with previously "
                  << "inserted curve " << std::endl;
        return Halfedge_handle();
      }
    }
  }    

  Halfedge_around_vertex_circulator prev2 = v2->incident_halfedges();
  after = prev2;
  infinite_loop = prev2;
  ++after;

  if (after != prev2) {
    while (!(traits->curve_is_between_cw(cv, prev2->curve(),
                                         after->curve(),v2->point())))
    {
      prev2 = after;
      ++after;
      if (prev2 == infinite_loop) // infinite loop indication
      {
        std::cerr << std::endl << "Planar_map_2::insert_at_vertices("
                  << "const X_curve_2 & cv, Vertex_const_handle v1,"
                  << "Vertex_const_handle v2) called with previously "
                  << "inserted curve " << std::endl;
        return Halfedge_handle();
      }
    }
  }    

  return insert_at_vertices(cv, prev1, prev2, en);
}

//-----------------------------------------------------------------------------
template < class Dcel, class Traits >
bool 
Planar_map_2< Dcel, Traits >::
prev1_inside_hole(
  typename Planar_map_2< Dcel, Traits >::Halfedge_const_handle prev1,
  typename Planar_map_2< Dcel, Traits >::Halfedge_const_handle prev2,
  const typename Planar_map_2< Dcel, Traits >::X_curve_2 & cv)
{
  
  // Defining geometrically whether there is a new face. If there is,
  // finds if prev1 is on the outside of the new face (send
  // prev1,prev2) or on the inside of the new face (send prev2,prev1)

  // The algorithm: 

  // 1. go over all the halfedges of the face which
  // will hold prev1 (since the new face is not constructed yet,
  // this is modeled by going from prev2->next to prev1 and
  // then over the new curve)

  // 2. find if the left-most-lower halfedge in the path (i.e, the one
  // with the leftmost down target and is the lowest to the right
  // among the incident edges of this vertex) is directed left (we are
  // on the outside) or right (we are inside ) (if not on same ccb
  // then it doesn't matter and return true)
  
  Ccb_halfedge_const_circulator left_edge(prev2);
  ++left_edge;
  Ccb_halfedge_const_circulator first(prev2),curr(left_edge),
    last(prev1);
  ++last; //we want the prev1 to be checked as well 

  Point_2 left = prev2->target()->point();
  bool b;

  do {
    //source
    b=false;
    if (traits->point_is_left( curr->source()->point(),left)) 
      b=true;
    else
      if (traits->point_is_same(curr->source()->point(),left)) {
        if (traits->curve_is_vertical(curr->curve()) &&
            traits->point_is_lower(curr->target()->point(),left) )
          b=true;
        else
          if (traits->curve_compare_at_x_right(curr->curve(),
                                               left_edge->curve(),
                                               left)==SMALLER ) 
            b=true;
      }

    if (b) {
      left=curr->source()->point();
      left_edge=curr;
    }

    //target
    b=false;
    if (traits->point_is_left( curr->target()->point(),left))
      b=true;
    if (traits->point_is_same(curr->target()->point(),left)) {
      if (traits->curve_is_vertical(curr->curve()) &&
          traits->point_is_lower(curr->source()->point(),left) )
        b=true;
      else
        if (traits->curve_compare_at_x_right(curr->curve(),
                                             left_edge->curve(),
                                             left)==SMALLER ) 
          b=true;

      //we want in the degenerate case to return the halfedge 
      //pointing _at_ the left point 
        else
          if ( (curr)==(left_edge->twin()) )
            b=true;
    }

    if (b) {
      left=curr->target()->point();
      left_edge=curr;
    }

    ++curr;
  } while ( (curr != first) && (curr != last) );

  //test the new curve against left_edge
  if (traits->point_is_same(traits->curve_target(cv),left)||
      traits->point_is_same(traits->curve_source(cv),left)) {
    if (traits->curve_is_vertical(cv)) {
      return (traits->point_is_lower(prev2->target()->point(),
                                     prev1->target()->point()));
    }
    else
      if (traits->curve_compare_at_x_right(cv,left_edge->curve(), 
                                           left)==SMALLER ) {  
        return (traits->point_is_left(prev1->target()->point(),
                                      prev2->target()->point()));
      }
  }

  //check if left_edge is from left to right
  if (traits->curve_is_vertical(left_edge->curve())) {
    if (traits->point_is_lower(left_edge->source()->point(),
                               left_edge->target()->point()))
      return false;
    else
      return true;
  }

  return (traits->point_is_left(left_edge->source()->point(),
                                left_edge->target()->point()));
}

/*!
 */
template < class Dcel, class Traits >
typename Planar_map_2< Dcel, Traits >::Halfedge_handle 
Planar_map_2< Dcel, Traits >::
insert(const typename Planar_map_2< Dcel, Traits >::X_curve_2 & cv, 
       Change_notification * en)
{
  CGAL_precondition_msg( ! traits->curve_is_degenerate(cv),
                         "Curve length should be greater than zero.");
  CGAL_precondition_msg(bb, "Bounding box should not be null.");

  bb->insert(cv);
    
  Locate_type lt1, lt2;
  Point_2 src = traits->curve_source(cv);
  Point_2 tgt = traits->curve_target(cv);
    
  // The point location may not change the bounding box.
  Halfedge_handle h1 = ((const Point_location_base*)pl)->locate(src, lt1);
  Halfedge_handle h2 = ((const Point_location_base*)pl)->locate(tgt, lt2);
   
  // In principal, the result of a locate should not be an edge, 
  // because the planar map does not accept proper intersections.
  // It is only possible in case a bounding box curve was hit.
  if (lt1 == EDGE || lt1 == UNBOUNDED_EDGE) 
  {
    // the curve intersects the bounding box.
    Halfedge_handle h = h1, h2;
    bb->split_boundary_edge(h, h1, h2, src);
    // make sure the intersection point is in the map, 
    // i.e. split the halfedge that contains its.
    lt1 = VERTEX; 
  }

  if (lt2 == EDGE || lt2 == UNBOUNDED_EDGE) 
  {
    Halfedge_handle h1, h = h2;
    bb->split_boundary_edge(h, h1, h2, tgt);
    // make sure the intersection point is in the map, 
    // i.e. split the halfedge that contains its.
    lt2 = VERTEX;
  }

  if (lt1 == VERTEX && lt2 == VERTEX) 
    return insert_at_vertices(cv, h1->target(), h2->target(), en); 
    
  if (lt1 == VERTEX && lt2 != VERTEX)
    return insert_from_vertex(cv, h1->target(), true, en); 

  if (lt1 != VERTEX && lt2 == VERTEX)
    return insert_from_vertex(cv, h2->target(), false, en)->twin();
    
  if (lt1 == UNBOUNDED_FACE)
    return insert_in_face_interior(cv, unbounded_face(), en);
    
  if (lt1 == FACE)
    return insert_in_face_interior(cv, h1->face(), en);
    
  CGAL_assertion_msg(lt1 == VERTEX || lt1 == UNBOUNDED_FACE || lt1 == FACE,
  "Endpoints should not coinside with an edge. No intersections allowed.");

  return Halfedge_handle();
}

//-----------------------------------------------------------------------------
template < class Dcel, class Traits >
typename Planar_map_2< Dcel, Traits >::Halfedge_handle 
Planar_map_2< Dcel, Traits >::
split_edge(typename Planar_map_2< Dcel, Traits >::Halfedge_handle   e, 
           const typename Planar_map_2< Dcel, Traits >::X_curve_2   & c1, 
           const typename Planar_map_2< Dcel, Traits >::X_curve_2   & c2,
           Change_notification                                    * en )
{
  CGAL_precondition(traits->point_is_same(traits->curve_source(c2),
                                          traits->curve_target(c1)));

  CGAL_precondition(traits->point_is_same(traits->curve_source(c1),
                                          e->source()->point()) &&
                    traits->point_is_same(traits->curve_target(c2),
                                          e->target()->point()) ||
                    traits->point_is_same(traits->curve_source(c1),
                                          e->target()->point()) &&
                    traits->point_is_same(traits->curve_target(c2),
                                          e->source()->point()));

  X_curve_2 cv(e->curve());

  Halfedge_handle h = Topological_map<Dcel>::split_edge(e);

  if (traits->point_is_same(traits->curve_source(c1),h->source()->point())) {
    h->set_curve(c1);
    h->twin()->set_curve(c1);
    h->next_halfedge()->set_curve(c2);
    h->next_halfedge()->twin()->set_curve(c2);
    h->target()->set_point(traits->curve_target(c1));
    pl->split_edge(cv,h,h->next_halfedge(),c1,c2);

    if (en != NULL) 
      en->split_edge(h, h->next_halfedge(), c1, c2);
  }
  else {
    h->set_curve(c2);
    h->twin()->set_curve(c2);
    h->next_halfedge()->set_curve(c1);
    h->next_halfedge()->twin()->set_curve(c1);
    h->target()->set_point(traits->curve_target(c1));
    pl->split_edge(cv,h,h->next_halfedge(),c2,c1);
      
    if (en != NULL) 
      en->split_edge(h, h->next_halfedge(), c2, c1);
  }

  //if (en != NULL) 
  //  en->split_edge(h, h->next_halfedge(), c1, c2);

  return h;
}
//-----------------------------------------------------------------------------
template < class Dcel, class Traits >
typename Planar_map_2< Dcel, Traits >::Halfedge_handle 
Planar_map_2< Dcel, Traits >::
merge_edge(typename Planar_map_2< Dcel, Traits >::Halfedge_handle   e1, 
           typename Planar_map_2< Dcel, Traits >::Halfedge_handle   e2, 
           const typename Planar_map_2< Dcel, Traits >::X_curve_2   & cv, 
           Change_notification                                    * en)
{
  CGAL_precondition((traits->point_is_same(traits->curve_source(cv),
                                           e1->source()->point()) &&
                     traits->point_is_same(traits->curve_target(cv),
                                           e2->target()->point())) || 
                    (traits->point_is_same(traits->curve_target(cv),
                                           e1->source()->point()) &&
                     traits->point_is_same(traits->curve_source(cv),
                                           e2->target()->point())));

  // problematic: since we assume e1 will be the new merged halfedge
  // after merging.  en->merge(e1,e2,cv);
    
  X_curve_2 c1(e1->curve()), c2(e2->curve());

  Halfedge_handle h = Topological_map<Dcel>::merge_edge(e1,e2); 
  h->set_curve(cv);
  h->twin()->set_curve(cv);

  //pl->merge_edge(c1,c2,h);
  //iddo - for arrangement
  pl->merge_edge(c1, c2, h, cv);
    
  // problematic: e2 does not exist anymore 
  //if (en != NULL) 
  //  en->merge_edge(h, e1, e2, cv);

  return h;
} 
//-----------------------------------------------------------------------------
template < class Dcel, class Traits >
typename Planar_map_2< Dcel, Traits >::Face_handle
Planar_map_2< Dcel, Traits >::
remove_edge(typename Planar_map_2< Dcel, Traits >::Halfedge_handle e )
{

  // en->remove_edge(e);
    
  pl->remove_edge(e);

  //if a new hole can be created define geometrically the 
  //halfedge (e or e->twin) that points at the new hole.
  //if the leftmost point in the path e...e->twin
  //is left of the leftmost point in the path e->twin ... e
  //then e->twin  points at the hole created.

  if (e->face() == e->twin()->face() ) {
    Ccb_halfedge_circulator ccb_e=e->ccb() ;
    Ccb_halfedge_circulator ccb_t=e->twin()->ccb();

    Point_2 e_left=e->target()->point();
    Point_2 t_left=ccb_t->target()->point();

    //find the leftmost point in the path from e to its twin
    Ccb_halfedge_circulator aux=ccb_e;
    do {
      if (traits->compare_x(aux->target()->point(),e_left)==SMALLER) {
        e_left=aux->target()->point();
      }
    } while (++aux!=ccb_t);

    //find the leftmost point in the path from the twin to e
    aux=ccb_t;
    do {
      if (traits->compare_x(aux->target()->point(),t_left)==SMALLER) {
        t_left=aux->target()->point();
      }        
    } while (++aux!=ccb_e);

    //compare the two left points
    if (traits->compare_x(t_left,e_left) == SMALLER) //e points at hole 
      return Topological_map<Dcel>::remove_edge(e);
    else
      return Topological_map<Dcel>::remove_edge(e->twin());
  }
  else {
    return Topological_map<Dcel>::remove_edge(e);
  }
}
//-----------------------------------------------------------------------------
template < class Dcel, class Traits >
bool Planar_map_2< Dcel, Traits >::
point_is_in( const Point_2           & p, 
             Halfedge_const_handle     ne,
             const X_curve_2         & ncv) const
{
  // count stores the number of curves that intersect the upward vertical 
  // ray shot from p (except for a degenerate case which is explained in 
  // the code)
  int count = 0;

  // 1. Find the first halfedge, whose curve is non-vertical, along
  // the ccb that includes input halfedge ne.
  Ccb_halfedge_const_circulator circ = ne;
  do {
    ++circ;
  } while ( circ != ne && traits->curve_is_vertical(circ->curve()) );

  // If the whole ccb is vertical then there is no face, so point p
  // cannot be in it
  if ( circ == ne && traits->curve_is_vertical(ncv) )
    return false; 

  // 2. Go over all curves of the ccb and count those which are above p.
  Ccb_halfedge_const_circulator last = circ;
  do {

    // Put curve of current halfedge in circv.
    X_curve_2 circv;
    // If not on the new halfedge circ definitely has a curve
    if (circ != ne) 
    { 
      circv=circ->curve();
    }
    // o/w, circ might not have a curve yet (e.g in arrangement)
    // so we take the input curve.
    else { 
      circv=ncv;
    }

    // If query point is vertex point on the outer ccb
    if (traits->point_is_same(circ->target()->point(), p)) 
      return false;

    // If current curve is not vertical
    if ( ! traits->curve_is_vertical(circv)) 
    {
      // If point is under current curve in the range (source,target] of it
      if ((traits->curve_get_point_status(circv,p) == 
           Traits::UNDER_CURVE) && 
          !(traits->point_is_same_x(circ->source()->point(), p))) 
      {  
        // If p is exactly under a vertex of the ccb 
        if (traits->point_is_same_x(circ->target()->point(), p)) 
        {
          // Put curve of next halfedge that is not vertical in nextcv
          Ccb_halfedge_const_circulator next = circ;
          ++next;
          X_curve_2 nextcv;
          if (next != ne) {
            nextcv = next->curve();
          }
          else {
            nextcv = ncv;
          }
          if (traits->curve_is_vertical(nextcv)) {
            //advance to non-vertical edge
            while (traits->curve_is_vertical(nextcv)) {
              if (next!=ne) {
                nextcv=next->curve();
              }
              else {
                nextcv=ncv;
              }
              ++next;
            }
          }
          // If nextcv is on the same side of the vertical line
          // from p as circv is 
          if ((traits->point_is_right(circ->source()->point(), p) &&
               traits->point_is_left(next->target()->point(), p)) ||
              (traits->point_is_left(circ->source()->point(), p) &&
               traits->point_is_right(next->target()->point(), p))) {
            // then we raise the count
            ++count;
          }
        }
        else 
        {
          // o/w, point p is under the interior of the current curve
          // so we raise the count
          ++count;
        }
      }
    } // If current curve is not vertical
  } while (++circ != last);
    
  return (count%2 != 0);  //if count is odd return true
}
//-----------------------------------------------------------------------------
template < class Dcel, class Traits >
Planar_map_2< Dcel, Traits >& 
Planar_map_2< Dcel, Traits >::
operator=(const Planar_map_2< Dcel, Traits >& pm)
{
  if( this != &pm ){
    clear();
    assign(pm);
               
    Halfedge_iterator h_iter;
    for( h_iter = halfedges_begin(); 
         h_iter != halfedges_end(); 
         h_iter++, h_iter++)
      pl->insert(h_iter, h_iter->curve());
                 
    for( Vertex_iterator v_iter = vertices_begin(); 
         v_iter != vertices_end(); 
         v_iter++)
      bb->insert(v_iter->point());
                 
    for( h_iter = halfedges_begin(); 
         h_iter !=  halfedges_end(); 
         h_iter++, h_iter++)
      bb->insert(h_iter->curve());
  }
  return *this;
}
//-----------------------------------------------------------------------------
template < class Dcel, class Traits >
void Planar_map_2< Dcel, Traits >:: clear()
{
  pl->clear();
  TPM::clear();
  // Halfedge_iterator it=halfedges_begin(),prev=it,it_e=halfedges_end();
  // while (it!=it_e) {++it;++it;remove_edge(prev);prev=it;}
  bb->clear();
}
//-----------------------------------------------------------------------------
template < class Dcel, class Traits >
void Planar_map_2< Dcel, Traits >::
x_curve_container(X_curve_2_container &l) const
{
  Halfedge_const_iterator it=halfedges_begin(),it_e=halfedges_end();
  while (it!=it_e){
    l.push_back(it->curve());
    ++it;
    ++it;
  }
}
//-----------------------------------------------------------------------------
  
CGAL_END_NAMESPACE

#endif // CGAL_PLANAR_MAP_2_H
// EOF
