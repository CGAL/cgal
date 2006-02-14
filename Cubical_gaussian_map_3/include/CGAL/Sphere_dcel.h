// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
//
// Author(s)     : Kapelushnik Lior <liorkape@post.tau.ac.il>

/*! \file
 * spherical arrangements of none intersecting arcs of great circles on a sphere
 */

#ifndef CGAL_SPHERE_DCEL_H
#define CGAL_SPHERE_DCEL_H

#include <CGAL/Spherical_cgm_arr_dcel.h>
#include <CGAL/Cubical_gaussian_map_3.h>
#include <CGAL/Sphere_arc.h>

#include <list>
//temporary:
#include <iostream>
// end temporary

CGAL_BEGIN_NAMESPACE

/*
 the class represents the spherical arrangement elements which are
 halfedge, vertex, face and circulators around vertex or connected compononent boundary

 _Kernel - the kernel to base uppon the spherical map and the internal cubical gaussian map
 _T_Dcel - the cubical gaussian map basic dcel structure template class
*/
template <class _Kernel, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
  template <class T>
#endif
  class _T_Dcel = Spherical_cgm_arr_dcel>
class SphereTopologicalMap { 
public:
  
  class Halfedge;
  class Vertex;
  class Face;
  class Halfedge_around_vertex_circulator;
  class Ccb_halfedge_circulator;

  // public definitions
  typedef CGAL::Cubical_gaussian_map_3<_Kernel, _T_Dcel>  CGM;
  typedef _Kernel                      Kernel;

  typedef std::list<Halfedge>               HalfedgeList;
  typedef typename HalfedgeList::iterator         Halfedge_iterator;
  typedef typename HalfedgeList::const_iterator      Halfedge_const_iterator;
  typedef Halfedge_iterator                 Halfedge_handle;
  typedef Halfedge_const_iterator             Halfedge_const_handle;

  typedef std::list<Vertex>                 VertexList;
  typedef typename VertexList::iterator           Vertex_iterator;
  typedef typename VertexList::const_iterator       Vertex_const_iterator;
  typedef Vertex_iterator                 Vertex_handle;
  typedef Vertex_const_iterator               Vertex_const_handle;

  typedef std::list<Face>                 FaceList;
  typedef typename FaceList::iterator           Face_iterator;
  typedef typename FaceList::const_iterator         Face_const_iterator;
  typedef Face_iterator                   Face_handle;
  typedef Face_const_iterator               Face_const_handle;

  typedef std::list<Ccb_halfedge_circulator>         HolesList;
  typedef typename HolesList::iterator          Holes_iterator;
     
  /*
   represents a circulator around vertex that allows for circulating
   the incoming halfedges of a vertex
  */
  class Halfedge_around_vertex_circulator: public Halfedge_handle {
  private:
    typedef Halfedge_handle                        Base;
    typedef Halfedge_around_vertex_circulator              Self;
    // basic arrangement circulator
    typedef typename CGM::Arr_halfedge_around_vertex_const_circulator  PM_Circulator;
    typedef typename CGM::Arr_halfedge_const_iterator      Int_Halfedge_const_iterator;
    // cubical gaussian map circulator
    typedef typename CGM::Halfedge_around_vertex_const_circulator  CGM_Circulator;
    // the corresponding cubical circulator around the vertex
    CGM_Circulator m_cubeCirculator; 
    CGM *m_cgm; // pointer to the cubical gaussian map representing the spherical map
  
  public:
    /*
     empty constructor
    */
    Halfedge_around_vertex_circulator() {}

    /*
     constructor from a cubical map circulator and a cubical map pointer

     inCir - a cubical map circulator
     cgm - a cubical gaussian map of the circulator
    */
    Halfedge_around_vertex_circulator(CGM_Circulator &inCir, CGM *cgm):
      m_cubeCirculator(inCir), m_cgm(cgm) {
      
      Base * base = this;
      Halfedge_handle spEdge;
      // spEdge will hold the spherical map halfedge pointed by the cubical map
      // circulator halfedge, since the target of the cubical map is a real vertex
      // it is guaranteed that the cubical halfedge will have a pointer to 
      // a spherical halfedge
      spEdge = *((Halfedge_handle *)(m_cubeCirculator->getSphereEdge()));
      // set circulator halfedge value to the spherical map halfedge
      *base = spEdge;
    }
    
    /*
     advance to next halfedge around vertex

     return value - the next halfedge around spherical vertex
    */
    Self & operator++() {
      Base * base = this;
      ++m_cubeCirculator; // advance internal cgm circulator
      Halfedge_handle spEdge;
      // get spherical halfedge pointed by cgm advanced circulator
      spEdge = *((Halfedge_handle *)(m_cubeCirculator->getSphereEdge()));
      *base = spEdge; // set value to next halfedge
      return *this;
    }

    /*
     advance to next halfedge around vertex and return current halfedge as if
     the circulator was not advanced

     return value - the current circulator
    */
    Self operator++(int) {
      Self tmp = *this; // save current circulator
      ++(*this); // advance circulator
      return tmp; // return circulator before advancing
    }

    /*
     test equality of two circulators

     circulator - circulator to compare object with
     return value - true if current object and circulator are the same
    */
    bool operator==(const Self &circulator) const {
      // compare the halfedge values with the circulator
      return Base::operator==(circulator);
    }

    /*
     test inequality of two circulators

     circulator - circulator to compare object with
     return value - true if current object and circulator are not the same
    */
    bool operator!=(const Self &circulator) const {
      // use composition of the equality operator
      return (!(*this == circulator));
    }  
  };


  /*
   represents a circulator around a connected component boudary that allows for 
   circulating a ccb of a face
  */  
  class Ccb_halfedge_circulator: public Halfedge_handle {
    typedef Halfedge_handle        Base;
    typedef Ccb_halfedge_circulator    Self;
  
  public:
    /*
     constructor from a regular halfedge to a ccb circulator

     i - a spherical map halfedge
    */
    Ccb_halfedge_circulator(Halfedge_handle i):
      Base(i) {}

    /*
     empty constructor
    */
    Ccb_halfedge_circulator(): Base() {};

    /*
     test for equality of two ccb circulators

     i - the ccb circulator to compare with
     return value - true if the circulators are equal
    */
    bool operator==(const Self &i) const {
      // compare the halfedge part of the two circulators
      return Base::operator==(i);
    }

    /*
     test for inequality of two ccb circulators

     i - the ccb circulator to compare with
     return value - true if the circulators are not equal
    */
    bool operator!=(const Self &i) const {
      // by composition check if the two circulators are not equal
      return (!((*this) == i));
    }

    /*
     advance the circulator

     return value - the next ccb circulator
    */
    Self & operator++() {
      Base *base = this;
      // advancing is moving to the next halfedge
      *base = (*base)->next();
      return *this;
    }

    /*
     advance the circulator and return the circulator before advancing

     return value - the circulator before advancing
    */
    Self operator++(int) {
      Self tmp = *this; // remember circulator before advance is done
      ++(*this); // advance current circulator
      return tmp; // return the circulator before advancing
    }
  };
 
  /*
   this class represents a spherical vertex
  */
  class Vertex {  
  public:
    // public definitions
    typedef typename CGM::Arr_vertex_iterator         Int_Vertex_iterator;
    typedef typename CGM::Arr_halfedge_iterator         Int_Halfedge_iterator;
    typedef typename _Kernel::Direction_3        Direction_3;
    
    /*
     constructor,
      
     inVer - the cubical vertex that represents the sphere vertex
     cgm - pointer to the cubical gaussian map holding the cgm vertex
    */
    Vertex(Int_Vertex_iterator inVer, CGM *cgm):
      m_cubeVertex(inVer), m_cgm(cgm) {};

    ~Vertex() {
    }

    /*
     get the sphere vertex handle pointing to the internal cgm vertex that represents
     the spherical vertex

     return value - a pointer to Vertex_handle as defined in the Spherical_map
      which points to the cubical vertex
    */
    inline void *sphereVertexPointer() {
      return m_cubeVertex->getSphereVertex();
    }

    /*
     test equality of vertices

     sphere_it - the other vertex to compare with
     return value - true if sphere_it and current vertex are the same
    */
    bool operator==(const Vertex & sphere_it) {
      // two sphere vertices are equal if they are represented by the sane
      // cubical vertex
      return (m_cubeVertex == sphere_it.m_cubeVertex);
    }

    /*
     test inequality of vertices

     sphere_it - the other vertex to compare with
     return value - true if sphere_it and current vertex are not the same
    */
    bool operator!=(const Vertex & sphere_it) {
      // use composition with the equality operator
      return (!(*this == sphere_it));
    }

    /*
     the direction of the vertex (the intersection of the direction
     and the sphere is the vertex point)
      
     return value - the vertex direction
    */
    inline Direction_3 direction() {
      // get direction stored in cubical map vertex
      return m_cubeVertex->get_direction();
    }

    /*
     the direction of the vertex (the intersection of the direction
     and the sphere is the vertex point)
      
     return value - the vertex direction
     a const version
    */
    inline const Direction_3 direction() const {
      // get direction stored in cubical map vertex
      return m_cubeVertex->get_direction();
    }

    /*
     returns a circulator that allows to traverse the halfedges that have the vertex
     as their target

     return value - a halfedge circulator around the vertex
    */
    Halfedge_around_vertex_circulator incident_halfedges() const {
      // build a cgm circulator around the corresponding cgm vertex
      typename CGM::Halfedge_around_vertex_const_circulator hec(m_cubeVertex->            incident_halfedges());
      // return a circulator based on the cgm circulator
      return Halfedge_around_vertex_circulator(hec, m_cgm);
    }
    
      /*
     returns true if e is incident to v (i.e., if v is the source or the 
     target of e). 

     e - the halfedge to check incidence to
     return value - true if e is incident to the vertex
    */
    bool is_incident_edge (Halfedge_const_handle e) {
      // the halfedge is incident if it's target or source is the vertex
      if ((e->target())->direction() == direction()) {
        return true;
      }
      if ((e->source())->direction() == direction()) {
        return true;
      }
      return false;
    }
    
    /*
     returns true if face f is incident to the vertex

     f - the face to check if the vertex is incident to
     return value - true if f is incident to the vertex
    */
    bool is_incident_face(Face_const_handle f) {
      /*
      IntFaceList nearFaces;      
      Int_face_handle nearFace;

      nearFace = f->m_cubeFace;
      connectedFaces(*m_cgm, nearFaces, nearFace);
      typename IntFaceList::iterator lit;
      */
      /*
      see if one the the cgm faces of the sphere face is incident to
      the representative vertex
      */
      /*
      for (lit = nearFaces.begin(); lit != nearFaces.end(); lit++) {
      if (m_cubeVertex->is_incident_face(*lit)) {
      return true;
      }
      }
      return false;
      */
      Halfedge_around_vertex_circulator adjEdges=incident_halfedges();
      Halfedge_around_vertex_circulator circEnd=adjEdges;
      // circulate the halfedges around the vertex, for each halfedge check if
      // it's pointed spherical face is f, if so f is incident, if all the surrounding
      // halfedges faces are not f, f is not incident to the vertex
      do {
        if (adjEdges->face() == f) {
          return true;
        }
        ++adjEdges;
      }  while (adjEdges != circEnd);
      return false;
    }

    /*
     find the number of halfedges pointing to the vertex

     return value - the number of incoming halfedges pointing to the vertex
    */
    unsigned int degree() const {
      Halfedge_around_vertex_circulator circ=incident_halfedges();
      Halfedge_around_vertex_circulator startCirc = circ;
      unsigned int counter;
      counter = 0;
      // circulate the halfedges around the vertex and count
      do {
        ++counter;
        ++circ;
      } while (circ != startCirc);
      //      return m_cgm->degree(m_cubeVertex);
      return counter;
    }
 
  private:
    // the cubical vertex projected by the sphere vertex
    Int_Vertex_iterator          m_cubeVertex;
    CGM *m_cgm; // the cubical gaussian map holding the cubical corresponding vertex
  };

  class Halfedge {
  public:
    // public definitions
    typedef typename CGM::Arr_halfedge_handle         Int_halfedge_handle;
    typedef typename CGM::Arr_vertex_handle           Int_vertex_handle;
    typedef Sphere_arc<_Kernel>              Curve;
    
    /*
     constructor, create a new halfedge
      
     inEdge - a halfedge on the cubical gaussian map that is part of the
       sphere halfedge
     cgm - the cubical gaussian map holding the internal halfedge
    */
    Halfedge(Int_halfedge_handle inEdge, CGM *cgm):
      m_cgm(cgm), m_cubeHalfedge(inEdge) {};

    ~Halfedge() {
    }
    
    /*
     get the sphere halfedge handle pointing to the internal cgm halfedge that represents
     the spherical halfedge

     return value - a pointer to Halfedge_handle as defined in the Spherical_map
      which points to the cubical halfedge
    */
    inline void *sphereHalfedgePointer() {
      return m_cubeHalfedge->getSphereEdge();
    }


    /*
     returns the destination vertex of the spherical halfedge

     return value - the spherical handle to the target vertex of the halfedge
    */
    Vertex_handle target() {
      Int_halfedge_handle lEdge;
      lEdge = lastEdge(); // get last cubical halfedge of the spherical halfedge
      Int_vertex_handle targ = lEdge->target(); // get the target of last halfedge
      // find the representative target vertex 
      // which may be the last vertex or one of it's adjacent vertices
      while (!targ->is_rep()) {
        void *tmp = targ->get_adjacent_vertex();
        targ = *((Int_vertex_handle *) (&tmp));
      }
      return *((Vertex_handle *)(targ->getSphereVertex()));
    }

    /*
     returns the curve represented by the halfedge

     return value - the curve represented by the halfedge
    */
    Curve curve() {
      // create a curve with the halfedge endpoints (short arc)
      return Curve(source()->direction(), target()->direction());
    }

    /*
     returns the destination vertex of the spherical halfedge
     a const version of the none const target
    */
    Vertex_const_handle target() const {
      Int_halfedge_handle lEdge;
      lEdge = lastEdge();
      Int_vertex_handle targ = lEdge->target();
      // find the representative target vertex
      while (!targ->is_rep()) {
        void *tmp = targ->get_adjacent_vertex();
        targ = *((Int_vertex_handle *) (&tmp));
      }
      return *((Vertex_const_handle *)(targ->getSphereVertex()));
    }

    /*
     returns the source vertex of the spherical halfedge

     return value - the spherical handle to the source vertex of the halfedge
    */
    Vertex_handle source() {
      Halfedge_handle theTwin;
      theTwin = twin();
      // source is the target of the twin
      return theTwin->target();
    }

    /*
     returns the source vertex of the spherical halfedge            
     a const version of the none const source
    */
    Vertex_const_handle source() const {
      Halfedge_const_handle theTwin;
      theTwin = twin();
      return theTwin->target();
    }

    /*
     get the face next to the halfedge

     return value - a spherical face handle to the face adjacent of the halfedge
    */
    Face_handle face() {
      // get the spherical face handle pointed by the adjacent face to the
      // cubical halfedge
      return *((Face_handle *)((m_cubeHalfedge->face())->getSphereFace()));
    }

    /*
     get the face next to the halfedge, a const version

     return value - a spherical face handle to the face adjacent of the halfedge
    */
    Face_const_handle face() const {
      // get the spherical face handle pointed by the adjacent face to the
      // cubical halfedge
      return *((Face_handle *)((m_cubeHalfedge->face())->getSphereFace()));
    }

    /* 
     the next halfedge around the face

     return value - next face halfedge
    */   
    Halfedge_handle next() {
      Int_halfedge_handle cur; // current halfedge of arc

      // find the last halfedge of the arc
      cur = lastEdge();

      // found last edge, now find it's next
      cur = cur->next();
      // if cube unreal halfedge, traverse adjacent until real halfedge
      while (!cur->get_is_real()) {
        cur = m_cgm->get_adjacent_halfedge_handle(cur);
        cur = cur->next();
      }
      return Face::outHandle(cur, m_cgm);
    }

    /*
     the twin halfedge

     return value - the spherical twin halfedge 
    */
    Halfedge_handle twin() {
      Int_halfedge_handle twin;
      twin = m_cubeHalfedge->twin(); // get cubical twin
      if ((twin->face())->is_unbounded()) {
        // if cubical twin in unbounded face, twin is the adjacent on another
        // cubical face
        twin  = m_cgm->get_adjacent_halfedge_handle(m_cubeHalfedge);
      }

      // find last cubical halfedge of the spherical halfedge
      while (!(twin->target()->getReal())) {
        twin = twin->next(); // advance cubical halfedge
  
        while (!twin->get_is_real()) {
          // traverse cubical unreal halfedges until a real halfedge
          twin = m_cgm->get_adjacent_halfedge_handle(twin);
          twin = twin->next();
        }
    
      }

      return *((Halfedge_handle *)(twin->getSphereEdge()));
    }
    
    /*
     the twin halfedge
     a const version of the none const twin
    */
    Halfedge_const_handle twin() const {
      Int_halfedge_handle twin;
      twin = m_cubeHalfedge->twin(); // get cubical twin
      if ((twin->face())->is_unbounded()) {
        // if cubical twin in unbounded face, twin is the adjacent on another
        // cubical face
        twin  = m_cgm->get_adjacent_halfedge_handle(m_cubeHalfedge);
      }

      // find last cubical halfedge of the spherical halfedge
      while (!(twin->target()->getReal())) {
        twin = twin->next(); // advance cubical halfedge
  
        while (!twin->get_is_real()) {
          // traverse cubical unreal halfedges until a real halfedge
          twin = m_cgm->get_adjacent_halfedge_handle(twin);
          twin = twin->next();
        }
    
      }

      return *((Halfedge_const_handle *)(twin->getSphereEdge()));
    }

    /*
     get a halfedge around ccb circulator starting at the halfedge

     return value - a ccb halfedge circulator starting from the halfedge
    */    
    Ccb_halfedge_circulator ccb() {
      Halfedge_handle thisHandle;
      thisHandle = *((Halfedge_handle *)(m_cubeHalfedge->getSphereEdge()));
      return Ccb_halfedge_circulator(thisHandle);
    }
  
    /*
     get the last cgm halfedge of the spherical halfedge

     return value - the last cgm halfedge of the spherical halfedge
    */
    inline const Int_halfedge_handle cubeHalfedge() const {
      return m_cubeHalfedge;
    }
  private:
      CGM * m_cgm; // the internal cubical gaussian map
    // a cubical gaussian map halfedge corresponding the spherical halfedge,
    // the target of this cgm halfedge is real and it is part of the sphere halfedge arc
    Int_halfedge_handle m_cubeHalfedge; 

    /*
     given a halfedge, finds the last edge of the arc
     that contains the halfedge

     return value - the last cgm halfedge that correspond the sphere halfedge,
       it's cgm target is a real vertex
    */ 
    Int_halfedge_handle lastEdge() const {
      Int_halfedge_handle cur;
      
      cur = m_cubeHalfedge;
      // move to next real cgm halfedge until target is real
      while (!(cur->target()->getReal())) {
        cur = cur->next();
        
        // traverse cubical unreal halfedges until a real halfedge  
        while (!cur->get_is_real()) {
          cur = m_cgm->get_adjacent_halfedge_handle(cur);
          cur = cur->next();
        }
    
      }
      return cur;
    }  
  };

  /*
   this class represents a spherical face
  */
  class Face {
  public:
    // allow the spherical map to access m_doesOuterExist and set it's value
    template <class SphericalMapDcel, class Sphere_traits>
    friend class Spherical_map;

    // public definitions
    typedef typename CGM::Arr_face_handle           Int_face_handle;
    typedef typename CGM::Arr_halfedge_handle         Int_halfedge_handle;
    typedef typename CGM::Arrangement::Holes_iterator  Int_holes_iterator;

    /*
     constructor from a cubical face

     inFace - an internal cubical face that is part of the spherical face
     cgm - the cubical gaussian map containing the cgm face
    */
    Face(Int_face_handle inFace, CGM *cgm): 
      m_cgm(cgm), m_cubeFace(inFace), m_dirty(true), m_doesOuterExist(true) {}
      
    /*
     get the sphere face handle pointing to the internal cgm face that represents
     the spherical face (occupying part of the spherical face area)

     return value - a pointer to Face_handle as defined in the Spherical_map
      which points to a cubical face that is part of the spherical face
    */
    inline void *sphereFacePointer() {
      return m_cubeFace->getSphereFace();
    }

    /*
     get an iterator to first hole in face

     return value - an iterator to the first hole of the face
    */
    Holes_iterator holes_begin() {
      update();
      return m_holes.begin();
    }

    /*
     get an iterator to past the end hole in face

     return value - an iterator to the past the end hole of the face
    */
    Holes_iterator holes_end() {
      update();
      return m_holes.end();
    }

    /*
     get a halfedge on the face outer connected component boundary

     return value - a halfedge on the face outer ccb
    */
    Halfedge_handle halfedge_on_outer_ccb() {
      update();
      return m_outerCcb;
    }

    /*
     get a ccb circulator on the face outer connected component boundary

     return value - a ccb circulator on the face outer ccb
    */
    Ccb_halfedge_circulator outer_ccb() {
      update();
      return m_outerCcb;
    }

    /*
     find if the face has an outer ccb (are there any halfedges in the face)

     return value - true if the face has an outer ccb
    */
    inline bool does_outer_ccb_exist() {
      return m_doesOuterExist;
    }
    
    /*
     find if a halfedge is on the face's outer ccb

     e - the halfedge to check if on outer ccb
     return value - true if e is a halfedge on the face outer ccb
    */
    bool is_halfedge_on_outer_ccb(Halfedge_const_handle e) {
      update();
      if (!m_doesOuterExist) {
        // no outer ccb, halfedge can't be on outer ccb
        return false;
      }
      Ccb_halfedge_circulator circ = m_outerCcb;
      // circulate on the halfedge ccb and check if e is one of the halfedges
      do {
        if (e==circ) {
          return true; // e found on the outer ccb
        }
        ++circ;
      } while (circ != m_outerCcb);
      // e not found on the outer ccb
      return false;
    }

    /*
     find if a halfedge is on a face's inner ccb

     e - the halfedge to check if on an inner ccb
     return value - true if e is a halfedge one of the face's inner ccb
    */
    bool is_halfedge_on_inner_ccb(Halfedge_const_handle e) {
      if (e->face() == *((Face_const_handle *)m_cubeFace->getSphereFace())) {      
        // if halfedge is adjacent to face,
        // using composition, if not on on outer ccb, on inner ccb
        return (!(this->is_halfedge_on_outer_ccb(e)));
      } else {
        // halfedge not adjacent to face
        return false;
      }
    }

    /*
     update the ccb structures of the face
    */
    void update() {    
      if (!m_dirty || !m_doesOuterExist) {
        // already updated or one face of entire cube
        return;
      }
      IntFaceList nearFaces;      
      // find all cgm faces that are part of the spherical face
      connectedFaces(*m_cgm, nearFaces, m_cubeFace);

      typename IntFaceList::iterator lit;      
      bool firstCcb = true; // seperate first ccb as outer ccb
            
      // loop over all cgm faces of the spherical face
      for (lit = nearFaces.begin(); lit != nearFaces.end(); lit++) {
        // process one internal cgm face

        // first process internal face outer ccb
        Int_halfedge_handle curIntHe = (*lit)->outer_ccb();
        Int_halfedge_handle baseIntHe = curIntHe;

        // search for a real halfedge on the cgm face ccb
        // if halfedge is not real but it is marked, it is known
        // that the ccb have already been processed
        while (!curIntHe->get_is_real() && !curIntHe->getMark()) {
          curIntHe = curIntHe->next();
          if (curIntHe == baseIntHe) {
            break;
          }
        }    

        if (curIntHe->getMark() || !curIntHe->get_is_real()) {
          // the ccb has already been processed or isn't real,
          // any ccb will have one cgm face with real halfedge
          // and will be processed at laest in one cgm face
        } else {
          // process new ccb

          // get spherical handle
          Halfedge_handle curHe = outHandle(curIntHe, m_cgm);
          Ccb_halfedge_circulator currentCcb(curHe);
          markCcb(curIntHe); // mark ccb not to be processed again
          if (firstCcb) {
            // spherical outer ccb
            m_outerCcb = currentCcb;
            firstCcb = false;
          } else {
            // a hole
            m_holes.push_back(currentCcb);
          }  
        }
  
        // process holes of internal face    
        Int_holes_iterator holesIt;
        for (holesIt=(*lit)->holes_begin(); 
          holesIt!= (*lit)->holes_end(); 
          holesIt++) {
          
          // process each hole similar to the outer ccb processing
          Int_halfedge_handle curIntHe = (*holesIt);
          Int_halfedge_handle baseIntHe = curIntHe;
          // get a real or marked halfedge on the ccb if any
          while (!curIntHe->get_is_real() && !curIntHe->getMark()) {
            curIntHe = curIntHe->next();
            if (curIntHe == baseIntHe) {
              break;
            }
          }    
    
          if (curIntHe->getMark() || !curIntHe->get_is_real()) {
            // the ccb has already been processed or isn't real,
            // any ccb will have one cgm face with real halfedge
            // and will be processed at laest in one cgm face
            continue;
          } else {
            // process new ccb
            Halfedge_handle curHe = outHandle(curIntHe, m_cgm);      
            Ccb_halfedge_circulator currentCcb(curHe);
            markCcb(curIntHe); // mark ccb not to be processed again

            if (firstCcb) {
              // spherical outer ccb
              m_outerCcb = currentCcb;
              firstCcb = false;
            } else {
              // a hole
              m_holes.push_back(currentCcb);
            }  
          }
        }
      }     
      m_dirty = false;
    }


    
    /*
     find the spherical halfedge from a cgm halfedge

     ihandle - a cubical halfedge
     cgm - the cgm holding the cubical halfedge

     return value - the spherical halfedge handle that ihandle is part of it
    */
    static
    Halfedge_handle outHandle(Int_halfedge_handle ihandle, CGM *cgm) {
      // search for halfedeg with real target
      while (!(ihandle->target()->getReal())) {
        // advance
        ihandle = ihandle->next();
        // if on cube unreal halfedge, traverse adjacents until real halfedge
        while (!ihandle->get_is_real()) {
          ihandle = cgm->get_adjacent_halfedge_handle(ihandle);
          ihandle = ihandle->next();
        }    
      }
      return *((Halfedge_handle *)(ihandle->getSphereEdge()));
    }        
  
  private:

    /*
     mark a connected component boundary starting at ccbHandle

     ccbHandle - a cubical halfedge handle on the ccb to be marked
    */
    void markCcb(Int_halfedge_handle ccbHandle) {
      Int_halfedge_handle cur;
      
      cur = ccbHandle;
      // traverse cgm ccb halfedges
      do {
        cur = cur->next(); // advance a halfedge
        // if on cube unreal halfedge, traverse adjacents until real halfedge
        while (!cur->get_is_real()) {
          cur = m_cgm->get_adjacent_halfedge_handle(cur);
          cur = cur->next();
        }
        cur->setMark(); // mark cgm halfedge
      } while (!(cur == ccbHandle));
    }

    CGM * m_cgm; // the internal cubical gaussian map
    // a cubical gaussian map face corresponding the spherical face, part of the
    // spherical face
    Int_face_handle m_cubeFace;
    // are ccb structures valid or needed to be updated
    bool m_dirty;
    bool m_doesOuterExist; // is there an outer ccb to the spherical face
    Ccb_halfedge_circulator m_outerCcb;
    HolesList  m_holes; // list of face holes
  };
  
private:
public:
  typedef typename CGM::Arr_ccb_halfedge_circulator    Int_Ccb_halfedge_circulator;
  typedef typename CGM::Arr_face_handle             Int_face_handle;
  typedef std::list<Int_face_handle>            IntFaceList;
  
public:
  /*
   find all faces on cube that represents the same sphere face
   
   cgm - the cubical gaussian map
   faces - will hold the cubical faces that represents the sphere face
   startFace - one cubical face of the sphere face
  */
  static void connectedFaces(CGM &cgm, IntFaceList &faces, Int_face_handle &startFace) {
    connectedFacesRec(cgm, faces, startFace); // find connected sub faces
    typename IntFaceList::iterator fit;
    // clear the sub-faces mark
    for (fit = faces.begin(); fit != faces.end(); fit++) {
      (*fit)->setMark(false);
    }
  }

  /*
   find all faces on cube that represents the same sphere face by recursivley 
   traversing all adjacent faces, leave faces marked

   cgm - the cubical gaussian map
   faces - will hold the cubical faces that represents the sphere face
   startFace - one cubical face of the sphere face
  */
  static void connectedFacesRec(CGM &cgm, IntFaceList &faces, 
             Int_face_handle &startFace) {
    if (startFace->getMark() || startFace->is_unbounded()) {
      // face has already been processed or outside the cube
      return;
    }
    startFace->setMark(true); // mark face as processed
    faces.push_back(startFace); // add to list of faces
    Int_Ccb_halfedge_circulator ccbIt, ccbEnd;
    if (startFace->is_unbounded()) {
      // the cgm face should be on the cube
      std::cerr << "unbounded, no boudary" << std::endl;
    }
    ccbIt = startFace->outer_ccb();
    ccbEnd = ccbIt;    
    do {
      if (!ccbIt->get_is_real()) {
        // an unreal boundary, continue in adjacent face
        Int_face_handle nextFace = (cgm.get_adjacent_halfedge_handle(ccbIt))->face();
        // process adjacent cgm face
        connectedFacesRec(cgm, faces, nextFace);
      }
      ccbIt++;
    } while (ccbIt != ccbEnd);
  };
};

CGAL_END_NAMESPACE

#endif
