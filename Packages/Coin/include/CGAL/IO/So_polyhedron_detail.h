// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// file          : So_polyhedron_detail.h
// package       : OpenInventor
// author(s)     : Radu Ursu<rursu@sophia.inria.fr>
// release       : 
// release_date  : 
//
// coordinator   : Andreas Fabri<afabri@sophia.inria.fr>
//
// ============================================================================


#ifndef CGAL_SO_POLYHEDRON_DETAIL_H
#define CGAL_SO_POLYHEDRON_DETAIL_H

#include <Inventor/details/SoSubDetail.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/SoPrimitiveVertex.h>
#include <Inventor/SbName.h>



template <class Polyhedron_3>
class SoPolyhedronDetail : public SoDetail{
  typedef SoDetail inherited;  
public:

  typedef typename Polyhedron_3                     Polyhedron;
  typedef typename Polyhedron_3::Vertex_handle      Vertex_handle;
  typedef typename Polyhedron_3::Vertex_iterator    Vertex_iterator;
  typedef typename Polyhedron_3::Facet_handle       Facet_handle;
  typedef typename Polyhedron_3::Facet_iterator     Facet_iterator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator
                                              Halfedge_around_facet_circulator;
  //typedef typename Polyhedron_3::Point_iterator          Point_iterator;

  SoPolyhedronDetail() : p(p_temp){};
  SoPolyhedronDetail(SoRayPickAction * action,
                        const SoPrimitiveVertex * v1,
                        const SoPrimitiveVertex * v2,
                        const SoPrimitiveVertex * v3,
                        SoPickedPoint *pp,
                        Polyhedron_3 &P) : p(P){
    vertex[0].setPoint(v1->getPoint());
    vertex[0].setNormal(v1->getNormal());
    vertex[1].setPoint(v2->getPoint());
    vertex[1].setNormal(v2->getNormal());
    vertex[2].setPoint(v3->getPoint());
    vertex[2].setNormal(v3->getNormal());
  }
  virtual ~SoPolyhedronDetail(){};

  static void initClass(void){
    SO_DETAIL_INIT_CLASS(SoPolyhedronDetail, SoDetail);
    do {
      // Make sure we only initialize once.
      assert(SoPolyhedronDetail::classTypeId == SoType::badType());
      // Make sure superclass get initialized before subclass.
      assert(SoDetail::getClassTypeId() != SoType::badType());
    
      SoPolyhedronDetail::classTypeId =
           SoType::createType(SoDetail::getClassTypeId(),
                              SO__QUOTE(SoPolyhedronDetail));
    } while (0)
  }
  virtual SoDetail * copy(void) const {
    SoPolyhedronDetail<Polyhedron_3> * copy = new SoPolyhedronDetail<Polyhedron_3>;
    return copy;
  }

  const SoPrimitiveVertex * get_vertex(int i){
    return &vertex[i];
  }

  const Facet_handle find_face(){    
    Facet_iterator fit = p.facets_begin();    
    while(fit != p.facets_end())
    {      
      int index = 0;
      Halfedge_around_facet_circulator haf = (*fit).facet_begin ();
      bool found = false; //all are equal at the beginning
      do{
        Vertex_handle vh = (*haf).vertex();
        const SbVec3f v = vertex[index].getPoint();
        if( static_cast<float>((*vh).point().x()) != v[0] ||
            static_cast<float>((*vh).point().y()) != v[1] ||
            static_cast<float>((*vh).point().z()) != v[2])
        {
          found = true; //if there is at least one different
        }

        if(index == 2 && !found){          
          return (&(*fit));
        }
        haf++;
        index++;
      }while (haf != (*fit).facet_begin() && !found);

      fit++;
    }
    
    return Facet_handle();
  }

  //defined from SoSubDetail.h because of template declaration of this class
  SoType getTypeId(void) const { return SoPolyhedronDetail::classTypeId; }
  static SoType getClassTypeId(void) { return SoPolyhedronDetail::classTypeId; }
private:
  static SoType classTypeId; //defined from SoSubDetail.h
  SoPrimitiveVertex vertex[3];
  Polyhedron_3 &p;
  Polyhedron_3 p_temp;
};

template <class Polyhedron_3>
SoType SoPolyhedronDetail<Polyhedron_3>::classTypeId = SoType::badType();

#endif