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

  typedef Polyhedron_3                              Polyhedron;
  typedef typename Polyhedron_3::Vertex_handle      Vertex_handle;
  typedef typename Polyhedron_3::Vertex_iterator    Vertex_iterator;
  typedef typename Polyhedron_3::Facet_handle       Facet_handle;
  typedef typename Polyhedron_3::Facet_iterator     Facet_iterator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator
                                              Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Traits               Traits;  
  typedef typename Traits::Point_3                  Point;
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
    } while (0);
  }
  virtual SoDetail * copy(void) const {
    SoPolyhedronDetail<Polyhedron_3> * copy = new SoPolyhedronDetail<Polyhedron_3>;
    return copy;
  }

  const SoPrimitiveVertex * get_vertex(int i){
    return &vertex[i];
  }

  Facet_handle find_face(){
    Facet_iterator fit = p.facets_begin();
    int number_of_matching_vertices;
    while(fit != p.facets_end())
    {      
      Halfedge_around_facet_circulator haf = (*fit).facet_begin ();
      bool found_one_equal = false;
      number_of_matching_vertices = 0;
      do{
        Vertex_handle vh = (*haf).vertex();
        found_one_equal = false;
        for(int i=0; i<3; i++)
        {
          const SbVec3f v = vertex[i].getPoint();
          float x = static_cast<float>((*vh).point().x());
          float y = static_cast<float>((*vh).point().y());
          float z = static_cast<float>((*vh).point().z());
          if( x == v[0])
            if( y == v[1])
              if( z == v[2]){
                found_one_equal = true;
                number_of_matching_vertices++;
              }
        }
        haf++;
      }while(found_one_equal && haf != (*fit).facet_begin());
      if(number_of_matching_vertices == 3) //all are equal
        return (fit);
      fit++;
    }    
    exit(1);
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
