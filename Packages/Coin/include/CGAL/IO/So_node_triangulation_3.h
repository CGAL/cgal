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
// file          : So_node_triangulation_3.h
// package       : OpenInventor
// author(s)     : Radu Ursu<rursu@sophia.inria.fr>
// release       : 
// release_date  : 
//
// coordinator   : Andreas Fabri<afabri@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_SO_NODE_TRIANGULATION_3_H
#define CGAL_SO_NODE_TRIANGULATION_3_H

#include <CGAL/Triangulation_3.h>

#if HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H
#if HAVE_WINDOWS_H
#include <windows.h> // gl.h needs types from this file on MSWindows.
#endif // HAVE_WINDOWS_H

#include <qgl.h>
#include <Inventor/SbBox.h>
#include <Inventor/caches/SoNormalCache.h>
#include <Inventor/caches/SoBoundingBoxCache.h>
#include <Inventor/SoPrimitiveVertex.h>
#include <Inventor/elements/SoGLCoordinateElement.h>
#include <Inventor/elements/SoNormalBindingElement.h>
#include <Inventor/elements/SoGLLazyElement.h>
#include <Inventor/elements/SoGLTextureCoordinateElement.h>
#include <Inventor/elements/SoGLTextureEnabledElement.h>
#include <Inventor/elements/SoMaterialBindingElement.h>
#include <Inventor/elements/SoModelMatrixElement.h>
#include <Inventor/elements/SoDrawStyleElement.h>
#include <Inventor/elements/SoGLLightModelElement.h>

#include <Inventor/misc/SoState.h>

#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoGetPrimitiveCountAction.h>

#include <Inventor/bundles/SoMaterialBundle.h>
#include <Inventor/bundles/SoTextureCoordinateBundle.h>
#include <Inventor/details/SoLineDetail.h>

#include <Inventor/nodes/SoNonIndexedShape.h>
#include <Inventor/fields/SoMFInt32.h>

template <class Triangulation_3>
class Node_triangulation_3 : public SoNonIndexedShape{
  
  //SO_NODE_HEADER(Node_triangulation_3<Triangulation_3>);  //defined in Inventor/nodes/SoSubNode.h

public:
  typedef typename Triangulation_3::Geom_traits     Traits;
  typedef typename Triangulation_3::Cell_handle		  Cell_handle;
  typedef typename Triangulation_3::Vertex_handle	  Vertex_handle;
  typedef typename Triangulation_3::Locate_type		  Locate_type;
  typedef typename Triangulation_3::Edge_iterator	  Edge_iterator;
  typedef typename Triangulation_3::Vertex_iterator	Vertex_iterator;
  typedef typename Triangulation_3::Finite_vertices_iterator
                                            Finite_vertices_iterator;
  typedef typename Triangulation_3::Finite_edges_iterator
                                            Finite_edges_iterator;
  typedef typename Triangulation_3::Facet_iterator   Facet_iterator;


  static void initClass(){
    do {
      const char * classname = SO__QUOTE(Node_triangulation_3);
      //PRIVATE_COMMON_INIT_CODE(_class_, classname, &_class_::createInstance, _parentclass_);
      do {
        // Make sure we only initialize once.
        assert(Node_triangulation_3::classTypeId == SoType::badType() && "don't init() twice!");
        // Make sure superclass gets initialized before subclass.
        assert(SoShape::getClassTypeId() != SoType::badType() && "you forgot init() on parentclass!");

        // Set up entry in the type system.
        Node_triangulation_3::classTypeId =
        SoType::createType(SoShape::getClassTypeId(),
                         classname,
                         &Node_triangulation_3::createInstance,
                         SoNode::getNextActionMethodIndex());
        SoNode::incNextActionMethodIndex();

        // Store parent's fielddata pointer for later use in the constructor.
        Node_triangulation_3::parentFieldData = SoShape::getFieldDataPtr();
      } while (0);
    } while (0);
  }// Initializes this class
  Node_triangulation_3() : t(t_temp) {
    do {
      Node_triangulation_3::classinstances++;
      // Catch attempts to use a node class which has not been initialized.
      assert(Node_triangulation_3::classTypeId != SoType::badType() && "you forgot init()!");
      // Initialize a fielddata container for the class only once. 
      if (!Node_triangulation_3::fieldData) {
        Node_triangulation_3::fieldData =
          new SoFieldData(Node_triangulation_3::parentFieldData ? \
                        *Node_triangulation_3::parentFieldData : NULL);
      }
      // Extension classes from the application programmers should not be
      // considered native. This is important to get the export code to do
      // the Right Thing. 
      this->isBuiltIn = FALSE;
    } while (0);
  }// The constructor
  Node_triangulation_3(Triangulation_3 &T) : t(T) {
    do {
      Node_triangulation_3::classinstances++;
      // Catch attempts to use a node class which has not been initialized.
      assert(Node_triangulation_3::classTypeId != SoType::badType() && "you forgot init()!");
      // Initialize a fielddata container for the class only once. 
      if (!Node_triangulation_3::fieldData) {
        Node_triangulation_3::fieldData =
          new SoFieldData(Node_triangulation_3::parentFieldData ? \
                        *Node_triangulation_3::parentFieldData : NULL);
      }
      // Extension classes from the application programmers should not be
      // considered native. This is important to get the export code to do
      // the Right Thing. 
      this->isBuiltIn = FALSE;
    } while (0);
  }// The constructor
  
  SoMFInt32 numVertices;

public:
  static SoType getClassTypeId(void){
    return Node_triangulation_3::classTypeId;
  }
  virtual SoType getTypeId(void) const{
    return Node_triangulation_3::classTypeId;
  }
private:
  static SoType classTypeId;
  static void * createInstance(void){
    return new Node_triangulation_3;
  }
protected:
  static const SoFieldData ** getFieldDataPtr(void){
    return (const SoFieldData **)(&Node_triangulation_3::fieldData);
  }
  virtual const SoFieldData * getFieldData(void) const{
    return Node_triangulation_3::fieldData;
  }
private:
  static const SoFieldData ** parentFieldData;
  static const SoFieldData * fieldData;
  static unsigned int classinstances;

protected:
  virtual ~Node_triangulation_3(){}          //The destructor
  
  virtual void  getPrimitiveCount(SoGetPrimitiveCountAction * action){
    if (!this->shouldPrimitiveCount(action)) return;
      int32_t dummyarray[1];
      const int32_t *ptr = this->numVertices.getValues(0);
      const int32_t *end = ptr + this->numVertices.getNum();
      this->fixNumVerticesPointers(action->getState(), ptr, end, dummyarray);

    if (action->canApproximateCount()) {
      action->addNumLines(end-ptr);
    }
    else {
      int cnt = 0;
      while (ptr < end) {
        cnt += *ptr++ - 1;
      }
      action->addNumLines(cnt);
    }
  }

  virtual void  GLRender(SoGLRenderAction *action){
    SoState * state = action->getState();
    SbBool didpush = FALSE;
    if (this->vertexProperty.getValue()) {
      state->push();
      didpush = TRUE;
      this->vertexProperty.getValue()->GLRender(action);
    }
    const SoCoordinateElement * tmp;
    const SbVec3f * normals;
    SbBool needNormals =
      (SoLightModelElement::get(state) !=
      SoLightModelElement::BASE_COLOR);
    SoVertexShape::getVertexData(state, tmp, normals,
                               needNormals);

    if (normals == NULL && needNormals) {
      needNormals = FALSE;
      if (!didpush) {
        state->push();
        didpush = TRUE;
      }
      SoLightModelElement::set(state, SoLightModelElement::BASE_COLOR);
    }

    // First see if the object is visible and should be rendered
    // now. This is a method on SoShape that checks for INVISIBLE
    // draw style, BOUNDING_BOX complexity, and delayed
    // transparency.
    if (!this->shouldGLRender(action)) {
      if (didpush)
        state->pop();
      return;
    }

    SbBool doTextures;
    const SoGLCoordinateElement * coords = (SoGLCoordinateElement *)tmp;

    SoTextureCoordinateBundle tb(action, TRUE, FALSE);
    doTextures = tb.needCoordinates();

    /*
    Binding mbind = findMaterialBinding(action->getState());

    Binding nbind;
    if (!needNormals) nbind = OVERALL;
    else nbind = findNormalBinding(action->getState());
    */
    SoMaterialBundle mb(action);
    mb.sendFirst(); // make sure we have the correct material
  
    glPushMatrix();  
      Finite_vertices_iterator vit;
      glBegin(GL_POINTS);
      for (vit = t.finite_vertices_begin(); vit != t.finite_vertices_end(); ++vit)
        glVertex3f(CGAL::to_double((*vit).point().x()), CGAL::to_double((*vit).point().y()), CGAL::to_double((*vit).point().z()));
      glEnd();

    Finite_edges_iterator eit;
    for (eit = t.finite_edges_begin(); eit != t.finite_edges_end(); ++eit) {
      Point p1( ((*eit).first)->vertex((*eit).second)->point().x(), 
		            ((*eit).first)->vertex((*eit).second)->point().y(), 
		            ((*eit).first)->vertex((*eit).second)->point().z());
      Point p2( ((*eit).first)->vertex((*eit).third)->point().x(), 
		            ((*eit).first)->vertex((*eit).third)->point().y(), 
		            ((*eit).first)->vertex((*eit).third)->point().z());
      glBegin(GL_LINES);
        glVertex3f(CGAL::to_double(p1.x()), CGAL::to_double(p1.y()), CGAL::to_double(p1.z()));
        glVertex3f(CGAL::to_double(p2.x()), CGAL::to_double(p2.y()), CGAL::to_double(p2.z()));
      glEnd();
    }    

    glPopMatrix();
  }
  
  virtual void  computeBBox(SoAction *action,
    SbBox3f &box, SbVec3f &center){
    Finite_vertices_iterator vit;
    typename Traits::FT xmin = 0, ymin = 0, zmin = 0, xmax = 0, ymax = 0, zmax = 0;    
    for (vit = t.finite_vertices_begin(); vit != t.finite_vertices_end(); ++vit) {
      if((*vit).point().x() < xmin)
        xmin = (*vit).point().x();
      if((*vit).point().y() < ymin)
        ymin = (*vit).point().y();
      if((*vit).point().z() < zmin)
        zmin = (*vit).point().z();
      if((*vit).point().x() > xmax)
        xmax = (*vit).point().x();
      if((*vit).point().y() > ymax)
        ymax = (*vit).point().y();
      if((*vit).point().z() > zmax)
        zmax = (*vit).point().z();
      vit++;
    }
    SbVec3f min, max;
    min.setValue(xmin, ymin, zmin);
    max.setValue(xmax, ymax, zmax);
    // Set the box to bound the two extreme points
    box.setBounds(min, max);
  }
  // Generates triangles representing the triangulation
  virtual void  generatePrimitives(SoAction *action){}

private:
  virtual SbBool generateDefaultNormals(SoState *, SoNormalCache * nc){
    // not possible to generate normals for LineSet
    nc->set(0, NULL);
    return TRUE;
  }
  virtual SbBool generateDefaultNormals(SoState * state,
    SoNormalBundle * bundle){return FALSE;}

  Triangulation_3 &t;
  Triangulation_3 t_temp;

};

template<class Triangulation_3>
const SoFieldData ** Node_triangulation_3<Triangulation_3>::parentFieldData = NULL;

template<class Triangulation_3>
unsigned int Node_triangulation_3<Triangulation_3>::classinstances = 0;

template<class Triangulation_3>
const SoFieldData * Node_triangulation_3<Triangulation_3>::fieldData = NULL;

template <class Triangulation_3>
SoType Node_triangulation_3<Triangulation_3>::classTypeId = SoType::badType();


#endif
