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
// file          : So_node_tetrahedron_3.h
// package       : OpenInventor
// author(s)     : Radu Ursu<rursu@sophia.inria.fr>
// release       : 
// release_date  : 
//
// coordinator   : Andreas Fabri<afabri@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_SO_NODE_TETRAHEDRON_3_H
#define CGAL_SO_NODE_TETRAHEDRON_3_H


#include <CGAL/Vector_3.h>
#include <CGAL/Tetrahedron_3.h>

#if HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H
#if HAVE_WINDOWS_H
#include <windows.h> // gl.h needs types from this file on MSWindows.
#endif // HAVE_WINDOWS_H

#include <qgl.h>
#include <Inventor/SbLinear.h>
#include <Inventor/SbBox.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/SoPrimitiveVertex.h>
#include <Inventor/elements/SoGLLazyElement.h>
#include <Inventor/elements/SoGLTextureCoordinateElement.h>
#include <Inventor/elements/SoGLTextureEnabledElement.h>
#include <Inventor/elements/SoMaterialBindingElement.h>
#include <Inventor/elements/SoModelMatrixElement.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/misc/SoState.h>

#include <Inventor/SbName.h>
#include <Inventor/SoType.h>
#include <Inventor/fields/SoFieldData.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoNonIndexedShape.h>

#include "SoPolyhedronDetail.h"


template <class Kernel>
class Node_tetrahedron_3 : public SoNonIndexedShape{
  
  //SO_NODE_HEADER(Node_tetrahedron_3);  //defined in Inventor/nodes/SoSubNode.h

public:
  typedef typename CGAL::Tetrahedron_3<Kernel>         Tetrahedron;
  typedef typename CGAL::Vector_3<Kernel>              Vecor_3;
  typedef typename CGAL::Point_3<Kernel>               Point_3;

  static void initClass(){
    do {
      const char * classname = SO__QUOTE(Node_tetrahedron_3);
      //PRIVATE_COMMON_INIT_CODE(_class_, classname, &_class_::createInstance, _parentclass_);
      do {
        // Make sure we only initialize once.
        assert(Node_tetrahedron_3::classTypeId == SoType::badType() && "don't init() twice!");
        // Make sure superclass gets initialized before subclass.
        assert(SoShape::getClassTypeId() != SoType::badType() && "you forgot init() on parentclass!");

        // Set up entry in the type system.
        Node_tetrahedron_3::classTypeId =
        SoType::createType(SoShape::getClassTypeId(),
                         classname,
                         &Node_tetrahedron_3::createInstance,
                         SoNode::getNextActionMethodIndex());
        SoNode::incNextActionMethodIndex();

        // Store parent's fielddata pointer for later use in the constructor.
        Node_tetrahedron_3::parentFieldData = SoShape::getFieldDataPtr();
      } while (0);
    } while (0);
  };            // Initializes this class
  Node_tetrahedron_3() : t(t_temp){
    do {
      Node_tetrahedron_3::classinstances++;
      // Catch attempts to use a node class which has not been initialized.
      assert(Node_tetrahedron_3::classTypeId != SoType::badType() && "you forgot init()!");
      // Initialize a fielddata container for the class only once. 
      if (!Node_tetrahedron_3::fieldData) {
        Node_tetrahedron_3::fieldData =
          new SoFieldData(Node_tetrahedron_3::parentFieldData ? \
                        *Node_tetrahedron_3::parentFieldData : NULL);
      }
      // Extension classes from the application programmers should not be
      // considered native. This is important to get the export code to do
      // the Right Thing. 
      this->isBuiltIn = FALSE;
    } while (0);
  };                // The constructor
  Node_tetrahedron_3(Tetrahedron &T) : t(T){
    do {
      Node_tetrahedron_3::classinstances++;
      // Catch attempts to use a node class which has not been initialized.
      assert(Node_tetrahedron_3::classTypeId != SoType::badType() && "you forgot init()!");
      // Initialize a fielddata container for the class only once. 
      if (!Node_tetrahedron_3::fieldData) {
        Node_tetrahedron_3::fieldData =
          new SoFieldData(Node_tetrahedron_3::parentFieldData ? \
                        *Node_tetrahedron_3::parentFieldData : NULL);
      }
      // Extension classes from the application programmers should not be
      // considered native. This is important to get the export code to do
      // the Right Thing. 
      this->isBuiltIn = FALSE;
    } while (0);
  };   // The constructor

public:
  static SoType getClassTypeId(void){
    return Node_tetrahedron_3::classTypeId;
  }
  virtual SoType getTypeId(void) const{
    return Node_tetrahedron_3::classTypeId;
  }
private:
  static SoType classTypeId;
  static void * createInstance(void){
    return new Node_tetrahedron_3;
  }
protected:
  static const SoFieldData ** getFieldDataPtr(void){
    return (const SoFieldData **)(&Node_tetrahedron_3::fieldData);
  }
  virtual const SoFieldData * getFieldData(void) const{
    return Node_tetrahedron_3::fieldData;
  }
private:
  static const SoFieldData ** parentFieldData;
  static const SoFieldData * fieldData;
  static unsigned int classinstances;

protected:
  virtual void  GLRender(SoGLRenderAction *action){
    SoState * state = action->getState();

    // First see if the object is visible and should be rendered
    // now. This is a method on SoShape that checks for INVISIBLE
    // draw style, BOUNDING_BOX complexity, and delayed
    // transparency.
    if (! shouldGLRender(action))
      return;
  
    // Determine if we need to send normals. Normals are
    // necessary if we are not doing BASE_COLOR lighting.
         
    // we use the Lazy element to get the light model.   
    SbBool sendNormals = (SoLazyElement::getLightModel(state) 
	    != SoLazyElement::BASE_COLOR);      

    // See if texturing is enabled. If so, we will have to
    // send explicit texture coordinates. The "doTextures" flag
    // will indicate if we care about textures at all.
   
    // Note this has changed slightly in Inventor version 2.1.
    // The texture coordinate type now is either FUNCTION or DEFAULT.
    // Texture coordinates are needed only for DEFAULT textures.

    SbBool doTextures =
      ( SoGLTextureEnabledElement::get(state) &&
        SoTextureCoordinateElement::getType(state) !=
        SoTextureCoordinateElement::FUNCTION);

    // Determine if there's a material bound per part
    SoMaterialBindingElement::Binding binding = 
        SoMaterialBindingElement::get(state);
    SbBool materialPerPart =
        (binding == SoMaterialBindingElement::PER_PART ||
        binding == SoMaterialBindingElement::PER_PART_INDEXED);

    // issue a lazy element send.
    // This send will ensure that all material state in GL is current. 
    SoGLLazyElement::sendAllMaterial(state);
  
    double sqnorm;
    Vector_3 v_n, normal;
    Point_3 p0, p1, p2, p3;

    
    if(t.orientation() == CGAL::NEGATIVE){
      p0 = t[0]; p1 = t[2]; p2 = t[1]; p3 = t[3];
    } else {
      p0 = t[0]; p1 = t[1]; p2 = t[2]; p3 = t[3];
    }

    glPushMatrix();
      glBegin(GL_TRIANGLES);
        normal = CGAL::cross_product(p2 - p0, p1 - p0);
        sqnorm = normal * normal;
        if(sqnorm != 0){
          v_n = normal / std::sqrt(sqnorm);
          glNormal3f(v_n.x(), v_n.y(), v_n.z());          
        }        
        glVertex3f(p0.x(),p0.y(),p0.z());
        glVertex3f(p1.x(),p1.y(),p1.z());
        glVertex3f(p2.x(),p2.y(),p2.z());

        normal = CGAL::cross_product(p3 - p0, p2 - p0);
        sqnorm = normal * normal;
        if(sqnorm != 0){
          v_n = normal / std::sqrt(sqnorm);          
          glNormal3f(v_n.x(), v_n.y(), v_n.z());
        }

        glVertex3f(p0.x(),p0.y(),p0.z());
        glVertex3f(p2.x(),p2.y(),p2.z());
        glVertex3f(p3.x(),p3.y(),p3.z());

        normal = CGAL::cross_product(p1 - p0,
          p3 - p0);
        sqnorm = normal * normal;
        if(sqnorm != 0){
          v_n = normal / std::sqrt(sqnorm);          
          glNormal3f(v_n.x(), v_n.y(), v_n.z());
        }        
        glVertex3f(p0.x(),p0.y(),p0.z());
        glVertex3f(p1.x(),p1.y(),p1.z());
        glVertex3f(p3.x(),p3.y(),p3.z());

        normal = CGAL::cross_product(p2 - p1,
          p3 - p1);
        sqnorm = normal * normal;
        if(sqnorm != 0){
          v_n = normal / std::sqrt(sqnorm);          
          glNormal3f(v_n.x(), v_n.y(), v_n.z());
        }        
        glVertex3f(p1.x(),p1.y(),p1.z());
        glVertex3f(p2.x(),p2.y(),p2.z());
        glVertex3f(p3.x(),p3.y(),p3.z());

      glEnd();
    glPopMatrix();
  };
  
  virtual void  computeBBox(SoAction *action,
    SbBox3f &box, SbVec3f &center){
    Kernel::FT xmin = 0, ymin = 0, zmin = 0, xmax = 0, ymax = 0, zmax = 0;
    for(int i=0; i<4; i++){
      if(t.vertex(i).x() < xmin)
        xmin = t.vertex(i).x();
      if(t.vertex(i).y() < ymin)
        ymin = t.vertex(i).y();
      if(t.vertex(i).z() < zmin)
        zmin = t.vertex(i).z();
      if(t.vertex(i).x() > xmax)
        xmax = t.vertex(i).x();
      if(t.vertex(i).y() > ymax)
        ymax = t.vertex(i).y();
      if(t.vertex(i).z() > zmax)
        zmax = t.vertex(i).z();
      
    }    
    SbVec3f min, max;
    min.setValue(xmin, ymin, zmin);
    max.setValue(xmax, ymax, zmax);
    // Set the box to bound the two extreme points
    box.setBounds(min, max);
  };
  
  
  virtual void getPrimitiveCount(SoGetPrimitiveCountAction *action)
  {
    if (!this->shouldPrimitiveCount(action)) return;
    else
      SoShape::getPrimitiveCount(action);
  }


  // Generates triangles representing a polyhedron
  virtual void  generatePrimitives(SoAction *action){
    SoPrimitiveVertex   pv;

    // Access the state from the action
    SoState  *state = action->getState();    
    
    // See if we have to use a texture coordinate function,
    // rather than generating explicit texture coordinates.
    SbBool useTexFunc = 
      (SoTextureCoordinateElement::getType(state) ==
       SoTextureCoordinateElement::FUNCTION);
    
    // If we need to generate texture coordinates with a
    // function, we'll need an SoGLTextureCoordinateElement.
    // Otherwise, we'll set up the coordinates directly.
    const SoTextureCoordinateElement *tce;
    SbVec4f texCoord;
    if (useTexFunc)
      tce = SoTextureCoordinateElement::getInstance(state);
    else {
      texCoord[2] = 0.0;
      texCoord[3] = 1.0;
    }

    SoMaterialBindingElement::Binding bind = SoMaterialBindingElement::get(action->getState());
    float complexity = this->getComplexityValue(action);

    // We'll use this macro to make the code easier. It uses the
    // "point" variable to store the primitive vertex's point.
    SbVec3f  sbpoint;

    #define GEN_VERTEX(pv, x, y, z, s, t, sbnormal) \
      sbpoint.setValue(x, y, z);                    \
      if (useTexFunc)                               \
         texCoord = tce->get(sbpoint, sbnormal);    \
      else {                                        \
         texCoord[0] = s;                           \
         texCoord[1] = t;                           \
          }                                         \
      pv.setPoint(sbpoint);                         \
      pv.setNormal(sbnormal);                       \
      pv.setTextureCoords(texCoord);                \
      shapeVertex(&pv)
    
/*
    Facet_iterator fit = p.facets_begin();
    while(fit != p.facets_end()){
      Halfedge_around_facet_circulator h = (*fit).facet_begin();
    
      Vector_3 normal = CGAL::cross_product(
        h->next()->vertex()->point() - h->vertex()->point(),
        h->next()->next()->vertex()->point() - h->next()->vertex()->point());

      beginShape(action, TRIANGLES);        
        double sqnorm = normal * normal;
        SbVec3f sbnormal;
        if(sqnorm != 0){
          Vector_3 v_n = normal / std::sqrt(sqnorm);
          Point pn = Point(0, 0, 0) + v_n;
          //glNormal3f(pn.x(), pn.y(), pn.z());
          sbnormal.setValue(pn.x(), pn.y(), pn.z());
        }
        do{      
          Point point = h->vertex()->point();
          //glVertex3f(point[0],point[1],point[2]);          
          GEN_VERTEX(pv, point[0], point[1],  point[2], .25,  0.0, sbnormal);
        }while(++h != (*fit).facet_begin());
      endShape();
      fit++;
    }//end while
*/
  };


private:
  virtual ~Node_tetrahedron_3(){};       //The destructor

  Tetrahedron &t;
  Tetrahedron t_temp;

};

template<class Tetrahedron_3>
const SoFieldData ** Node_tetrahedron_3<Tetrahedron_3>::parentFieldData = NULL;

template<class Tetrahedron_3>
unsigned int Node_tetrahedron_3<Tetrahedron_3>::classinstances = 0;

template<class Tetrahedron_3>
const SoFieldData * Node_tetrahedron_3<Tetrahedron_3>::fieldData = NULL;

template <class Tetrahedron_3>
SoType Node_tetrahedron_3<Tetrahedron_3>::classTypeId = SoType::badType();

#endif