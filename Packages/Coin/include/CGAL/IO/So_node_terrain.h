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
// file          : So_node_terrain.h
// package       : OpenInventor
// author(s)     : Radu Ursu<rursu@sophia.inria.fr>
// release       : 
// release_date  : 
//
// coordinator   : Andreas Fabri<afabri@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_SO_NODE_TERRAIN_H
#define CGAL_SO_NODE_TERRAIN_H

#include <CGAL/Triangulation_2.h>

#if HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H
#if HAVE_WINDOWS_H
#include <windows.h> // gl.h needs types from this file on MSWindows.
#endif // HAVE_WINDOWS_H

#include <qgl.h>

#include <CGAL/Cartesian.h>

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

template <class Triangulation_2>
class Node_terrain : public SoNonIndexedShape{
  
  //SO_NODE_HEADER(Node_terrain<Triangulation_2>);  //defined in Inventor/nodes/SoSubNode.h

public:
  
  typedef typename Triangulation_2::Geom_traits     Traits;
  typedef typename Triangulation_2::Finite_vertices_iterator
                                            Finite_vertices_iterator;
  typedef typename Triangulation_2::Finite_faces_iterator
                                            Finite_faces_iterator;  


  typedef CGAL::Cartesian<double> FT;
  typedef CGAL::Point_3<FT>       CPoint3;
  typedef CGAL::Vector_3<FT>      CVector3;

  static void initClass(){
    do {
      const char * classname = SO__QUOTE(Node_terrain);
      //PRIVATE_COMMON_INIT_CODE(_class_, classname, &_class_::createInstance, _parentclass_);
      do {
        // Make sure we only initialize once.
        assert(Node_terrain::classTypeId == SoType::badType() && "don't init() twice!");
        // Make sure superclass gets initialized before subclass.
        assert(SoShape::getClassTypeId() != SoType::badType() && "you forgot init() on parentclass!");

        // Set up entry in the type system.
        Node_terrain::classTypeId =
        SoType::createType(SoShape::getClassTypeId(),
                         classname,
                         &Node_terrain::createInstance,
                         SoNode::getNextActionMethodIndex());
        SoNode::incNextActionMethodIndex();

        // Store parent's fielddata pointer for later use in the constructor.
        Node_terrain::parentFieldData = SoShape::getFieldDataPtr();
      } while (0);
    } while (0);
  }// Initializes this class
  Node_terrain() : t(t_temp) {
    do {
      Node_terrain::classinstances++;
      // Catch attempts to use a node class which has not been initialized.
      assert(Node_terrain::classTypeId != SoType::badType() && "you forgot init()!");
      // Initialize a fielddata container for the class only once. 
      if (!Node_terrain::fieldData) {
        Node_terrain::fieldData =
          new SoFieldData(Node_terrain::parentFieldData ? \
                        *Node_terrain::parentFieldData : NULL);
      }
      // Extension classes from the application programmers should not be
      // considered native. This is important to get the export code to do
      // the Right Thing. 
      this->isBuiltIn = FALSE;
    } while (0);
  }// The constructor
  Node_terrain(Triangulation_2 &T) : t(T) {
    do {
      Node_terrain::classinstances++;
      // Catch attempts to use a node class which has not been initialized.
      assert(Node_terrain::classTypeId != SoType::badType() && "you forgot init()!");
      // Initialize a fielddata container for the class only once. 
      if (!Node_terrain::fieldData) {
        Node_terrain::fieldData =
          new SoFieldData(Node_terrain::parentFieldData ? \
                        *Node_terrain::parentFieldData : NULL);
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
    return Node_terrain::classTypeId;
  }
  virtual SoType getTypeId(void) const{
    return Node_terrain::classTypeId;
  }
private:
  static SoType classTypeId;
  static void * createInstance(void){
    return new Node_terrain;
  }
protected:
  static const SoFieldData ** getFieldDataPtr(void){
    return (const SoFieldData **)(&Node_terrain::fieldData);
  }
  virtual const SoFieldData * getFieldData(void) const{
    return Node_terrain::fieldData;
  }
private:
  static const SoFieldData ** parentFieldData;
  static const SoFieldData * fieldData;
  static unsigned int classinstances;

protected:
  virtual ~Node_terrain(){}          //The destructor
  
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

    // First see if the object is visible and should be rendered
    // now. This is a method on SoShape that checks for INVISIBLE
    // draw style, BOUNDING_BOX complexity, and delayed
    // transparency.
    if (!this->shouldGLRender(action)) {
      if (didpush)
        state->pop();
      return;
    }

    // Determine if there's a material bound per part
    SoMaterialBindingElement::Binding material_binding = 
        SoMaterialBindingElement::get(state);
    SbBool materialPerPart =
        (material_binding == SoMaterialBindingElement::PER_PART ||
         material_binding == SoMaterialBindingElement::PER_PART_INDEXED);

    // issue a lazy element send.
    // This send will ensure that all material state in GL is current. 
    SoGLLazyElement::sendAllMaterial(state);
  


    glPushMatrix();  
      Finite_vertices_iterator vit;
      glBegin(GL_POINTS);
      for (vit = t.finite_vertices_begin(); vit != t.finite_vertices_end(); ++vit)
        glVertex3f(CGAL::to_double((*vit).point().x()), CGAL::to_double((*vit).point().y()), CGAL::to_double((*vit).point().z()));
      glEnd();

      Finite_faces_iterator fit;
      glBegin(GL_TRIANGLES);      
      for (fit = t.finite_faces_begin(); fit != t.finite_faces_end(); ++fit){
        CPoint3 p1(CGAL::to_double((*(*fit).vertex(0)).point().x()),
                   CGAL::to_double((*(*fit).vertex(0)).point().y()),
                   CGAL::to_double((*(*fit).vertex(0)).point().z()));
        CPoint3 p2(CGAL::to_double((*(*fit).vertex(1)).point().x()),
                   CGAL::to_double((*(*fit).vertex(1)).point().y()),
                   CGAL::to_double((*(*fit).vertex(1)).point().z()));
        CPoint3 p3(CGAL::to_double((*(*fit).vertex(2)).point().x()),
                   CGAL::to_double((*(*fit).vertex(2)).point().y()),
                   CGAL::to_double((*(*fit).vertex(2)).point().z()));

        CVector3 normal = CGAL::cross_product(p2 - p1, p3 - p1);          
        double sqnorm = normal * normal;
        if(sqnorm != 0){
          CVector3 v_n = normal / std::sqrt(sqnorm);
          CPoint3 pn = CPoint3(0, 0, 0) + v_n;
          glNormal3f(pn.x(), pn.y(), pn.z());
        }

        glVertex3f(CGAL::to_double((*(*fit).vertex(0)).point().x()), CGAL::to_double((*(*fit).vertex(0)).point().y()), CGAL::to_double((*(*fit).vertex(0)).point().z()));
        glVertex3f(CGAL::to_double((*(*fit).vertex(1)).point().x()), CGAL::to_double((*(*fit).vertex(1)).point().y()), CGAL::to_double((*(*fit).vertex(1)).point().z()));
        glVertex3f(CGAL::to_double((*(*fit).vertex(2)).point().x()), CGAL::to_double((*(*fit).vertex(2)).point().y()), CGAL::to_double((*(*fit).vertex(2)).point().z()));      
      }
      glEnd();

    glPopMatrix();
  }
  
  virtual void  computeBBox(SoAction *action,
    SbBox3f &box, SbVec3f &center){
    Finite_vertices_iterator vit;
    double xmin = 0, ymin = 0, zmin = 0, xmax = 0, ymax = 0, zmax = 0;    
    for (vit = t.finite_vertices_begin(); vit != t.finite_vertices_end(); ++vit) {
      if(CGAL::to_double((*vit).point().x()) < xmin)
        xmin = CGAL::to_double((*vit).point().x());
      if(CGAL::to_double((*vit).point().y()) < ymin)
        ymin = CGAL::to_double((*vit).point().y());
      if(CGAL::to_double((*vit).point().z()) < zmin)
        zmin = CGAL::to_double((*vit).point().z());
      if(CGAL::to_double((*vit).point().x()) > xmax)
        xmax = CGAL::to_double((*vit).point().x());
      if(CGAL::to_double((*vit).point().y()) > ymax)
        ymax = CGAL::to_double((*vit).point().y());
      if(CGAL::to_double((*vit).point().z()) > zmax)
        zmax = CGAL::to_double((*vit).point().z());
      vit++;
    }
    SbVec3f min, max;
    // Set the box to bound the two extreme points
    min.setValue(xmin, ymin, zmin);
    max.setValue(xmax, ymax, zmax);
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

  Triangulation_2 &t;
  Triangulation_2 t_temp;

};

template<class Triangulation_2>
const SoFieldData ** Node_terrain<Triangulation_2>::parentFieldData = NULL;

template<class Triangulation_2>
unsigned int Node_terrain<Triangulation_2>::classinstances = 0;

template<class Triangulation_2>
const SoFieldData * Node_terrain<Triangulation_2>::fieldData = NULL;

template <class Triangulation_2>
SoType Node_terrain<Triangulation_2>::classTypeId = SoType::badType();


#endif
