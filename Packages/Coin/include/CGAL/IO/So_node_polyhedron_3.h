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
// file          : So_node_polyhedron_3.h
// package       : OpenInventor
// author(s)     : Radu Ursu<rursu@sophia.inria.fr>
// release       : 
// release_date  : 
//
// coordinator   : Andreas Fabri<afabri@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_CGAL_NODE_POLYHEDRON_3_H
#define CGAL_CGAL_NODE_POLYHEDRON_3_H


#include <CGAL/Vector_3.h>
#include <CGAL/Polyhedron_3.h>

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
#include <Inventor/elements/SoTextureCoordinateBindingElement.h>
#include <Inventor/elements/SoModelMatrixElement.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/misc/SoState.h>

#include <Inventor/SbName.h>
#include <Inventor/SoType.h>
#include <Inventor/fields/SoFieldData.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoNonIndexedShape.h>

#include "So_polyhedron_detail.h"
#include <qprogressdialog.h>

#include <map>


template<class Handle>
class my_less{
public:
  bool operator()(const Handle& f1, const Handle& f2) const { return (&(*f1))<(&(*f2));}
};

template <class Polyhedron_3>
class Node_polyhedron_3 : public SoNonIndexedShape{
  
  //SO_NODE_HEADER(Node_polyhedron_3);  
  //defined in Inventor/nodes/SoSubNode.h

public:
  typedef Polyhedron_3                                  Polyhedron;
  typedef typename Polyhedron::Traits                   Traits;  
  typedef typename Traits::Vector_3                     Vector_3;
  typedef typename Polyhedron::Halfedge_handle          Halfedge_handle;
  typedef typename Polyhedron::Facet                    Facet;
  typedef typename Polyhedron::Facet_handle             Facet_handle;
  typedef typename Polyhedron::Facet_const_handle       Facet_const_handle;
  typedef typename Polyhedron::Vertex                   Vertex;
  typedef typename Polyhedron::Vertex_handle            Vertex_handle;
  typedef typename Polyhedron::Vertex_const_handle      Vertex_const_handle;
  typedef typename Polyhedron::Facet_iterator           Facet_iterator;
  typedef typename Polyhedron::Facet_const_iterator     Facet_const_iterator;
  typedef typename Polyhedron::Vertex_iterator          Vertex_iterator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator
                                Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator
                                Halfedge_around_vertex_circulator;
  typedef typename Traits::Point_3                      Point;
  typedef my_less<Vertex_handle>                     Vertex_handle_less;
  typedef my_less<Facet_handle>                      Facet_handle_less;

  //typedef typename Polyhedron::Size                    Size;
  
  //friend bool less(Facet_handle&, Facet_handle&);

  static void initClass(){
    do {
      const char * classname = SO__QUOTE(Node_polyhedron_3);      
      do {
        // Make sure we only initialize once.
        assert(Node_polyhedron_3::classTypeId == SoType::badType() && 
                "don't init() twice!");
        // Make sure superclass gets initialized before subclass.
        assert(SoShape::getClassTypeId() != SoType::badType() && 
                "you forgot init() on parentclass!");

        // Set up entry in the type system.
        Node_polyhedron_3::classTypeId =
        SoType::createType(SoShape::getClassTypeId(),
                         classname,
                         &Node_polyhedron_3::createInstance,
                         SoNode::getNextActionMethodIndex());
        SoNode::incNextActionMethodIndex();

        // Store parent's fielddata pointer for later use in the constructor.
        Node_polyhedron_3::parentFieldData = SoShape::getFieldDataPtr();
      } while (0);
    } while (0);
  };            // Initializes this class
  Node_polyhedron_3() : p(p_temp), LOCK(0){
    do {
      Node_polyhedron_3::classinstances++;
      // Catch attempts to use a node class which has not been initialized.
      assert(Node_polyhedron_3::classTypeId != SoType::badType() &&
                "you forgot init()!");
      // Initialize a fielddata container for the class only once. 
      if (!Node_polyhedron_3::fieldData) {
        Node_polyhedron_3::fieldData =
          new SoFieldData(Node_polyhedron_3::parentFieldData ? \
                        *Node_polyhedron_3::parentFieldData : NULL);
      }
      // Extension classes from the application programmers should not be
      // considered native. This is important to get the export code to do
      // the Right Thing. 
      this->isBuiltIn = FALSE;
    } while (0);
  };                // The constructor
  Node_polyhedron_3(Polyhedron &P) : p(P), LOCK(0){
    do {
      Node_polyhedron_3::classinstances++;
      // Catch attempts to use a node class which has not been initialized.
      assert(Node_polyhedron_3::classTypeId != SoType::badType() && 
                "you forgot init()!");

      compute_normals_for_faces();
      compute_normals_for_vertices();
      // Initialize a fielddata container for the class only once. 
      if (!Node_polyhedron_3::fieldData) {
        Node_polyhedron_3::fieldData =
          new SoFieldData(Node_polyhedron_3::parentFieldData ? \
                        *Node_polyhedron_3::parentFieldData : NULL);
      }
      // Extension classes from the application programmers should not be
      // considered native. This is important to get the export code to do
      // the Right Thing. 
      this->isBuiltIn = FALSE;
    } while (0);
  };   // The constructor
  
  void compute_normals_for_faces(){
    if (LOCK)
      return;
    else
      lock();
    faces_normals.erase(faces_normals.begin(), faces_normals.end());
    QProgressDialog progress( "Computing normals for faces...", 
      "Cancel computing", p.size_of_facets(), NULL, "progress", true );
    progress.setMinimumDuration(0);
    int faces_count = 0;
    Facet_iterator fit;
    for(fit = p.facets_begin(); fit != p.facets_end(); fit++){
      progress.setProgress( faces_count );
      Halfedge_around_facet_circulator h = (*fit).facet_begin();
      Vector_3 normal = CGAL::cross_product(
          h->next()->vertex()->point() - h->vertex()->point(),
          h->next()->next()->vertex()->point() - 
          h->next()->vertex()->point());
      double sqnorm = normal * normal;
      if(sqnorm != 0){
        Vector_3 v_n = normal / std::sqrt(sqnorm);        
        faces_normals[fit] = v_n;
      }      
      faces_count++;
    }
    
    progress.setProgress( faces_count );
    unlock();
  }
  void compute_normals_for_vertices(){
    if (LOCK)
      return;
    else
      lock();    
    //vertices_normals.erase(vertices_normals.begin(), vertices_normals.end());
    QProgressDialog progress( "Computing normals for vertices...", 
      "Cancel computing", p.size_of_vertices(), NULL, "progress", true );
    progress.setMinimumDuration(0);
    int vertices_count = 0;

    Vertex_iterator vit;
    for(vit = p.vertices_begin(); vit != p.vertices_end(); vit++){
      progress.setProgress(vertices_count);
      Halfedge_around_vertex_circulator vh = (*vit).vertex_begin();
      unsigned int normals_count = 0;
      Vector_3 normals_sum(0, 0, 0);
      do{
          Vector_3 normal = CGAL::cross_product(
              vh->next()->vertex()->point() - vh->vertex()->point(),
              vh->next()->next()->vertex()->point() - 
              vh->next()->vertex()->point());
          normals_sum = normals_sum + normal;
          normals_count++;
        }while(++vh != (*vit).vertex_begin());
      normals_sum = normals_sum/normals_count;
      double sqnorm = normals_sum * normals_sum;
      if(sqnorm != 0){
        Vector_3 v_n = normals_sum / std::sqrt(sqnorm);
        vertices_normals[vit] = v_n;
      }
      vertices_count++;
    }    
    unlock();
  }

  void lock(){LOCK++;}
  void unlock(){
    if(LOCK>0)
      LOCK--;
    else
      assert( LOCK != 0 && "lock is already 0. Be sure you have the same \
      number of locks as the number of unlocks");
  }

public:
  static SoType getClassTypeId(void){
    return Node_polyhedron_3::classTypeId;
  }
  virtual SoType getTypeId(void) const{
    return Node_polyhedron_3::classTypeId;
  }
private:
  static SoType classTypeId;
  static void * createInstance(void){
    return new Node_polyhedron_3;
  }
protected:
  static const SoFieldData ** getFieldDataPtr(void){
    return (const SoFieldData **)(&Node_polyhedron_3::fieldData);
  }
  virtual const SoFieldData * getFieldData(void) const{
    return Node_polyhedron_3::fieldData;
  }
private:
  static const SoFieldData ** parentFieldData;
  static const SoFieldData * fieldData;
  static unsigned int classinstances;

protected:
  virtual void  GLRender(SoGLRenderAction *action){
    if (LOCK)
      return;
    else
      lock();
    //std::cout << "called GLRENDER";
    SoState * state = action->getState();

    // First see if the object is visible and should be rendered
    // now. This is a method on SoShape that checks for INVISIBLE
    // draw style, BOUNDING_BOX complexity, and delayed
    // transparency.
    if (! shouldGLRender(action))
      return;
  
    // Determine if we need to send normals. Normals are
    // necessary if we are not doing BASE_COLOR lighting.
    /*     
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
    if(doTextures)
      // Determine if there's a texture bound per vertex or per indexed vertex
      SoTextureCoordinateBindingElement::Binding texture_binding = 
        SoTextureCoordinateBindingElement::get(state);
    

    // Determine if there's a material bound per part
    SoMaterialBindingElement::Binding material_binding = 
        SoMaterialBindingElement::get(state);
    
    SbBool materialPerPart =
        (material_binding == SoMaterialBindingElement::PER_PART ||
         material_binding == SoMaterialBindingElement::PER_PART_INDEXED);
    */
    // issue a lazy element send.
    // This send will ensure that all material state in GL is current. 
    SoGLLazyElement::sendAllMaterial(state);
  
    float complexity = SbClamp(this->getComplexityValue(action), 0.0f, 1.0f);
    //Complexity value, valid settings range 
    //from 0.0 (worst appearance, best perfomance)
    //to 1.0 (optimal appearance, lowest rendering speed)    

    glPushMatrix();

    if(complexity == 0){//render the bounding box
      SbVec3f min = polyhedron_bounding_box.getMin();
      SbVec3f max = polyhedron_bounding_box.getMax();      

      glBegin(GL_QUADS);
      {
        Vector_3 normal = CGAL::cross_product(
          Point(min[0], min[1], max[2]) - Point(min[0], min[1], min[2]),
          Point(min[0], max[1], max[2]) - Point(min[0], min[1], max[2]));          
        double sqnorm = normal * normal;
        if(sqnorm != 0){
          Vector_3 v_n = normal / std::sqrt(sqnorm);
          Point pn = Point(0, 0, 0) + v_n;
          glNormal3f(pn.x(), pn.y(), pn.z());
        }
      }
        glVertex3f(min[0], min[1], min[2]);
        glVertex3f(min[0], min[1], max[2]);
        glVertex3f(min[0], max[1], max[2]);
        glVertex3f(min[0], max[1], min[2]);
      {
        Vector_3 normal = CGAL::cross_product(
          Point(max[0], max[1], min[2]) - Point(max[0], min[1], min[2]),
          Point(max[0], min[1], min[2]) - Point(min[0], min[1], min[2]));          
        double sqnorm = normal * normal;
        if(sqnorm != 0){
          Vector_3 v_n = normal / std::sqrt(sqnorm);
          Point pn = Point(0, 0, 0) + v_n;
          glNormal3f(pn.x(), pn.y(), pn.z());
        }
      }

        glVertex3f(min[0], min[1], min[2]);
        glVertex3f(max[0], min[1], min[2]);
        glVertex3f(max[0], max[1], min[2]);
        glVertex3f(min[0], max[1], min[2]);

      {
        Vector_3 normal = CGAL::cross_product(
          Point(max[0], min[1], max[2]) - Point(min[0], min[1], max[2]),
          Point(max[0], max[1], max[2]) - Point(max[0], min[1], max[2]));
        double sqnorm = normal * normal;
        if(sqnorm != 0){
          Vector_3 v_n = normal / std::sqrt(sqnorm);
          Point pn = Point(0, 0, 0) + v_n;
          glNormal3f(pn.x(), pn.y(), pn.z());
        }
      }
        glVertex3f(min[0], min[1], max[2]);
        glVertex3f(max[0], min[1], max[2]);
        glVertex3f(max[0], max[1], max[2]);
        glVertex3f(min[0], max[1], max[2]);

      {
        Vector_3 normal = CGAL::cross_product(
          Point(max[0], min[1], max[0]) - Point(min[0], min[1], max[0]),
          Point(min[0], min[1], max[0]) - Point(min[0], min[1], min[2]));          
        double sqnorm = normal * normal;
        if(sqnorm != 0){
          Vector_3 v_n = normal / std::sqrt(sqnorm);
          Point pn = Point(0, 0, 0) + v_n;
          glNormal3f(pn.x(), pn.y(), pn.z());
        }
      }

        glVertex3f(min[0], min[1], min[2]);
        glVertex3f(min[0], min[1], max[2]);
        glVertex3f(max[0], min[1], max[2]);
        glVertex3f(max[0], min[1], min[2]);

      {
        Vector_3 normal = CGAL::cross_product(
          Point(max[0], max[1], max[2]) - Point(max[0], min[1], max[2]),
          Point(max[0], min[1], max[2]) - Point(max[0], min[1], min[2]));          
        double sqnorm = normal * normal;
        if(sqnorm != 0){
          Vector_3 v_n = normal / std::sqrt(sqnorm);
          Point pn = Point(0, 0, 0) + v_n;
          glNormal3f(pn.x(), pn.y(), pn.z());
        }
      }

        glVertex3f(max[0], min[1], min[2]);
        glVertex3f(max[0], min[1], max[2]);
        glVertex3f(max[0], max[1], max[2]);
        glVertex3f(max[0], max[1], min[2]);
       
      {
        Vector_3 normal = CGAL::cross_product(
          Point(min[0], max[1], min[2]) - Point(min[0], max[1], max[2]),
          Point(min[0], max[1], max[2]) - Point(max[0], max[1], max[2]));
        double sqnorm = normal * normal;
        if(sqnorm != 0){
          Vector_3 v_n = normal / std::sqrt(sqnorm);
          Point pn = Point(0, 0, 0) + v_n;
          glNormal3f(pn.x(), pn.y(), pn.z());
        }
      }

        glVertex3f(max[0], max[1], max[2]);
        glVertex3f(min[0], max[1], max[2]);
        glVertex3f(min[0], max[1], min[2]);
        glVertex3f(max[0], max[1], min[2]);

      glEnd();
    } else if(complexity==1){
      //render in smooth (specify the normals for the vertices)
      Facet_iterator fit = p.facets_begin();
      while(fit != p.facets_end()){
        Halfedge_around_facet_circulator h = (*fit).facet_begin();
        glBegin(GL_POLYGON);
          do{
            Point pn = Point(0, 0, 0) + vertices_normals[h->vertex()];
            glNormal3f(pn.x(), pn.y(), pn.z());
            Point point = h->vertex()->point();
            glVertex3f(point[0],point[1],point[2]);
          }while(++h != (*fit).facet_begin());
        glEnd();        
        fit++;          
      }//end while
    } else { //render specifying vertices only for normals
      double nr_of_facets = p.size_of_facets();
      Facet_iterator fit = p.facets_begin();
      unsigned int complexity_count = 0;
      unsigned int complexity_step = (unsigned int)(1/complexity);
      while(fit != p.facets_end() && complexity_count < nr_of_facets){
        glBegin(GL_POLYGON);        
          Point pn = Point(0, 0, 0) + faces_normals[fit];
          glNormal3f(pn.x(), pn.y(), pn.z());
          Halfedge_around_facet_circulator h = (*fit).facet_begin();
          do{      
            Point point = h->vertex()->point();
            glVertex3f(point[0],point[1],point[2]);
          }while(++h != (*fit).facet_begin());
        glEnd();
        for(unsigned int step_i = 0; step_i < complexity_step; step_i++)
          fit++;          
        complexity_count += complexity_step;
      }//end while
    }
    glPopMatrix();
    unlock();
  };
  
  virtual void  computeBBox(SoAction *action, SbBox3f &box, SbVec3f &center){
    Vertex_iterator vit = p.vertices_begin();
    Kernel::FT xmin = 0, ymin = 0, zmin = 0, xmax = 0, ymax = 0, zmax = 0;
    while(vit != p.vertices_end()){
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
    polyhedron_bounding_box.setBounds(min, max);
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

    //SoMaterialBindingElement::Binding bind = 
    //  SoMaterialBindingElement::get(action->getState());
    //float complexity = this->getComplexityValue(action);

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
    


    //Generate POLYGON primitives
    Facet_iterator fit = p.facets_begin();
    while(fit != p.facets_end()){
      Halfedge_around_facet_circulator h = (*fit).facet_begin();
    
      Vector_3 normal = CGAL::cross_product(
        h->next()->vertex()->point() - h->vertex()->point(),
        h->next()->next()->vertex()->point() - h->next()->vertex()->point());

      beginShape(action, POLYGON);        
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
          GEN_VERTEX(pv, static_cast<float>(point[0]), static_cast<float>(point[1]),  static_cast<float>(point[2]), .25,  0.0, sbnormal);
        }while(++h != (*fit).facet_begin());
      endShape();
      fit++;
    }//end while
/*
    //Generate POINTS primitives
    Vertex_iterator vit = p.vertices_begin();
    while(vit != p.vertices_end()){
      //Halfedge_around_facet_circulator h = (*fit).facet_begin();
    
      //Vector_3 normal = CGAL::cross_product(
        //h->next()->vertex()->point() - h->vertex()->point(),
        //h->next()->next()->vertex()->point() - h->next()->vertex()->point());

      beginShape(action, POINTS);
        //double sqnorm = normal * normal;
        SbVec3f sbnormal;        
        //if(sqnorm != 0){
        //  Vector_3 v_n = normal / std::sqrt(sqnorm);
        //  Point pn = Point(0, 0, 0) + v_n;
        //  //glNormal3f(pn.x(), pn.y(), pn.z());
        //  sbnormal.setValue(pn.x(), pn.y(), pn.z());
        //}
          Point point = (*vit).point();
          //glVertex3f(point[0],point[1],point[2]);          
          GEN_VERTEX(pv, point[0], point[1],  point[2], .25,  0.0, sbnormal);
      endShape();
      vit++;
    }//end while
  */


  };


  //    The following method is used to create triangle detail
  //    Ex. :
  //    SoRayPickAction rp(viewer->getViewportRegion());
  //    rp.setPoint(mbe->getPosition());
  //    rp.apply(viewer->getSceneManager()->getSceneGraph());
  //    Now the detail instance was constructed and stored in the *rp* object
  //    You don't have to delete it, it will be automatically deleted
  //    SoPickedPoint * point = rp.getPickedPoint();

  virtual SoDetail* 
  createTriangleDetail( SoRayPickAction * action,
                        const SoPrimitiveVertex * v1,
                        const SoPrimitiveVertex * v2,
                        const SoPrimitiveVertex * v3,
                        SoPickedPoint *pp){
    
    SoPolyhedronDetail<Polyhedron_3> 
      *copy = new SoPolyhedronDetail<Polyhedron_3>(action, v1, v2, v3, pp, p);
    return copy;
  }

  //    The following method is used to create point detail
  virtual SoDetail *
  createPointDetail(SoRayPickAction * action ,
                           const SoPrimitiveVertex * v,
                           SoPickedPoint * pp)
  {
    SoPolyhedronDetail<Polyhedron_3> 
      *copy = new SoPolyhedronDetail<Polyhedron_3>(action, v, v, v, pp, p);
    return copy;  
  }

private:
  virtual ~Node_polyhedron_3(){};       //The destructor
  
  std::map<Vertex_handle, Vector_3, Vertex_handle_less>
                                      vertices_normals;  
  std::map<Facet_handle, Vector_3, Facet_handle_less>
                                      faces_normals;
  Polyhedron &p;
  Polyhedron p_temp;
  int LOCK;   //used to secure rendering
              //the node is rendered only when the data is ready

  SbBox3f    polyhedron_bounding_box;
};

template<class Polyhedron_3>
const SoFieldData ** Node_polyhedron_3<Polyhedron_3>::parentFieldData = NULL;

template<class Polyhedron_3>
unsigned int Node_polyhedron_3<Polyhedron_3>::classinstances = 0;

template<class Polyhedron_3>
const SoFieldData * Node_polyhedron_3<Polyhedron_3>::fieldData = NULL;

template <class Polyhedron_3>
SoType Node_polyhedron_3<Polyhedron_3>::classTypeId = SoType::badType();

#endif
