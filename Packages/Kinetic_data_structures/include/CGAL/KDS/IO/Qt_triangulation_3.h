#ifndef CGAL_COIN_KINETIC_TRIANGULATION_3_H
#define CGAL_COIN_KINETIC_TRIANGULATION_3_H

#include <CGAL/KDS/basic.h>

// Inventor parts
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodekits/SoShapeKit.h>
#include <Inventor/nodekits/SoAppearanceKit.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoShapeHints.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#include <Inventor/nodes/SoIndexedLineSet.h>
#include <Inventor/nodes/SoMaterialBinding.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/SoOutput.h>
#include <Inventor/actions/SoWriteAction.h>

#include <Inventor/events/SoEvent.h>
#include <Inventor/events/SoButtonEvent.h>
#include <Inventor/events/SoKeyboardEvent.h>
#include <Inventor/nodes/SoEventCallback.h>

#include <CGAL/KDS/IO/Coin_pointer.h>

#include <CGAL/KDS/internal/triangulation_helpers_3.h>
#include <CGAL/KDS/Ref_counted.h>

CGAL_KDS_BEGIN_NAMESPACE


//! A class to display a 3D triangulation.
/*!
  
*/
template <class KDel, class Qtgui, class Qtmpt> 
class Qt_triangulation_3: public Ref_counted<Qt_triangulation_3<KDel, Qtgui, Qtmpt> > {
protected:

  typedef Qt_triangulation_3<KDel, Qtgui, Qtmpt> This;
  typedef typename KDel::Triangulation::Facet Facet;
  typedef typename KDel::Triangulation::Edge Edge;
  typedef typename KDel::Triangulation::Facet_circulator Facet_circulator;
  typedef typename KDel::Triangulation::Cell_circulator Cell_circulator;
  typedef typename KDel::Triangulation::Finite_edges_iterator Finite_edges_iterator;
  typedef typename KDel::Triangulation::Finite_facets_iterator Finite_facets_iterator;
  typedef typename KDel::Triangulation::Edge_iterator Edge_iterator;
  typedef typename KDel::Triangulation::Facet_iterator Facet_iterator;  
  typedef typename KDel::Triangulation::Vertex_handle Vertex_handle;
  typedef typename KDel::Triangulation::Cell_handle Cell_handle;
  typedef typename KDel::Triangulation::All_cells_iterator All_cells_iterator;
  typedef typename KDel::Triangulation::All_edges_iterator All_edges_iterator;

  typedef typename KDel::Triangulation::Geom_traits::Point_3 Object_key;

   
  // I just want the root()
  class Guil: public Qtgui::Listener {
  public:
    Guil(typename Qtgui::Pointer& h): Qtgui::Listener(h){}
    void new_notification(typename Qtgui::Listener::Notification_type ){
    }
  };

  class KDell: public KDel::Listener {
  public:
    KDell(typename KDel::Pointer& h, This *t): KDel::Listener(h), t_(t){}
    void new_notification(typename KDel::Listener::Notification_type ){
      t_->generate_geometry();
    }
  protected:
    This *t_;
  };
  typedef enum {NO_CERT=0, UNFAILING_CERT=1, CERT=2, NEXT_CERT=3, HIDE=-1} Color_id;

  typename KDel::Triangulation::Vertex_handle facet_vertex(typename KDel::Triangulation::Facet f, int i) {
    return vertex_of_facet(f, i);
  }

  typename KDel::Triangulation::Vertex_handle edge_vertex(typename KDel::Triangulation::Edge f, int i) {
    return vertex_of_facet(f, i);
  }


public:

  
  //! The type for draw styles.
  /*!
    The currently supported styles are HIDDEN, SHOWN and TRANSPARENT.
  */
  typedef enum {HIDDEN, SHOWN, TRANSPARENT} Draw_style;
  
  //! Initialize everything.
  /*!
    Sorry about all the arguments, but they are needed, I think.
  */
  Qt_triangulation_3(typename KDel::Pointer kdel, 
		     typename Qtgui::Pointer qtgui,
		     typename Qtmpt::Pointer mps): coordinates_(mps->coordinate_node()),
						   convex_hull_(SHOWN), 
						   facets_style_(SHOWN),
						   kdell_(kdel, this), 
						   guil_(qtgui) {
    set_scene_graph_parent(guil_.root());
    //if (0) kk.orientation_3_object(); 
    // for some reason it does not parse if I remove kk above, I want to get rid of the variable not used warning
  }
 
  //! The field for how the convex hull is drawn.
  Draw_style convex_hull_draw_style() {
    return convex_hull_;
  }

  //! Setting the convex hull draw style field.
  void set_convex_hull_draw_style(Draw_style ds) {
    convex_hull_=ds;
    generate_geometry();
  }
  
  //! The draw style for non-convex-hull faces.
  Draw_style facets_draw_style() const {
	return facets_style_;
  }
  //! Set the draw style for non0convex hull faces
  void set_facets_draw_style(Draw_style ds) {
    facets_style_= ds;
    generate_geometry();
  }
  
protected:
  CGAL::Coin_pointer<SoSeparator> parent_;
  CGAL::Coin_pointer<SoShapeKit> facets_kit_;
  CGAL::Coin_pointer<SoIndexedFaceSet> facets_;
  CGAL::Coin_pointer<SoShapeKit> edges_kit_;
  CGAL::Coin_pointer<SoIndexedLineSet> edges_;
  CGAL::Coin_pointer<SoCoordinate3> coordinates_;
  Draw_style convex_hull_;
  Draw_style facets_style_;
  KDell kdell_;
  Guil guil_;

  const typename KDel::Triangulation& triangulation() const {
    return kdell_.notifier()->triangulation();
  }

  void set_scene_graph_parent(SoSeparator* sep);

  static void keyboard_callback(void *data, SoEventCallback *eventCB){
    This *th = reinterpret_cast<This*>(data);
    const SoEvent *event= eventCB->getEvent();
    CGAL_assertion(event->isOfType(SoKeyboardEvent::getClassTypeId()));
    const SoKeyboardEvent *kbe= reinterpret_cast<const SoKeyboardEvent*>(event);
    //std::cout << "Pressed " << kbe->getPrintableCharacter() << std::endl;
    if (kbe->getKey()== SoKeyboardEvent::H && kbe->getState()== SoButtonEvent::UP) {
      Draw_style ds= th->convex_hull_draw_style();
      switch (ds){
      case HIDDEN: ds= SHOWN; break;
      case SHOWN: ds= HIDDEN; break;
      default: // disable error message
	; // this is needed to compile
      }
      eventCB->setHandled();
      //std::cout << "H pushed.\n";
      th->set_convex_hull_draw_style(ds);
    } else if (kbe->getKey() == SoKeyboardEvent::F && kbe->getState()==SoButtonEvent::UP){
	 Draw_style ds= th->facets_draw_style();
      switch (ds){
      case HIDDEN: ds= SHOWN; break;
      case SHOWN: ds= HIDDEN; break;
      default: // disable error message
        ; // this is needed to compile
      }
      eventCB->setHandled();
      //std::cout << "H pushed.\n";
      th->set_facets_draw_style(ds);
    }
  }

  Color_id color(const Facet &f) const {
    if (facets_draw_style()== HIDDEN) return HIDE;
    bool hinf=false;
    for (int i=0; i<4; ++i) {
      if (!f.first->vertex(i)->point()){
	hinf=true;
	break;
      }
    }
    if (!triangulation().mirror_vertex(f.first, f.second)->point()) hinf=true;

    if (hinf && convex_hull_== HIDDEN) return HIDE;

    if (internal::has_degree_3(triangulation(), f)){
      return NO_CERT;
    } else if (internal::facet_label(f)){
      return CERT;
    } else {
      for (unsigned int i=0; i< 3; ++i){
	CGAL_assertion(!internal::edge_label(internal::facet_edge(f, i)));
      }
      return UNFAILING_CERT;
    }
  }

  Color_id color(const Edge &f) const {
    if (!internal::has_degree_3(triangulation(), f)) return HIDE;
    if (!internal::edge_label(f)){
      return CERT;
    } else {
      return UNFAILING_CERT;
    }
  }
  
  void generate_geometry();
};


template <class K, class G, class M>
void Qt_triangulation_3<K, G, M>::set_scene_graph_parent(SoSeparator* sep) {
  parent_=sep;
  
  CGAL::Coin_pointer<SoEventCallback> myevcb= new SoEventCallback;
  myevcb->addEventCallback(SoKeyboardEvent::getClassTypeId(),keyboard_callback, this);
  parent_->addChild(myevcb.get());
  
  {
    facets_kit_ = new SoShapeKit;
    facets_kit_->setName("Delaunay_facets_kit");
    {
      CGAL::Coin_pointer<SoAppearanceKit> ap= new SoAppearanceKit;
      {
	CGAL::Coin_pointer<SoMaterial> mat= new SoMaterial;
	mat->setName("Facet_material");
	mat->ambientColor.setNum(4);
	mat->diffuseColor.setNum(4);
	mat->specularColor.setNum(4);
	mat->emissiveColor.setNum(4);
	mat->shininess.setNum(4);
	mat->transparency.setNum(4);
	for (unsigned int i=0; i< 4; ++i){
	  mat->ambientColor.set1Value(i, SbColor(0.2,0.2,0.2));
	  mat->specularColor.set1Value(i, SbColor(0.0, 0.0, 0.0));
	  mat->emissiveColor.set1Value(i, SbColor(0.0,0.0,0.0));
	  mat->shininess.set1Value(i, .2);
	  mat->transparency.set1Value(i, 0.0);
	}
	mat->diffuseColor.set1Value(0, SbColor(0.3, 0.3, 0.3));
	mat->diffuseColor.set1Value(1, SbColor(0.1, 0.4, 0.1));
	mat->diffuseColor.set1Value(2, SbColor(0.1, 0.9, 0.1));
	mat->diffuseColor.set1Value(3, SbColor(1.0, 0.0, 0.0));
	ap->setPart("material", mat.get());
      }
      facets_kit_->setPart("appearance", ap.get());
    }
    facets_kit_->setPart("coordinate3", coordinates_.get());
    {
      CGAL::Coin_pointer< SoShapeHints> hint= new SoShapeHints;
      hint->vertexOrdering.setValue(SoShapeHints::CLOCKWISE);
      hint->shapeType.setValue(SoShapeHints::UNKNOWN_SHAPE_TYPE);
      hint->faceType.setValue(SoShapeHints::CONVEX);
      hint->creaseAngle.setValue(0);
      facets_kit_->setPart("shapeHints", hint.get());
    }
    {
      CGAL::Coin_pointer<SoMaterialBinding> bind= new SoMaterialBinding;
      bind->value.setValue(SoMaterialBinding::PER_FACE_INDEXED);
      facets_kit_->setPart("materialBinding", bind.get());
    }
    {
      facets_= new SoIndexedFaceSet;
      facets_->coordIndex.setNum(0);
      facets_kit_->setPart("shape",facets_.get());
    }
    parent_->addChild(facets_kit_.get());
  }
  
  {
    edges_kit_= new SoShapeKit;
    edges_kit_->setName("Delaunay_edges_kit");
    {
      CGAL::Coin_pointer<SoAppearanceKit> ap= new SoAppearanceKit;
      {
	CGAL::Coin_pointer<SoMaterial> mat= new SoMaterial;
	mat->setName("Edges_material");
	mat->ambientColor.setNum(4);
	mat->diffuseColor.setNum(4);
	mat->specularColor.setNum(4);
	mat->emissiveColor.setNum(4);
	mat->shininess.setNum(4);
	mat->transparency.setNum(4);
	for (unsigned int i=0; i< 4; ++i){
	  mat->specularColor.set1Value(i, SbColor(0.0, 0.0, 0.0));
	  mat->emissiveColor.set1Value(i, SbColor(1.0,0.0,0.0));
	  mat->ambientColor.set1Value(i, SbColor(0.0,1.0,0.0));
	  mat->diffuseColor.set1Value(i, SbColor(0.0,0.0,1.0));
	  mat->shininess.set1Value(i, .2);
	  mat->transparency.set1Value(i, 0.0);
	}
	mat->diffuseColor.set1Value(0, SbColor(0.3, 0.3, 0.3));
	mat->diffuseColor.set1Value(1, SbColor(0.1, 0.4, 0.1));
	mat->diffuseColor.set1Value(2, SbColor(0.1, 0.9, 0.1));
	mat->diffuseColor.set1Value(3, SbColor(1.0, 0.0, 0.0));
	   
	  
	ap->setPart("material", mat.get());
      }
      {
	CGAL::Coin_pointer<SoDrawStyle> ds= new SoDrawStyle;
	ds->lineWidth.setValue(2);
	ap->setPart("drawStyle", ds.get());
      }
      edges_kit_->setPart("appearance", ap.get());
    }
    edges_kit_->setPart("coordinate3", coordinates_.get());
    edges_kit_->setPart("shapeHints", facets_kit_->getPart("shapeHints", true));
    edges_kit_->setPart("materialBinding", facets_kit_->getPart("materialBinding", true));
    {
      edges_= new SoIndexedLineSet;
      edges_->coordIndex.setNum(0);
      edges_kit_->setPart("shape", edges_.get());
    }
    parent_->addChild(edges_kit_.get());
  }
}

template <class K, class G, class M>
void Qt_triangulation_3<K,G,M>::generate_geometry(){
  if (parent_==NULL) return;
  if (triangulation().dimension() != 3) return;
  //print();
  unsigned int facet_count=0;
  unsigned int edge_count=0;
  for (Finite_facets_iterator ffi= triangulation().finite_facets_begin(); 
       ffi != triangulation().finite_facets_end(); ++ffi){
    if (color(*ffi) != HIDE) ++facet_count;
  }
  for (Finite_edges_iterator fei= triangulation().finite_edges_begin();
       fei != triangulation().finite_edges_end(); ++fei){
    if (color(*fei) != HIDE) ++edge_count;
  }
  facets_->coordIndex.setNum(4*facet_count);
  facets_->materialIndex.setNum(facet_count);
  edges_->coordIndex.setNum(3*edge_count);
  edges_->materialIndex.setNum(edge_count);
  {
    int *coords= facets_->coordIndex.startEditing();
    int *mat = facets_->materialIndex.startEditing();
    unsigned int index=0, matindex=0;
    for (Finite_facets_iterator ffi= triangulation().finite_facets_begin(); 
	 ffi != triangulation().finite_facets_end(); ++ffi){
      Color_id id= color(*ffi);
      if (id == HIDE){
	  
      } else {
	coords[index]= facet_vertex(*ffi, 0)->point().index();
	//coords[index]= ffi->vertex(0)->point().index();
	++index;;
	coords[index]= facet_vertex(*ffi, 1)->point().index();
	++index;
	coords[index]= facet_vertex(*ffi, 2)->point().index();
	++index;
	coords[index]= -1;
	++index;

	mat[matindex]= id;
	++matindex;
      }
    }
    facets_->coordIndex.finishEditing();
    facets_->materialIndex.finishEditing();
  }

  {
    int *coords= edges_->coordIndex.startEditing();
    int *mat = edges_->materialIndex.startEditing();
    unsigned int index=0, matindex=0;
    for (Finite_edges_iterator fei= triangulation().finite_edges_begin();
	 fei != triangulation().finite_edges_end(); ++fei){
      Color_id id= color(*fei);
      if (id == HIDE){
	  
      } else {
	coords[index]= edge_vertex(*fei, 0)->point().index();
	++index;
	coords[index]= edge_vertex(*fei, 1)->point().index();
	++index;
	coords[index]= -1;
	++index;

	mat[matindex]= id;
	++matindex;
      }
    }
    edges_->coordIndex.finishEditing();
    edges_->materialIndex.finishEditing();
  }

  SoOutput out;
  out.openFile("output.iv");
  SoWriteAction wa(&out);
  wa.apply(parent_.get());
};

CGAL_KDS_END_NAMESPACE

#endif
