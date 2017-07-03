#include "config_mesh_3.h"
#include "Mesh_3_plugin_cgal_code.h"

#include <CGAL/Real_timer.h>
#include <C3t3_type.h>
#include <CGAL/Mesh_3/polylines_to_protect.h>
#include <CGAL/Bbox_3.h>
#include "CGAL/boost/graph/Graph_with_descriptor_with_graph.h"
#include <Polyhedron_type.h>
#include <CGAL/Mesh_3/properties_Surface_mesh.h>

#include <Scene_c3t3_item.h>

#include <vector>

#include "Mesh_function.h"
#include "Facet_extra_criterion.h"


using namespace CGAL::Three;

typedef Tr::Bare_point Bare_point;

template<typename Mesh>
struct Polyhedral_mesh_domain_selector
{

  typedef Polyhedral_mesh_domain type;
  Polyhedral_mesh_domain_selector(const Mesh&)
  {}
};

template<>
struct Polyhedral_mesh_domain_selector<SMwgd>
{

  typedef Polyhedral_mesh_domain_sm type;
  Polyhedral_mesh_domain_selector(const SMwgd&)
  {}
};
template<class Mesh>
Meshing_thread* cgal_code_mesh_3_templated(const Mesh* pMesh,
                                 const Polylines_container& polylines,
                                 const Mesh* pBoundingMesh,
                                 QString filename,
                                 const double facet_angle,
                                 const double facet_sizing,
                                 const double facet_approx,
                                 const double tet_sizing,
                                 const double edge_size,
                                 const double tet_shape,
                                 bool protect_features,
                                 bool protect_borders,
                                 const double sharp_edges_angle,
                                 const int manifold,
                                 const bool surface_only,
                                 CGAL::Three::Scene_interface* scene)
{
  if(!pMesh) return 0;

  std::cerr << "Meshing file \"" << qPrintable(filename) << "\"\n";
  std::cerr << "  angle: " << facet_angle << std::endl
            << "  edge size bound: " << edge_size << std::endl
            << "  facets size bound: " << facet_sizing << std::endl
            << "  approximation bound: " << facet_approx << std::endl;
  if (is_closed(*pMesh))
    std::cerr << "  tetrahedra size bound: " << tet_sizing << std::endl;

  std::cerr << "Build AABB tree...";
  CGAL::Real_timer timer;
  timer.start();

  typedef typename Polyhedral_mesh_domain_selector<Mesh>::type Polyhedral_mesh_domain;
  // Create domain
  Polyhedral_mesh_domain* p_domain = NULL;
  if (!surface_only && is_closed(*pMesh))
    p_domain = new Polyhedral_mesh_domain(*pMesh);
  else if (!surface_only && pBoundingMesh != NULL && is_closed(*pBoundingMesh))
    p_domain = new Polyhedral_mesh_domain(*pMesh, *pBoundingMesh);
  else
  {
    std::vector<const Mesh*> poly_ptrs_vector(1, pMesh);
    p_domain = new Polyhedral_mesh_domain(poly_ptrs_vector.begin(), poly_ptrs_vector.end());
  }

  // Features
  if(polylines.empty() && protect_features) {
      //includes detection of borders in the surface case
      p_domain->detect_features(sharp_edges_angle);
  }
  else if (polylines.empty() && protect_borders)
  {
    p_domain->detect_borders();
  }
  if(! polylines.empty()){
    p_domain->add_features(polylines.begin(), polylines.end());
    protect_features = true; // so that it will be passed in make_mesh_3
  }

  std::cerr << " done (" << timer.time() * 1000 << " ms)" << std::endl;

  Scene_c3t3_item* p_new_item = new Scene_c3t3_item;
  p_new_item->setScene(scene);

  QString tooltip = QString("<div>From \"") + filename +
    QString("\" with the following mesh parameters"
            "<ul>"
            "<li>Angle: %1</li>"
            "<li>Edge size bound: %2</li>"
            "<li>Facets size bound: %3</li>"
            "<li>Approximation bound: %4</li>")
    .arg(facet_angle)
    .arg(edge_size)
    .arg(facet_sizing)
    .arg(facet_approx);
  if (is_closed(*pMesh))
    tooltip += QString("<li>Tetrahedra size bound: %1</li>" )
        .arg(tet_sizing);
  tooltip += "</ul></div>";

  p_new_item->setProperty("toolTip",tooltip);
  Mesh_parameters param;
  param.facet_angle = facet_angle;
  param.facet_sizing = facet_sizing;
  param.facet_approx = facet_approx;
  param.tet_sizing = tet_sizing;
  param.tet_shape = tet_shape;
  param.edge_sizing =  edge_size;
  param.manifold = manifold;
  param.protect_features = protect_features || protect_borders;
  param.use_sizing_field_with_aabb_tree = polylines.empty() && protect_features;

  typedef ::Mesh_function<Polyhedral_mesh_domain,
                          Mesh_fnt::Polyhedral_domain_tag> Mesh_function;
  Mesh_function* p_mesh_function = new Mesh_function(p_new_item->c3t3(),
                                                     p_domain, param);
  return new Meshing_thread(p_mesh_function, p_new_item);
}


Meshing_thread* cgal_code_mesh_3(const Polyhedron* pMesh,
                                 const Polylines_container& polylines,
                                 const Polyhedron* pBoundingMesh,
                                 QString filename,
                                 const double facet_angle,
                                 const double facet_sizing,
                                 const double facet_approx,
                                 const double tet_sizing,
                                 const double edge_size,
                                 const double tet_shape,
                                 bool protect_features,
                                 bool protect_borders,
                                 const double sharp_edges_angle,
                                 const int manifold,
                                 const bool surface_only,
                                 CGAL::Three::Scene_interface* scene)
{
  return cgal_code_mesh_3_templated(pMesh,
                          polylines,
                          pBoundingMesh,
                          filename,
                          facet_angle,
                          facet_sizing,
                          facet_approx,
                          tet_sizing,
                          edge_size,
                          tet_shape,
                          protect_features,
                          protect_borders,
                          sharp_edges_angle,
                          manifold,
                          surface_only,
                          scene);
}

Meshing_thread* cgal_code_mesh_3(const SMwgd* pMesh,
                                 const Polylines_container& polylines,
                                 const SMwgd* pBoundingMesh,
                                 QString filename,
                                 const double facet_angle,
                                 const double facet_sizing,
                                 const double facet_approx,
                                 const double tet_sizing,
                                 const double edge_size,
                                 const double tet_shape,
                                 bool protect_features,
                                 bool protect_borders,
                                 const double sharp_edges_angle,
                                 const int manifold,
                                 const bool surface_only,
                                 CGAL::Three::Scene_interface* scene)
{
  return cgal_code_mesh_3_templated(pMesh,
                          polylines,
                          pBoundingMesh,
                          filename,
                          facet_angle,
                          facet_sizing,
                          facet_approx,
                          tet_sizing,
                          edge_size,
                          tet_shape,
                          protect_features,
                          protect_borders,
                          sharp_edges_angle,
                          manifold,
                          surface_only,
                          scene);
}
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS

Meshing_thread* cgal_code_mesh_3(const Implicit_function_interface* pfunction,
                                 const double facet_angle,
                                 const double facet_sizing,
                                 const double facet_approx,
                                 const double tet_sizing,
                                 const double edge_size,
                                 const double tet_shape,
                                 const int manifold,
                                 CGAL::Three::Scene_interface* scene)
{
  if (pfunction == NULL) { return NULL; }

  CGAL::Bbox_3 domain_bbox(pfunction->bbox().xmin(),
                           pfunction->bbox().ymin(),
                           pfunction->bbox().zmin(),
                           pfunction->bbox().xmax(),
                           pfunction->bbox().ymax(),
                           pfunction->bbox().zmax());

  Function_mesh_domain* p_domain =
    new Function_mesh_domain(Function_wrapper(*pfunction), domain_bbox, 1e-7);

  Scene_c3t3_item* p_new_item = new Scene_c3t3_item;
  p_new_item->setScene(scene);

  Mesh_parameters param;
  param.protect_features = false;
  param.facet_angle = facet_angle;
  param.facet_sizing = facet_sizing;
  param.facet_approx = facet_approx;
  param.tet_sizing = tet_sizing;
  param.tet_shape = tet_shape;
  param.edge_sizing = edge_size;
  param.manifold = manifold;

  typedef ::Mesh_function<Function_mesh_domain,
                          Mesh_fnt::Implicit_domain_tag> Mesh_function;
  Mesh_function* p_mesh_function = new Mesh_function(p_new_item->c3t3(),
                                                     p_domain, param);
  return new Meshing_thread(p_mesh_function, p_new_item);
}
#endif // CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS


#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
#include <CGAL/Gray_image_mesh_domain_3.h>

typedef CGAL::Gray_image_mesh_domain_3<CGAL::Image_3, Kernel, float, CGAL::Compare_to_isovalue > Gray_image_domain1;
typedef CGAL::Polyhedron_demo_labeled_mesh_domain_3<Gray_image_domain1> Gray_image_domain;
typedef CGAL::Mesh_domain_with_polyline_features_3<Gray_image_domain> Gray_Image_mesh_domain;

Meshing_thread* cgal_code_mesh_3(const Image* pImage,
                                 const Polylines_container& polylines,
                                 const double facet_angle,
                                 const double facet_sizing,
                                 const double facet_approx,
                                 const double tet_sizing,
                                 const double edge_size,
                                 const double tet_shape,
                                 bool protect_features,
                                 const int manifold,
                                 CGAL::Three::Scene_interface* scene,
                                 bool detect_connected_components,
                                 bool is_gray,
                                 float iso_value,
                                 float value_outside,
                                 bool inside_is_less)
{
  if (NULL == pImage) { return NULL; }

  Image_mesh_domain* p_domain = new Image_mesh_domain(*pImage, 1e-6);

  if(protect_features && polylines.empty())
  {
    std::vector<std::vector<Bare_point> > polylines_on_bbox;
    if(pImage->image()->wordKind == WK_FLOAT)
      CGAL::polylines_to_protect<Bare_point, float>(*pImage, polylines_on_bbox);
    else //WK_FIXED
      CGAL::polylines_to_protect<Bare_point, unsigned char>(*pImage, polylines_on_bbox);

    p_domain->add_features(polylines_on_bbox.begin(), polylines_on_bbox.end());
  }
  if(! polylines.empty()){
    // Insert edge in domain
    p_domain->add_features(polylines.begin(), polylines.end());
    protect_features = true; // so that it will be passed in make_mesh_3
  }
  Mesh_parameters param;
  param.protect_features = protect_features;
  param.detect_connected_components = detect_connected_components;
  param.facet_angle = facet_angle;
  param.facet_sizing = facet_sizing;
  param.facet_approx = facet_approx;
  param.tet_sizing = tet_sizing;
  param.edge_sizing = edge_size;
  param.tet_shape = tet_shape;
  param.manifold = manifold;
  param.image_3_ptr = pImage;
  Scene_c3t3_item* p_new_item = new Scene_c3t3_item;
  p_new_item->setScene(scene);
  if(!is_gray)
  {
    Image_mesh_domain* p_domain = new Image_mesh_domain(*pImage, 1e-6);

    if(protect_features && polylines.empty()){
      std::vector<std::vector<Bare_point> > polylines_on_bbox;
      CGAL::polylines_to_protect<
        Bare_point,
        Image_mesh_domain::Image_word_type>(*pImage, polylines_on_bbox);
      p_domain->add_features(polylines_on_bbox.begin(), polylines_on_bbox.end());
    }
    if(! polylines.empty()){
      // Insert edge in domain
      p_domain->add_features(polylines.begin(), polylines.end());
    }
    typedef ::Mesh_function<Image_mesh_domain,
                            Mesh_fnt::Labeled_image_domain_tag> Mesh_function;
    Mesh_function* p_mesh_function = new Mesh_function(p_new_item->c3t3(),
                                                       p_domain, param);
    return new Meshing_thread(p_mesh_function, p_new_item);
  }
  else
  {
    if(CGAL::Compare_to_isovalue(iso_value, inside_is_less)(value_outside) !=0)
    {
      std::cerr << "Warning : \"Value inside is less than iso value\"'s value has been inverted to avoid crash.  "
                << " " << std::endl;
      inside_is_less = !inside_is_less;
    }
    Gray_Image_mesh_domain* p_domain = new Gray_Image_mesh_domain(*pImage, CGAL::Compare_to_isovalue(iso_value, inside_is_less), value_outside);
    if(protect_features && polylines.empty())
    {
      std::cerr << "Warning : Automatic detection of features"
                << " in Gray images is not implemented yet" << std::endl;
      //std::vector<std::vector<Bare_point> > polylines_on_bbox;
      //CGAL::polylines_to_protect<Bare_point>(*pImage, polylines_on_bbox);
      //p_domain->add_features(polylines_on_bbox.begin(), polylines_on_bbox.end());
    }
    if(! polylines.empty()){
      // Insert edge in domain
      p_domain->add_features(polylines.begin(), polylines.end());
    }
    typedef ::Mesh_function<Gray_Image_mesh_domain,
                            Mesh_fnt::Gray_image_domain_tag> Mesh_function;
    Mesh_function* p_mesh_function = new Mesh_function(p_new_item->c3t3(),
                                                       p_domain, param);
    return new Meshing_thread(p_mesh_function, p_new_item);
  }
}

#endif //CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES


//#include "Polyhedron_demo_mesh_3_plugin_cgal_code.moc"
//#include "Scene_c3t3_item.moc" //Check this one, it's strange moc include.

