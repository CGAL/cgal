#include "config_mesh_3.h"
#include "Mesh_3_plugin_cgal_code.h"

#include <CGAL/Real_timer.h>
#include <C3t3_type.h>
#include <CGAL/Mesh_3/polylines_to_protect.h>
#include <CGAL/Bbox_3.h>

#include <Scene_c3t3_item.h>

#include <vector>

#include "Mesh_function.h"
#include "Facet_extra_criterion.h"


using namespace CGAL::Three;

typedef Tr::Bare_point Bare_point;

struct Compare_to_isovalue {
  double iso_value;
  bool less;
  typedef bool result_type;

  Compare_to_isovalue(double iso_value, bool less)
    : iso_value(iso_value), less(less) {}

  bool operator()(double x) const {
    return (x < iso_value) == less;
  }
};

Meshing_thread* cgal_code_mesh_3(QList<const SMesh*> pMeshes,
                                 const Polylines_container& polylines,
                                 const SMesh* pBoundingMesh,
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
                                 const bool surface_only)
{
  if(pMeshes.empty() && nullptr == pBoundingMesh) return 0;

  std::cerr << "Meshing file \"" << qPrintable(filename) << "\"\n";
  std::cerr << "  angle: " << facet_angle << std::endl
            << "  edge size bound: " << edge_size << std::endl
            << "  facets size bound: " << facet_sizing << std::endl
            << "  approximation bound: " << facet_approx << std::endl;
  if (!surface_only)
    std::cerr << "  tetrahedra size bound: " << tet_sizing << std::endl;

  std::cerr << "Build AABB tree...";
  CGAL::Real_timer timer;
  timer.start();

  // Create domain
  Polyhedral_mesh_domain* p_domain = nullptr;
  if (surface_only || nullptr == pBoundingMesh)
    p_domain = new Polyhedral_mesh_domain(pMeshes.begin(), pMeshes.end());
  else if(pMeshes.empty())
    p_domain = new Polyhedral_mesh_domain(*pBoundingMesh);
  else
    p_domain = new Polyhedral_mesh_domain(pMeshes.begin(), pMeshes.end(),
                                          *pBoundingMesh);
  // Features
  if(polylines.empty()) {
    if(protect_features) {
      //includes detection of borders in the surface case
      p_domain->detect_features(sharp_edges_angle);
    }
    else if (protect_borders) {
      p_domain->detect_borders();
    }
  } else {
    p_domain->add_features(polylines.begin(), polylines.end());
    protect_features = true; // so that it will be passed in make_mesh_3
  }

  std::cerr << " done (" << timer.time() * 1000 << " ms)" << std::endl;

  Scene_c3t3_item* p_new_item = new Scene_c3t3_item(surface_only);
  if(polylines.empty()) {
    if(protect_features) {
      p_new_item->set_sharp_edges_angle(sharp_edges_angle);
    }
    else if (protect_borders) {
      p_new_item->set_detect_borders(true);
    }
  }

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
  if (!surface_only)
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

#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS

Meshing_thread* cgal_code_mesh_3(const Implicit_function_interface* pfunction,
                                 const double facet_angle,
                                 const double facet_sizing,
                                 const double facet_approx,
                                 const double tet_sizing,
                                 const double edge_size,
                                 const double tet_shape,
                                 const int manifold,
                                 const bool surface_only)
{
  if (pfunction == nullptr) { return nullptr; }

  CGAL::Bbox_3 domain_bbox(pfunction->bbox().xmin(),
                           pfunction->bbox().ymin(),
                           pfunction->bbox().zmin(),
                           pfunction->bbox().xmax(),
                           pfunction->bbox().ymax(),
                           pfunction->bbox().zmax());
  namespace p = CGAL::parameters;
  Function_mesh_domain* p_domain =
    new Function_mesh_domain(Function_wrapper(*pfunction), domain_bbox, 1e-7,
                             p::construct_surface_patch_index =
                               [](int i, int j) { return (i * 1000 + j); }
                             );

  Scene_c3t3_item* p_new_item = new Scene_c3t3_item(surface_only);

  Mesh_parameters param;
  param.protect_features = false;
  param.facet_angle = facet_angle;
  param.facet_sizing = facet_sizing;
  param.facet_approx = facet_approx;
  param.tet_sizing = tet_sizing;
  param.tet_shape = tet_shape;
  param.edge_sizing = edge_size;
  param.manifold = manifold;
  param.detect_connected_components = false; // to avoid random values
                                             // in the debug displays


  typedef ::Mesh_function<Function_mesh_domain,
                          Mesh_fnt::Implicit_domain_tag> Mesh_function;
  Mesh_function* p_mesh_function = new Mesh_function(p_new_item->c3t3(),
                                                     p_domain, param);
  return new Meshing_thread(p_mesh_function, p_new_item);
}
#endif // CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS


#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
#include <CGAL/Labeled_mesh_domain_3.h>

typedef CGAL::Labeled_mesh_domain_3<Kernel, int, int> Gray_image_domain;
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
                                 const bool surface_only,
                                 bool detect_connected_components,
                                 bool is_gray,
                                 float iso_value,
                                 float value_outside,
                                 bool inside_is_less)
{
  if (nullptr == pImage) { return nullptr; }

  if(! polylines.empty()){
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
  Scene_c3t3_item* p_new_item = new Scene_c3t3_item(surface_only);
  if(!is_gray)
  {
    namespace p = CGAL::parameters;

    Image_mesh_domain* p_domain = new Image_mesh_domain
      (Image_mesh_domain::create_labeled_image_mesh_domain
       (p::image = *pImage,
        p::relative_error_bound = 1e-6,
        p::construct_surface_patch_index =
          [](int i, int j) { return (i * 1000 + j); }
        )
       );

    if(protect_features && polylines.empty()){
      std::vector<std::vector<Bare_point> > polylines_on_bbox;

      CGAL_IMAGE_IO_CASE(pImage->image(),
                         {
                           typedef Word Image_word_type;
                           (CGAL::polylines_to_protect<
                              Bare_point,
                              Image_word_type>(*pImage, polylines_on_bbox));
                           p_domain->add_features(polylines_on_bbox.begin(),
                                                  polylines_on_bbox.end());
                         }
                         );
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
    if(Compare_to_isovalue(iso_value, inside_is_less)(value_outside) !=0)
    {
      std::cerr << "Warning : \"Value inside is less than iso value\"'s value has been inverted to avoid crash.  "
                << " " << std::endl;
      inside_is_less = !inside_is_less;
    }
    Compare_to_isovalue domain_functor(iso_value, inside_is_less);
    namespace p = CGAL::parameters;
    Gray_Image_mesh_domain* p_domain = new Gray_Image_mesh_domain
      (Gray_Image_mesh_domain::create_gray_image_mesh_domain
       (p::image = *pImage,
        p::relative_error_bound = 1e-6,
        p::image_values_to_subdomain_indices = domain_functor,
        p::value_outside = value_outside,
        p::construct_surface_patch_index =
          [](int i, int j) { return (i * 1000 + j); }
        )
       );
    if(protect_features && polylines.empty())
    {
      using CGAL::Mesh_3::internal::Linear_interpolator;
      std::vector<std::vector<Bare_point> > polylines_on_bbox;

      CGAL_IMAGE_IO_CASE(pImage->image(),

      typedef Word Image_word_type;
      (CGAL::polylines_to_protect<Bare_point, Image_word_type>
        ( *pImage,
          pImage->vx(),
          pImage->vy(),
          pImage->vz(),
          polylines_on_bbox,
          (Image_word_type*)0,
          CGAL::Null_subdomain_index(),
          domain_functor,
          Linear_interpolator<Kernel, Image_word_type>(),
          0,
          0,
          Image_word_type(iso_value),
          (std::max)(10, int(10 * (std::min)((std::min)(pImage->vx(),
                                                        pImage->vy()),
                                             pImage->vz())
                             / edge_size))));

       );

      p_domain->add_features(polylines_on_bbox.begin(), polylines_on_bbox.end());
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

