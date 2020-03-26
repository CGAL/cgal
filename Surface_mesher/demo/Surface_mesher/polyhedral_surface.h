#ifndef _POLYHEDRAL_SURFACE_H
#define _POLYHEDRAL_SURFACE_H

#define CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION

#include "surface.h"
#include "viewer.h"

#include <map>
#include <utility>

typedef std::pair<const char*, const char*> Signal_slot_pair;
typedef std::map<const char*, Signal_slot_pair> Connection_map;

// CGAL
#include <CGAL/basic.h>

// kernel
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// surface
#include <CGAL/Polyhedral_surface_3.h>

// piece-wise loop subdivision
#include <CGAL/pws_loop_subdivision.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Poly_kernel;
typedef Poly_kernel::FT FT;
typedef Poly_kernel::Point_3 Poly_point;
typedef Poly_kernel::Sphere_3 Sphere;
typedef Poly_kernel::Vector_3 Vector;
typedef Poly_kernel::Triangle_3 Triangle_3;

typedef CGAL::Polyhedral_surface_3<Poly_kernel,
              CGAL::Surface_mesher::Has_edges> CGAL_polyhedral_surface;


class Polyhedral_surface : public Surface
{
  Q_OBJECT;

public:
  Polyhedral_surface(QObject* parent,
                     double sharp_edges_angle_lower_bound = 90,
                     double sharp_edges_angle_upper_bound = 120);

  ~Polyhedral_surface();
public Q_SLOTS:
  void clear();
  void connect_actions();
  void display_nb_elements_in_status_bar() const;
  void set_dirty();
  void busy() const;
  void not_busy() const;
  void set_sharp_edges_angle_bounds(const double lower_bound,
                                    const double upper_bound);
  void update_data_structures();
  void construct_octree();
  void toggle_display_octree(bool b);
  void toggle_display_edges_octree(bool b);
  void toggle_display_surface(bool b);
  void toggle_display_all_edges(bool b);
  void toggle_display_control_edges(bool b);
  void make_one_subdivision_step();
  void on_action_Options_triggered();

public:
  bool open(const QString& filename);
  void close();
  void draw();
  void drawWithNames();
  void postSelection(const QPoint&);
  void draw(bool with_names);
  void get_bbox(float& xmin, float& ymin, float& zmin,
                float& xmax, float& ymax, float& zmax);

public Q_SLOTS:
  void set_inverse_normals(const bool b);

public:
  bool inverse_normals() const;

protected:
  bool m_inverse_normals;
  CGAL_polyhedral_surface* surface_ptr;
  QObject* parent;
  bool display_octree;
  bool display_edges_octree;
  bool display_surface;
  bool display_all_edges;
  bool display_control_edges;
  double sharp_edges_angle_lower_bound;
  double sharp_edges_angle_upper_bound;
  bool is_octree_initialized;
  int selected_edge;
  int selected_facet;
  bool is_dirty;
  GLint list_id;
  Connection_map connection_map;
};

#endif // _POLYHEDRAL_SURFACE_H
