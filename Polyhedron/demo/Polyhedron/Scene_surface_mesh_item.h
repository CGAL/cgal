#ifndef CGAL_SCENE_SURFACE_MESH_ITEM_H
#define CGAL_SCENE_SURFACE_MESH_ITEM_H
//Defines the precision of the positions (for performance/precision sake)
#define CGAL_GL_DATA GL_FLOAT
#define cgal_gl_data float
#define CGAL_IS_FLOAT 1

#include "Scene_surface_mesh_item_config.h"
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Viewer_interface.h>
#include <vector>

#include <boost/scoped_ptr.hpp>
#include <boost/array.hpp>

#include <CGAL/Surface_mesh/Surface_mesh_fwd.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>




class SCENE_SURFACE_MESH_ITEM_EXPORT Scene_surface_mesh_item
  : public CGAL::Three::Scene_item
{
  Q_OBJECT
public:
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef Kernel::Point_3 Point;
  typedef CGAL::Surface_mesh<Point> SMesh;
  typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;

  // Takes ownership of the argument.
  Scene_surface_mesh_item(SMesh*);

  Scene_surface_mesh_item(const Scene_surface_mesh_item& other);

  ~Scene_surface_mesh_item(){delete smesh_;}

  Scene_surface_mesh_item* clone() const;
  void draw(CGAL::Three::Viewer_interface *) const;
  void drawEdges(CGAL::Three::Viewer_interface *) const;
  void drawPoints(CGAL::Three::Viewer_interface *) const;

  bool supportsRenderingMode(RenderingMode m) const;
  bool isFinite() const { return true; }
  bool isEmpty() const;
  Bbox bbox() const;

  enum VAOs {
   Flat_facets = 0,
   Smooth_facets,
   Edges,
   NbOfVaos
  };
  enum VBOs {
    Flat_vertices = 0,
    Smooth_vertices,
    Flat_normals,
    Smooth_normals,
    VColors,
    FColors,
    NbOfVbos
  };

  QString toolTip() const;

  SMesh* polyhedron() { return smesh_; }
  const SMesh* polyhedron() const { return smesh_; }
public Q_SLOTS:
  virtual void selection_changed(bool);
private:
  mutable bool floated;
  mutable bool has_vcolors;
  mutable bool has_fcolors;
  SMesh* smesh_;
  void initializeBuffers(CGAL::Three::Viewer_interface *) const;
  void addFlatData(Point, Kernel::Vector_3, CGAL::Color *) const;
  mutable std::vector<unsigned int> idx_data_;
  std::vector<unsigned int> idx_edge_data_;
  mutable std::vector<cgal_gl_data> smooth_vertices;
  mutable std::vector<cgal_gl_data> smooth_normals;
  mutable std::vector<cgal_gl_data> flat_vertices;
  mutable std::vector<cgal_gl_data> flat_normals;
  mutable std::vector<cgal_gl_data> f_colors;
  mutable std::vector<cgal_gl_data> v_colors;
  mutable QOpenGLShaderProgram *program;

  //! \param fd a face_descriptor of the facet that needs to be triangulated.
  //! \param fnormals a property_map containing the normals of the mesh.
  //! \param p_cdt a reference to an empty CDT that will be filled by this function.
  //! \param v2v a reference to an empty flat_map that will be filled by this function.
  //!

  //!
  //! \brief triangulate_facet Triangulates a facet.
  //! \param fd a face_descriptor of the facet that needs to be triangulated.
  //! \param fnormals a property_map containing the normals of the mesh.
  //! \param fcolors a property_map containing the colors of the mesh
  //! \param im a property_map containing the indices of the vertices of the mesh
  //! \param index if true, the function will fill the index vector. If false, the function will
  //! fill the flat data vectors.
  void
  triangulate_facet(face_descriptor fd,
                    SMesh::Property_map<face_descriptor, Kernel::Vector_3 > *fnormals,
                    SMesh::Property_map<face_descriptor, CGAL::Color> *fcolors,
                    boost::property_map< SMesh, boost::vertex_index_t >::type* im,
                    bool index) const;
  mutable bool is_filled;
  void checkFloat() const;
  void compute_elements();
  mutable bool isinit;
};


#endif /* CGAL_SCENE_SURFACE_MESH_ITEM_H */
