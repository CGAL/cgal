/**
 * @file   ui/gl/MainOpenGLWindow.h
 * @author Gernot Walzl
 * @date   2011-12-20
 */

#ifndef UI_GL_MAINOPENGLWINDOW_H
#define UI_GL_MAINOPENGLWINDOW_H

#include "typedefs_thread.h"
#include "algo/ptrs.h"
#include "algo/Controller.h"
#include "data/2d/ptrs.h"
#include "data/2d/skel/ptrs.h"
#include "data/2d/mesh/ptrs.h"
#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "ui/gl/ptrs.h"
#include "ui/gl/OpenGLWindow.h"
#include "ui/ps/ptrs.h"
#include "util/Configuration.h"
#include <cmath>
#include <list>

namespace ui { namespace gl {

using algo::ControllerSPtr;
using data::_2d::PolygonSPtr;
using data::_3d::PolyhedronSPtr;
using data::_3d::SphericalPolygonSPtr;
using util::Configuration;
using util::ConfigurationSPtr;

class MainOpenGLWindow : public OpenGLWindow {
public:
    virtual ~MainOpenGLWindow();
    static MainOpenGLWindowSPtr create(int argc, const char* argv[],
            unsigned int width, unsigned int height,
            ControllerSPtr controller);

    static const int MODE_2D = 2;
    static const int MODE_3D = 3;
    static const int MODE_SPHERICAL = 4;

    void setTranslate(const vec3f translate);
    void setScale(float scale);

    void setPolygon(PolygonSPtr polygon);
    void setSkel2d(data::_2d::skel::StraightSkeletonSPtr skel_2d);
    void setMesh2d(data::_2d::mesh::MeshSPtr mesh_2d);
    void setPolyhedron(PolyhedronSPtr polyhedron);
    void setSkel3d(data::_3d::skel::StraightSkeletonSPtr skel_3d);
    void setSphericalPolygon(SphericalPolygonSPtr sphericalpolygon);
    void setSphericalSkel(data::_3d::skel::SphericalSkeletonSPtr sphericalskel);

    void toggleRoof();
    void togglePoly();
    void toggleSkel();
    void toggleExperimental();
    void incThickness();
    void decThickness();

    void saveSkel();
    void saveLastPoly();
    void dumpWin();
    void printScreen();
    void printCutPattern();

    void run();
    ThreadSPtr startThread();

protected:
    MainOpenGLWindow(int argc, const char* argv[],
            unsigned int width, unsigned int height,
            ControllerSPtr controller);

    void convert(data::_2d::Vector2SPtr in, vec3f& out);
    void convert(data::_2d::Point2SPtr in, vec3f& out);
    void convert(data::_3d::Vector3SPtr in, vec3f& out);
    void convert(data::_3d::Point3SPtr in, vec3f& out);

    void drawCoordAxes(float size);
    void drawPolygon(PolygonSPtr polygon, float height, bool bold, bool vertices_only);
    void drawSkel2d(data::_2d::skel::StraightSkeletonSPtr skel_2d);
    void drawMesh2d(data::_2d::mesh::MeshSPtr mesh_2d);
    void drawPolyhedron(PolyhedronSPtr polyhedron, bool bold, bool vertices_only);
    void drawSkel3d(data::_3d::skel::StraightSkeletonSPtr skel_3d);
    void drawColoredNodes(data::_3d::skel::StraightSkeletonSPtr skel_3d);
    void drawSphericalPolygon(SphericalPolygonSPtr sphericalpolygon, bool bold, bool vertices_only);
    void drawSphericalSkel(data::_3d::skel::SphericalSkeletonSPtr sphericalskel);

    virtual void drawContent();
    virtual void handleKeyPressed(int key, int x, int y);
    virtual void handleKeyReleased(int key, int x, int y);
    virtual void handleMousePressed(int button, int x, int y);
    virtual void handleMouseReleased(int button, int x, int y);
    virtual void handleMotion(int x, int y);

    ControllerSPtr controller_;
    KeyboardAdapterSPtr keyboard_adapter_;
    MouseAdapterSPtr mouse_adapter_;

    int mode_;
    vec3f translate_;
    float scale_;
    PolygonSPtr polygon_;
    data::_2d::skel::StraightSkeletonSPtr skel_2d_;
    data::_2d::mesh::MeshSPtr mesh_2d_;
    PolyhedronSPtr polyhedron_;
    data::_3d::skel::StraightSkeletonSPtr skel_3d_;
    SphericalPolygonSPtr sphericalpolygon_;
    data::_3d::skel::SphericalSkeletonSPtr sphericalskel_;
    bool toggle_roof_;
    unsigned int toggle_poly_;  // 0 = off, 1 = first, 2 = first, 3 = all, 4 = all, 5 = last
    bool toggle_skel_;
    bool toggle_experimental_;
    float thickness_;
    float crosshair_size_;
    float coord_axes_size_;
    bool highlight_;
    bool draw_dirs_;
};

} }

#endif /* UI_GL_MAINOPENGLWINDOW_H */
