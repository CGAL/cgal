/**
 * @file   ui/gl/MainOpenGLWindow.cpp
 * @author Gernot Walzl
 * @date   2011-12-20
 */

#include "ui/gl/MainOpenGLWindow.h"

#include "algo/3d/KernelWrapper.h"

#include "data/2d/Polygon.h"
#include "data/2d/Vertex.h"
#include "data/2d/Edge.h"
#include "data/2d/VertexData.h"
#include "data/2d/EdgeData.h"
#include "data/2d/skel/AbstractEvent.h"
#include "data/2d/skel/Node.h"
#include "data/2d/skel/Arc.h"
#include "data/2d/skel/StraightSkeleton.h"
#include "data/2d/mesh/Mesh.h"
#include "data/2d/mesh/MeshCell.h"
#include "data/2d/mesh/MeshVertex.h"

#include "data/3d/Polyhedron.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"
#include "data/3d/Triangle.h"
#include "data/3d/VertexData.h"
#include "data/3d/EdgeData.h"
#include "data/3d/SphericalPolygon.h"
#include "data/3d/CircularVertex.h"
#include "data/3d/CircularEdge.h"
#include "data/3d/CircularVertexData.h"
#include "data/3d/CircularEdgeData.h"
#include "data/3d/skel/StraightSkeleton.h"
#include "data/3d/skel/AbstractEvent.h"
#include "data/3d/skel/Arc.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/Sheet.h"
#include "data/3d/skel/SphericalSkeleton.h"
#include "data/3d/skel/SphericalSkelVertexData.h"
#include "data/3d/skel/SkelEdgeData.h"
#include "data/3d/skel/SphericalAbstractEvent.h"
#include "data/3d/skel/CircularNode.h"
#include "data/3d/skel/CircularArc.h"
#include "data/3d/skel/ConstOffsetEvent.h"
#include "data/3d/skel/EdgeEvent.h"
#include "data/3d/skel/EdgeMergeEvent.h"
#include "data/3d/skel/TriangleEvent.h"
#include "data/3d/skel/DblEdgeMergeEvent.h"
#include "data/3d/skel/DblTriangleEvent.h"
#include "data/3d/skel/TetrahedronEvent.h"
#include "data/3d/skel/VertexEvent.h"
#include "data/3d/skel/FlipVertexEvent.h"
#include "data/3d/skel/SurfaceEvent.h"
#include "data/3d/skel/PolyhedronSplitEvent.h"
#include "data/3d/skel/SplitMergeEvent.h"
#include "data/3d/skel/EdgeSplitEvent.h"
#include "data/3d/skel/PierceEvent.h"

#include "db/2d/DAOFactory.h"
#include "db/2d/PolygonDAO.h"
#include "db/2d/StraightSkeletonDAO.h"
#include "db/3d/DAOFactory.h"
#include "db/3d/PolyhedronDAO.h"
#include "db/3d/StraightSkeletonDAO.h"
#include "db/3d/OBJFile.h"

#include "ui/gl/Camera.h"
#include "ui/gl/KeyboardAdapter.h"
#include "ui/gl/MouseAdapter.h"
#include "ui/ps/PlanePSPrinter.h"
#include "ui/ps/SpacePSPrinter.h"
#include "ui/ps/CutPatternPrinter.h"

#include "util/StringFactory.h"

#include <fstream>
#include <sstream>

namespace ui { namespace gl {

MainOpenGLWindow::~MainOpenGLWindow() {
    // intentionally does nothing
}

MainOpenGLWindow::MainOpenGLWindow(int argc, const char* argv[],
        unsigned int width, unsigned int height,
        ControllerSPtr controller) :
        OpenGLWindow(argc, argv, width, height, "Straight Skeleton") {
    // delegating constructor
    this->controller_ = controller;
    this->keyboard_adapter_ = KeyboardAdapter::create(this->camera_, controller);
    this->mouse_adapter_ = MouseAdapter::create(this->camera_);
    this->mode_ = 0;
    this->toggle_roof_ = true;
    this->toggle_poly_ = 3;
    this->toggle_skel_ = true;
    this->toggle_experimental_ = false;
    this->thickness_ = 1.0f;
    this->crosshair_size_ = 0.0f;
    this->coord_axes_size_ = 1.0f;
    this->highlight_ = true;
    this->draw_dirs_ = true;
    for (unsigned int i = 0; i < 3; i++) {
        this->translate_[i] = 0.0f;
    }
    this->scale_ = 1.0f;
    ConfigurationSPtr config = Configuration::getInstance();
    if (config->isLoaded()) {
        std::string section("ui_gl_MainOpenGLWindow");
        crosshair_size_ = (float)config->getDouble(section, "crosshair_size");
        coord_axes_size_ = (float)config->getDouble(section, "coord_axes_size");
        highlight_ = config->getBool(section, "highlight");
        draw_dirs_ = config->getBool(section, "draw_dirs");
    }
}

MainOpenGLWindowSPtr MainOpenGLWindow::create(int argc, const char* argv[],
        unsigned int width, unsigned int height,
        ControllerSPtr controller) {
    MainOpenGLWindowSPtr result = MainOpenGLWindowSPtr(
            new MainOpenGLWindow(argc, argv, width, height, controller));
    result->keyboard_adapter_->setWindow(result);
    return result;
}


void MainOpenGLWindow::convert(data::_2d::Vector2SPtr in, vec3f& out) {
    for (unsigned int i = 0; i < 2; i++) {
        out[i] = (float)(*in)[i];
    }
    out[2] = 0.0f;
}

void MainOpenGLWindow::convert(data::_2d::Point2SPtr in, vec3f& out) {
    for (unsigned int i = 0; i < 2; i++) {
        out[i] = (float)(*in)[i];
    }
    out[2] = 0.0f;
}

void MainOpenGLWindow::convert(data::_3d::Vector3SPtr in, vec3f& out) {
    for (unsigned int i = 0; i < 3; i++) {
        out[i] = (float)(*in)[i];
    }
}

void MainOpenGLWindow::convert(data::_3d::Point3SPtr in, vec3f& out) {
    for (unsigned int i = 0; i < 3; i++) {
        out[i] = (float)(*in)[i];
    }
}


void MainOpenGLWindow::toggleRoof() {
    Lock l(mutex_);
    toggle_roof_ = !toggle_roof_;
}

void MainOpenGLWindow::togglePoly() {
    Lock l(mutex_);
    unsigned int num_toggles = 6;
    toggle_poly_ = (toggle_poly_+1)%num_toggles;
}

void MainOpenGLWindow::toggleSkel() {
    Lock l(mutex_);
    toggle_skel_ = !toggle_skel_;
}

void MainOpenGLWindow::toggleExperimental() {
    Lock l(mutex_);
    toggle_experimental_ = !toggle_experimental_;
}


void MainOpenGLWindow::incThickness() {
    Lock l(mutex_);
//    if (thickness_ < 1.95f) {
//        thickness_ += 0.1f;
//    }
    thickness_ *= 2.0f;
}

void MainOpenGLWindow::decThickness() {
    Lock l(mutex_);
//    if (thickness_ > 0.15f) {
//        thickness_ -= 0.1f;
//    }
    thickness_ /= 2.0f;
}


void MainOpenGLWindow::saveSkel() {
    if (mode_ == MODE_2D) {
        if (polygon_) {
            db::_2d::PolygonDAOSPtr polygon_dao =
                    db::_2d::DAOFactory::getPolygonDAO();
            int polygon_id = polygon_dao->insert(skel_2d_->getPolygon());
            if (polygon_id > 0 && skel_2d_) {
                db::_2d::StraightSkeletonDAOSPtr skel_dao =
                        db::_2d::DAOFactory::getStraightSkeletonDAO();
                int skelid = skel_dao->insert(skel_2d_);
                std::cout << "SkelID=" << skelid << std::endl;
            }
        }
    } else if (mode_ == MODE_3D) {
        if (polyhedron_) {
            db::_3d::PolyhedronDAOSPtr polyhedron_dao =
                    db::_3d::DAOFactory::getPolyhedronDAO();
            int polyhedron_id = polyhedron_dao->insert(skel_3d_->getPolyhedron());
            if (polyhedron_id > 0 && skel_3d_) {
                db::_3d::StraightSkeletonDAOSPtr skel_dao =
                        db::_3d::DAOFactory::getStraightSkeletonDAO();
                int skelid = skel_dao->insert(skel_3d_);
                std::cout << "SkelID=" << skelid << std::endl;
            }
        }
    }
}

void MainOpenGLWindow::saveLastPoly() {
    if (mode_ == MODE_2D) {
        PolygonSPtr polygon = polygon_;
        if (skel_2d_) {
            ReadLock l(skel_2d_->mutex());
            if (skel_2d_->events().size() > 0) {
                data::_2d::skel::AbstractEventSPtr event = skel_2d_->events().back();
                if (event->getPolygonResult()) {
                    polygon = event->getPolygonResult();
                }
            }
        }
        db::_2d::DAOFactory::getPolygonDAO()->insert(polygon);
        std::cout << "PolygonID=" << polygon->getID() << std::endl;
    } else if (mode_ == MODE_3D) {
        PolyhedronSPtr polyhedron = polyhedron_;
        if (skel_3d_) {
            ReadLock l(skel_3d_->mutex());
            if (skel_3d_->events().size() > 0) {
                data::_3d::skel::AbstractEventSPtr event = skel_3d_->events().back();
                if (event->getPolyhedronResult()) {
                    polyhedron = event->getPolyhedronResult();
                }
            }
        }
        db::_3d::DAOFactory::getPolyhedronDAO()->insert(polyhedron);
        std::cout << "PolyhedronID=" << polyhedron->getID() << std::endl;

        std::string now = util::StringFactory::now(util::StringFactory::DATE_FORMAT);
        std::string filename_obj = now + ".obj";
        if (db::_3d::OBJFile::save(filename_obj, polyhedron)) {
            std::cout << "filename_obj=" << filename_obj << std::endl;
        }
    }
}

void MainOpenGLWindow::dumpWin() {
    std::string now = util::StringFactory::now("%Y-%m-%d_%H%M%S%f");
    std::string filename_bmp = now + ".bmp";
    DEBUG_VAR(filename_bmp);
    dumpWindow(filename_bmp.c_str());
}

void MainOpenGLWindow::printScreen() {
    if (mode_ == MODE_2D) {
        ui::ps::PlanePSPrinterSPtr printer2 = ui::ps::PlanePSPrinter::create();
        printer2->setScale(printer2->getScale() * scale_);
        printer2->initBoundingBox(polygon_);
        std::string now = util::StringFactory::now(util::StringFactory::DATE_FORMAT);
        std::string filename_eps = now + ".eps";
        DEBUG_VAR(filename_eps);
        std::ofstream of;
        of.open(filename_eps.c_str());
        if (of.is_open()) {
            printer2->printHead(of);
            printer2->printComment(polygon_->getDescription(), of);
            if (polygon_ && (1 <= toggle_poly_ && toggle_poly_ <= 3)) {
                printer2->setLineWidth(1.0, of);
                printer2->printPolygon(polygon_, of);
            }
            if (skel_2d_) {
                if (toggle_skel_) {
                    printer2->setLineWidth(0.5, of);
                    printer2->printSkel(skel_2d_, of);
                }
                if (toggle_poly_ >= 3) {
                    ReadLock l(skel_2d_->mutex());
                    std::list<data::_2d::skel::AbstractEventSPtr>::reverse_iterator it_e =
                            skel_2d_->events().rbegin();
                    while (it_e != skel_2d_->events().rend()) {
                        data::_2d::skel::AbstractEventSPtr event = *it_e++;
                        std::stringstream stream_offset;
                        stream_offset << "offset=" << event->getOffset();
                        printer2->printComment(stream_offset.str(), of);
                        printer2->setLineWidth(1.0, of);
                        printer2->printPolygon(event->getPolygonResult(), of);
                        if (toggle_poly_ >= 4) {
                            break;
                        }
                    }
                }
            }
            if (mesh_2d_ && toggle_experimental_) {
                printer2->setLineWidth(0.1, of);
                printer2->printMesh(mesh_2d_, of);
            }
            of.close();
        }
    } else if (mode_ == MODE_3D || mode_ == MODE_SPHERICAL) {
        ui::ps::SpacePSPrinterSPtr printer3 = ui::ps::SpacePSPrinter::create();
        float modelview[16];
        float projection[16];
        int viewport[4];
        vec3f cam_eye;
        vec3f cam_center;
        for (unsigned int i = 0; i < 3; i++) {
            cam_eye[i] = (float)camera_->eye(i) / scale_ - translate_[i];
            cam_center[i] = (float)camera_->center(i) / scale_ - translate_[i];
        }
        glGetFloatv(GL_MODELVIEW_MATRIX, modelview);
        glGetFloatv(GL_PROJECTION_MATRIX, projection);
        glGetIntegerv(GL_VIEWPORT, viewport);
        printer3->setModelviewMatrix(modelview);
        printer3->setProjectionMatrix(projection);
        printer3->setViewport(viewport);
        printer3->setCamEye(cam_eye);
        printer3->setCamCenter(cam_center);
        std::string now = util::StringFactory::now(util::StringFactory::DATE_FORMAT);
        std::string filename_eps = now + ".eps";
        DEBUG_VAR(filename_eps);
        std::ofstream of;
        of.open(filename_eps.c_str());
        if (of.is_open()) {
            printer3->printHead(of);
            if (polyhedron_->getDescription().length() > 0) {
                printer3->printComment(polyhedron_->getDescription(), of);
            }
            printer3->printCommentCamera(of);
            if (mode_ == MODE_3D) {
                if (polyhedron_) {
                    if (toggle_poly_ == 1) {
                        printer3->setLineWidth(2.0, of);
                        printer3->printPolyhedronShade(polyhedron_,
                                0.25f, 0.95f, true, of);
                        printer3->setGray(0.0, of);
                    }
                    if (1 <= toggle_poly_ && toggle_poly_ <= 3) {
                        printer3->setLineWidth(1.0, of);
                        printer3->printPolyhedron(polyhedron_, of);
                    }
                }
                if (skel_3d_) {
                    if (toggle_skel_) {
                        printer3->setLineWidth(0.5, of);
                        printer3->printSkel(skel_3d_, of);
                    }
                    if (toggle_poly_ >= 3) {
                        ReadLock l(skel_3d_->mutex());
                        std::list<data::_3d::skel::AbstractEventSPtr>::reverse_iterator it_e =
                                skel_3d_->events().rbegin();
                        while (it_e != skel_3d_->events().rend()) {
                            data::_3d::skel::AbstractEventSPtr event = *it_e++;
                            std::stringstream stream_offset;
                            stream_offset << "offset=" << event->getOffset();
                            printer3->printComment(stream_offset.str(), of);
                            printer3->setLineWidth(1.0, of);
                            printer3->printPolyhedron(event->getPolyhedronResult(), of);
                            if (toggle_poly_ >= 4) {
                                break;
                            }
                        }
                    }
                }
            } else if (mode_ == MODE_SPHERICAL) {
                if (sphericalpolygon_) {
                    data::_3d::Point3SPtr p_center =
                            data::_3d::KernelFactory::createPoint3(sphericalpolygon_->getSphere());
                    vec3f center;
                    convert(p_center, center);
                    double radius = sphericalpolygon_->getRadius();
                    printer3->setLineWidth(0.25, of);
                    printer3->printSphere(center, radius, of);
                    printer3->setLineWidth(1.0, of);
                    printer3->printSphericalPolygon(sphericalpolygon_, of);
                    printer3->setLineWidth(0.25, of);
                    printer3->printSphericalIntersections(sphericalpolygon_, of);
                }
                if (sphericalskel_) {
                    if (toggle_skel_) {
                        printer3->setLineWidth(0.5, of);
                        printer3->printSphericalSkel(sphericalskel_, of);
                    }
                }
            }
            of.close();
        }
    }
}

void MainOpenGLWindow::printCutPattern() {
    if (mode_ == MODE_3D) {
        PolyhedronSPtr polyhedron = polyhedron_;
        if (skel_3d_) {
            ReadLock l(skel_3d_->mutex());
            if (skel_3d_->events().size() > 0) {
                data::_3d::skel::AbstractEventSPtr event = skel_3d_->events().back();
                if (event->getPolyhedronResult()) {
                    polyhedron = event->getPolyhedronResult();
                }
            }
        }
        ui::ps::CutPatternPrinter::printCutPattern(polyhedron, "cut_pattern");
    }
}

void MainOpenGLWindow::run() {
    this->createWindow();
    OpenGLWindow::mainLoop();
}


ThreadSPtr MainOpenGLWindow::startThread() {
    return ThreadSPtr(new std::thread(
            std::bind(&MainOpenGLWindow::run, this)));
}


void MainOpenGLWindow::setTranslate(const vec3f translate) {
    for (unsigned int i = 0; i < 3; i++) {
        translate_[i] = translate[i];
    }
}

void MainOpenGLWindow::setScale(float scale) {
    scale_ = scale;
}


void MainOpenGLWindow::setPolygon(PolygonSPtr polygon) {
    Lock l(mutex_);
    this->polygon_ = polygon;
    this->mode_ = MODE_2D;
}

void MainOpenGLWindow::setSkel2d(data::_2d::skel::StraightSkeletonSPtr skel_2d) {
    Lock l(mutex_);
    this->skel_2d_ = skel_2d;
    this->mode_ = MODE_2D;
}

void MainOpenGLWindow::setMesh2d(data::_2d::mesh::MeshSPtr mesh_2d) {
    Lock l(mutex_);
    this->mesh_2d_ = mesh_2d;
    this->mode_ = MODE_2D;
}

void MainOpenGLWindow::setPolyhedron(PolyhedronSPtr polyhedron) {
    Lock l(mutex_);
    this->polyhedron_ = polyhedron;
    this->mode_ = MODE_3D;
}

void MainOpenGLWindow::setSkel3d(data::_3d::skel::StraightSkeletonSPtr skel_3d) {
    Lock l(mutex_);
    this->skel_3d_ = skel_3d;
    this->mode_ = MODE_3D;
}

void MainOpenGLWindow::setSphericalPolygon(SphericalPolygonSPtr sphericalpolygon) {
    Lock l(mutex_);
    this->sphericalpolygon_ = sphericalpolygon;
    this->mode_ = MODE_SPHERICAL;
}

void MainOpenGLWindow::setSphericalSkel(data::_3d::skel::SphericalSkeletonSPtr sphericalskel) {
    Lock l(mutex_);
    this->sphericalskel_ = sphericalskel;
    this->mode_ = MODE_SPHERICAL;
}


void MainOpenGLWindow::drawCoordAxes(float size) {
    vec3f origin = {0.0, 0.0, 0.0};
    vec3f dir_x = {size, 0.0, 0.0};
    vec3f dir_y = {0.0, size, 0.0};
    vec3f dir_z = {0.0, 0.0, size};
    setColor(c_red);
    drawArrow(origin, dir_x);
    setColor(c_green);
    drawArrow(origin, dir_y);
    setColor(c_blue);
    drawArrow(origin, dir_z);
    setColor(c_white);
}


void MainOpenGLWindow::drawPolygon(PolygonSPtr polygon, float height, bool bold, bool vertices_only) {
    if (!polygon) {
        return;
    }
    ReadLock l(polygon->mutex());
    vec4f color_begin;
    getColor(color_begin);
    std::list<data::_2d::VertexSPtr>::iterator it_v = polygon->vertices().begin();
    while (it_v != polygon->vertices().end()) {
        data::_2d::VertexSPtr vertex = *it_v++;
        if (highlight_ && vertex->hasData()) {
            if (vertex->getData()->isHighlight()) {
                setColor(c_yellow);
            }
        }
        vec3f p;
        convert(vertex->getPoint(), p);
        p[2] = height;
        if (bold) {
            drawSphere(p, 0.2f * thickness_/scale_);
        } else {
            drawSphere(p, 0.1f * thickness_/scale_);
        }
        setColor(color_begin);
    }
    if (vertices_only) {
        return;
    }
    std::list<data::_2d::EdgeSPtr>::iterator it_e = polygon->edges().begin();
    while (it_e != polygon->edges().end()) {
        data::_2d::EdgeSPtr edge = *it_e++;
        data::_2d::VertexSPtr vertex_src = edge->getVertexSrc();
        data::_2d::VertexSPtr vertex_dst = edge->getVertexDst();
        if (vertex_src->getPoint() == vertex_dst->getPoint()) {
            continue;
        }
        if (highlight_ && edge->hasData()) {
            if (edge->getData()->isHighlight()) {
                setColor(c_yellow);
            }
        }
        vec3f src;
        vec3f dst;
        convert(vertex_src->getPoint(), src);
        convert(vertex_dst->getPoint(), dst);
        src[2] = height;
        dst[2] = height;
        if (bold) {
            drawCylinder(src, dst, 0.1f * thickness_/scale_);
        } else {
            drawCylinder(src, dst, 0.05f * thickness_/scale_);
        }
        setColor(color_begin);
    }
}


void MainOpenGLWindow::drawSkel2d(data::_2d::skel::StraightSkeletonSPtr skel_2d) {
    if (!skel_2d) {
        return;
    }
    ReadLock l(skel_2d->mutex());
    vec4f color_begin;
    getColor(color_begin);
    if (toggle_skel_) {
        std::list<data::_2d::skel::NodeSPtr>::iterator it_n = skel_2d->nodes().begin();
        while (it_n != skel_2d->nodes().end()) {
            data::_2d::skel::NodeSPtr node = *it_n++;
            vec3f p;
            convert(node->getPoint(), p);
            if (toggle_roof_) {
                p[2] = (float)node->getHeight();
            }
            drawSphere(p, 0.1f * thickness_/scale_);
        }
        std::list<data::_2d::skel::ArcSPtr>::iterator it_a = skel_2d->arcs().begin();
        while (it_a != skel_2d->arcs().end()) {
            data::_2d::skel::ArcSPtr arc = *it_a++;
            data::_2d::skel::NodeSPtr node = arc->getNodeSrc();
            vec3f p;
            convert(node->getPoint(), p);
            if (toggle_roof_) {
                p[2] = (float)node->getHeight();
            }
            if (arc->hasNodeDst()) {
                node = arc->getNodeDst();
                vec3f q;
                convert(node->getPoint(), q);
                if (toggle_roof_) {
                    q[2] = (float)node->getHeight();
                }
                drawCylinder(p, q, 0.05f * thickness_/scale_);
            } else {
                data::_2d::Vector2SPtr v_dir = arc->getDirection();
                if (v_dir) {
                    vec3f dir;
                    convert(v_dir, dir);
                    normalize(dir);
                    scale(dir, 0.5f * thickness_/scale_);
                    drawArrow(p, dir);
                }
            }
        }
    }
    if (toggle_poly_ >= 3) {
        setColor(c_white);
        bool vertices_only = false;
        std::list<data::_2d::skel::AbstractEventSPtr>::reverse_iterator it_e = skel_2d->events().rbegin();
        while (it_e != skel_2d->events().rend()) {
            data::_2d::skel::AbstractEventSPtr event = *it_e++;
            float height = 0.0f;
            if (toggle_roof_) {
                height = event->getOffset();
            }
            drawPolygon(event->getPolygonResult(), height, false, vertices_only);
            if (toggle_poly_ == 4) {
                vertices_only = true;
            } else if (toggle_poly_ == 5) {
                break;
            }
        }
    }
    setColor(color_begin);
}

void MainOpenGLWindow::drawMesh2d(data::_2d::mesh::MeshSPtr mesh_2d) {
    if (!mesh_2d) {
        return;
    }
    ReadLock l(mesh_2d->mutex());
    vec4f color_begin;
    getColor(color_begin);
    std::list<data::_2d::mesh::MeshCellSPtr>::iterator it_c = mesh_2d->cells().begin();
    while (it_c != mesh_2d->cells().end()) {
        data::_2d::mesh::MeshCellSPtr cell = *it_c++;
        if (cell->vertices().size() != 4) {
            setColor(c_red);
        }
        data::_2d::mesh::MeshVertexSPtr vertex_prev = cell->vertices().back();
        std::list<data::_2d::mesh::MeshVertexSPtr>::iterator it_v = cell->vertices().begin();
        while (it_v != cell->vertices().end()) {
            data::_2d::mesh::MeshVertexSPtr vertex = *it_v++;
            vec3f p_src;
            vec3f p_dst;
            convert(vertex_prev->getPoint(), p_src);
            convert(vertex->getPoint(), p_dst);
            drawSphere(p_src, 0.1f * thickness_/scale_);
            drawCylinder(p_src, p_dst, 0.05f * thickness_/scale_);
            vertex_prev = vertex;
        }
        setColor(color_begin);
    }
}

void MainOpenGLWindow::drawPolyhedron(PolyhedronSPtr polyhedron, bool bold, bool vertices_only) {
    if (!polyhedron) {
        return;
    }
    ReadLock l(polyhedron->mutex());
    vec4f color_begin;
    getColor(color_begin);
    std::list<data::_3d::VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        data::_3d::VertexSPtr vertex = *it_v++;
        if (highlight_ && vertex->hasData()) {
            if (vertex->getData()->isHighlight()) {
                setColor(c_yellow);
            }
        }
        vec3f p;
        convert(vertex->getPoint(), p);
        if (bold) {
            drawSphere(p, 0.2f * thickness_/scale_);
        } else {
            drawSphere(p, 0.1f * thickness_/scale_);
        }
        setColor(color_begin);
    }
    if (vertices_only) {
        return;
    }
    std::list<data::_3d::EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        data::_3d::EdgeSPtr edge = *it_e++;
        data::_3d::VertexSPtr vertex_src = edge->getVertexSrc();
        data::_3d::VertexSPtr vertex_dst = edge->getVertexDst();
        if (vertex_src->getPoint() == vertex_dst->getPoint()) {
            continue;
        }
        if (highlight_ && edge->hasData()) {
            if (edge->getData()->isHighlight()) {
                setColor(c_yellow);
            }
        }
        vec3f src;
        vec3f dst;
        convert(vertex_src->getPoint(), src);
        convert(vertex_dst->getPoint(), dst);
        if (bold) {
            drawCylinder(src, dst, 0.1f * thickness_/scale_);
        } else {
            drawCylinder(src, dst, 0.05f * thickness_/scale_);
        }
        if (draw_dirs_ && edge->hasData()) {
            data::_3d::skel::SkelEdgeDataSPtr data =
                    std::dynamic_pointer_cast<data::_3d::skel::SkelEdgeData>(edge->getData());
            data::_3d::skel::SheetSPtr sheet = data->getSheet();
            if (sheet) {
                data::_3d::Vector3SPtr normal =
                        data::_3d::KernelFactory::createVector3(sheet->getPlane());
                data::_3d::Vector3SPtr line =
                        data::_3d::KernelFactory::createVector3(edge->line());
                vec3f v_normal;
                convert(normal, v_normal);
                vec3f v_line;
                convert(line, v_line);
                vec3f v_dir;
                cross(v_normal, v_line, v_dir);
                normalize(v_dir);
                scale(v_dir, 0.5 * thickness_/scale_);
                vec3f pos = {
                    (float)((src[0]+dst[0])/2.0),
                    (float)((src[1]+dst[1])/2.0),
                    (float)((src[2]+dst[2])/2.0)};
                drawArrow(pos, v_dir);
            }
        }
        setColor(color_begin);
    }
    if (toggle_poly_ >= 2) {
        return;
    }
    setColor(c_trans_grey);
    std::list<data::_3d::FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        data::_3d::FacetSPtr facet = *it_f++;
        std::list<data::_3d::TriangleSPtr>::iterator it_t = facet->triangles().begin();
        while (it_t != facet->triangles().end()) {
            data::_3d::TriangleSPtr triangle = *it_t++;
            vec3f a;
            vec3f b;
            vec3f c;
            convert(triangle->getVertex(0)->getPoint(), a);
            convert(triangle->getVertex(1)->getPoint(), b);
            convert(triangle->getVertex(2)->getPoint(), c);
            drawTriangle(a, b, c);

            data::_3d::Vector3SPtr normal =
                    data::_3d::KernelFactory::createVector3(triangle->plane());
            vec3f pos_norm;
            vec3f dir_norm;
            double length_norm = 0.0f;
            for (unsigned int i = 0; i < 3; i++) {
                pos_norm[i] = (a[i] + b[i] + c[i]) / 3;
                length_norm += (*normal)[i] * (*normal)[i];
            }
            length_norm = sqrt(length_norm);
            for (unsigned int i = 0; i < 3; i++) {
                dir_norm[i] = ((*normal)[i] / length_norm) * 0.5 * thickness_/scale_;
            }
            drawArrow(pos_norm, dir_norm);
        }
    }
    setColor(color_begin);
}

void MainOpenGLWindow::drawSkel3d(data::_3d::skel::StraightSkeletonSPtr skel_3d) {
    if (!skel_3d) {
        return;
    }
    ReadLock l(skel_3d->mutex());
    vec4f color_begin;
    getColor(color_begin);
    if (toggle_skel_) {
        std::list<data::_3d::skel::NodeSPtr>::iterator it_n = skel_3d->nodes().begin();
        while (it_n != skel_3d->nodes().end()) {
            data::_3d::skel::NodeSPtr node = *it_n++;
            vec3f p;
            convert(node->getPoint(), p);
            drawSphere(p, 0.1f * thickness_/scale_);
        }
        std::list<data::_3d::skel::ArcSPtr>::iterator it_a = skel_3d->arcs().begin();
        while (it_a != skel_3d->arcs().end()) {
            data::_3d::skel::ArcSPtr arc = *it_a++;
            data::_3d::skel::NodeSPtr node = arc->getNodeSrc();
            vec3f p;
            convert(node->getPoint(), p);
            if (arc->hasNodeDst()) {
                node = arc->getNodeDst();
                vec3f q;
                convert(node->getPoint(), q);
                drawCylinder(p, q, 0.05f * thickness_/scale_);
            } else {
                data::_3d::Vector3SPtr v_dir = arc->getDirection();
                if (v_dir) {
                    vec3f dir;
                    convert(v_dir, dir);
                    normalize(dir);
                    scale(dir, 0.5f * thickness_/scale_);
                    drawArrow(p, dir);
                }
            }
        }
    }
    if (toggle_poly_ >= 3) {
        setColor(c_white);
        bool vertices_only = false;
        std::list<data::_3d::skel::AbstractEventSPtr>::reverse_iterator it_e = skel_3d->events().rbegin();
        while (it_e != skel_3d->events().rend()) {
            data::_3d::skel::AbstractEventSPtr event = *it_e++;
            drawPolyhedron(event->getPolyhedronResult(), false, vertices_only);
            if (toggle_poly_ == 4) {
                vertices_only = true;
            } else if (toggle_poly_ == 5) {
                break;
            }
        }
    }
    setColor(color_begin);
}


void MainOpenGLWindow::drawColoredNodes(data::_3d::skel::StraightSkeletonSPtr skel_3d) {
    if (!skel_3d) {
        return;
    }
    ReadLock l(skel_3d->mutex());
    vec4f color_begin;
    getColor(color_begin);
    std::list<data::_3d::skel::AbstractEventSPtr>::iterator it_e = skel_3d->events().begin();
    while (it_e != skel_3d->events().end()) {
        data::_3d::skel::AbstractEventSPtr event = *it_e++;
        data::_3d::skel::NodeSPtr node;
        if (event->getType() == data::_3d::skel::AbstractEvent::EDGE_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::EdgeEvent>(event)->getNode();
            setColor(c_blue);
        } else if (event->getType() == data::_3d::skel::AbstractEvent::EDGE_MERGE_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::EdgeMergeEvent>(event)->getNode();
            setColor(c_blue);
        } else if (event->getType() == data::_3d::skel::AbstractEvent::TRIANGLE_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::TriangleEvent>(event)->getNode();
            setColor(c_blue);
        } else if (event->getType() == data::_3d::skel::AbstractEvent::DBL_EDGE_MERGE_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::DblEdgeMergeEvent>(event)->getNode();
            setColor(c_blue);
        } else if (event->getType() == data::_3d::skel::AbstractEvent::DBL_TRIANGLE_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::DblTriangleEvent>(event)->getNode();
            setColor(c_blue);
        } else if (event->getType() == data::_3d::skel::AbstractEvent::TETRAHEDRON_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::TetrahedronEvent>(event)->getNode();
            setColor(c_blue);
        } else if (event->getType() == data::_3d::skel::AbstractEvent::VERTEX_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::VertexEvent>(event)->getNode();
            setColor(c_red);
        } else if (event->getType() == data::_3d::skel::AbstractEvent::FLIP_VERTEX_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::FlipVertexEvent>(event)->getNode();
            setColor(c_red);
        } else if (event->getType() == data::_3d::skel::AbstractEvent::SURFACE_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::SurfaceEvent>(event)->getNode();
            setColor(c_red);
        } else if (event->getType() == data::_3d::skel::AbstractEvent::POLYHEDRON_SPLIT_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::PolyhedronSplitEvent>(event)->getNode();
            setColor(c_red);
        } else if (event->getType() == data::_3d::skel::AbstractEvent::SPLIT_MERGE_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::SplitMergeEvent>(event)->getNode();
            setColor(c_red);
        } else if (event->getType() == data::_3d::skel::AbstractEvent::EDGE_SPLIT_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::EdgeSplitEvent>(event)->getNode();
            setColor(c_red);
        } else if (event->getType() == data::_3d::skel::AbstractEvent::PIERCE_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::PierceEvent>(event)->getNode();
            setColor(c_red);
        } else {
            continue;  // CONST_OFFSET_EVENT
        }
        vec3f p;
        convert(node->getPoint(), p);
        drawSphere(p, 0.11f * thickness_/scale_);
    }
    setColor(color_begin);
}


void MainOpenGLWindow::drawSphericalPolygon(SphericalPolygonSPtr sphericalpolygon, bool bold, bool vertices_only) {
    if (!sphericalpolygon) {
        return;
    }
    ReadLock l(sphericalpolygon->mutex());
    vec4f color_begin;
    getColor(color_begin);
    double radius = sphericalpolygon->getRadius();
    data::_3d::Sphere3SPtr sphere = sphericalpolygon->getSphere();
    data::_3d::Point3SPtr p_center =
            data::_3d::KernelFactory::createPoint3(sphere);
    vec3f center;
    convert(p_center, center);
    std::list<data::_3d::CircularVertexSPtr>::iterator it_v = sphericalpolygon->vertices().begin();
    while (it_v != sphericalpolygon->vertices().end()) {
        data::_3d::CircularVertexSPtr vertex = *it_v++;
        if (!vertex->isPointValid()) {
            continue;
        }
        if (highlight_ && vertex->hasData()) {
            if (vertex->getData()->isHighlight()) {
                setColor(c_yellow);
            }
        }
        vec3f p;
        convert(vertex->getPoint(), p);
        if (bold) {
            drawSphere(p, 0.2f * thickness_/scale_);
        } else {
            drawSphere(p, 0.1f * thickness_/scale_);
        }
        setColor(color_begin);
    }
    if (vertices_only) {
        return;
    }
    std::list<data::_3d::CircularEdgeSPtr>::iterator it_e = sphericalpolygon->edges().begin();
    while (it_e != sphericalpolygon->edges().end()) {
        data::_3d::CircularEdgeSPtr edge = *it_e++;
        data::_3d::CircularVertexSPtr vertex_src = edge->getVertexSrc();
        data::_3d::CircularVertexSPtr vertex_dst = edge->getVertexDst();
        if (vertex_src->getPoint() == vertex_dst->getPoint()) {
            continue;
        }
        if (!(vertex_src->isPointValid() && vertex_dst->isPointValid())) {
            // continue;
            vec4f c_darker;
            for (unsigned int i = 0; i < 3; i++) {
                c_darker[i] = color_begin[i] * 0.5f;
            }
            c_darker[3] = color_begin[3];
            setColor(c_darker);
        }
        if (highlight_ && edge->hasData()) {
            if (edge->getData()->isHighlight()) {
                setColor(c_yellow);
            }
        }
        data::_3d::Plane3SPtr plane_edge = edge->supportingPlane();
        data::_3d::Vector3SPtr normal_edge =
                algo::_3d::KernelWrapper::normalize(
                data::_3d::KernelFactory::createVector3(plane_edge));
        data::_3d::Point3SPtr p_normal_edge =
                data::_3d::KernelFactory::createPoint3(
                (*p_center) + (*normal_edge) * radius);
        vec3f axis;
        convert(normal_edge, axis);

        data::_3d::Point3SPtr p_src;
        data::_3d::Point3SPtr p_dst;
        if (vertex_src->isPointValid()) {
            p_src = vertex_src->getPoint();
        } else {
            data::_3d::CircularEdgeSPtr edge_in = vertex_src->getEdgeIn();
            data::_3d::Plane3SPtr plane_in = edge_in->supportingPlane();
            data::_3d::Vector3SPtr normal_in =
                    algo::_3d::KernelWrapper::normalize(
                    data::_3d::KernelFactory::createVector3(plane_in));
            data::_3d::Point3SPtr p_normal_in =
                    data::_3d::KernelFactory::createPoint3(
                    (*p_center) + (*normal_in) * radius);
            data::_3d::Plane3SPtr plane_jump =
                    data::_3d::KernelFactory::createPlane3(
                    p_center, p_normal_edge, p_normal_in);
            data::_3d::Plane3SPtr plane = plane_edge;
            if (vertex_src->hasData()) {
                data::_3d::skel::SphericalSkelVertexDataSPtr data_src =
                        std::dynamic_pointer_cast<data::_3d::skel::SphericalSkelVertexData>(
                        vertex_src->getData());
                data::_3d::skel::CircularArcSPtr arc_src = data_src->getArc();
                if (arc_src) {
                    data::_3d::CircularEdgeSPtr edge_origin = arc_src->getEdgeRight();
                    if (algo::_3d::KernelWrapper::angle(
                            edge_origin->getSupportingPlane(), plane_edge) > M_PI/2.0) {
                        // inversion happened
                        plane = algo::_3d::KernelWrapper::opposite(plane_edge);
                    }
                }
            }
            data::_3d::Line3SPtr line_in =
                    algo::_3d::KernelWrapper::intersection(plane_jump, plane);
            p_src = algo::_3d::KernelWrapper::intersection(sphere, line_in);
        }
        if (vertex_dst->isPointValid()) {
            p_dst = vertex_dst->getPoint();
        } else {
            data::_3d::CircularEdgeSPtr edge_out = vertex_dst->getEdgeOut();
            data::_3d::Plane3SPtr plane_out = edge_out->supportingPlane();
            data::_3d::Vector3SPtr normal_out =
                    algo::_3d::KernelWrapper::normalize(
                    data::_3d::KernelFactory::createVector3(plane_out));
            data::_3d::Point3SPtr p_normal_out =
                    data::_3d::KernelFactory::createPoint3(
                    (*p_center) + (*normal_out) * radius);
            data::_3d::Plane3SPtr plane_jump =
                    data::_3d::KernelFactory::createPlane3(
                    p_center, p_normal_edge, p_normal_out);
            data::_3d::Plane3SPtr plane = plane_edge;
            if (vertex_dst->hasData()) {
                data::_3d::skel::SphericalSkelVertexDataSPtr data_dst =
                        std::dynamic_pointer_cast<data::_3d::skel::SphericalSkelVertexData>(
                        vertex_dst->getData());
                data::_3d::skel::CircularArcSPtr arc_dst = data_dst->getArc();
                if (arc_dst) {
                    data::_3d::CircularEdgeSPtr edge_origin = arc_dst->getEdgeLeft();
                    if (algo::_3d::KernelWrapper::angle(
                            edge_origin->getSupportingPlane(), plane_edge) > M_PI/2.0) {
                        // inversion happened
                        plane = algo::_3d::KernelWrapper::opposite(plane_edge);
                    }
                }
            }
            data::_3d::Line3SPtr line_out =
                    algo::_3d::KernelWrapper::intersection(plane_jump, plane);
            p_dst = algo::_3d::KernelWrapper::intersection(sphere, line_out);
        }

        if (p_src && p_dst) {
            vec3f src;
            vec3f dst;
            convert(p_src, src);
            convert(p_dst, dst);
            if (bold) {
                drawCircularCylinder(center, axis, src, dst, 0.1f * thickness_/scale_);
            } else {
                drawCircularCylinder(center, axis, src, dst, 0.05f * thickness_/scale_);
            }
        }
        setColor(color_begin);
    }
}

void MainOpenGLWindow::drawSphericalSkel(data::_3d::skel::SphericalSkeletonSPtr sphericalskel) {
    if (!sphericalskel) {
        return;
    }
    ReadLock l(sphericalskel->mutex());
    vec4f color_begin;
    getColor(color_begin);
    if (toggle_skel_) {
        data::_3d::Point3SPtr p_center =
                data::_3d::KernelFactory::createPoint3(sphericalskel->getSphere());
        vec3f center;
        convert(p_center, center);
        std::list<data::_3d::skel::CircularNodeSPtr>::iterator it_n = sphericalskel->nodes().begin();
        while (it_n != sphericalskel->nodes().end()) {
            data::_3d::skel::CircularNodeSPtr node = *it_n++;
            vec3f p;
            convert(node->getPoint(), p);
            drawSphere(p, 0.1f * thickness_/scale_);
        }
        std::list<data::_3d::skel::CircularArcSPtr>::iterator it_a = sphericalskel->arcs().begin();
        while (it_a != sphericalskel->arcs().end()) {
            data::_3d::skel::CircularArcSPtr arc = *it_a++;
            data::_3d::skel::CircularNodeSPtr node = arc->getNodeSrc();
            vec3f src;
            convert(node->getPoint(), src);
            if (arc->hasNodeDst()) {
                node = arc->getNodeDst();
                vec3f dst;
                convert(node->getPoint(), dst);
                data::_3d::Vector3SPtr normal =
                        data::_3d::KernelFactory::createVector3(arc->getSupportingPlane());
                vec3f axis;
                convert(normal, axis);
                drawCircularCylinder(center, axis, src, dst, 0.05f * thickness_/scale_);
            } else {
                data::_3d::Vector3SPtr v_dir = arc->getDirection();
                if (v_dir) {
                    vec3f dir;
                    convert(v_dir, dir);
                    normalize(dir);
                    scale(dir, 0.5f * thickness_/scale_);
                    drawArrow(src, dir);
                }
            }
        }
    }
    if (toggle_poly_ >= 3) {
        setColor(c_white);
        bool vertices_only = false;
        std::list<data::_3d::skel::SphericalAbstractEventSPtr>::reverse_iterator it_e = sphericalskel->events().rbegin();
        while (it_e != sphericalskel->events().rend()) {
            data::_3d::skel::SphericalAbstractEventSPtr event = *it_e++;
            drawSphericalPolygon(event->getPolygonResult(), false, vertices_only);
            if (toggle_poly_ == 4) {
                vertices_only = true;
            } else if (toggle_poly_ == 5) {
                break;
            }
        }
    }
    setColor(color_begin);
}


void MainOpenGLWindow::drawContent() {
    bool screenshot = false;
    if (controller_) {
        PolygonSPtr polygon = controller_->getDispPolygon();
        if (polygon) setPolygon(polygon);
        data::_2d::skel::StraightSkeletonSPtr skel_2d = controller_->getDispSkel2d();
        if (skel_2d) setSkel2d(skel_2d);
        PolyhedronSPtr polyhedron = controller_->getDispPolyhedron();
        if (polyhedron) setPolyhedron(polyhedron);
        data::_3d::skel::StraightSkeletonSPtr skel_3d = controller_->getDispSkel3d();
        if (skel_3d) setSkel3d(skel_3d);
        SphericalPolygonSPtr sphericalpolygon = controller_->getDispSphericalPolygon();
        if (sphericalpolygon) setSphericalPolygon(sphericalpolygon);
        data::_3d::skel::SphericalSkeletonSPtr sphericalskel = controller_->getDispSphericalSkel();
        if (sphericalskel) setSphericalSkel(sphericalskel);
        screenshot = controller_->getScreenshot();
    }

    Lock l(mutex_);
    if (coord_axes_size_ > 0.0f) {
        drawCoordAxes(coord_axes_size_);
    }
    glScalef(scale_, scale_, scale_);
    glTranslatef(translate_[0], translate_[1], translate_[2]);
    if (mode_ == MODE_2D) {
        if (skel_2d_) {
            setColor(c_grey);
            drawSkel2d(skel_2d_);
        }
        if (mesh_2d_ && toggle_experimental_) {
            setColor(c_blue);
            drawMesh2d(mesh_2d_);
        }
        if (polygon_ && (1 <= toggle_poly_ && toggle_poly_ <= 3)) {
            setColor(c_white);
            drawPolygon(polygon_, 0.0f, true, false);
        }
    } else if (mode_ == MODE_3D) {
        if (skel_3d_) {
            setColor(c_grey);
            drawSkel3d(skel_3d_);
        }
        if (skel_3d_ && toggle_experimental_) {
            drawColoredNodes(skel_3d_);
        }
        if (polyhedron_ && (1 <= toggle_poly_ && toggle_poly_ <= 3)) {
            setColor(c_white);
            drawPolyhedron(polyhedron_, true, false);
        }
    } else if (mode_ == MODE_SPHERICAL) {
        if (sphericalskel_) {
            setColor(c_grey);
            drawSphericalSkel(sphericalskel_);
        }
        if (sphericalpolygon_ && (1 <= toggle_poly_ && toggle_poly_ <= 3)) {
            setColor(c_white);
            drawSphericalPolygon(sphericalpolygon_, true, false);
        }
        if (polyhedron_ && (1 <= toggle_poly_ && toggle_poly_ <= 3)) {
            setColor(c_white);
            drawPolyhedron(polyhedron_, true, false);
        }
        if (sphericalpolygon_) {
            data::_3d::Point3SPtr p_center =
                    data::_3d::KernelFactory::createPoint3(sphericalpolygon_->getSphere());
            vec3f center;
            convert(p_center, center);
            double radius = sphericalpolygon_->getRadius();
            setColor(c_trans_grey);
            drawSphere(center, radius);
        }
    }
    if (crosshair_size_ > 0.0f) {
        drawCrosshair(crosshair_size_);
    }
    if (screenshot) {
        dumpWin();
    }
}

void MainOpenGLWindow::handleKeyPressed(int key, int x, int y) {
    this->keyboard_adapter_->pressed(key, x, y);
}

void MainOpenGLWindow::handleKeyReleased(int key, int x, int y) {
    this->keyboard_adapter_->released(key, x, y);
}

void MainOpenGLWindow::handleMousePressed(int button, int x, int y) {
    this->mouse_adapter_->pressed(button, x, y);
}

void MainOpenGLWindow::handleMouseReleased(int button, int x, int y) {
    this->mouse_adapter_->released(button, x, y);
}

void MainOpenGLWindow::handleMotion(int x, int y) {
    this->mouse_adapter_->dragged(x, y);
}

} }
