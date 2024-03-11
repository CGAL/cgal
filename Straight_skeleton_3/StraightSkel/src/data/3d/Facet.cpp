/**
 * @file   data/3d/Facet.cpp
 * @author Gernot Walzl
 * @date   2011-11-26
 */

#include "data/3d/Facet.h"

#include "data/2d/Vertex.h"
#include "data/2d/Edge.h"
#include "data/2d/Polygon.h"
#include "data/2d/KernelFactory.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Triangle.h"
#include "data/3d/Polyhedron.h"
#include "data/3d/KernelFactory.h"
#include "util/StringFactory.h"
#include <cmath>
#include <sstream>

namespace data { namespace _3d {

Facet::Facet() {
    this->id_ = -1;
}

Facet::~Facet() {
    triangles_.clear();
    edges_.clear();
    vertices_.clear();
}

FacetSPtr Facet::create() {
    FacetSPtr result = FacetSPtr(new Facet());
    return result;
}

FacetSPtr Facet::create(unsigned int num_vertices, VertexSPtr vertices[]) {
    FacetSPtr result = FacetSPtr(new Facet());
    for (unsigned int i = 0; i < num_vertices; i++) {
        result->addVertex(vertices[i]);
    }
    for (unsigned int i = 0; i < num_vertices; i++) {
        EdgeSPtr edge = vertices[i]->findEdge(vertices[(i+1)%num_vertices]);
        if (!edge) {
            edge = Edge::create(vertices[i], vertices[(i+1)%num_vertices]);
        }
        result->addEdge(edge);
    }
    return result;
}

FacetSPtr Facet::create(unsigned int num_edges, EdgeSPtr edges[]) {
    FacetSPtr result = FacetSPtr(new Facet());
    for (unsigned int i = 0; i < num_edges; i++) {
        result->addEdge(edges[i]);
    }
    return result;
}

FacetSPtr Facet::clone() const {
    FacetSPtr result = Facet::create();
    std::list<VertexSPtr>::const_iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
        VertexSPtr vertex = *it_v++;
        VertexSPtr vertex_c = vertex->clone();
        result->addVertex(vertex_c);
    }
    std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
        Point3SPtr p_src = edge->getVertexSrc()->getPoint();
        Point3SPtr p_dst = edge->getVertexDst()->getPoint();
        VertexSPtr src;
        VertexSPtr dst;
        std::list<VertexSPtr>::const_iterator it_v = vertices_.begin();
        while (it_v != vertices_.end()) {
            VertexSPtr vertex = *it_v++;
            if (vertex->getPoint() == p_src) {
                src = vertex;
            }
            if (vertex->getPoint() == p_dst) {
                dst = vertex;
            }
        }
        EdgeSPtr edge_c = Edge::create(src, dst);
        if (edge->getFacetL() == shared_from_this()) {
            edge_c->setFacetL(result);
        }
        if (edge->getFacetR() == shared_from_this()) {
            edge_c->setFacetR(result);
        }
        result->addEdge(edge_c);
    }
    return result;
}

void Facet::addVertex(VertexSPtr vertex) {
    vertices_.insert(vertices_.end(), vertex);
    vertex->addFacet(shared_from_this());
}

bool Facet::removeVertex(VertexSPtr vertex) {
    bool result = false;
    std::list<VertexSPtr>::iterator it_v =
            std::find(vertices_.begin(), vertices_.end(), vertex);
    if (it_v != vertices_.end()) {
        vertices_.erase(it_v);
        vertex->removeFacet(shared_from_this());
        std::list<TriangleSPtr>::iterator it_t = triangles_.begin();
        while (it_t != triangles_.end()) {
            TriangleSPtr triangle = *it_t++;
            for (unsigned int i = 0; i < 3; i++) {
                if (triangle->getVertex(i) == vertex) {
                    removeTriangle(triangle);
                    break;
                }
            }
        }
        result = true;
    }
    return result;
}

void Facet::addEdge(EdgeSPtr edge) {
    std::list<EdgeSPtr>::iterator it = edges_.insert(edges_.end(), edge);
    FacetSPtr facet_l = edge->getFacetL();
    FacetSPtr facet_r = edge->getFacetR();
    if (facet_l == shared_from_this() && facet_r == shared_from_this()) {
        std::list<EdgeSPtr>::iterator it_e =
                std::find(edges_.begin(), edges_.end(), edge);
        if (it_e == edge->getFacetLListIt()) {
            edge->setFacetRListIt(it);
        } else {
            edge->setFacetLListIt(it);
        }
    } else if (facet_l == shared_from_this()) {
        edge->setFacetLListIt(it);
    } else if (facet_r == shared_from_this()) {
        edge->setFacetRListIt(it);
    } else if (!facet_l) {
        edge->setFacetL(shared_from_this());
        edge->setFacetLListIt(it);
    } else if (!facet_r) {
        edge->setFacetR(shared_from_this());
        edge->setFacetRListIt(it);
    } else {
        DEBUG_VAR(edge->toString());
        throw std::runtime_error("The given edge already has a left and a right facet.");
    }
    VertexSPtr vertex_src = edge->src(shared_from_this());
    if (!containsVertex(vertex_src)) {
        addVertex(vertex_src);
    }
    VertexSPtr vertex_dst = edge->dst(shared_from_this());
    if (!containsVertex(vertex_dst)) {
        addVertex(vertex_dst);
    }
}

bool Facet::removeEdge(EdgeSPtr edge) {
    bool result = false;
    if (edge->getFacetL() == shared_from_this()) {
        edges_.erase(edge->getFacetLListIt());
        edge->setFacetL(FacetSPtr());
        edge->setFacetLListIt(std::list<EdgeSPtr>::iterator());
        result = true;
    } else if (edge->getFacetR() == shared_from_this()) {
        edges_.erase(edge->getFacetRListIt());
        edge->setFacetR(FacetSPtr());
        edge->setFacetRListIt(std::list<EdgeSPtr>::iterator());
        result = true;
    }
    return result;
}

EdgeSPtr Facet::findEdge(FacetSPtr facet) const {
    EdgeSPtr result = EdgeSPtr();
    std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
        if (edge->getFacetL() == facet || edge->getFacetR() == facet) {
            result = edge;
            break;
        }
    }
    return result;
}

std::list<EdgeSPtr> Facet::findEdges(FacetSPtr facet) const {
    std::list<EdgeSPtr> result;
    std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
        if (edge->getFacetL() == facet || edge->getFacetR() == facet) {
            result.push_back(edge);
        }
    }
    return result;
}

bool Facet::containsVertex(VertexSPtr vertex) const {
    bool result = (vertices_.end() !=
            std::find(vertices_.begin(), vertices_.end(), vertex));
    return result;
}

bool Facet::containsEdge(EdgeSPtr edge) const {
    bool result = (edges_.end() !=
            std::find(edges_.begin(), edges_.end(), edge));
    return result;
}

void Facet::sortVertices() {
    FacetSPtr self(shared_from_this());
    std::list<VertexSPtr> tmp;
    std::list<VertexSPtr>::iterator it_v_tmp;
    while (vertices_.size() > 0) {
        VertexSPtr first = VertexSPtr();
        VertexSPtr vertex = vertices_.front();
        EdgeSPtr edge = EdgeSPtr();
        std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
        while (it_e != edges_.end()) {
            EdgeSPtr my_edge = *it_e++;
            if (vertex == my_edge->src(self)) {
                edge = my_edge;
                break;
            }
        }
        while (vertex != first) {
            if (!first) {
                first = vertex;
                it_v_tmp = tmp.insert(tmp.end(), vertex);
            } else {
                tmp.insert(tmp.end(), vertex);
            }
            vertex = edge->dst(self);
            edge = edge->next(self);
        }
        while (it_v_tmp != tmp.end()) {
            VertexSPtr vertex = *it_v_tmp++;
            std::list<VertexSPtr>::iterator it_v =
                    std::find(vertices_.begin(), vertices_.end(), vertex);
            if (it_v != vertices_.end()) {
                vertices_.erase(it_v);
            }
        }
    }
    vertices_.clear();
    it_v_tmp = tmp.begin();
    while (it_v_tmp != tmp.end()) {
        VertexSPtr vertex = *it_v_tmp++;
        vertices_.insert(vertices_.end(), vertex);
    }
}

void Facet::sortEdges() {
    std::list<EdgeSPtr> tmp;
    std::list<EdgeSPtr>::iterator it_e_tmp;
    while (edges_.size() > 0) {
        EdgeSPtr first = EdgeSPtr();
        EdgeSPtr edge = edges_.front();
        while (edge != first) {
            if (!first) {
                first = edge;
                it_e_tmp = tmp.insert(tmp.end(), edge);
            } else {
                tmp.insert(tmp.end(), edge);
            }
            edge = edge->next(shared_from_this());
        }
        while (it_e_tmp != tmp.end()) {
            EdgeSPtr edge = *it_e_tmp++;
            removeEdge(edge);
        }
    }
    edges_.clear();
    it_e_tmp = tmp.begin();
    while (it_e_tmp != tmp.end()) {
        EdgeSPtr edge = *it_e_tmp++;
        addEdge(edge);
    }
}

void Facet::addTriangle(TriangleSPtr triangle) {
    std::list<TriangleSPtr>::iterator it = triangles_.insert(triangles_.end(), triangle);
    triangle->setFacet(shared_from_this());
    triangle->setFacetListIt(it);
    for (unsigned int i = 0; i < 3; i++) {
        VertexSPtr vertex = triangle->getVertex(i);
        if (!containsVertex(vertex)) {
            addVertex(vertex);
        }
    }
}

bool Facet::removeTriangle(TriangleSPtr triangle) {
    bool result = false;
    if (triangle->getFacet() == shared_from_this()) {
        triangles_.erase(triangle->getFacetListIt());
        triangle->setFacet(FacetSPtr());
        triangle->setFacetListIt(std::list<TriangleSPtr>::iterator());
        result = true;
    }
    return result;
}

PolyhedronSPtr Facet::getPolyhedron() const {
    DEBUG_WPTR(this->polyhedron_);
    if (this->polyhedron_.expired())
        return PolyhedronSPtr();
    else
        return PolyhedronSPtr(this->polyhedron_);
}

void Facet::setPolyhedron(PolyhedronSPtr polyhedron) {
    this->polyhedron_ = polyhedron;
}

std::list<FacetSPtr>::iterator Facet::getPolyhedronListIt() const {
    return this->polyhedron_list_it_;
}

void Facet::setPolyhedronListIt(std::list<FacetSPtr>::iterator list_it) {
    this->polyhedron_list_it_ = list_it;
}

FacetDataSPtr Facet::getData() const {
    DEBUG_SPTR(this->data_);
    return this->data_;
}

void Facet::setData(FacetDataSPtr data) {
    this->data_ = data;
}

bool Facet::hasData() const {
    bool result = false;
    if (data_) {
        result = true;
    }
    return result;
}

std::list<VertexSPtr>& Facet::vertices() {
    return this->vertices_;
}

std::list<EdgeSPtr>& Facet::edges() {
    return this->edges_;
}

std::list<TriangleSPtr>& Facet::triangles() {
    return this->triangles_;
}

FacetSPtr Facet::next(VertexSPtr vertex) const {
    FacetSPtr result = FacetSPtr();
    std::list<FacetWPtr>::const_iterator it_f = vertex->facets().begin();
    while (it_f != vertex->facets().end()) {
        FacetWPtr facet_wptr = *it_f;
        if (!facet_wptr.expired()) {
            FacetSPtr facet(facet_wptr);
            if (facet == shared_from_this()) {
                if (vertex->degree() == 1) {
                    result = facet;
                }
                break;
            }
        }
        it_f++;
    }
    if (it_f != vertex->facets().end()) {
        std::list<FacetWPtr>::const_iterator it_f_begin = it_f++;
        if (it_f == vertex->facets().end()) {
            it_f = vertex->facets().begin();
        }
        while (it_f != it_f_begin) {
            FacetWPtr facet_wptr = *it_f++;
            if (it_f == vertex->facets().end()) {
                it_f = vertex->facets().begin();
            }
            if (!facet_wptr.expired()) {
                FacetSPtr facet(facet_wptr);
                std::list<EdgeWPtr>::const_iterator it_e = vertex->edges().begin();
                while (it_e != vertex->edges().end()) {
                    EdgeWPtr edge_wptr = *it_e++;
                    if (!edge_wptr.expired()) {
                        EdgeSPtr edge(edge_wptr);
                        FacetSPtr facet_l = edge->getFacetL();
                        FacetSPtr facet_r = edge->getFacetR();
                        if ((facet_l == shared_from_this() &&
                                    facet_r == facet &&
                                    edge->getVertexDst() == vertex) ||
                                (facet_r == shared_from_this() &&
                                    facet_l == facet &&
                                    edge->getVertexSrc() == vertex)) {
                            result = facet;
                            break;
                        }
                    }
                }
            }
        }
    }
    DEBUG_SPTR(result);
    return result;
}

FacetSPtr Facet::prev(VertexSPtr vertex) const {
    FacetSPtr result = FacetSPtr();
    std::list<FacetWPtr>::const_reverse_iterator it_f = vertex->facets().rbegin();
    while (it_f != vertex->facets().rend()) {
        FacetWPtr facet_wptr = *it_f;
        if (!facet_wptr.expired()) {
            FacetSPtr facet(facet_wptr);
            if (facet == shared_from_this()) {
                if (vertex->degree() == 1) {
                    result = facet;
                }
                break;
            }
        }
        it_f++;
    }
    if (it_f != vertex->facets().rend()) {
        std::list<FacetWPtr>::const_reverse_iterator it_f_begin = it_f++;
        if (it_f == vertex->facets().rend()) {
            it_f = vertex->facets().rbegin();
        }
        while (it_f != it_f_begin) {
            FacetWPtr facet_wptr = *it_f++;
            if (it_f == vertex->facets().rend()) {
                it_f = vertex->facets().rbegin();
            }
            if (!facet_wptr.expired()) {
                FacetSPtr facet(facet_wptr);
                std::list<EdgeWPtr>::const_iterator it_e = vertex->edges().begin();
                while (it_e != vertex->edges().end()) {
                    EdgeWPtr edge_wptr = *it_e++;
                    if (!edge_wptr.expired()) {
                        EdgeSPtr edge(edge_wptr);
                        FacetSPtr facet_l = edge->getFacetL();
                        FacetSPtr facet_r = edge->getFacetR();
                        if ((facet_l == shared_from_this() &&
                                    facet_r == facet &&
                                    edge->getVertexSrc() == vertex) ||
                                (facet_r == shared_from_this() &&
                                    facet_l == facet &&
                                    edge->getVertexDst() == vertex)) {
                            result = facet;
                            break;
                        }
                    }
                }
            }
        }
    }
    DEBUG_SPTR(result);
    return result;
}

void Facet::merge(FacetSPtr facet) {
    std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
    while (it_v != facet->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        facet->removeVertex(vertex);
        if (!containsVertex(vertex)) {
            addVertex(vertex);
        }
    }
    std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
    while (it_e != facet->edges().end()) {
        EdgeSPtr edge = *it_e++;
        FacetSPtr facet_l = edge->getFacetL();
        FacetSPtr facet_r = edge->getFacetR();
        if ((facet_l == facet && facet_r == shared_from_this()) ||
                (facet_r == facet && facet_l == shared_from_this())) {
            facet->removeEdge(edge);
            removeEdge(edge);
            edge->getPolyhedron()->removeEdge(edge);
        } else {
            if (facet_l == facet) {
                facet->removeEdge(edge);
                edge->setFacetL(shared_from_this());
                addEdge(edge);
            }
            if (facet_r == facet) {
                facet->removeEdge(edge);
                edge->setFacetR(shared_from_this());
                addEdge(edge);
            }
        }
    }
    std::list<TriangleSPtr>::iterator it_t = facet->triangles_.begin();
    while (it_t != facet->triangles().end()) {
        TriangleSPtr triangle = *it_t++;
        facet->removeTriangle(triangle);
        addTriangle(triangle);
    }
}

int Facet::getID() const {
  return this->id_;
}

void Facet::setID(int id) {
    this->id_ = id;
}

Plane3SPtr Facet::getPlane() const {
    DEBUG_SPTR(this->plane_);
    return this->plane_;
}

void Facet::setPlane(Plane3SPtr plane) {
    this->plane_ = plane;
}

bool Facet::initPlane() {
    bool result = false;
    std::list<TriangleSPtr>::iterator it_t = triangles_.begin();
    while (it_t != triangles_.end()) {
        TriangleSPtr triangle = *it_t++;
        Plane3SPtr plane = triangle->plane();
        if (plane) {
            plane_ = plane;
            result = true;
            break;
        }
    }
    if (!result) {
        Point3SPtr points[3];
        unsigned int i = 0;
        Point3SPtr point_prev;
        std::list<VertexSPtr>::iterator it_v = vertices_.begin();
        while (i < 3 && it_v != vertices_.end()) {
            VertexSPtr vertex = *it_v++;
            Point3SPtr point = vertex->getPoint();
            if (point_prev != point) {
                points[i] = point;
                i++;
            }
            point_prev = point;
        }
        if (i >= 3) {
            Plane3SPtr plane = KernelFactory::createPlane3(
                    points[0], points[1], points[2]);
            if (plane) {
                plane_ = plane;
                result = true;
            }
        }
    }
    return result;
}

Plane3SPtr Facet::plane() {
    if (!this->plane_) {
        this->initPlane();
    }
    DEBUG_SPTR(this->plane_);
    return this->plane_;
}

bool Facet::makeFirstConvex() {
    bool result = false;
    if (!plane_) {
        return false;
    }
    EdgeSPtr edge_begin;
    FacetSPtr self(shared_from_this());
    Vector3SPtr normal = KernelFactory::createVector3(plane_);
    EdgeSPtr edge = edges_.front();
    EdgeSPtr first = EdgeSPtr();
    while (edge != first) {
        if (!first) {
            first = edge;
        }
        EdgeSPtr edge_next = edge->next(self);
        Point3SPtr points[3];
        points[0] = edge->src(self)->getPoint();
        points[1] = edge->dst(self)->getPoint();
        points[2] = edge_next->dst(self)->getPoint();
        if (points[0] != points[1] &&
                points[1] != points[2] &&
                points[2] != points[0]) {
            Plane3SPtr plane_current = KernelFactory::createPlane3(
                    points[0], points[1], points[2]);
            Vector3SPtr normal_current = KernelFactory::createVector3(plane_current);
            double angle = 0.0;
            double arg = 0.0;
#ifdef USE_CGAL
            arg = ((*normal)*(*normal_current)) /
                    CGAL::sqrt(normal->squared_length() * normal_current->squared_length());
#else
            arg = ((*normal)*(*normal_current)) /
                    sqrt(normal->squared_length() * normal_current->squared_length());
#endif
            // fixes issues with floating point precision
            if (arg <= -1.0) {
                angle = M_PI;
            } else if (arg >= 1.0) {
                angle = 0.0;
            } else {
                angle = acos(arg);
            }
            if (angle < M_PI/4.0) {
                edge_begin = edge;
                result = true;
                break;
            }
        }
        edge = edge_next;
    }
    if (edge_begin) {
        std::list<EdgeSPtr>::iterator it_e = edges_.insert(edges_.begin(), edge_begin);
        if (edge_begin->getFacetL() == shared_from_this()) {
            edges_.erase(edge->getFacetLListIt());
            edge_begin->setFacetLListIt(it_e);
        } else if (edge_begin->getFacetR() == shared_from_this()) {
            edges_.erase(edge->getFacetRListIt());
            edge_begin->setFacetRListIt(it_e);
        }
        sortEdges();
        VertexSPtr vertex_begin = edge_begin->src(shared_from_this());
        std::list<VertexSPtr>::iterator it_v =
                std::find(vertices_.begin(), vertices_.end(), vertex_begin);
        vertices_.erase(it_v);
        vertices_.insert(vertices_.begin(), vertex_begin);
        sortVertices();
    }
    if (!result) {
        DEBUG_VAL("Unable to make first 3 vertices convex.");
        DEBUG_VAR(toString());
    }
    return result;
}

/**
 * http://de.wikipedia.org/wiki/Drehmatrix
 *
 *              [ n_x^2 (1 - cos(alpha)) + cos(alpha)         n_x n_y (1 - cos(alpha)) - n_z sin(alpha)   n_x n_z (1 - cos(alpha)) + n_y sin(alpha) ]
 * R_n(alpha) = [ n_y n_x (1 - cos(alpha)) + n_z sin(alpha)   n_y^2 (1 - cos(alpha)) + cos(alpha)         n_y n_z (1 - cos(alpha)) - n_x sin(alpha) ]
 *              [ n_z n_x (1 - cos(alpha)) - n_y sin(alpha)   n_z n_y (1 - cos(alpha)) + n_x sin(alpha)   n_z^2 (1 - cos(alpha)) + cos(alpha)       ]
 *
 * n_z = 0
 */
data::_2d::PolygonSPtr Facet::toPolygon() {
    data::_2d::PolygonSPtr result = data::_2d::Polygon::create();
    Vector3SPtr normal = KernelFactory::createVector3(plane());
    double n[3];   // n = normal.cross([0 0 1])
    n[0] = (*normal)[1];
    n[1] = -(*normal)[0];
    n[2] = 0.0;
    double length = 0.0;
    if ((*normal)[0] != 0.0 || (*normal)[1] != 0.0) {
        // normalize n
        for (unsigned int i = 0; i < 3; i++) {
            length += n[i]*n[i];
        }
        length = sqrt(length);
        for (unsigned int i = 0; i < 3; i++) {
            n[i] /= length;
        }
    }
    length = 0.0;
    for (unsigned int i = 0; i < 3; i++) {
        length += (*normal)[i] * (*normal)[i];
    }
    length = sqrt(length);
    // alpha = acos((normal * [0 0 1]) / normal.length())
    double alpha = acos((*normal)[2] / length);
    Vector3SPtr r_x = KernelFactory::createVector3(
            n[0] * n[0] * (1 - cos(alpha)) + cos(alpha),
            n[0] * n[1] * (1 - cos(alpha)) - n[2] * sin(alpha),
            n[0] * n[2] * (1 - cos(alpha)) + n[1] * sin(alpha));
    Vector3SPtr r_y = KernelFactory::createVector3(
            n[1] * n[0] * (1 - cos(alpha)) + n[2] * sin(alpha),
            n[1] * n[1] * (1 - cos(alpha)) + cos(alpha),
            n[1] * n[2] * (1 - cos(alpha)) - n[0] * sin(alpha));
    std::list<VertexSPtr>::const_iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
        VertexSPtr vertex = *it_v++;
        Vector3SPtr v_p = KernelFactory::createVector3(vertex->getPoint());
        double x2 = (*r_x) * (*v_p);
        double y2 = (*r_y) * (*v_p);
        data::_2d::Point2SPtr p2 = data::_2d::KernelFactory::createPoint2(x2, y2);
        data::_2d::VertexSPtr vertex2 = data::_2d::Vertex::create(p2);
        result->addVertex(vertex2);
    }
    std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
        VertexSPtr src = edge->src(shared_from_this());
        VertexSPtr dst = edge->dst(shared_from_this());
        data::_2d::VertexSPtr src2;
        data::_2d::VertexSPtr dst2;
        std::list<VertexSPtr>::const_iterator it_v = vertices_.begin();
        std::list<data::_2d::VertexSPtr>::const_iterator it_v2 = result->vertices().begin();
        while (it_v != vertices_.end() && it_v2 != result->vertices().end()) {
            VertexSPtr vertex = *it_v++;
            data::_2d::VertexSPtr vertex2 = *it_v2++;
            if (src == vertex) {
                src2 = vertex2;
            }
            if (dst == vertex) {
                dst2 = vertex2;
            }
        }
        data::_2d::EdgeSPtr edge2 = data::_2d::Edge::create(src2, dst2);
        result->addEdge(edge2);
    }
    return result;
}

std::string Facet::toString() const {
    std::stringstream sstr;
    sstr << "Facet(";
    if (id_ != -1) {
        sstr << "id=" << util::StringFactory::fromInteger(id_) << ", ";
    } else {
        sstr << util::StringFactory::fromPointer(this) << ", ";
    }
    if (plane_) {
#ifdef USE_CGAL
        sstr << "<" << util::StringFactory::fromDouble(plane_->a()) << ", "
             << util::StringFactory::fromDouble(plane_->b()) << ", "
             << util::StringFactory::fromDouble(plane_->c()) << ", "
             << util::StringFactory::fromDouble(plane_->d()) << ">, ";
#else
        sstr << "<" << util::StringFactory::fromDouble(plane_->getA()) << ", "
             << util::StringFactory::fromDouble(plane_->getB()) << ", "
             << util::StringFactory::fromDouble(plane_->getC()) << ", "
             << util::StringFactory::fromDouble(plane_->getD()) << ">, ";
#endif
    }
    sstr << "Vertices:" + util::StringFactory::fromInteger(vertices_.size()) + ", ";
    sstr << "Edges:" + util::StringFactory::fromInteger(edges_.size()) + ",";
    if (vertices_.size() > 0) {
        sstr << "\n";
        std::list<VertexSPtr>::const_iterator it_v = vertices_.begin();
        while (it_v != vertices_.end()) {
            VertexSPtr vertex = *it_v++;
            sstr << "\t" << vertex->toString() << "\n";
        }
    }
    if (edges_.size() > 0) {
        sstr << "\n";
        std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
            while (it_e != edges_.end()) {
            EdgeSPtr edge = *it_e++;
            sstr << "\t" << edge->toString() << "\n";
        }
    }
    if (triangles_.size() > 0) {
        sstr << "\n";
        std::list<TriangleSPtr>::const_iterator it_t = triangles_.begin();
        while (it_t != triangles_.end()) {
            TriangleSPtr triangle = *it_t++;
            sstr << "\t" << triangle->toString() << "\n";
        }
    }
    sstr << ")\n";

    return sstr.str();
}

} }
