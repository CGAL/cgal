#include "config.h"

#include "Scene_polygon_soup.h"
#include <CGAL/IO/Polyhedron_iostream.h>

#include <QObject>
#include <QtDebug>

#include <set>
#include <stack>
#include <algorithm>
#include <boost/array.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/File_scanner_OFF.h>
#include <CGAL/IO/File_writer_OFF.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;

struct Polygon_soup {
  typedef std::vector<Point_3> Points;
  typedef std::vector<std::size_t> Polygon_3;
  typedef std::map<std::pair<std::size_t, std::size_t>, std::set<std::size_t> > Edges_map;
  typedef boost::array<std::size_t, 2> Edge;
  typedef std::vector<Polygon_3> Polygons;
  typedef std::vector<Edge> Edges;
  typedef Polygons::size_type size_type;
  Points points;
  Polygons polygons;
  Edges_map edges;
  Edges non_manifold_edges;
  bool display_non_manifold_edges;

  Polygon_soup* clone() const {
    Polygon_soup* result = new Polygon_soup();
    result->points = points;
    result->polygons = polygons;
    result->edges = edges;
    result->non_manifold_edges = non_manifold_edges;
    result->display_non_manifold_edges = display_non_manifold_edges;
    return result;
  }

  void clear() {
    points.clear();
    polygons.clear();
    edges.clear();
    non_manifold_edges.clear();
  }

  void fill_edges() {
    // Fill edges
    edges.clear();
    for(size_type i = 0; i < polygons.size(); ++i)
    {
      const size_type size = polygons[i].size();
      for(size_type j = 0; j < size; ++j) {
        const std::size_t& i0 = polygons[i][j];
        const std::size_t& i1 = polygons[i][ j+1 < size ? j+1: 0];
        edges[std::make_pair(i0, i1)].insert(i);
//         qDebug() << tr("edges[std::make_pair(%1, %2)].insert(%3). Size=%4")
//           .arg(i0).arg(i1).arg(i).arg(edges[std::make_pair(i0, i1)].size());
      }
    }

    // Fill non-manifold edges
    non_manifold_edges.clear();
    for(size_type i = 0; i < polygons.size(); ++i)
    {
      const size_type size = polygons[i].size();
      for(size_type j = 0; j < size; ++j) {
        const std::size_t& i0 = polygons[i][j];
        const std::size_t& i1 = polygons[i][ j+1 < size ? j+1: 0];
        if( (i0 < i1) && 
            (edges[std::make_pair(i0, i1)].size() +
             edges[std::make_pair(i1, i0)].size() > 2) )
        {
          Edge edge;
          edge[0] = i0;
          edge[1] = i1;
          non_manifold_edges.push_back(edge);
        }
      }
    }
  }

  void inverse_orientation(const std::size_t index) {
    std::reverse(polygons[index].begin(), polygons[index].end());
  }
};

Scene_polygon_soup::Scene_polygon_soup()
  : Scene_item(),
    soup(0),
    oriented(false)
{
    compile_shaders();
    are_buffers_initialized = false;
}

Scene_polygon_soup::~Scene_polygon_soup()
{
    for(int i=0; i<vboSize; i++)
        buffers[i].destroy();
    for(int i=0; i<vaoSize; i++)
        vao[i].destroy();

  delete soup;
}



/**************************************************
****************SHADER FUNCTIONS******************/

void Scene_polygon_soup::compile_shaders()
{

    for(int i=0; i< vboSize; i++)
        buffers[i].create();
    for(int i=0; i< vaoSize; i++)
        vao[i].create();

    //Vertex source code
    const char vertex_source[] =
    {
        "#version 120 \n"
        "attribute highp vec4 vertex;\n"
        "attribute highp vec3 normal;\n"

        "uniform highp mat4 mvp_matrix;\n"
        "uniform highp mat4 mv_matrix; \n"
        "varying highp vec4 fP; \n"
        "varying highp vec3 fN; \n"
        "void main(void)\n"
        "{\n"
        "   fP = mv_matrix * vertex; \n"
        "   fN = mat3(mv_matrix)* normal; \n"
        "   gl_Position = mvp_matrix * vertex; \n"
        "}"
    };
    //Fragment source code
    const char fragment_source[] =
    {
        "#version 120 \n"
        "varying highp vec4 fP; \n"
        "varying highp vec3 fN; \n"
        "uniform vec4 color; \n"
        "uniform bool is_two_side; \n"
        "uniform highp vec4 light_pos;  \n"
        "uniform highp vec4 light_diff; \n"
        "uniform highp vec4 light_spec; \n"
        "uniform highp vec4 light_amb;  \n"
        "uniform float spec_power ; \n"

        "void main(void) { \n"

        "   vec3 L = light_pos.xyz - fP.xyz; \n"
        "   vec3 V = -fP.xyz; \n"

        "   vec3 N; \n"
        "   if(fN == vec3(0.0,0.0,0.0)) \n"
        "       N = vec3(0.0,0.0,0.0); \n"
        "   else \n"
        "       N = normalize(fN); \n"
        "   L = normalize(L); \n"
        "   V = normalize(V); \n"

        "   vec3 R = reflect(-L, N); \n"
        "   vec4 diffuse; \n"
        "   if(!is_two_side) \n"
        "       diffuse = max(dot(N,L),0) * light_diff*color; \n"
        "   else \n"
        "       diffuse = max(abs(dot(N,L)),0) * light_diff*color; \n"
        "   vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec; \n"

         "gl_FragColor = color*light_amb + diffuse + specular; \n"
        "} \n"
        "\n"
    };
    QOpenGLShader *vertex_shader = new QOpenGLShader(QOpenGLShader::Vertex);
    if(!vertex_shader->compileSourceCode(vertex_source))
    {
        std::cerr<<"Compiling vertex source FAILED"<<std::endl;
    }

    QOpenGLShader *fragment_shader= new QOpenGLShader(QOpenGLShader::Fragment);
    if(!fragment_shader->compileSourceCode(fragment_source))
    {
        std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
    }

    if(!rendering_program.addShader(vertex_shader))
    {
        std::cerr<<"adding vertex shader FAILED"<<std::endl;
    }
    if(!rendering_program.addShader(fragment_shader))
    {
        std::cerr<<"adding fragment shader FAILED"<<std::endl;
    }
    if(!rendering_program.link())
    {
        std::cerr<<"linking Program FAILED"<<std::endl;
    }
}

void Scene_polygon_soup::compute_elements()
{
    normal.resize(0);
    v_poly.resize(0);
    v_edge.resize(0);
    vertex_nm.resize(0);

    typedef Polygon_soup::Polygons::const_iterator Polygons_iterator;
    typedef Polygon_soup::Polygons::size_type size_type;
    for(Polygons_iterator it = soup->polygons.begin();
        it != soup->polygons.end(); ++it)
    {
      const Point_3& pa = soup->points[it->at(0)];
      const Point_3& pb = soup->points[it->at(1)];
      const Point_3& pc = soup->points[it->at(2)];

      Kernel::Vector_3 n = CGAL::cross_product(pb-pa, pc -pa);
      n = n / std::sqrt(n * n);

      for(size_type i = 0; i < it->size(); ++i) {
        const Point_3& p = soup->points[it->at(i)];
        normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());
        v_poly.push_back(p.x()); v_poly.push_back(p.y()); v_poly.push_back(p.z());
      }
    }

    //edges
    for(size_type i = 0; i < soup->polygons.size(); ++i)
    {
      const size_type size = soup->polygons[i].size();
      for(size_type j = 0; j < size; ++j) {
        const std::size_t& i0 = soup->polygons[i][j];
        const std::size_t& i1 = soup->polygons[i][ j+1 < size ? j+1: 0];
        Polygon_soup::Edge edge;
        edge[0] = i0;
        edge[1] = i1;
        const Point_3& a = soup->points[edge[0]];
        const Point_3& b = soup->points[edge[1]];
        v_edge.push_back(a.x()); v_edge.push_back(a.y());  v_edge.push_back(a.z());
        v_edge.push_back(b.x()); v_edge.push_back(b.y()); v_edge.push_back(b.z());
      }
    }

    //non_manifold_edges

    for(Polygon_soup::size_type
        i = 0,
        end = soup->non_manifold_edges.size();
        i < end; ++i)
    {
        const Polygon_soup::Edge& edge = soup->non_manifold_edges[i];
        const Point_3& a = soup->points[edge[0]];
        const Point_3& b = soup->points[edge[1]];

        vertex_nm.push_back(a.x()); vertex_nm.push_back(a.y());  vertex_nm.push_back(a.z());
        vertex_nm.push_back(b.x()); vertex_nm.push_back(b.y()); vertex_nm.push_back(b.z());
    }
}

void Scene_polygon_soup::initialize_buffers() const
{
    rendering_program.bind();

        vao[0].bind();
        buffers[0].bind();
        buffers[0].allocate(v_poly.data(), static_cast<int>(v_poly.size()*sizeof(float)));
        poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(poly_vertexLocation[0]);
        rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
        buffers[0].release();

        buffers[1].bind();
        buffers[1].allocate(normal.data(), static_cast<int>(normal.size()*sizeof(float)));
        normalsLocation[0] = rendering_program.attributeLocation("normal");
        rendering_program.enableAttributeArray(normalsLocation[0]);
        rendering_program.setAttributeBuffer(normalsLocation[0],GL_FLOAT,0,3);
        buffers[1].release();

        vao[0].release();

        vao[1].bind();
        buffers[2].bind();
        buffers[2].allocate(v_edge.data(), static_cast<int>(v_edge.size()*sizeof(float)));
        poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(poly_vertexLocation[0]);
        rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
        buffers[2].release();
        std::vector<float> empty_array;
        empty_array.resize(0);
        for(std::size_t i=0; i<v_edge.size(); i++)
            empty_array.push_back(0.0);
        buffers[3].bind();
        buffers[3].allocate(empty_array.data(), static_cast<int>(empty_array.size()*sizeof(float)));
        normalsLocation[0] = rendering_program.attributeLocation("normal");
        rendering_program.enableAttributeArray(normalsLocation[0]);
        rendering_program.setAttributeBuffer(normalsLocation[0],GL_FLOAT,0,3);
        buffers[3].release();

        vao[1].release();

        vao[2].bind();
        buffers[4].bind();
        buffers[4].allocate(v_edge.data(), static_cast<int>(v_edge.size()*sizeof(float)));
        poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(poly_vertexLocation[0]);
        rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
        buffers[4].release();
        empty_array.resize(0);
        for(std::size_t i=0; i<v_edge.size(); i++)
            empty_array.push_back(0.0);
        buffers[5].bind();
        buffers[5].allocate(empty_array.data(), static_cast<int>(empty_array.size()*sizeof(float)));
        normalsLocation[0] = rendering_program.attributeLocation("normal");
        rendering_program.enableAttributeArray(normalsLocation[0]);
        rendering_program.setAttributeBuffer(normalsLocation[0],GL_FLOAT,0,3);
        buffers[5].release();

        vao[2].release();

    rendering_program.release();
    are_buffers_initialized = true;

}

void Scene_polygon_soup::attrib_buffers(Viewer* viewer) const
{
    QMatrix4x4 mvpMatrix;
    QMatrix4x4 mvMatrix;
    double mat[16];
    viewer->camera()->getModelViewProjectionMatrix(mat);
    for(int i=0; i < 16; i++)
    {
        mvpMatrix.data()[i] = (float)mat[i];
    }
    viewer->camera()->getModelViewMatrix(mat);
    for(int i=0; i < 16; i++)
    {
        mvMatrix.data()[i] = (float)mat[i];
    }
    QVector4D	position(0.0f,0.0f,1.0f,1.0f );
    GLboolean isTwoSide;
    viewer->glGetBooleanv(GL_LIGHT_MODEL_TWO_SIDE,&isTwoSide);
    // define material
     QVector4D	ambient;
     QVector4D	diffuse;
     QVector4D	specular;
     GLfloat      shininess ;
    // Ambient
    ambient[0] = 0.29225f;
    ambient[1] = 0.29225f;
    ambient[2] = 0.29225f;
    ambient[3] = 1.0f;
    // Diffuse
    diffuse[0] = 0.50754f;
    diffuse[1] = 0.50754f;
    diffuse[2] = 0.50754f;
    diffuse[3] = 1.0f;
    // Specular
    specular[0] = 0.0f;
    specular[1] = 0.0f;
    specular[2] = 0.0f;
    specular[3] = 0.0f;
    // Shininess
    shininess = 51.2f;


    rendering_program.bind();
    colorLocation[0] = rendering_program.uniformLocation("color");
    twosideLocation = rendering_program.uniformLocation("is_two_side");
    mvpLocation[0] = rendering_program.uniformLocation("mvp_matrix");
    mvLocation[0] = rendering_program.uniformLocation("mv_matrix");
    lightLocation[0] = rendering_program.uniformLocation("light_pos");
    lightLocation[1] = rendering_program.uniformLocation("light_diff");
    lightLocation[2] = rendering_program.uniformLocation("light_spec");
    lightLocation[3] = rendering_program.uniformLocation("light_amb");
    lightLocation[4] = rendering_program.uniformLocation("spec_power");

    rendering_program.setUniformValue(lightLocation[0], position);
    rendering_program.setUniformValue(twosideLocation, isTwoSide);
    rendering_program.setUniformValue(mvpLocation[0], mvpMatrix);
    rendering_program.setUniformValue(mvLocation[0], mvMatrix);
    rendering_program.setUniformValue(lightLocation[1], diffuse);
    rendering_program.setUniformValue(lightLocation[2], specular);
    rendering_program.setUniformValue(lightLocation[3], ambient);
    rendering_program.setUniformValue(lightLocation[4], shininess);

    rendering_program.release();
}

void
Scene_polygon_soup::draw(Viewer* viewer) const  {
    if(!are_buffers_initialized)
        initialize_buffers();

    QColor color;
    vao[0].bind();
    float current_color[4];
    viewer->glGetFloatv(GL_CURRENT_COLOR, current_color);
    color.setRgbF(current_color[0],current_color[1],current_color[2],current_color[3]);
    attrib_buffers(viewer);
    rendering_program.bind();
    rendering_program.setUniformValue(colorLocation[0], color);
    viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(v_poly.size()/3));
    rendering_program.release();
    vao[0].release();
}

void
Scene_polygon_soup::draw_edges(Viewer* viewer) const  {
    if(!are_buffers_initialized)
        initialize_buffers();
    QColor color;
    vao[1].bind();
    float current_color[4];
    glGetFloatv(GL_CURRENT_COLOR, current_color);
    color.setRgbF(current_color[0],current_color[1],current_color[2],current_color[3]);
    attrib_buffers(viewer);
    rendering_program.bind();
    rendering_program.setUniformValue(colorLocation[0], color);
    glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(v_edge.size()/3));
    rendering_program.release();
    vao[1].release();

    if(soup->display_non_manifold_edges)
    {
        vao[2].bind();
        color.setRgbF(1.0,0.0,0.0);
        attrib_buffers(viewer);
        rendering_program.bind();
        rendering_program.setUniformValue(colorLocation[0], color);
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(v_edge.size()/3));
        rendering_program.release();
        vao[2].release();
    }
}

void
Scene_polygon_soup::draw_points(Viewer* viewer) const  {

    if(!are_buffers_initialized)
        initialize_buffers();
    QColor color;
    vao[0].bind();
    float current_color[4];
    viewer->glGetFloatv(GL_CURRENT_COLOR, current_color);
    color.setRgbF(current_color[0],current_color[1],current_color[2],current_color[3]);
    attrib_buffers(viewer);
    rendering_program.bind();
    rendering_program.setUniformValue(colorLocation[0], color);
    viewer->glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(v_poly.size()/3));
    rendering_program.release();
    vao[0].release();

}


Scene_polygon_soup* 
Scene_polygon_soup::clone() const {
  Scene_polygon_soup* new_soup = new Scene_polygon_soup();
  new_soup->soup = soup->clone();
  new_soup->oriented = oriented;
  return new_soup;
}

bool
Scene_polygon_soup::load(std::istream& in)
{
  if(!soup)
    soup = new Polygon_soup;
  CGAL::File_scanner_OFF scanner(in);
  soup->clear();
  soup->points.resize(scanner.size_of_vertices());
  soup->polygons.resize(scanner.size_of_facets());
  for (std::size_t i = 0; i < scanner.size_of_vertices(); ++i) {
    double x, y, z, w;
    scanner.scan_vertex( x, y, z, w);
    soup->points[i] = Point_3(x, y, z, w);
    scanner.skip_to_next_vertex( i);
  }
  if(!in)
    return false;

  for (std::size_t i = 0; i < scanner.size_of_facets(); ++i) {
    std::size_t no;
    scanner.scan_facet( no, i);
    soup->polygons[i].resize(no);
    for(std::size_t j = 0; j < no; ++j) {
      std::size_t id;
      scanner.scan_facet_vertex_index(id, i);
      if(id < scanner.size_of_vertices())
      {
        soup->polygons[i][j] = id;
      }
      else
        return false;
    }
  }
  soup->fill_edges();
  oriented = false;
  changed();
  return ! in.fail();
}

void
Scene_polygon_soup::setDisplayNonManifoldEdges(const bool b)
{
  soup->display_non_manifold_edges = b;
  changed();
}

bool
Scene_polygon_soup::displayNonManifoldEdges() const {
  return soup->display_non_manifold_edges;
}

void Scene_polygon_soup::shuffle_orientations()
{
  for(Polygon_soup::size_type i = 0, end = soup->polygons.size();
      i < end; ++i)
  {
    if(std::rand() % 2 == 0) soup->inverse_orientation(i);
  }
  soup->fill_edges();
  changed();
}

void Scene_polygon_soup::inside_out()
{
  for(Polygon_soup::size_type i = 0, end = soup->polygons.size();
      i < end; ++i)
  {
    soup->inverse_orientation(i);
  }
  soup->fill_edges();
  changed();
}

bool 
Scene_polygon_soup::orient()
{
  typedef Polygon_soup::Polygons::size_type size_type;
  typedef Polygon_soup::Edges_map Edges;

  if(isEmpty() || this->oriented)
    return true; // nothing to do

  Polygon_soup::Polygons& polygons = soup->polygons;
  Polygon_soup::Edges_map& edges = soup->edges;

  std::vector<bool> oriented;
  std::stack<std::size_t> stack;
  using std::make_pair;

  // no polygon is oriented
  oriented.resize(polygons.size());

  size_type polygon_index = 0;
  bool success = true;

  while (polygon_index != polygons.size()) 
  {
    while ( polygon_index != polygons.size() && oriented[polygon_index] ) {
      ++polygon_index;
    }
    if(polygon_index == polygons.size()) break;

//     qDebug() << tr("Seed %1...\n").arg(polygon_index);
    oriented[polygon_index] = true;
    stack.push(polygon_index);
    while(! stack.empty() )
    {
      const size_type to_be_oriented_index = stack.top();
//       qDebug() << tr("polygon #%1").arg(to_be_oriented_index);
      stack.pop();
      const size_type size = polygons[to_be_oriented_index].size();
      for(size_type ih = 0 ; ih < size ; ++ih) {
        size_type ihp1 = ih+1;
        if(ihp1>=size) ihp1 = 0;
        const std::size_t& i1 = polygons[to_be_oriented_index][ih];
        const std::size_t& i2 = polygons[to_be_oriented_index][ihp1];

//         qDebug() << tr("edge %3-%4 (%1,%2)").arg(i1).arg(i2).arg(ih).arg(ihp1);
        // edge (i1,i2)
        Edges::iterator it_same_orient = edges.find(make_pair(i1, i2));
        // edges (i2,i1)
        Edges::iterator it_other_orient = edges.find(make_pair(i2, i1));

        CGAL_assertion(it_same_orient != edges.end());
        if(it_same_orient->second.size() > 1) {
          if((it_other_orient != edges.end() && it_other_orient->second.size() > 0) ||
             it_same_orient->second.size() > 2) {
            // three polygons at the edge
//             qDebug() << "three polygons at the edge";
            success = false; // non-orientable
          }
          {
            // one neighbor polyhedron, opposite orientation
            size_type index = *(it_same_orient->second.begin());
            if(index == to_be_oriented_index)
              index = *(++it_same_orient->second.begin());
            if(oriented[index]) {
//               qDebug() << tr("neighbor polygon #%1 is already oriented, but in opposite orientation").arg(index);
              success = false; // non-orientable
              continue; // next edge
            }
   
            // reverse the orientation
            const size_type size = polygons[index].size();
            for(size_type j = 0; j < size; ++j) {
              const std::size_t& i0 = polygons[index][j];
              const std::size_t& i1 = polygons[index][ j+1 < size ? j+1: 0];
              CGAL_assertion_code(const bool r = )
                edges[std::make_pair(i0, i1)].erase(index);
              CGAL_assertion(r);
            }
            soup->inverse_orientation(index);
            for(size_type j = 0; j < size; ++j) {
              const std::size_t& i0 = polygons[index][j];
              const std::size_t& i1 = polygons[index][ j+1 < size ? j+1: 0];
              edges[std::make_pair(i0, i1)].insert(index);
            }
//             qDebug() << tr("inverse the orientation of polygon #%1\n").arg(index);
            oriented[index] = true;
            stack.push(index);
          }
        }
        else if(it_other_orient != edges.end() && it_other_orient->second.size() == 1) {
          // one polygon, same orientation
          const size_type index = *(it_other_orient->second.begin());
          if(oriented[index])
            continue;
          oriented[index] = true;
//           qDebug() << tr("keep the orientation of polygon #%1\n").arg(index);
          stack.push(index);
        }
        else {
//           qDebug() << "else" << it_same_orient->second.size() << it_other_orient->second.size();
          success = false; // non-orientable
        }
      } // end for on all edges of one 
    } // end while loop on the polygons of the connected component
  } // end while loop on all non-oriented polygons remaining 
  return success;
}


bool 
Scene_polygon_soup::save(std::ostream& out) const
{
  typedef Polygon_soup::size_type size_type;
  CGAL::File_writer_OFF writer(true); // verbose
  writer.write_header(out,
                      soup->points.size(),
                      0,
                      soup->polygons.size());
  for(size_type i = 0, end = soup->points.size();
      i < end; ++i)
  {
    const Point_3& p = soup->points[i];
    writer.write_vertex( p.x(), p.y(), p.z() );
  }
  writer.write_facet_header();
  for(size_type i = 0, end = soup->polygons.size();
      i < end; ++i)
  {
    const Polygon_soup::Polygon_3& polygon = soup->polygons[i]; 
    const size_type size = polygon.size();
    writer.write_facet_begin(size);
    for(size_type j = 0; j < size; ++j) {
      writer.write_facet_vertex_index(polygon[j]);
    }
    writer.write_facet_end();
  }
  writer.write_footer();

  return ! out.fail();
}

QString 
Scene_polygon_soup::toolTip() const
{
  if(!soup)
    return QString();

  return QObject::tr("<p><b>%1</b> (mode: %5, color: %6)<br />"
                     "<i>Polygons soup</i></p>"
                     "<p>Number of vertices: %2<br />"
                     "Number of polygons: %3</p>")
    .arg(this->name())
    .arg(soup->points.size())
    .arg(soup->polygons.size())
    .arg(this->renderingModeName())
    .arg(this->color().name());
}

bool
Scene_polygon_soup::isEmpty() const {
  return (soup == 0 || soup->points.empty());
}

Scene_polygon_soup::Bbox
Scene_polygon_soup::bbox() const {
  const Point_3& p = *(soup->points.begin());
  CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
  for(Polygon_soup::Points::const_iterator it = soup->points.begin();
      it != soup->points.end();
      ++it) {
    bbox = bbox + it->bbox();
  }
  return Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
              bbox.xmax(),bbox.ymax(),bbox.zmax());
}

void Scene_polygon_soup::changed()
{
   // Scene_item::changed();
    compute_elements();
    are_buffers_initialized = false;
}
