#ifndef CGAL_GOCAD_IO_H
#define CGAL_GOCAD_IO_H

#include <deque>
#include <iostream>
#include <string>
#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <boost/tuple/tuple.hpp>
#include <iostream>

template <class Facegraph, class P>
class Build_from_gocad_BGL
{
  typedef P Point_3;
  typedef std::deque<Point_3> Points_3;
  typedef boost::tuple<int,int,int> Facet;
  typedef std::deque<Facet> Surface;

  std::istream& is;
  int counter;
  Points_3 meshPoints;
  Surface mesh;

public:

  std::string name, color;


  Build_from_gocad_BGL(std::istream& is_)
    : is(is_), counter(0)
  {}

  void operator()( Facegraph& graph) {
    read(is, meshPoints, mesh);

//    graph.begin_surface( meshPoints.size(), mesh.size());
    typedef typename Points_3::size_type size_type;

    for(size_type i=0; i < mesh.size(); i++){
      std::array<typename boost::graph_traits<Facegraph>::vertex_descriptor, 3> face;
      face[0] = add_vertex( meshPoints[mesh[i].template get<0>()],graph);
      face[1] = add_vertex( meshPoints[mesh[i].template get<1>()],graph);
      face[2] = add_vertex( meshPoints[mesh[i].template get<2>()],graph);

      CGAL::Euler::add_face(face, graph);
    }
  }

  void
  read(std::istream& input, Points_3& points, Surface& surface, int offset = 0)
  {
    char c;
    std::string s, tface("TFACE");
    int i,j,k;
    Point_3 p;
    bool vertices_read = false;
    while(input >> s){
      if(s == tface){
        break;
      }
      std::string::size_type idx;

      if((idx = s.find("name")) != std::string::npos){
        std::istringstream str(s.substr(idx+5));
        str >> name;
      }
      if((idx = s.find("color")) != std::string::npos){
        std::istringstream str(s.substr(idx+6));
        str >> color;
      }
    }
    std::getline(input, s);

    while(input.get(c)){
      if((c == 'V')||(c == 'P')){
        input >> s >> i >> p;
        if(! vertices_read){
          vertices_read = true;
          offset -= i; // Some files start with index 0 others with 1
        }

        points.push_back(p);

      } else if(vertices_read && (c == 'T')){
        input >> c >> c >> c >>  i >> j >> k;
        surface.push_back(boost::make_tuple(offset+i, offset+j, offset+k));
      } else if(c == 'E'){
        break;
      }
      std::getline(input, s);
    }
  }

};

template <typename FaceGraph>
bool
read_gocad(FaceGraph& polyhedron, std::istream& in, std::string& name, std::string& color)
{
  //typedef typename Polyhedron::HalfedgeDS HDS;
  typedef typename boost::property_traits<typename boost::property_map<FaceGraph, CGAL::vertex_point_t>::type>::value_type Point_3;

  Build_from_gocad_BGL<FaceGraph, Point_3> builder(in);
  builder(polyhedron);
  name=builder.name;
  color=builder.color;

  return in.good() && polyhedron.is_valid();
}

template <typename FaceGraph>
bool
write_gocad(FaceGraph& polyhedron, std::ostream& os, const std::string& name)
{
  os << "GOCAD TSurf 1\n"
    "HEADER {\n"
    "name:";
  os << name << std::endl;
  os << "*border:on\n"
    "*border*bstone:on\n"
    "}\n"
    "GOCAD_ORIGINAL_COORDINATE_SYSTEM\n"
    "NAME Default\n"
    "AXIS_NAME \"X\" \"Y\" \"Z\"\n"
    "AXIS_UNIT \"m\" \"m\" \"m\"\n"
    "ZPOSITIVE Elevation\n"
    "END_ORIGINAL_COORDINATE_SYSTEM\n"
    "TFACE\n";

  os.precision(16);
  typedef typename boost::property_map<FaceGraph, CGAL::vertex_point_t>::type VPMap;
  VPMap vpmap = get(CGAL::vertex_point, polyhedron);
  std::map<typename boost::graph_traits<FaceGraph>::vertex_descriptor, int> id_map;
  {
    typename boost::graph_traits<FaceGraph>::vertex_iterator it, end;
    it = vertices(polyhedron).begin();
    end = vertices(polyhedron).end();
    int i=0;
    for(; it != end; ++it){
      id_map[*it] = i;
      os << "VRTX " << i << " " << get(vpmap, *it) << "\n";
      ++i;
    }
  }

  {
    typename boost::graph_traits<FaceGraph>::face_iterator it, end;
    it = faces(polyhedron).begin();
    end = faces(polyhedron).end();
    for(; it != end; ++it){
      os << "TRGL " << id_map[target(prev(halfedge(*it, polyhedron), polyhedron), polyhedron)] << " "
         << id_map[target(halfedge(*it, polyhedron), polyhedron)] << " "
         << id_map[target(next(halfedge(*it, polyhedron), polyhedron), polyhedron)] << "\n";
    }
  }

  os << "END" << std::endl;

  return true;
}

#endif // CGAL_GOCAD_IO_H
