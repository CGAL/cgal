#ifndef CGAL_GOCAD_IO_H
#define CGAL_GOCAD_IO_H

#include <deque>
#include <iostream>
#include <string>
#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <boost/tuple/tuple.hpp>
#include <iostream>


template <class HDS, class P>
class Build_from_gocad : public CGAL::Modifier_base<HDS> {


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


  Build_from_gocad(std::istream& is_)
    : is(is_), counter(0)
  {}

  void operator()( HDS& hds) {
    read(is, meshPoints, mesh);

    CGAL::Polyhedron_incremental_builder_3<Polyhedron::HalfedgeDS> B(hds, true);
    B.begin_surface( meshPoints.size(), mesh.size());
    typedef typename Points_3::size_type size_type;
    
    for(size_type i=0; i < meshPoints.size(); i++){
      B.add_vertex( meshPoints[i]);
    }
    for(size_type i=0; i < mesh.size(); i++){
      B.begin_facet(); 
      B.add_vertex_to_facet( mesh[i].template get<0>());
      B.add_vertex_to_facet( mesh[i].template get<1>());
      B.add_vertex_to_facet( mesh[i].template get<2>());
      B.end_facet();
    }
    if(B.error())
      {
        std::cerr << "An error occured while creating a Polyhedron" << std::endl;
        B.rollback();
      }
    
    B.end_surface();
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


template <typename Polyhedron>
bool
read_gocad(Polyhedron& polyhedron, std::istream& in, std::string& name, std::string& color)
{
  typedef typename Polyhedron::HalfedgeDS HDS;
  typedef typename Polyhedron::Point_3 Point_3;
  Build_from_gocad<HDS,Point_3> builder(in);

  polyhedron.delegate(builder);
  name=builder.name;
  color=builder.color;

  return in.good() && polyhedron.is_valid();
}

template <typename Polyhedron>
bool
write_gocad(Polyhedron& polyhedron, std::ostream& os, const std::string& name)
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
  {
    typename Polyhedron::Vertex_iterator it, end;
    it = polyhedron.vertices_begin();
    end = polyhedron.vertices_end();
    int i = 0;
    for(; it != end; ++it){
      it->id() = i;
      os << "VRTX " << i << " " << it->point() << "\n";
      ++i;
    }
  }

  {
    typename Polyhedron::Facet_iterator it, end;
    it = polyhedron.facets_begin();
    end = polyhedron.facets_end();
    for(; it != end; ++it){
      os << "TRGL " << it->halfedge()->prev()->vertex()->id() << " " 
         << it->halfedge()->vertex()->id() << " "
         << it->halfedge()->next()->vertex()->id()<< "\n";
    }
  }

  os << "END" << std::endl;

  return true;
}

#endif // CGAL_GOCAD_IO_H
