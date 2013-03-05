#include "Scene_polyhedron_item.h"
#include "Scene_polygon_soup_item.h"
#include "Kernel_type.h"
#include "Polyhedron_type.h"

#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include "Polyhedron_demo_io_plugin_interface.h"
#include <fstream>

#include <boost/tuple/tuple.hpp>

#include <QColor>

typedef Kernel::Point_3 Point_3;
typedef std::vector<Point_3> Points_3;
typedef boost::tuple<int,int,int> Facet;
typedef std::vector<Facet> Surface;

 

template <class HDS>
class Build_from_gocad : public CGAL::Modifier_base<HDS> {
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

class Polyhedron_demo_gocad_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_io_plugin_interface)

public:
  QString nameFilters() const;
  QString name() const { return "gocad_plugin"; }
  bool canLoad() const;
  Scene_item* load(QFileInfo fileinfo);

  bool canSave(const Scene_item*);
  bool save(const Scene_item*, QFileInfo fileinfo);
};

QString Polyhedron_demo_gocad_plugin::nameFilters() const {
  return "GOCAD files (*.ts *.xyz)";
}

bool Polyhedron_demo_gocad_plugin::canLoad() const {
  return true;
}


Scene_item* 
Polyhedron_demo_gocad_plugin::load(QFileInfo fileinfo) {

  // Open file
  std::ifstream in(fileinfo.filePath().toUtf8());
  if(!in) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    return NULL;
  }
    
  // Try to read GOCAD file in a polyhedron
  Polyhedron P;
  //Scene_polyhedron_item* item = new Scene_polyhedron_item(P);
  //item->setName(fileinfo.baseName());
  
  Build_from_gocad<Polyhedron::HalfedgeDS> builder(in);
  //  item->polyhedron()->

  P.delegate(builder);

  if((!  in.good()) || (! P.is_valid())){
    // std::cerr << "Error: Invalid polyhedron" << std::endl;
    return 0;
  }   

  Scene_polyhedron_item* item = new Scene_polyhedron_item(P);
  if(builder.name.size() == 0){
    item->setName(fileinfo.baseName());
  } else {
    item->setName(builder.name.c_str());
  }
  QColor color(builder.color.c_str());
  if(color.isValid()) 
  {  
    item->setColor(color);
    item->changed();
  }
  

  return item;
}

bool Polyhedron_demo_gocad_plugin::canSave(const Scene_item*)
{
  return false;
}

bool Polyhedron_demo_gocad_plugin::save(const Scene_item*, QFileInfo)
{
  return false;
}

#include <QtPlugin>
Q_EXPORT_PLUGIN2(Polyhedron_demo_gocad_plugin, Polyhedron_demo_gocad_plugin)
#include "Polyhedron_demo_gocad_plugin.moc"
