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
class Build_from_stl : public CGAL::Modifier_base<HDS> {
  std::istream& is;
  int counter;
  Points_3 meshPoints;
  Surface mesh;

public:

  std::string name, color;


  Build_from_stl(std::istream& is_)
    : is(is_), counter(0)
  {}

  void operator()( HDS& hds) {
    if(!read(is, meshPoints, mesh)) return;

    CGAL::Polyhedron_incremental_builder_3<Polyhedron::HalfedgeDS> B(hds);
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

bool
read(std::istream& input, Points_3& points, Surface& surface, int /*offset*/ = 0)
{
  std::string s, solid("solid"), facet("facet"), outer("outer"), loop("loop"), vertex("vertex"), endloop("endloop"), endsolid("endsolid");

  std::map<Point_3, int> vertex_index;
  int index = 0;
  int ijk[3];
  Point_3 p;

  input >> s;
  if(s == solid){
    std::getline(input, s);
  } else {
    std::cerr << "We expect keyword 'solid'" << std::endl;
    return false;
  }

  while(input >> s){
    if(s == endsolid){
      //std::cerr << "found endsolid" << std::endl;
    } else if(s == facet){
      //std::cerr << "found facet" << std::endl;
      std::getline(input, s); // ignore the normal
      input >> s;
      if(s != outer){
        std::cerr << "Expect 'outer' and got " << s << std::endl;
        return false;
      }
      input >> s;
      if(s != loop){
        std::cerr << "Expect 'loop' and got " << s << std::endl;
        return false;
     }
      int count = 0;
      do {
        input >> s;
        if(s == vertex){
          //      std::cerr << "found vertex" << std::endl;
          if(count < 3){
            input >> p;
            if(vertex_index.find(p) == vertex_index.end()){
              ijk[count] = index;
              vertex_index[p] = index++;
              points.push_back(p);
            } else {
              ijk[count] = vertex_index[p];
            }
            ++count;
          } else {
            std::cerr << "We can only read triangulated surfaces" << std::endl;
            return false;
          }
        }
      }while(s != endloop);

      surface.push_back(boost::make_tuple(ijk[0], ijk[1], ijk[2]));
    }
  }
  return true;
}

};

class Polyhedron_demo_stl_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_io_plugin_interface)

public:
  QString nameFilters() const;
  QString name() const { return "stl_plugin"; }
  bool canLoad() const;
  Scene_item* load(QFileInfo fileinfo);

  bool canSave(const Scene_item*);
  bool save(const Scene_item*, QFileInfo fileinfo);
};

QString Polyhedron_demo_stl_plugin::nameFilters() const {
  return "STL files (*.stl)";
}

bool Polyhedron_demo_stl_plugin::canLoad() const {
  return true;
}


Scene_item* 
Polyhedron_demo_stl_plugin::load(QFileInfo fileinfo) {

  // Open file
  std::ifstream in(fileinfo.filePath().toUtf8());
  if(!in) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    return NULL;
  }
    
  // Try to read STL file in a polyhedron
  Polyhedron P;
  //Scene_polyhedron_item* item = new Scene_polyhedron_item(P);
  //item->setName(fileinfo.baseName());
  
  Build_from_stl<Polyhedron::HalfedgeDS> builder(in);
  //  item->polyhedron()->

  P.delegate(builder);
 
  if(! P.is_valid()){
    std::cerr << "Error: Invalid polyhedron" << std::endl;
    return 0;
  }  

  if(P.empty()){
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

bool Polyhedron_demo_stl_plugin::canSave(const Scene_item*)
{
  return false;
}

bool Polyhedron_demo_stl_plugin::save(const Scene_item*, QFileInfo)
{
  return false;
}

#include <QtPlugin>
Q_EXPORT_PLUGIN2(Polyhedron_demo_stl_plugin, Polyhedron_demo_stl_plugin)
#include "Polyhedron_demo_stl_plugin.moc"
