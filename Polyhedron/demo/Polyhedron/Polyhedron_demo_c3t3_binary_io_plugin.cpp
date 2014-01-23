#include <CGAL/Mesh_3/io_signature.h>
#include "Scene_c3t3_item.h"

#include "Polyhedron_demo_io_plugin_interface.h"
#include <iostream>
#include <fstream>

class Polyhedron_demo_c3t3_binary_io_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_io_plugin_interface)

public:
  QString name() const { return "Polyhedron_demo_c3t3_binary_io_plugin"; }
  QString nameFilters() const { return "binary files (*.cgal)"; }

  bool canLoad() const;
  Scene_item* load(QFileInfo fileinfo);

  bool canSave(const Scene_item*);
  bool save(const Scene_item*, QFileInfo fileinfo);

private:
  bool try_load_other_binary_format(std::istream& in, C3t3& c3t3);
};


bool Polyhedron_demo_c3t3_binary_io_plugin::canLoad() const {
  return true;
}


Scene_item*
Polyhedron_demo_c3t3_binary_io_plugin::load(QFileInfo fileinfo) {

  // Open file
  std::ifstream in(fileinfo.filePath().toUtf8(),
                   std::ios_base::in|std::ios_base::binary);
  if(!in) {
    std::cerr << "Error! Cannot open file "
              << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    return NULL;
  }

  Scene_c3t3_item* item = new Scene_c3t3_item();
  item->setName(fileinfo.baseName());
  if(!item->load_binary(in))
  {
    item->c3t3().clear();
    in.seekg(0);
    if(try_load_other_binary_format(in, item->c3t3())) {
      item->changed();
      return item;
    }
    delete item;

    return NULL;
  }

  return item;
}

bool Polyhedron_demo_c3t3_binary_io_plugin::canSave(const Scene_item* item)
{
  // This plugin supports c3t3 items.
  return qobject_cast<const Scene_c3t3_item*>(item);
}

bool
Polyhedron_demo_c3t3_binary_io_plugin::
save(const Scene_item* item, QFileInfo fileinfo)
{
  if(!fileinfo.filePath().endsWith(".cgal"))
    return false;
  // This plugin supports polyhedrons and polygon soups
  const Scene_c3t3_item* c3t3_item =
    qobject_cast<const Scene_c3t3_item*>(item);

  if(!c3t3_item)
    return false;
  std::ofstream out(fileinfo.filePath().toUtf8(),
                   std::ios_base::out|std::ios_base::binary);

  return out && c3t3_item->save_binary(out);
}

struct Fake_mesh_domain {
  typedef CGAL::Tag_true Has_features;
  typedef int Subdomain_index;
  typedef std::pair<int, int> Surface_patch_index;
  typedef int Curve_segment_index;
  typedef int Corner_index;
  typedef boost::variant<Subdomain_index,Surface_patch_index> Index;
};

typedef c3t3_type_h::Geom_traits Fake_gt;
typedef CGAL::Mesh_vertex_base_3<Fake_gt, Fake_mesh_domain> Fake_vertex_base;
typedef CGAL::Compact_mesh_cell_base_3<Fake_gt, Fake_mesh_domain> Fake_cell_base;
typedef CGAL::Triangulation_data_structure_3<Fake_vertex_base,Fake_cell_base> Fake_tds;
typedef CGAL::Regular_triangulation_3<Fake_gt, Fake_tds> Fake_tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<
  Fake_tr,
  Fake_mesh_domain::Corner_index,
  Fake_mesh_domain::Curve_segment_index> Fake_c3t3;

#ifdef CGAL_MESH_3_IO_SIGNATURE_H
namespace CGAL {
template <>
struct Get_io_signature<Fake_mesh_domain::Surface_patch_index> {
  std::string operator()() const
  {
    return std::string("std::pair<i,i>");
  }
}; // end Get_io_signature<Fake_mesh_domain::Surface_patch_index>
} // end namespace CGAL
#endif

typedef Fake_mesh_domain::Surface_patch_index Patch_id;

namespace std {
  std::ostream& operator<<(std::ostream& out, const Patch_id& id) {
    return out << id.first << " " << id.second;
  }
  std::istream& operator>>(std::istream& in, Patch_id& id) {
    return in >> id.first >> id.second;
  }
} // end namespace std

namespace CGAL {
template <>
class Output_rep<Patch_id> {
  typedef Patch_id T;
  const T& t;
public:
  //! initialize with a const reference to \a t.
  Output_rep( const T& tt) : t(tt) {}
  //! perform the output, calls \c operator\<\< by default.
  std::ostream& operator()( std::ostream& out) const {
    if(is_ascii(out)) {
      out << t.first << " " << t.second;
    } else {
      CGAL::write(out, t.first);
      CGAL::write(out, t.second);
    }
    return out;
  }
};

template <>
class Input_rep<Patch_id> {
  typedef Patch_id T;
  T& t;
public:
  //! initialize with a const reference to \a t.
  Input_rep( T& tt) : t(tt) {}
  //! perform the output, calls \c operator\<\< by default.
  std::istream& operator()( std::istream& in) const {
    if(is_ascii(in)) {
      in >> t.first >> t.second;
    } else {
      CGAL::read(in, t.first);
      CGAL::read(in, t.second);
    }
    return in;
  }
};
} // end namespace CGAL

struct Update_vertex {
  typedef Fake_mesh_domain::Surface_patch_index Sp_index;
  template <typename V1, typename V2>
  bool operator()(const V1& v1, V2& v2) {
    v2.set_point(v1.point());
    v2.set_dimension(v1.in_dimension());
    v2.set_special(v1.is_special());
    switch(v1.in_dimension()) {
    case 0:
    case 1:
    case 3:
      v2.set_index(boost::get<int>(v1.index()));
      break;
    default: // case 2
      const typename V1::Index& index = v1.index();
      const Sp_index sp_index =
        boost::get<Sp_index>(index);
      const int i = (static_cast<int>(sp_index.first)
                     << (std::numeric_limits<int>::digits >> 1))
        + sp_index.second;
      v2.set_index(i);
    }
    return true;
  }
}; // end struct Update_vertex

struct Update_cell {
  typedef Fake_mesh_domain::Surface_patch_index Sp_index;
  template <typename C1, typename C2>
  bool operator()(const C1& c1, C2& c2) {
    c2.set_subdomain_index(c1.subdomain_index());
    for(int i = 0; i < 4; ++i) {
      const Sp_index sp_index = c1.surface_patch_index(i);
      const int new_index = (static_cast<int>(sp_index.first)
                             << (std::numeric_limits<int>::digits >> 1))
        + sp_index.second;
      c2.set_surface_patch_index(i, new_index);
    }
    return true;
  }
}; // end struct Update_cell

#include <CGAL/Triangulation_file_input.h>

bool
Polyhedron_demo_c3t3_binary_io_plugin::
try_load_other_binary_format(std::istream& is, C3t3& c3t3)
{
  CGAL::set_ascii_mode(is);
  std::string s;
  is >> s;
  if (s != "binary" ||
      !(is >> s) ||
      s != "CGAL" ||
      !(is >> s) ||
      s != "c3t3")
  {
    return false;
  }
  std::getline(is, s);
  if(s != "") {
    if(s != std::string(" ") + CGAL::Get_io_signature<Fake_c3t3>()()) {
      std::cerr << "Polyhedron_demo_c3t3_binary_io_plugin::try_load_other_binary_format:"
                << "\n  expected format: " << CGAL::Get_io_signature<Fake_c3t3>()()
                << "\n       got format:" << s << std::endl;
      return false;
    }
  }
  CGAL::set_binary_mode(is);
  Fake_c3t3 fake_c3t3;
  return CGAL::file_input<
    Fake_c3t3::Triangulation,
    C3t3::Triangulation,
    Update_vertex,
    Update_cell>(is, c3t3.triangulation());
}

#include <QtPlugin>
#include <QtPlugin>
Q_EXPORT_PLUGIN2(Polyhedron_demo_c3t3_binary_io_plugin, Polyhedron_demo_c3t3_binary_io_plugin)
#include "Polyhedron_demo_c3t3_binary_io_plugin.moc"
