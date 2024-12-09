#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/SMDS_3/io_signature.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/IO/File_binary_mesh_3.h>

#include <variant>

#include <string>

#include <CGAL/disable_warnings.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

//
// Define two fake Mesh_domain, with index types
//

// One with all types equal to int.
struct MD_homogeneous_types {
  typedef CGAL::Tag_false Has_features;
  typedef int Subdomain_index;
  typedef int Surface_patch_index;
  typedef int Curve_index;
  typedef int Corner_index;
  typedef int Index;

  static Subdomain_index get_sub_domain_index_1() { return 1; }
  static Subdomain_index get_sub_domain_index_2() { return 2; }
  static Surface_patch_index get_surface_patch_index_1() { return 3; }
  static Surface_patch_index get_surface_patch_index_2() { return 4; }
  static Curve_index get_curve_index_1() { return 5; }
  static Curve_index get_curve_index_2() { return 6; }
  static Corner_index get_corner_index_1() { return 7; }
  static Corner_index get_corner_index_2() { return 8; }

  static std::string reference_format_string()
  {
    return "Triangulation_3(Weighted_point<Point_3>,Vb(Tvb_3+i+i),Cb(i+RTcb_3+(i)[4]))";
  }
};

// One with heterogeneous index types
struct MD_heterogeneous_types {
  typedef CGAL::Tag_true Has_features;
  enum Subdomain_index_enum { Z = 0, A, B, C, D};
  typedef Subdomain_index_enum Subdomain_index;
  typedef std::pair<int, int> Surface_patch_index;
  typedef int Curve_index;
  typedef double Corner_index;
  typedef std::variant<Subdomain_index,
                         Surface_patch_index,
                         Curve_index,
                         Corner_index> Index;
  static Subdomain_index get_sub_domain_index_1() { return A; }
  static Subdomain_index get_sub_domain_index_2() { return B; }
  static Surface_patch_index get_surface_patch_index_1() { return std::make_pair(1, 2); }
  static Surface_patch_index get_surface_patch_index_2() { return std::make_pair(3, 4); }
  static Curve_index get_curve_index_1() { return 5; }
  static Curve_index get_curve_index_2() { return 6; }
  static Corner_index get_corner_index_1() { return 7.; }
  static Corner_index get_corner_index_2() { return 8.; }

  static std::string reference_format_string()
  {
    return "Triangulation_3(Weighted_point<Point_3>,Vb(Tvb_3+i+std::variant<enum,std::pair<i,i>,i,d>),Cb(enum+RTcb_3+(std::pair<i,i>)[4]))";
  }
};

//
// Define I/O for MD_heterogeneous_types::Subdomain_index (then enum)
//

// First, new technique: specialization of Output_rep and Input_rep (from
// CGAL/IO/io.h). That works for CGAL::read and CGAL::write, and ease the
// overload of << and >>
namespace CGAL {

template <>
class Output_rep<MD_heterogeneous_types::Subdomain_index> {
  typedef MD_heterogeneous_types::Subdomain_index T;
  const T& t;
public:
  //! initialize with a const reference to \a t.
  Output_rep( const T& tt) : t(tt) {}
  //! perform the output, calls \c operator\<\< by default.
  std::ostream& operator()( std::ostream& out) const {
    if(IO::is_ascii(out)) {
      out << static_cast<int>(t);
    } else {
      CGAL::write(out, static_cast<int>(t));
    }
    return out;
  }
};

template <>
class Input_rep<MD_heterogeneous_types::Subdomain_index> {
  typedef MD_heterogeneous_types::Subdomain_index T;
  T& t;
public:
  //! initialize with a const reference to \a t.
  Input_rep( T& tt) : t(tt) {}
  //! perform the output, calls \c operator\<\< by default.
  std::istream& operator()( std::istream& in) const {
    int i;
    if(IO::is_ascii(in)) {
      in >> i;
    } else {
      CGAL::read(in, i);
    }
    t = static_cast<T>(i);
    return in;
  }
};
} // end namespace CGAL

// Second: operator<< and >>
namespace std {
std::ostream& operator<<(std::ostream& out,
                         MD_heterogeneous_types::Subdomain_index index) {
  return out << CGAL::IO::oformat(index);
}
std::istream& operator>>(std::istream& in,
                         MD_heterogeneous_types::Subdomain_index& index) {
  return in >> CGAL::IO::iformat(index);
}
} // end namespace std

// Second: the specialization of CGAL::Get_io_signature (in
// <CGAL/Mesh_3/io_signature.h>).
namespace CGAL {
template <>
struct Get_io_signature<MD_heterogeneous_types::Subdomain_index> {
  std::string operator()() const {
    return "enum";
  }
};

//
// Define I/O for MD_heterogeneous_types::Surface_patch_index (pair of int)
//

// First, new technique: specialization of Output_rep and Input_rep (from
// CGAL/IO/io.h). That works for CGAL::read and CGAL::write, and ease the
// overload of << and >>
template <>
class Output_rep<MD_heterogeneous_types::Surface_patch_index> {
  typedef MD_heterogeneous_types::Surface_patch_index T;
  const T& t;
public:
  //! initialize with a const reference to \a t.
  Output_rep( const T& tt) : t(tt) {}
  //! perform the output, calls \c operator\<\< by default.
  std::ostream& operator()( std::ostream& out) const {
    if(IO::is_ascii(out)) {
      out << t.first << " " << t.second;
    } else {
      CGAL::write(out, t.first);
      CGAL::write(out, t.second);
    }
    return out;
  }
};

template <>
class Input_rep<MD_heterogeneous_types::Surface_patch_index> {
  typedef MD_heterogeneous_types::Surface_patch_index T;
  T& t;
public:
  //! initialize with a const reference to \a t.
  Input_rep( T& tt) : t(tt) {}
  //! perform the output, calls \c operator\<\< by default.
  std::istream& operator()( std::istream& in) const {
    if(IO::is_ascii(in)) {
      in >> t.first >> t.second;
    } else {
      CGAL::read(in, t.first);
      CGAL::read(in, t.second);
    }
    return in;
  }
};
} // end namespace CGAL

// Second: operator<< and >>
namespace std {
std::ostream& operator<<(std::ostream& out,
                         MD_heterogeneous_types::Surface_patch_index index) {
  return out << CGAL::IO::oformat(index);
}
std::istream& operator>>(std::istream& in,
                         MD_heterogeneous_types::Surface_patch_index& index) {
  return in >> CGAL::IO::iformat(index);
}
} // end namespace std


//
// Class to test I/O of Mesh_complex_3_in_triangulation_3<Mesh_domain,K>
//
template <typename Mesh_domain>
struct Test_c3t3_io {
  typedef typename Mesh_domain::Subdomain_index Subdomain_index;
  typedef typename Mesh_domain::Surface_patch_index Surface_patch_index;
  typedef typename Mesh_domain::Curve_index Curve_index;
  typedef typename Mesh_domain::Corner_index Corner_index;
  typedef typename Mesh_domain::Index Index;

  typedef typename CGAL::Mesh_triangulation_3<Mesh_domain, K>::Type Tr;
  typedef CGAL::Mesh_complex_3_in_triangulation_3<
    Tr,
    Corner_index,
    Curve_index
    > C3t3;

  typedef typename Tr::Point Point;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Edge Edge;
  typedef typename Tr::Facet Facet;
  typedef typename Tr::Cell_handle Cell_handle;

  static bool check_equality(const C3t3& c1, const C3t3& c2) {
    const Tr& t1 = c1.triangulation();
    const Tr& t2 = c2.triangulation();
    if(t1.number_of_vertices() != t2.number_of_vertices()) {
      assert(false);
      return false;
    }
    if(t1.number_of_vertices() != t2.number_of_vertices()) {
      assert(false);
      return false;
    }
    if(c1.number_of_cells_in_complex() != c2.number_of_cells_in_complex()) {
      assert(false);
      return false;
    }
    if(c1.number_of_facets_in_complex() != c2.number_of_facets_in_complex()) {
      assert(false);
      return false;
    }
    for(typename Tr::Finite_vertices_iterator
          vit1 = t1.finite_vertices_begin(),
          vit2 = t2.finite_vertices_begin(),
          end1 = t1.finite_vertices_end();
        vit1 != end1; ++vit1, ++vit2)
    {
      if(!( vit1->in_dimension() == vit2->in_dimension() &&
            vit1->index() == vit2->index()) )
      {
        std::cerr << "Error: vertices are different:\n";
        std::cerr << *vit1 << "\n"
                  << *vit2 << std::endl;
        assert(false);
        return false;
      }
    }
#if 0
    // Note: The Triangulation_3 facets iterator order changes after a reload
    for(typename Tr::Finite_facets_iterator
          fit1 = t1.finite_facets_begin(),
          fit2 = t2.finite_facets_begin(),
          end1 = t1.finite_facets_end();
        fit1 != end1; ++fit1, ++fit2)
    {
      typename Tr::Cell_handle c1 = fit1->first;
      typename Tr::Cell_handle c2 = fit2->first;
      int i1 = fit1->second;
      int i2 = fit2->second;
      // CJ: this may cause an assertion because the 2 C3T3s may have
      // facets in different orders.
      // This is because the Finite_facets_iterator compares the
      // addresses of cells to ensure parsing unique facets. The
      // facets are stored in a Compact_container, which doesn't
      /// guarantee any order in memory (and even in the container itself)
      assert(i1 == i2);
      if( c1->surface_patch_index(i1) !=
          c2->surface_patch_index(i2) )
      {
        std::cerr << "Error: facets #" << i1
                  << "of the following cells are different:\n";
        std::cerr << *c1 << "\n"
                  << *c2 << std::endl;
        assert(false);
        return false;
      }
    }
#endif // not WIN32
    for(typename Tr::Finite_cells_iterator
          cit1 = t1.finite_cells_begin(),
          cit2 = t2.finite_cells_begin(),
          end1 = t1.finite_cells_end();
        cit1 != end1; ++cit1, ++cit2)
    {
      if(cit1->subdomain_index() != cit2->subdomain_index() )
      {
        std::cerr << "Error: cells are different:\n";
        std::cerr << *cit1 << "\n"
                  << *cit2 << std::endl;
        assert(false);
        return false;
      }
      for(int i = 0; i < 4; ++i) {
        if( cit1->surface_patch_index(i) !=
            cit2->surface_patch_index(i) )
        {
          std::cerr << "Error: cells are different:\n";
          std::cerr << *cit1 << "\n"
                    << *cit2 << "\n"
                    << "surface_patch_index(" << i << "):\n"
                    << cit1->surface_patch_index(i) << "\n"
                    << cit2->surface_patch_index(i) << "\n";
          assert(false);
          return false;
        }
      }
    }
    return true;
  }

  static bool test_io(const C3t3& c3t3, const char* prefix, bool binary) {
    std::string filename(prefix);
    std::cout << "test_io";
    std::ios_base::openmode mode = std::ios_base::in | std::ios_base::out;
    if(binary) {
      filename += ".binary";
      std::cout << " in binary mode";
      mode |= std::ios_base::binary;
    }
    filename += ".cgal";
    std::cout << std::endl;
    std::cout << "IO format: " << CGAL::Get_io_signature<C3t3>()() << std::endl;
    std::stringstream stream(mode);
    if(binary) {
      CGAL::IO::set_binary_mode(stream);
    }
    stream << c3t3;
    if(!binary) {
      std::cout << "Content of the stream:\n"
                << "*****begin*****\n"
                << stream.str()
                << "\n******end******" << std::endl;
    }
    stream.seekg(0);
    assert(stream);
    C3t3 c3t3_bis;
    stream >> c3t3_bis;
    std::cout << "Content of c3t3_bis:\n"
              << "*****begin*****\n"
              << c3t3_bis
              << "\n******end******" << std::endl;
    assert(stream);

    {
      std::ostringstream ss_c3t3, ss_c3t3_bis;
      ss_c3t3 << c3t3;
      ss_c3t3_bis << c3t3_bis;
      assert(ss_c3t3);
      assert(ss_c3t3_bis);
      assert(ss_c3t3.str() == ss_c3t3_bis.str());
      if(!check_equality(c3t3, c3t3_bis)) return false;
    }

    c3t3_bis.clear();
    {
      std::stringstream ss(std::ios_base::out |
                           ( binary ?
                             (std::ios_base::in | std::ios_base::binary)
                             : std::ios_base::in));
      CGAL::IO::save_binary_file(ss, c3t3, binary);
      assert(ss);
      CGAL::IO::load_binary_file(ss, c3t3_bis);
      assert(ss);
    }
    if(!check_equality(c3t3, c3t3_bis)) return false;

#ifndef CGAL_LITTLE_ENDIAN
    // skip binary I/O with the existing file for big endian
    return true;
#endif

    if (sizeof(void*) != 8) {
      // skip 32bits architectures as well
      return true;
    }

    {
      std::string filename(prefix);
      filename += "_new";
      if(binary) filename += ".binary";
      filename += ".cgal";
      std::ostringstream output(binary
                           ? (std::ios_base::out | std::ios_base::binary)
                           : std::ios_base::out);
      CGAL::IO::save_binary_file(output, c3t3_bis, binary);
    }

    c3t3_bis.clear();
    {
      std::ifstream input(filename.c_str(),
                          binary ? (std::ios_base::in | std::ios_base::binary)
                          : std::ios_base::in);
      assert(input);
      CGAL::IO::load_binary_file(input, c3t3_bis);
      assert(input);
    }
    if(!check_equality(c3t3, c3t3_bis)) return false;

    return true;
  }

  bool operator()(const char* prefix) const {
    CGAL::Get_io_signature<C3t3> get_io_signature; CGAL_USE(get_io_signature);
    assert(get_io_signature() == Mesh_domain::reference_format_string());
    // Create a C3t3
    C3t3 c3t3;
    Tr& tr = c3t3.triangulation();
    Vertex_handle v1 = tr.insert(Point(10,11,12));
    Vertex_handle v2 = tr.insert(Point(11,13,10));
    Vertex_handle v3 = tr.insert(Point(7,4,6));
    Vertex_handle v4 = tr.insert(Point(5,2,14));
    Vertex_handle v5 = tr.insert(Point(1,2,3));
    Vertex_handle v6 = tr.insert(Point(3,9,13));

    Edge e1 = *(c3t3.triangulation().finite_edges_begin());
    Edge e2 = *(++c3t3.triangulation().finite_edges_begin());

    Facet f1 = *(c3t3.triangulation().finite_facets_begin());
    Facet f2 = *(++c3t3.triangulation().finite_facets_begin());

    Cell_handle c1 = c3t3.triangulation().finite_cells_begin();
    Cell_handle c2 = ++c3t3.triangulation().finite_cells_begin();

    c3t3.add_to_complex(c1, Mesh_domain::get_sub_domain_index_1());
    c3t3.add_to_complex(c2, Mesh_domain::get_sub_domain_index_2());
    c3t3.add_to_complex(f1, Mesh_domain::get_surface_patch_index_1());
    c3t3.add_to_complex(f2, Mesh_domain::get_surface_patch_index_2());
    c3t3.add_to_complex(e1, Mesh_domain::get_curve_index_1());
    c3t3.add_to_complex(e2, Mesh_domain::get_curve_index_2());
    c3t3.add_to_complex(v1, Mesh_domain::get_corner_index_1());
    c3t3.add_to_complex(v2, Mesh_domain::get_corner_index_2());

    // Fill indices in various faces (with different dimensions)
    v1->set_dimension(0);
    v1->set_index(Mesh_domain::get_corner_index_1());
    v2->set_dimension(0);
    v2->set_index(Mesh_domain::get_corner_index_2());
    v3->set_dimension(1);
    v3->set_index(Mesh_domain::get_curve_index_1());
    v4->set_dimension(1);
    v4->set_index(Mesh_domain::get_curve_index_2());
    v5->set_dimension(2);
    v5->set_index(Mesh_domain::get_surface_patch_index_1());
    v6->set_dimension(3);
    v6->set_index(Mesh_domain::get_sub_domain_index_1());

    // then test I/O
    return test_io(c3t3, prefix, false) && test_io(c3t3, prefix, true);
  }
};

int main()
{
  std::cout << "First test I/O when all indices are integers" << std::endl;
  bool ok = Test_c3t3_io<MD_homogeneous_types>()("data/c3t3_io-homo");
  if(!ok) {
    std::cerr << "Error\n";
    return -1;
  }
  std::cout << "Then test I/O when all indices are different types" << std::endl;
  ok = Test_c3t3_io<MD_heterogeneous_types>()("data/c3t3_io-hetero");
  if(!ok) {
    std::cerr << "Error\n";
    return -1;
  }
}
