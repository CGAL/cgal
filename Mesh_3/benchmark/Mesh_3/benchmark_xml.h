#ifndef CGAL_MESH_3_BENCHMARK_BENCHMARK_XML_H
#define CGAL_MESH_3_BENCHMARK_BENCHMARK_XML_H

#include "../../test/Mesh_3/XML_exporter.h"

#define CGAL_MESH_3_EXPORT_PERFORMANCE_DATA
#define CGAL_MESH_3_SET_PERFORMANCE_DATA(value_name, value) \
        XML_perf_data::set(value_name, value);

class XML_perf_data
{
public:
  typedef Streaming_XML_exporter<std::string> XML_exporter;

private:
  XML_exporter m_xml;
  XML_exporter::Element_with_map m_current_element;

public:
  static std::string default_filename;

  XML_perf_data(const std::string& filename)
    : m_xml(filename, "ContainerPerformance", "Perf", construct_subelements_names())
  {}

  virtual ~XML_perf_data()
  {}

  static XML_perf_data& get()
  {
    static XML_perf_data singleton(build_filename());
    return singleton;
  }

  template <typename Value_type>
  static void set(const std::string& name, Value_type value)
  {
    get().set_data(name, value);
  }

  static void commit()
  {
    get().commit_current_element();
  }

protected:
  static std::string build_filename()
  {
    if(default_filename != "")
      return default_filename;

    std::stringstream sstr;
    sstr << "Performance_log_" << time(0) << ".xml";
    return sstr.str();
  }

  static std::vector<std::string> construct_subelements_names()
  {
    std::vector<std::string> subelements;
    subelements.push_back("Domain");

    subelements.push_back("Facet_angle");
    subelements.push_back("Facet_size");
    subelements.push_back("Facet_approx");
    subelements.push_back("Cell_size");
    subelements.push_back("Cell_shape");

    subelements.push_back("Technique");
    subelements.push_back("Num_threads");
    subelements.push_back("Lockgrid_size");
    subelements.push_back("Lock_radius");
    subelements.push_back("Statgrid_size");
    subelements.push_back("Num_work_items_per_batch");

    subelements.push_back("V");
    subelements.push_back("F");
    subelements.push_back("C");

    subelements.push_back("Facets_scan_time");
    subelements.push_back("Facets_refine_time");
    subelements.push_back("Cells_scan_time");
    subelements.push_back("Cells_refine_time");
    subelements.push_back("Lloyd_optim_time");
    subelements.push_back("Odt_optim_time");
    subelements.push_back("Perturber_optim_time");
    subelements.push_back("Exuder_optim_time");
    subelements.push_back("Total_time");

    subelements.push_back("Mem");

    subelements.push_back("Minimum_edge_length");
    subelements.push_back("Mean_edge_length");
    subelements.push_back("Maximum_edge_length");

    subelements.push_back("Minimum_facet_area");
    subelements.push_back("Mean_facet_area");
    subelements.push_back("Maximum_facet_area");
    subelements.push_back("Total_area");

    subelements.push_back("Minimum_facet_angle");
    subelements.push_back("Maximum_facet_angle");

    subelements.push_back("Minimum_cell_volume");
    subelements.push_back("Mean_cell_volume");
    subelements.push_back("Maximum_cell_volume");
    subelements.push_back("Total_volume");

    subelements.push_back("Minimum_dihedral_angle");
    subelements.push_back("Mean_dihedral_angle");
    subelements.push_back("Maximum_dihedral_angle");

    subelements.push_back("Smallest_edge_radius_ratio");
    subelements.push_back("Smallest_radius_radius_ratio");
    subelements.push_back("Biggest_V_SMA");

    return subelements;
  }

  void set_data(const std::string& name, const std::string& value)
  {
    m_current_element[name] = value;
  }

  template <typename Value_type>
  void set_data(const std::string& name, Value_type value)
  {
    std::stringstream sstr;
    sstr << value;
    set_data(name, sstr.str());
  }

  void commit_current_element()
  {
    m_xml.add_element(m_current_element);
    m_current_element.clear();
  }
};

#endif // CGAL_MESH_3_BENCHMARK_BENCHMARK_XML_H
