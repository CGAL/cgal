// Copyright (c) 2022, 2020 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel <efifogel@gmail.com>

#include <string>
#include <vector>

#include <QApplication>

// #include <FileGDBAPI.h>
#include <nlohmann/json.hpp>

#include "Globus_window.h"

enum Error_id {
  FILE_NOT_FOUND,
  ILLEGAL_EXTENSION,
  UNABLE_TO_OPEN,
  FILE_IS_EMPTY,
  INVALID_INITIAL_POLYGON,
  UNSUPPORTED,
  INVALID_OUTPUT_FILE,
  ILLEGAL_SOLUTION_FILE
};

struct Illegal_input : public std::logic_error {
  Illegal_input(Error_id /* err */, const std::string &msg,
    const std::string &filename) :
    std::logic_error(std::string(msg).append(" (").append(filename).
      append(")!"))
  {}

  Illegal_input(Error_id /* err */, const std::string &msg) :
    std::logic_error(std::string(msg).append("!"))
  {}
};

struct Input_file_missing_error : public std::logic_error {
  Input_file_missing_error(std::string &str) : std::logic_error(str) {}
};

//
struct country {
  //! Constructor
  country(std::string& name) :
    m_name(std::move(name))
  {}

  std::string m_name;
};

// namespace FGA = FileGDBAPI;

#if 0
/*!
 */
int load_countries(std::vector<country>& countries) {
  // Open the geodatabase and the table.
  fgdbError hr;
  FGA::Geodatabase geodatabase;
  hr = FGA::OpenGeodatabase(L"/home/efif/tmp/esri/19113835-45b5-40b9-989b-fd04aad324da.gdb", geodatabase);
  if (hr != S_OK) return hr;

  // std::vector<std::wstring> types;
  // hr = geodatabase.GetDatasetTypes(types);
  // if (hr != S_OK) return hr;
  // std::cout << "Types:\n";
  // for (const auto& type : types)
  //   std::wcout << type << std::endl;

  FGA::Table table;
  hr = geodatabase.OpenTable(L"World_Countries", table);
  if (hr != S_OK) {
    std::cerr << "Cannot open table!\n";
    return hr;
  }

  // std::string doc;
  // hr = table.GetDocumentation(doc);
  // if (hr != S_OK) return hr;
  // std::cout << "Table Information:\n";
  // std::cout << doc << std::endl;

  // FGA::fieldDefs field_defs;
  // hr = table.GetFields(fields_defs);
  // if (hr != S_OK) return hr;
  // std::cout << "Fields:\n";
  // for (const auto& field : Fields)
  //   std::wcout << field << std::endl;
  FGA::FieldInfo field_info;
  hr = table.GetFieldInformation(field_info);
  if (hr != S_OK) return hr;
  int count;
  hr = field_info.GetFieldCount(count);
  std::cout << "Number of fileds: " << count << std::endl;
  for (auto i = 0; i < count; ++i) {
    std::wstring name;
    hr = field_info.GetFieldName(i, name);
    if (hr != S_OK) return hr;
    std::wcout << "Name: " << name << std::endl;
  }

  FGA::Envelope envelope;
  FGA::EnumRows rows;
  hr = table.Search(L"SHAPE, COUNTRY", L"", envelope, true, rows);
  if (hr != S_OK) {
    std::cerr << "Cannot find rows!\n";
    return hr;
  }

  countries.clear();
  FGA::Row row;
  while (rows.Next(row) == S_OK) {
    std::wstring name;
    row.GetString(L"COUNTRY", name);
    std::string simple_name;
    simple_name.assign(name.begin(), name.end());
    countries.emplace_back(simple_name);
  }

  // Close the table
  hr = geodatabase.CloseTable(table);
  if (hr != S_OK) {
    std::wcout << "An error occurred while closing Cities." << endl;
    std::wcout << "Error code: " << hr << endl;
    return -1;
  }

  // Close the geodatabase
  hr = FGA::CloseGeodatabase(geodatabase);
  if (hr != S_OK) {
    std::wcout << "An error occurred while closing the geodatabase." << endl;
    std::wcout << "Error code: " << hr << endl;
    return -1;
  }

  return 0;
}
#endif

/*! Read a json file.
 */
bool read_json(const std::string& filename, nlohmann::json& data) {
  using json = nlohmann::json;
  std::ifstream infile(filename);
  if (! infile.is_open()) {
    throw Illegal_input(UNABLE_TO_OPEN, "Cannot open file", filename);
    return false;
  }
  data = json::parse(infile);
  infile.close();
  if (data.empty()) {
    throw Illegal_input(FILE_IS_EMPTY, "File is empty", filename);
    return false;
  }
  return true;
}

bool read_arrangement(const std::string& filename) {
  using json = nlohmann::json;
  json data;
  auto rc = read_json(filename, data);
  if (! rc) return false;

  // points
  auto it = data.find("points");
  if (it == data.end()) {
    std::cerr << "The points item is missing " << " (" << filename << ")\n";
    return false;
  }
  const auto& points = it.value();
  std::cout << "# points: " << points.size() << std::endl;

  // curves
  it = data.find("curves");
  if (it == data.end()) {
    std::cerr << "The curves item is missing " << " (" << filename << ")\n";
    return false;
  }
  const auto& curves = it.value();
  std::cout << "# curves: " << curves.size() << std::endl;

  // vertices
  it = data.find("vertices");
  if (it == data.end()) {
    std::cerr << "The vertices item is missing " << " (" << filename << ")\n";
    return false;
  }
  const auto& vertices = it.value();
  std::cout << "# vertices: " << vertices.size() << std::endl;

  // halfedges
  it = data.find("halfedges");
  if (it == data.end()) {
    std::cerr << "The halfedges item is missing " << " (" << filename << ")\n";
    return false;
  }
  const auto& halfedges = it.value();
  std::cout << "# halfedges: " << halfedges.size() << std::endl;

  // faces
  it = data.find("faces");
  if (it == data.end()) {
    std::cerr << "The faces item is missing " << " (" << filename << ")\n";
    return false;
  }
  const auto& faces = it.value();
  std::cout << "# faces: " << faces.size() << std::endl;

  return true;
}

//
int main(int argc, char* argv[]) {
  const char* filename = (argc > 1) ? argv[1] : "arr.json";

  auto rc = read_arrangement(filename);
  if (! rc) {
    std::cerr << "Failed to load database!\n";
    return -1;
  }

  QApplication app(argc, argv);
  QCoreApplication::setOrganizationName("CGAL");
  QCoreApplication::setApplicationName("globus");

  // Import resources from libCGAL (Qt5).
  CGAL_QT_INIT_RESOURCES;

  Globus_window window;
  window.show();
  return app.exec();
}
