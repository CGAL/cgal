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

#include <FileGDBAPI.h>

#include "Globus_window.h"

//
struct country {
  //! Constructor
  country(std::string& name) :
    m_name(std::move(name))
  {}

  std::string m_name;
};

namespace FGA = FileGDBAPI;

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

  return S_OK;
}

//
int main(int argc, char* argv[]) {
  // Load the countries based on a search extent.
  std::vector<country> countries;
  size_t hr = load_countries(countries);
  if (hr != S_OK) {
    std::cerr << "Failed to load database!\n";
    return -1;
  }
  for (const auto& country : countries) {
    std::cout << country.m_name << std::endl;
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
