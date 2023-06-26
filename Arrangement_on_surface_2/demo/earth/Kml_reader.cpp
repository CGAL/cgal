
#include "Kml_reader.h"

#include <iostream>
#include <unordered_map>

#include <qdebug.h>
#include <qfile.h>
#include <qxmlstream.h>


QVector3D Kml::Node::get_coords_3d(const double r) const
{
  const auto phi = qDegreesToRadians(lat);
  const auto theta = qDegreesToRadians(lon);
  const auto z = r * std::sin(phi);
  const auto rxy = r * std::cos(phi);
  const auto x = rxy * std::cos(theta);
  const auto y = rxy * std::sin(theta);

  return QVector3D(x, y, z);
}


Kml::Placemarks  Kml::read(const std::string& file_name)
{
  Placemarks    placemarks;
  Placemark     placemark;
  LinearRing    lring;

  QFile file(file_name.c_str());
  if (file.open(QIODevice::ReadOnly))
  {
    QXmlStreamReader xmlReader;
    xmlReader.setDevice(&file);

    xmlReader.readNext();

    // Reading from the file
    while (!xmlReader.isEndDocument())
    {
      QString name = xmlReader.name().toString();

      if (xmlReader.isStartElement())
      {
        if (name == "Placemark")
        {
          placemark = Placemark{};
        }
        else if (name == "Polygon")
        {
          lring = LinearRing{};
        }
        else if (name == "LinearRing")
        {
        }
        else if (name == "coordinates")
        {
          xmlReader.readNext();
          auto str = xmlReader.text().toString();
          auto node_strs = str.split(" ");
          for (const auto& node_str : node_strs)
          {
            if (node_str.isEmpty())
              continue;

            auto coord_strs = node_str.split(",");
            const auto lon = coord_strs[0].toDouble();
            const auto lat = coord_strs[1].toDouble();
            lring.nodes.push_back(Node{ lon, lat });
          }
        }
        else if (name == "SimpleData")
        {
          auto attributes = xmlReader.attributes();
          auto attr_name = attributes[0].name().toString();
          auto attr_value = attributes[0].value().toString();
          if ((attr_name == "name") && (attr_value == "name"))
          {
            xmlReader.readNext();
            placemark.name = xmlReader.text().toString().toStdString();;
          }
        }
      }
      else if (xmlReader.isEndElement())
      {
        if (name == "Placemark")
        {
          placemarks.push_back(std::move(placemark));
        }
        else if (name == "Polygon")
        {
          placemark.polygons.push_back(std::move(lring));
        }
        else if (name == "LinearRing")
        {
        }
        else if (name == "coordinates")
        {
          // no need to do anything here: the coordinates are read above!
        }
      }

      xmlReader.readNext();
    }

    if (xmlReader.hasError())
    {
      std::cout << "XML error: " << xmlReader.errorString().data() << std::endl;
    }
  }

  return placemarks;
}

Kml::Nodes Kml::get_duplicates(const Placemarks& placemarks)
{
  // collect all nodes into a single vector
  int polygon_count = 0;
  std::vector<Kml::Node> nodes;
  for (const auto& pm : placemarks)
  {
    for (const auto& pg : pm.polygons)
    {
      polygon_count++;
      for (const auto& node : pg.nodes)
        nodes.push_back(node);
    }
  }
  qDebug() << "polygon count = " << polygon_count;

  int count = nodes.size();
  std::vector<int> num_duplicates(count, 0);
  qDebug() << "node count (with duplicates) = " << count;
  int dup_count = 0;

  // this keeps track of how many nodes there are with certain dup-count
  std::unordered_map<int, int> dup_count_map;

  Nodes duplicate_nodes;
  for (int i = 0; i < count; ++i)
  {
    // if the current node has been detected as duplicate skip it
    if (num_duplicates[i] > 0)
      continue;

    const auto& curr_node = nodes[i];
    std::vector<int> curr_dup; // current set of duplicates
    for (int j = i + 1; j < count; ++j)
    {
      if (curr_node == nodes[j])
      {
        curr_dup.push_back(j);
      }
    }

    // if duplicates found
    if (!curr_dup.empty())
    {
      dup_count++;
      int num_dup = curr_dup.size() + 1; // +1 for the i'th node
      num_duplicates[i] = num_dup;
      for (const auto di : curr_dup)
        num_duplicates[di] = num_dup;

      duplicate_nodes.push_back(curr_node);
      dup_count_map[num_dup]++;
    }
  }
  qDebug() << "dup count = " << dup_count;
  for (const auto& p : dup_count_map)
  {
    const auto dup_count = p.first;
    const auto num_nodes_with_this_dup_count = p.second;
    qDebug() << dup_count << ": " << num_nodes_with_this_dup_count;
  }

  return duplicate_nodes;
}
