
#include "Kml_reader.h"

#include <iostream>

#include <qfile.h>
#include <qxmlstream.h>


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
