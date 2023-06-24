
#include "Kml_reader.h"

#include <iostream>

#include <qdebug.h>
#include <qfile.h>
#include <qxmlstream.h>

namespace {
 
  std::vector<std::string> split(const std::string& str, const char *delim)
  {
    std::string sc = str;
    char* token = strtok(sc.data(), delim);
    char* str_end = token + str.length();

    // Keep printing tokens while one of the delimiters present in str[].
    std::vector<std::string> results;
    while (token != NULL)
    {
      const auto first = token;
      //printf("%s\n", token);
      token = strtok(NULL, " ");
      results.push_back(std::string(first, token==nullptr ? str_end : token));
    }

    return results;
  }
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
          auto qstr = xmlReader.text().toString();
          auto ptr = qstr.data();
          auto str = qstr.toUtf8().toStdString();
          auto node_strs = split(str, " ");
          
          for (const auto& node_str : node_strs)
          {
            if (node_str.empty())
              continue;

            auto coord_strs = split(node_str, ",");
            const auto lon = std::stod(coord_strs[0]);
            const auto lat = std::stod(coord_strs[1]);
            lring.nodes.push_back(Node{ lon, lat });
          }


          //qDebug() << "---------------------------------------";
          //for (const auto& node_str : node_strs)
          //  std::cout << node_str << std::endl;

          //qDebug() << qstr;
          //auto node_qstrs = qstr.split(" ");
          //qDebug() << node_qstrs.size();
          //for (const auto& node_str : node_strs)
          //{
          //  if (node_str.isEmpty())
          //    continue;

          //  auto coord_strs = node_str.split(",");
          //  const auto lon = coord_strs[0].toDouble();
          //  const auto lat = coord_strs[1].toDouble();
          //  lring.nodes.push_back(Node{ lon, lat });
          //}
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
