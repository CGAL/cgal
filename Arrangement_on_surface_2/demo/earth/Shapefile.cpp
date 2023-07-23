//
//#include "Shapefile.h"
//
//#include <iostream>
//
//#include <shapefil.h>
//
//
//void Shapefile::read(const std::string& filename) 
//{
//  // Open the shapefile
//  SHPHandle shp = SHPOpen(filename.c_str(), "rb");
//  if (shp == nullptr) {
//    std::cerr << "Failed to open shapefile: " << filename << std::endl;
//    return;
//  }
//
//  // Get shapefile information
//  int numEntities, shapeType;
//  double minBounds[4], maxBounds[4];
//  SHPGetInfo(shp, &numEntities, &shapeType, minBounds, maxBounds);
//  std::cout << "Number of entities: " << numEntities << std::endl;
//  std::cout << "Shape type: " << shapeType << std::endl;
//  std::cout << "Bounds: (" << minBounds[0] << ", " << minBounds[1] << "), ("
//    << maxBounds[0] << ", " << maxBounds[1] << ")" << std::endl;
//  //SHPT_POLYGON
//  // Read individual shapes
//  for (int i = 0; i < numEntities; ++i) {
//    SHPObject* shape = SHPReadObject(shp, i);
//
//    // Process the shape data
//    // Example: Print the shape's type and number of points
//    std::cout << "Shape " << i << ": Type " << shape->nSHPType
//      << ", Number of parts: " << shape->nParts
//      << ", Number of points: " << shape->nVertices << std::endl;
//
//    // Clean up the shape object
//    SHPDestroyObject(shape);
//  }
//
//  // Close the shapefile
//  SHPClose(shp);
//}
