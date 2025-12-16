#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/GLTF.h>

#include <vector>
#include <iostream>

int main() {
    using Kernel = CGAL::Simple_cartesian<double>;
    using Point = Kernel::Point_3;

    std::string filename = "./data//Cube.gltf";

    std::cout << "Attempting to load GLTF file: " << filename << std::endl;

    tinygltf::Model model;
    tinygltf::TinyGLTF loader;
    std::string err;
    std::string warn;

    bool ret = loader.LoadASCIIFromFile(&model, &err, &warn, filename);

    if (!warn.empty()) {
        std::cout << "GLTF warning: " << warn << std::endl;
    }
    if (!err.empty()) {
        std::cerr << "GLTF error: " << err << std::endl;
    }
    if (!ret) {
        std::cerr << "Failed to parse GLTF file: " << filename << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Successfully parsed GLTF file." << std::endl;
    std::cout << "Number of meshes: " << model.meshes.size() << std::endl;
    std::cout << "Number of accessors: " << model.accessors.size() << std::endl;
    std::cout << "Number of bufferViews: " << model.bufferViews.size() << std::endl;
    std::cout << "Number of buffers: " << model.buffers.size() << std::endl;

    // Now use the CGAL reader
    std::vector<Point> points;
    std::vector<std::vector<std::size_t>> polygons;

    bool success = CGAL::IO::internal::read_GLTF(filename, points, polygons, true);

    if (!success) {
        std::cerr << "Failed to read GLTF file: " << filename << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Successfully loaded GLTF model.\n";
    std::cout << "→ Number of vertices: " << points.size() << std::endl;
    std::cout << "→ Number of faces: " << polygons.size() << std::endl;

    std::cout << "\nSample vertices:\n";
    for (size_t i = 0; i < std::min(points.size(), size_t(5)); ++i) {
        std::cout << "  " << points[i] << std::endl;
    }

    std::cout << "\nSample faces:\n";
    for (size_t i = 0; i < std::min(polygons.size(), size_t(5)); ++i) {
        std::cout << "  ";
        for (auto idx : polygons[i])
            std::cout << idx << " ";
        std::cout << std::endl;
    }

    return 0;
}
