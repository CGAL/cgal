#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/GLTF/read_GLTF.h>
#include <vector>
#include <iostream>

int main() {
    using Kernel = CGAL::Simple_cartesian<double>;
    using Point = Kernel::Point_3;

    std::vector<Point> points;
    std::vector<std::vector<std::size_t>> polygons;

    std::string filename = "Box.gltf";

    std::cout << "Attempting to load GLTF file: " << filename << std::endl;

    bool success = CGAL::IO::read_GLTF(filename, points, polygons);

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

    return EXIT_SUCCESS;
}
