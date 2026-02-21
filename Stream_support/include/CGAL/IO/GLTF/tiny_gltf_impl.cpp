// Compile tinygltf + stb implementation in exactly one translation unit.

#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <CGAL/IO/GLTF/tiny_gltf.h>
