# Try to find TinyGLTF
# Once done, this will define:
#
#  TINYGLTF_FOUND       - TRUE if TinyGLTF was found
#  TINYGLTF_INCLUDE_DIR - include directory containing tiny_gltf.h

# If CGAL was found via vcpkg, derive the vcpkg include dir from CGAL_DIR
# (e.g. .../vcpkg/installed/x64-windows/share/cgal -> .../vcpkg/installed/x64-windows/include)
if(CGAL_DIR)
  get_filename_component(_tinygltf_vcpkg_hint "${CGAL_DIR}/../../include" ABSOLUTE)
endif()

find_path(TINYGLTF_INCLUDE_DIR
  NAMES tiny_gltf.h
  HINTS
    ${_tinygltf_vcpkg_hint}
  PATHS
    /usr/include
    /usr/local/include
    ENV TINYGLTF_INC_DIR
)

if(TINYGLTF_INCLUDE_DIR)
  set(TINYGLTF_FOUND TRUE)
endif()