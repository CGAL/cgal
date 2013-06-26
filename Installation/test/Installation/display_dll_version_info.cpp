// This program is for MSVC only.
#ifndef WIN32
int main() {
  return 0;
}
#else

#include <iostream>
#include <cassert>
#include <windows.h>

bool get_version_info(const char* fileName, 
                      int& major,
                      int& minor,
                      int& patch,
                      int& build)
{
  std::cout << "Query FileVersion of \"" << fileName << "\"\n";
  DWORD handle = 0;
  DWORD size = GetFileVersionInfoSize(fileName, &handle);
  BYTE* versionInfo = new BYTE[size];
  if (!GetFileVersionInfo(fileName, handle, size, versionInfo))
  {
    delete[] versionInfo;
    return false;
  }
  // we have version information
  UINT len = 0;
  VS_FIXEDFILEINFO*   vsfi = NULL;
  VerQueryValue(versionInfo, "\\", (void**)&vsfi, &len);
  major = HIWORD(vsfi->dwFileVersionMS);
  minor = LOWORD(vsfi->dwFileVersionMS);
  patch = HIWORD(vsfi->dwFileVersionLS);
  build = LOWORD(vsfi->dwFileVersionLS);
  delete[] versionInfo;
  return true;
}

int main(int argc, char** argv) {
  if(argc == 2) {
    int major, minor, patch, build;
    bool result = get_version_info(argv[1], major, minor, patch, build);
    assert(result);
    std::cout << "version "
              << major << "."
              << minor << "."
              << patch << "."
              << build << "\n";
  } else {
    std::cerr << "Usage:\n"
              << "  display_dll_version_info /path/to/a.dll\n";
    return 0;
  }
}
#endif
