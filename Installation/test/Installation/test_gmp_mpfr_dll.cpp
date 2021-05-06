// This test is for MSVC only.
#ifndef _MSC_VER
int main() {
  return 0;
}
#else

#define GMP_SONAME "libgmp-10"
#define MPFR_SONAME "libmpfr-4"
#define GMP_SONAME_BACKUP "gmp"
#define MPFR_SONAME_BACKUP "mpfr-6"
#define GMP_MAJOR 5
#define MPFR_MAJOR 3


#include <iostream>
#include <cassert>
#include <windows.h>

#pragma warning(disable:4244 4146)
// conversion with loss of data
// warning on - applied on unsigned number

#include "gmp.h"
#include <mpfr.h>

bool get_version_info(const LPCTSTR name,
                      int& major,
                      int& minor,
                      int& patch,
                      int& build)
{
  HMODULE g_dllHandle = GetModuleHandle(name);
  if(!g_dllHandle) {
    std::cerr << name << " is not loaded!\n";
    return false;
  }
  char fileName[_MAX_PATH];
  DWORD size = GetModuleFileName(g_dllHandle, fileName, _MAX_PATH);
  fileName[size] = NULL;
  std::cerr << "Query FileVersion of \"" << fileName << "\"\n";
  DWORD handle = 0;
  size = GetFileVersionInfoSize(fileName, &handle);
  BYTE* versionInfo = new BYTE[size];
  if (!GetFileVersionInfo(fileName, handle, size, versionInfo))
  {
    delete[] versionInfo;
    std::cerr << name << " has no VersionInfo!\n";
    return true;
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

int main() {
  std::cout << "Hello GMP version " << gmp_version << std::endl;
  std::cout << "Hello MPFR version " << mpfr_get_version() << std::endl;
  int major, minor, patch, build;
  if(!get_version_info(GMP_SONAME, major, minor, patch, build)) {
    if(!get_version_info(GMP_SONAME_BACKUP, major, minor, patch, build)) {
      return 1;
    }
  }

  std::cout << "GMP version "
            << major << "."
            << minor << "."
            << patch << "."
            << build << "\n";
  major = 0;
  if(!get_version_info(MPFR_SONAME, major, minor, patch, build)) {
    if(!get_version_info(MPFR_SONAME_BACKUP, major, minor, patch, build)) {
      return 1;
    }
  }
  std::cout << "MPFR version "
            << major << "."
            << minor << "."
            << patch << "."
            << build << "\n";
}
#endif
