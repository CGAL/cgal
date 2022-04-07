@echo off

:go
${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/pkglist_filter.py %1

@echo on
