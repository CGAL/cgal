#pragma once

#include <iostream>

// helper macro

#define CENTER_NOP do { } while (0)

// printing macros

#ifdef NDEBUG
#define DEBUG(x) CENTER_NOP
#else
#define DEBUG(x) do { std::cout << x << std::endl; } while (0)
#endif

#if defined(NVERBOSE) && defined(NDEBUG)
#define PRINT(x) CENTER_NOP
#else
#define PRINT(x) do { std::cout << x << std::endl; } while (0)
#endif

#define ERROR(x) do { std::cerr << "Error: " << x << std::endl;\
                      std::exit(EXIT_FAILURE); } while (0)

#include <cassert>
