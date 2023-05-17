#include <test_equations.h>
#include <domains.h>
#include <boundaries.h>
#include <container.h>
#include <collocation_solver_full.h>

#ifndef COLLOCATION_2D_TEST_H
#define COLLOCATION_2D_TEST_H

void collocation_2D_container_test(const std::string &test_folder, int argc, char const *argv[]);

void collocation_2D_collocation_fill_test(const std::string &test_folder, int argc, char const *argv[]);

// 3D tests

void collocation_3D_container_test(const std::string &test_folder, int argc, char const *argv[]);

void collocation_3D_collocation_fill_test(const std::string &test_folder, int argc, char const *argv[]);

#endif
