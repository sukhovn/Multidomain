#include <domains.h>
#include <boundaries.h>
#include <spectral.h>

#ifndef DOMAINS_TEST_H
#define DOMAINS_TEST_H

void test_2D_basic_partitioning(const std::string &test_folder, int argc, char const *argv[]);
void test_3D_basic_partitioning(const std::string &test_folder, int argc, char const *argv[]);

void test_2D_boundary_partitioning(const std::string &test_folder, int argc, char const *argv[]);
void test_3D_boundary_partitioning(const std::string &test_folder, int argc, char const *argv[]);

void test_2D_regridding(const std::string &test_folder, int argc, char const *argv[]);

#endif