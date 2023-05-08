#include <domains_tests.h>
#include <collocation_tests.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

int main(int argc, char const *argv[]){
	std::string test_folder = "test_results";

	struct stat sttus = {0};
	if (stat(test_folder.c_str(), &sttus) == -1) {
		mkdir(test_folder.c_str(), 0700);
	}

	// Domain partitioning tests

	// test_2D_basic_partitioning(test_folder, argc, argv);
	// test_3D_basic_partitioning(test_folder, argc, argv);
	// test_2D_boundary_partitioning(test_folder, argc, argv);
	// test_3D_boundary_partitioning(test_folder, argc, argv);

	// test_2D_regridding(test_folder, argc, argv);
	
	// 2D Collocation tests

	collocation_2D_container_test(test_folder, argc, argv);
	
	return 0;
}