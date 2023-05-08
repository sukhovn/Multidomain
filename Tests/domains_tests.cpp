#include <domains_tests.h>

void test_2D_basic_partitioning(const std::string &test_folder, int argc, char const *argv[]){
	std::vector<double> posit = {0.0, 0.0};
	std::vector<double> lengths = {1.0, 1.0};

	Generic_Domain_Partitioning part(posit, lengths);

	part.domain_split({0.3, 0.5}, 0);
	part.domain_split({0.15, 0.3}, 1);
	part.domain_split({0.65, 0.5}, 1);

	part.domain_split({0.1, 0.15}, 0);
	part.domain_split({0.2, 0.15}, 0);
	part.domain_split({0.15, 0.2}, 1);

	std::cout << "Partitioning 1\n" << std::endl;
	print_partitioning(part);
	
	part.verify_consistency();
	part.save_domains(test_folder + "/domains2D.dat");
	
	Generic_Domain_Partitioning part2("Tests/domains_2D_test.dat");

	part.print_domain_2D(test_folder + "/domains_2D.pbm");
	part2.print_domain_2D(test_folder + "/domains_2D_test.pbm");

	std::cout << "Partitioning 2\n" << std::endl;
	print_partitioning(part2);

	part = part2;
	std::cout << "Copy test" << std::endl;
	if(part == part2) std::cout << "Copy test pass" << std::endl;
	else std::cout << "Copy test fail" << std::endl;

	return;
}

void test_3D_basic_partitioning(const std::string &test_folder, int argc, char const *argv[]){
	int dim = 3;
	std::vector<double> posit = {0.0, 0.0, 0.0};
	std::vector<double> lengths = {1.0, 1.0, 1.0};

	Generic_Domain_Partitioning part(posit, lengths);

	part.domain_split({0.3, 0.5, 0.5}, 0);
	part.domain_split({0.15, 0.3, 0.5}, 1);
	part.domain_split({0.15, 0.15, 0.3}, 2);
	part.domain_split({0.65, 0.5, 0.5}, 1);

	part.save_domains(test_folder + "/domains3D.dat");
	part.print_domain_3D(test_folder + "/domains3D");
	part.verify_consistency();

	return;
}

void test_2D_boundary_partitioning(const std::string &test_folder, int argc, char const *argv[]){
	std::vector<double> posit = {0.0, 0.0};
	std::vector<double> lengths = {1.0, 1.0};
	std::vector<int> dim_nums = {0, 0};

	// Generic_Domain_Partitioning part(posit, lengths);
	Domain_Partitioning bpart(posit, lengths, dim_nums);

	bpart.domain_split({0.3, 0.5}, 0);
	bpart.domain_split({0.15, 0.3}, 1);
	bpart.domain_split({0.65, 0.5}, 1);

	bpart.domain_split({0.1, 0.15}, 0);
	bpart.domain_split({0.2, 0.15}, 0);
	bpart.domain_split({0.15, 0.2}, 1);

	// Domain_Partitioning bpart(part);
	bpart.insert_boundary_domains();

	bpart.print_domain_2D(test_folder + "/domains_2D.pbm");
	std::cout << "Partitioning 1\n" << std::endl;
	print_partitioning(bpart);
	bpart.save_domains(test_folder + "/domains2D.dat");

	Domain_Partitioning bpart2("Tests/domains_2D_test.dat");
	bpart2.insert_boundary_domains();

	bpart2.print_domain_2D(test_folder + "/domains_2D_test.pbm");
	std::cout << "Partitioning 2\n" << std::endl;
	print_partitioning(bpart2);
	bpart2.save_domains(test_folder + "/domains2D_test.dat");

	bpart2.clear_boundary_domains();
	bpart2.insert_boundary_domains();
	bpart2.fill_boundary_point_range(ch_x_inv_low);
	std::cout << "Copy test" << std::endl;
	Domain_Partitioning bpart_copy(bpart2);
	if(bpart_copy.compare_bulk_domains(bpart2)) std::cout << "Incomplete copy test pass" << std::endl;
	if(bpart_copy == bpart2) std::cout << "Copy test pass" << std::endl;
	else std::cout << "Copy test fail" << std::endl;

	std::cout << "\nBulk comparison test" << std::endl;
	bpart_copy.clear_boundary_domains();
	if(bpart_copy.compare_bulk_domains(bpart2)) std::cout << "Incomplete copy test pass" << std::endl;
	if(bpart_copy == bpart2) std::cout << "Copy test pass" << std::endl;
	else std::cout << "Copy test fail" << std::endl;

	// print_partitioning(bpart_copy);

	return;
}

void test_3D_boundary_partitioning(const std::string &test_folder, int argc, char const *argv[]){
	std::vector<double> posit = {0.0, 0.0, 0.0};
	std::vector<double> lengths = {1.0, 1.0, 1.0};
	std::vector<int> dim_nums = {0, 0, 0};

	// Generic_Domain_Partitioning part(posit, lengths);
	Domain_Partitioning bpart(posit, lengths, dim_nums);

	bpart.domain_split({0.3, 0.5, 0.5}, 0);
	bpart.domain_split({0.15, 0.3, 0.5}, 1);
	bpart.domain_split({0.15, 0.15, 0.3}, 2);
	bpart.domain_split({0.65, 0.5, 0.5}, 1);

	// Domain_Partitioning bpart(part);
	bpart.insert_boundary_domains();

	bpart.save_domains(test_folder + "/domains3D.dat");
	bpart.print_domain_3D(test_folder + "/domains3D");
	bpart.verify_consistency();
	print_partitioning(bpart);

	std::cout << "Copy test" << std::endl;
	Domain_Partitioning bpart_copy(bpart);
	if(bpart_copy == bpart) std::cout << "Copy test pass" << std::endl;
	else std::cout << "Copy test fail" << std::endl;

	return;
}

void test_2D_regridding(const std::string &test_folder, int argc, char const *argv[]){
	Domain_Partitioning bpart("Tests/domains_2D_test.dat");

	bpart.print_domain_2D(test_folder + "/domains_2D_test_before_regridding.pbm");

	bpart.domain_split({0.5, 0.6}, 1);
	bpart.move_domain_wall({0.5, 0.6}, 0.7);

	bpart.move_domain_wall({0.5, 0.8}, 0.9);
	bpart.move_domain_wall({0.8, 0.4}, 0.9);

	std::vector<int> domain_num = bpart({0.5, 0.6})->dim_num();
	domain_num[0] += 5;
	bpart.change_domain_resolution({0.5, 0.6}, domain_num);

	bpart.print_domain_2D(test_folder + "/domains_2D_test_after_regridding.pbm");

	print_partitioning(bpart);

	return;
}