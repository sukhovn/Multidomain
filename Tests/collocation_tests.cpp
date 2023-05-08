#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

#include <collocation_tests.h>

double vector_max(const std::vector<double> &vec){
	double max = vec[0];
	for(int i = 1; i < vec.size(); ++i)
		if(vec[i] > max) max = vec[i];
	return max;
}

double sin_func(int fi, double *x){
	if(fi == 0)
		return std::sin(M_PI*x[0])*std::sin(M_PI*x[1]);
	else
		return std::cos(2.0*M_PI*x[0])*std::cos(2.0*M_PI*x[1]);
}

class Trig_Func{
	std::vector<unsigned int> node_x;
	std::vector<unsigned int> node_y;
	int ord_x;
	int ord_y;
	double sign;
public:
	Trig_Func(std::vector<unsigned int> node_xi, std::vector<unsigned int> node_yi): 
		node_x(node_xi), node_y(node_yi), ord_x(0), ord_y(0), sign(1.0)
	{
		if(ord_x%4 > 1) sign *= -1.0;
		if(ord_y%4 > 1) sign *= -1.0;
	}
	Trig_Func(const Trig_Func &trg, int ord_x_i, int ord_y_i): 
		node_x(trg.node_x), node_y(trg.node_y), ord_x(trg.ord_x+ord_x_i), ord_y(trg.ord_y+ord_y_i), sign(1.0)
	{
		if(ord_x%4 > 1) sign *= -1.0;
		if(ord_y%4 > 1) sign *= -1.0;
	}

	double operator()(int fn, double *x){
		double nx = node_x[fn];
		double ny = node_y[fn];
		double tmp = sign*std::pow(nx*M_PI, ord_x)*std::pow(ny*M_PI, ord_y);
		if(ord_x%2) tmp *= std::cos(nx*M_PI*x[0]);
		else tmp *= std::sin(nx*M_PI*x[0]);
		if(ord_y%2) tmp *= std::cos(ny*M_PI*x[1]);
		else tmp *= std::sin(ny*M_PI*x[1]);
		return tmp;
	}
};

void collocation_2D_container_test(const std::string &test_folder, int argc, char const *argv[]){
	Domain_Partitioning part("Tests/simple_domains_2D.dat");
	part.print_domain_2D(test_folder + "/simple_domains_2D.pbm");
	// Domain_Partitioning part("Tests/domains_2D_test.dat");
	
	Container tcont(2, part);

	// std::cout << "Partitioning\n" << std::endl;
	// print_partitioning(tcont);

	tcont.fill_function(sin_func);

	tcont.save(test_folder + "/cont_test");

	Evaluator eval(tcont);
	eval.test(tcont);

	Domain_Partitioning part_new("Tests/domains_2D_test.dat");
	part_new.print_domain_2D(test_folder + "/domains_2D_test.pbm");

	tcont.regrid(part_new);
	std::cout << "The difference is: " << tcont.compare_function(sin_func) << std::endl;
	tcont.save(test_folder + "/cont_test_regrid");

	part_new.domain_split({0.6, 0.6}, 0);
	part_new.print_domain_2D(test_folder + "/domains_2D_test_split.pbm");

	tcont.regrid(part_new);
	std::cout << "The difference is: " << tcont.compare_function(sin_func) << std::endl;
	tcont.save(test_folder + "/cont_test_regrid_2");

	part_new.move_domain_wall({0.8, 0.5}, 0.7);
	part_new.print_domain_2D(test_folder + "/domains_2D_test_split_move.pbm");

	tcont.regrid(part_new);
	std::cout << "The difference is: " << tcont.compare_function(sin_func) << std::endl;
	tcont.save(test_folder + "/cont_test_regrid_3");

	return;
}