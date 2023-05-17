#include <solver_tests.h>

static double vector_max(const std::vector<double> &vec){
	double max = vec[0];
	for(int i = 1; i < vec.size(); ++i)
		if(vec[i] > max) max = vec[i];
	return max;
}

void collocation_2D_collocation_test_solution(const std::string &test_folder, int argc, char const *argv[]){
	Domain_Partitioning bpart("Tests/simple_domains_2D.dat");
	// Domain_Partitioning bpart("Tests/domains_2D_test.dat");
	
	Container cont(2, bpart);
	cont.set_doubling_begin(0, {false, false});
	cont.set_doubling_begin(1, {false, false});

	double a_param = 50.0;
	Full_Compactification compact(cont, true);
	std::cout << std::endl;

	Cylindric_Test eqn_system(a_param);
	Collocation_Solver_Full<Full_Compactification, Cylindric_Test> solver(compact, &eqn_system, &cont);

	std::cout << "Number of points: " << cont.size() << std::endl;

	solver.solve();
	std::cout << "The equation difference is: " << vector_max(solver.find_function_maximum(&Cylindric_Test::equations, 2)) << std::endl;
	std::cout << "The solution difference is: " << vector_max(solver.find_function_maximum(&Cylindric_Test::solution_check, 2, false)) << std::endl;
	cont.save(test_folder + "/colloc_compact");
}