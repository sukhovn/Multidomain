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
	bpart.print_domain_2D(test_folder + "/domains_2D_solver.pbm");
	
	Container cont(2, bpart);
	cont.set_doubling_begin(0, {false, false});
	cont.set_doubling_begin(1, {false, false});

	double a_param = 50.0;
	Full_Compactification compact(bpart, true);
	std::cout << std::endl;

	Cylindric_Test eqn_system(a_param);
	Collocation_Solver_Full<Full_Compactification, Cylindric_Test> solver(compact, &eqn_system, &cont);

	std::cout << "Number of points: " << cont.size() << std::endl;

	solver.solve();
	std::cout << "The equation difference is: " << vector_max(solver.find_function_maximum(&Cylindric_Test::equations, 2)) << std::endl;
	std::cout << "The solution difference is: " << vector_max(solver.find_function_maximum(&Cylindric_Test::solution_check, 2, false)) << std::endl;
	cont.save(test_folder + "/colloc_compact");

	std::cout << "\nRegridding" << std::endl;

	bpart.move_domain_wall({0.2, 0.6}, 0.65);
	bpart.domain_split({0.25, 0.25}, 0);
	bpart.change_domain_resolution({0.7, 0.7}, {25, 35});
	bpart.print_domain_2D(test_folder + "/domains_2D_regrid.pbm");

	Full_Compactification comp_old(compact);
	compact.adjust_compactification(bpart);
	cont.regrid<Full_Compactification>(bpart, comp_old, compact);
	std::cout << "Number of points: " << cont.size() << std::endl;
	std::cout << "The equation difference is: " << vector_max(solver.find_function_maximum(&Cylindric_Test::equations, 2)) << std::endl;
	std::cout << "The solution difference is: " << vector_max(solver.find_function_maximum(&Cylindric_Test::solution_check, 2, false)) << '\n' << std::endl;

	std::cout << "Resolving" << std::endl;

	solver.solve();
	std::cout << "The equation difference is: " << vector_max(solver.find_function_maximum(&Cylindric_Test::equations, 2)) << std::endl;
	std::cout << "The solution difference is: " << vector_max(solver.find_function_maximum(&Cylindric_Test::solution_check, 2, false)) << std::endl;
	cont.save(test_folder + "/colloc_compact_regrid");

	auto square_function = [](int fn, double *x, double *phi, double **dphi, double **ddphi){return phi[0]*phi[0];};
	Evaluator eval_square(solver.evaluate_function(square_function, 1));
	std::cout << "The integral is: " << eval_square.integrate()[0] << std::endl;
	std::cout << "The true integral value is: " << 3.0*M_PI/(16.0*a_param*a_param*a_param) << std::endl;

	struct{
		double xmax;
		double ymax;
		double phi_max = 0.0;

		double operator()(int fn, double *x, double *phi, double **dphi, double **ddphi){
			if(phi[0] > phi_max){
				xmax = x[0];
				ymax = x[1];
				phi_max = phi[0];
			}
			return -phi[0];
		}
	} inv_function;

	Container func_min = solver.evaluate_function(inv_function, 1);
	func_min.set_doubling_begin(0, {false});
	func_min.set_doubling_begin(1, {false});
	std::cout << "Evaluating minimum" << std::endl;
	Evaluator eval_min(func_min);
	Collocation_Nelder_Mead minimizer(eval_min);
	minimizer.set_simplex({inv_function.xmax, inv_function.ymax});
	minimizer.minimize_track_csv(test_folder + "/colloc_compact_minimization.csv", true);
	std::cout << "Maximal function value: " << -minimizer.min_val << std::endl;
	std::cout << "Maximal value occured at: (" << minimizer.min_x[0] << ", " << minimizer.min_x[1] << ")" << std::endl;
}