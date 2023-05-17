#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

#include <collocation_tests.h>

struct{
	double operator()(int fi, double *x) const{
		if(fi == 0)
			return std::sin(M_PI*x[0])*std::sin(M_PI*x[1]);
		else
			return std::cos(2.0*M_PI*x[0])*std::cos(2.0*M_PI*x[1]);
	}
	double d_x(int fi, double *x) const{
		if(fi == 0)
			return M_PI*std::cos(M_PI*x[0])*std::sin(M_PI*x[1]);
		else
			return -2.0*M_PI*std::sin(2.0*M_PI*x[0])*std::cos(2.0*M_PI*x[1]);
	}
	double d_y(int fi, double *x) const{
		if(fi == 0)
			return M_PI*std::sin(M_PI*x[0])*std::cos(M_PI*x[1]);
		else
			return -2.0*M_PI*std::cos(2.0*M_PI*x[0])*std::sin(2.0*M_PI*x[1]);
	}
	double d_xy(int fi, double *x) const{
		if(fi == 0)
			return M_PI*M_PI*std::cos(M_PI*x[0])*std::cos(M_PI*x[1]);
		else
			return 4.0*M_PI*M_PI*std::sin(2.0*M_PI*x[0])*std::sin(2.0*M_PI*x[1]);
	}
} sin_func;

struct{
	double a = 50.0;

	double operator()(int fi, double *x) const{
		if(fi == 0)
			return x[1]*x[1]*std::exp(-0.5*a*(x[0]*x[0] + x[1]*x[1]));
		else
			return std::exp(-a*(x[0]*x[0] + x[1]*x[1]));
	}
	double d_x(int fi, double *x) const{
		if(fi == 0)
			return -a*x[0]*x[1]*x[1]*std::exp(-0.5*a*(x[0]*x[0] + x[1]*x[1]));
		else
			return -2.0*a*x[0]*std::exp(-a*(x[0]*x[0] + x[1]*x[1]));
	}
	double d_y(int fi, double *x) const{
		if(fi == 0)
			return x[1]*(2.0 - a*x[1]*x[1])*std::exp(-0.5*a*(x[0]*x[0] + x[1]*x[1]));
		else
			return -2.0*a*x[1]*std::exp(-a*(x[0]*x[0] + x[1]*x[1]));
	}
	double d_xy(int fi, double *x) const{
		if(fi == 0)
			return a*x[0]*x[1]*(a*x[1]*x[1] - 2.0)*std::exp(-0.5*a*(x[0]*x[0] + x[1]*x[1]));
		else
			return 4*a*a*x[0]*x[1]*std::exp(-a*(x[0]*x[0] + x[1]*x[1]));
	}
} gauss_profile;

struct Differentiation_Matrix_Test{
	Container &cont;

	Differentiation_Matrix_Test(Container &cont_): cont(cont_) {};

	template<class Callable>
	void check_differentiation(const Callable &func, const std::string &save_name){
		std::cout << "Differentiation matrix size: " << cont.size() << std::endl;
		cont.fill_function(func);

		Container cont_check(cont);
		Collocation_Matrix col_matr(cont);

		col_matr.fill_Dn_matrix(-1, 0);
		cont.copy_to_vector(col_matr, col_matr.rhs);
		col_matr.delta = col_matr.matr*col_matr.rhs;
		cont_check.fill_function([func](int fn, double *x){return func.d_x(fn, x);});
		cont_check.subtract_delta(col_matr);
		std::cout << "Maximal deviation x derivative: " << cont_check.max_abs() << std::endl;

		col_matr.fill_Dn_matrix(-1, 1);
		cont.copy_to_vector(col_matr, col_matr.rhs);
		col_matr.delta = col_matr.matr*col_matr.rhs;
		cont_check.fill_function([func](int fn, double *x){return func.d_y(fn, x);});
		cont_check.subtract_delta(col_matr);
		std::cout << "Maximal deviation y derivative: " << cont_check.max_abs() << std::endl;

		col_matr.fill_Dn_matrix(0, 1);
		cont.copy_to_vector(col_matr, col_matr.rhs);
		col_matr.delta = col_matr.matr*col_matr.rhs;
		cont_check.fill_function([func](int fn, double *x){return func.d_xy(fn, x);});
		cont_check.subtract_delta(col_matr);
		std::cout << "Maximal deviation xy derivative: " << cont_check.max_abs() << std::endl;

		if(!save_name.empty())
			cont_check.save(save_name);
	}
};

void collocation_doubling_test(int ndim){
	// On each dimension 5 possible combinations:
	// 0 - no doubling
	// 1 - even doubling from the left
	// 2 - odd doubling from the left
	// 3 - even doubling from the right
	// 4 - odd doubling from the right
	struct Doubling_Check{
		std::vector<int> doubling_condition;
		Doubling_Check(const std::vector<int> &condition): doubling_condition(condition) {};

		double operator()(int fn, double *x) const{
			double res = 1.0;
			for(int dim = 0; dim < doubling_condition.size(); ++dim){
				switch(doubling_condition[dim]){
					case 0:
						res *= std::exp(x[dim]);
						break;
					case 1:
						res *= x[dim]*std::sin(M_PI*x[dim]);
						break;
					case 2:
						res *= x[dim]*std::cos(M_PI*x[dim]);
						break;
					case 3:
						res *= (1.0-x[dim])*std::sin(M_PI*x[dim]);
						break;
					case 4:
						res *= (1.0-x[dim])*std::cos(M_PI*x[dim]);
						break;
				}
			}
			return res;
		}
	};

	std::cout << "Running comprehensive domain doubling check" << std::endl;
	std::vector<double> start(ndim, 0.0), lengths(ndim, 1.0);
	std::vector<int> sizes(ndim, 20), condition(ndim, 0), max_discrepancy;
	double tmp, max = 0.0;
	while(1){
		std::cout << "Testing condition ";
		for(int el : condition)
			std::cout << el;
		std::cout << std::endl;
		Domain_Partitioning part_test(start, lengths, sizes);
		Container tcont(1, part_test);
		for(int dim = 0; dim < condition.size(); ++dim){
			switch(condition[dim]){
				case 1:
					tcont.set_doubling_begin(dim, {false});
					break;
				case 2:
					tcont.set_doubling_begin(dim, {true});
					break;
				case 3:
					tcont.set_doubling_end(dim, {false});
					break;
				case 4:
					tcont.set_doubling_end(dim, {true});
					break;
			}
		}
		Doubling_Check check_function(condition);
		tcont.fill_function(check_function);
		Evaluator eval(tcont);
		tmp = eval.test(tcont);
		if(tmp > max){
			max = tmp;
			max_discrepancy = condition;
		}
		std::cout << std::endl;

		int dim;
		for(dim = condition.size()-1; dim >= 0; --dim){
			++condition[dim];
			if(condition[dim] < 5) break;
			condition[dim] = 0;
		}
		if(dim < 0)
			break;
	}

	std::cout << "Maximum discrepancy " << max << " found for configuration ";
	for(int el : max_discrepancy)
		std::cout << el;
	std::cout << std::endl;
}

template<class Callable>
void point_class_test(Container &cont, const Callable &func){
	cont.fill_function(func);

	Point pnt(&cont);

	double tmp;
	std::vector<double> max(4, 0.0);

	for(const Bulk_Domain *dom = cont.partition.first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
		std::vector<unsigned int> ind(cont.partition.n_dim, 0);
		for(int indx = 0; indx < dom->size; indx++){
			for(int i = cont.partition.n_dim-1, tmpi = indx; i >= 0; i--){
				if((ind[i] = tmpi % dom->dim_num(i)) > 1) break;
				tmpi /= dom->dim_num(i);
			}

			pnt.get_point(dom->bulk_domain_number, ind, true);

			for(int fn = 0; fn < cont.number_of_functions(); ++fn){
				tmp = std::fabs(func(fn, pnt.x) - pnt.phi[fn]);
				if(tmp > max[0]) max[0] = tmp;
				tmp = std::fabs(func.d_x(fn, pnt.x) - pnt.dphi[0][fn]);
				if(tmp > max[1]) max[1] = tmp;
				tmp = std::fabs(func.d_y(fn, pnt.x) - pnt.dphi[1][fn]);
				if(tmp > max[2]) max[2] = tmp;
				tmp = std::fabs(func.d_xy(fn, pnt.x) - pnt.ddphi[2][fn]);
				if(tmp > max[3]) max[3] = tmp;
			}
		}
	}

	std::cout << "Maximal function deviation: " << max[0] << std::endl;
	std::cout << "Maximal x derivative deviation: " << max[1] << std::endl;
	std::cout << "Maximal y derivative deviation: " << max[2] << std::endl;
	std::cout << "Maximal xy derivative deviation: " << max[3] << std::endl;
}

void collocation_2D_container_test(const std::string &test_folder, int argc, char const *argv[]){
	Domain_Partitioning part_sin("Tests/diff_domains_2D.dat");
	// Domain_Partitioning part_sin("Tests/domains_2D_test.dat");
	part_sin.print_domain_2D(test_folder + "/diff_domains_sin.pbm");
	
	Container tcont(2, part_sin);
	Differentiation_Matrix_Test dmat_test(tcont);

	std::cout << "Testing differentiation function on sin function" << std::endl;
	tcont.set_doubling_begin(0, {true, false});
	tcont.set_doubling_begin(1, {true, false});
	tcont.set_doubling_end(0, {true, false});
	tcont.set_doubling_end(1, {true, false});
	dmat_test.check_differentiation(sin_func, test_folder + "/diff_test_sin");
	std::cout << std::endl;

	tcont.save(test_folder + "/cont_test");

	std::cout << "Testing differentiation function on Gaussian profile function" << std::endl;
	Domain_Partitioning part_func("Tests/simple_domains_2D.dat");
	tcont.adjust_partitioning(part_func);
	tcont.reset_doubling();
	tcont.set_doubling_begin(0, {false, false});
	tcont.set_doubling_begin(1, {false, false});
	dmat_test.check_differentiation(gauss_profile, test_folder + "/diff_test_gauss");
	std::cout << std::endl;

	std::cout << "Evaluator test" << std::endl;
	collocation_doubling_test(2);

	std::cout << "\nRunning evaluator test on gaussian profile function" << std::endl;
	Evaluator eval(tcont);
	eval.test(tcont);
	std::cout << std::endl;

	std::cout << "Testing point class on Gaussian profile" << std::endl;
	point_class_test(tcont, gauss_profile);
	std::cout << std::endl;

	std::cout << "Testing point class on sin function" << std::endl;
	Domain_Partitioning part_new("Tests/domains_2D_test.dat");
	tcont.adjust_partitioning(part_new);
	tcont.reset_doubling();
	tcont.set_doubling_begin(0, {true, false});
	tcont.set_doubling_begin(1, {true, false});
	tcont.set_doubling_end(0, {true, false});
	tcont.set_doubling_end(1, {true, false});
	point_class_test(tcont, sin_func);
	std::cout << std::endl;

	std::cout << "Testing regridding" << std::endl;
	No_Compactification no_compact;

	tcont.adjust_partitioning(part_sin);
	tcont.reset_doubling();
	tcont.set_doubling_begin(0, {true, false});
	tcont.set_doubling_begin(1, {true, false});
	tcont.set_doubling_end(0, {true, false});
	tcont.set_doubling_end(1, {true, false});
	tcont.fill_function(sin_func);

	part_new.print_domain_2D(test_folder + "/domains_2D_test.pbm");
	tcont.regrid(part_new, no_compact, no_compact);
	std::cout << "The difference is: " << tcont.compare_function(sin_func) << std::endl;
	tcont.save(test_folder + "/cont_test_regrid");

	part_new.domain_split({0.6, 0.6}, 0);
	part_new.print_domain_2D(test_folder + "/domains_2D_test_split.pbm");

	tcont.regrid(part_new, no_compact, no_compact);
	std::cout << "The difference is: " << tcont.compare_function(sin_func) << std::endl;
	tcont.save(test_folder + "/cont_test_regrid_2");

	part_new.move_domain_wall({0.8, 0.5}, 0.7);
	part_new.print_domain_2D(test_folder + "/domains_2D_test_split_move.pbm");

	tcont.regrid(part_new, no_compact, no_compact);
	std::cout << "The difference is: " << tcont.compare_function(sin_func) << std::endl;
	tcont.save(test_folder + "/cont_test_regrid_3");

	return;
}

static double vector_max(const std::vector<double> &vec){
	double max = vec[0];
	for(int i = 1; i < vec.size(); ++i)
		if(vec[i] > max) max = vec[i];
	return max;
}

void collocation_2D_collocation_fill_test(const std::string &test_folder, int argc, char const *argv[]){
	// Domain_Partitioning bpart("Tests/simple_domains_2D.dat");
	Domain_Partitioning bpart("Tests/domains_2D_test.dat");

	Container cont(2, bpart);
	cont.set_doubling_begin(0, {true, false});
	cont.set_doubling_begin(1, {true, false});
	cont.set_doubling_end(0, {true, false});
	cont.set_doubling_end(1, {true, false});

	No_Compactification no_compact;
	Test_Equations test_eqn;
	Collocation_Solver_Full<No_Compactification, Test_Equations> solver(no_compact, &test_eqn, &cont);

	std::cout << "Filling collocation matrix" << std::endl;
	solver.collocation_fill<true, true>();
	std::cout << "Filled" << std::endl;
	std::cout << "Number of points: " << cont.size() << std::endl;

	std::cout << "Computing delta vector" << std::endl;
	solver.colm.compute_delta();
	cont.subtract_delta(solver.colm);
	std::cout << "Solver step taken" << std::endl;

	solver.collocation_fill<true, true>();
	solver.colm.compute_delta();
	cont.subtract_delta(solver.colm);

	cont.save(test_folder + "/colloc_simple");

	std::cout << "The difference is: " << vector_max(solver.find_function_maximum(&Test_Equations::solution_check, 2)) << std::endl;

	// std::cout << "The integral is: " << solver.integrate_function(pure_field_square) << std::endl;
	// std::cout << "The true integral value is: 0.25" << std::endl;

	// solver.fill_seed(sin_func);
	// Evaluator eval = solver.get_function_evaluator([](double *x, double *phi, double **dphi, double **ddphi){return -phi[0];});

	// Container cont_test(1, bpart);
	// cont_test.fill_function([](int fi, double *x){return -sin_func(fi, x);});
	// eval.test(cont_test);

	// Collocation_Nelder_Mead minimizer(eval);
	// minimizer.safe_minimization();
	// std::cout << "Maximal function value: " << -minimizer.min_val << std::endl;
	// std::cout << "Maximal value occured at: (" << minimizer.min_x[0] << ", " << minimizer.min_x[1] << ")" << std::endl;

	// int n_test = 30;
	// Evaluator_1D_Spec eval1d(n_test);

	// std::cout << "Incomplete integral test along dimension 0" << std::endl;
	// fill_integral<No_Compactification>(
	// 	eval1d,
	// 	solver,
	// 	[](double *x, double *phi, double **dphi, double **ddphi){return phi[0];},
	// 	0
	// );

	// eval1d.compare_test([](double x){return 2.0/M_PI/M_PI*(1.0 - std::cos(M_PI*x));});

	// std::cout << "Incomplete integral test along dimension 1" << std::endl;
	// fill_integral<No_Compactification>(
	// 	eval1d,
	// 	solver,
	// 	[](double *x, double *phi, double **dphi, double **ddphi){return phi[0];},
	// 	1
	// );

	// eval1d.compare_test([](double x){return 2.0/M_PI/M_PI*(1.0 - std::cos(M_PI*x));});
}


#include <test_equations_3D.h>

double sin_func_3D(int fi, double *x){
	if(fi == 0)
		return std::sin(M_PI*x[0])*std::sin(M_PI*x[1])*std::sin(M_PI*x[2]);
	else
		throw("Wrong function number");
}

void collocation_3D_container_test(const std::string &test_folder, int argc, char const *argv[]){
	Domain_Partitioning bpart("Tests/simple_domains_3D.dat");

	Container cont(1, bpart);
	cont.set_doubling_begin(0, {true});
	cont.set_doubling_begin(1, {true});
	cont.set_doubling_begin(2, {true});
	cont.set_doubling_end(0, {true});
	// cont.set_doubling_end(1, {true});
	// cont.set_doubling_end(2, {true});

	cont.partition.print_domain_3D(test_folder + "/domains_3D_test");

	cont.fill_function(sin_func_3D);
	cont.save(test_folder + "/cont_test_3D");

	std::cout << "Number of points: " << cont.size() << std::endl;
	std::cout << "Check value, should be zero: " << cont.compare_function(sin_func_3D) << std::endl;

	std::cout << "Evaluator test" << std::endl;
	collocation_doubling_test(3);
	std::cout << std::endl;

	Evaluator eval(cont);
	eval.test(cont);
	std::cout << std::endl;

	std::cout << "Testing regridding" << std::endl;
	No_Compactification no_compact;
	bpart.domain_split({0.2, 0.2, 0.2}, 0);

	cont.regrid(bpart, no_compact, no_compact);
	std::cout << "The difference is: " << cont.compare_function(sin_func_3D) << std::endl;

	return;
}

void collocation_3D_collocation_fill_test(const std::string &test_folder, int argc, char const *argv[]){
	// Domain_Partitioning bpart({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, {10, 10, 10});
	Domain_Partitioning bpart("Tests/simple_domains_3D.dat");
	bpart.save_domains(test_folder + "/domains.dat");

	Container cont(1, bpart);
	cont.set_doubling_begin(0, {false});
	cont.set_doubling_begin(1, {true});
	cont.set_doubling_begin(2, {false});
	cont.set_doubling_end(0, {false});
	// cont.set_doubling_end(1, {true});
	// cont.set_doubling_end(2, {false});

	No_Compactification no_compact;
	Test_Equations_3D test_eqn;
	Collocation_Solver_Full<No_Compactification, Test_Equations_3D> solver(no_compact, &test_eqn, &cont);

	std::cout << "Filling collocation matrix" << std::endl;
	solver.collocation_fill<true, true>();
	std::cout << "Filled" << std::endl;
	std::cout << "Number of points: " << cont.size() << std::endl;

	std::cout << "Computing delta vector" << std::endl;
	solver.colm.compute_delta();
	cont.subtract_delta(solver.colm);
	std::cout << "Solver step taken" << std::endl;

	solver.collocation_fill<true, true>();
	solver.colm.compute_delta();
	cont.subtract_delta(solver.colm);

	cont.save(test_folder + "/cont_test_3D");

	std::cout << "The difference is: " << vector_max(solver.find_function_maximum(&Test_Equations_3D::solution_check)) << std::endl;
	// std::cout << "The integral is: " << solver.integrate_function(pure_field_square) << std::endl;
	// std::cout << "The true integral value is: 0.125" << std::endl;

	// // solver.fill_seed([](int fn, double *x){
	// // 	return std::cos(M_PI*x[0])*std::sin(M_PI*x[1])*std::cos(M_PI*x[2]);
	// // });

	// Evaluator eval = solver.get_function_evaluator([](double *x, double *phi, double **dphi, double **ddphi){return -phi[0];});
	// Collocation_Nelder_Mead minimizer(eval);
	// minimizer.safe_minimization();
	// std::cout << "Maximal function value: " << -minimizer.min_val << std::endl;
	// std::cout << "Maximal value occured at: (" << minimizer.min_x[0] << ", " << minimizer.min_x[1] << ", " << minimizer.min_x[2] << ")" << std::endl;

	// int int_dim;
	// int n_test = 30;

	// Evaluator_1D_Spec eval1d(n_test);

	// //Dimension 0
	// fill_integral<No_Compactification>(
	// 	eval1d,
	// 	solver,
	// 	[](double *x, double *phi, double **dphi, double **ddphi){return phi[0]*phi[0];},
	// 	0
	// );

	// std::cout << "Incomplete integral test along dimension 0" << std::endl;
	// eval1d.compare_test([](double x){return 0.125*(x + std::sin(2.0*M_PI*x)/(2.0*M_PI));});

	// //Dimension 1
	// fill_integral<No_Compactification>(
	// 	eval1d,
	// 	solver,
	// 	[](double *x, double *phi, double **dphi, double **ddphi){return phi[0]*phi[0];},
	// 	1
	// );

	// std::cout << "Incomplete integral test along dimension 1" << std::endl;
	// eval1d.compare_test([](double x){return 0.125*(x - std::sin(2.0*M_PI*x)/(2.0*M_PI));});

	// //Dimension 2
	// fill_integral<No_Compactification>(
	// 	eval1d,
	// 	solver,
	// 	[](double *x, double *phi, double **dphi, double **ddphi){return phi[0]*phi[0];},
	// 	2
	// );

	// std::cout << "Incomplete integral test along dimension 2" << std::endl;
	// eval1d.compare_test([](double x){return 0.125*(x + std::sin(2.0*M_PI*x)/(2.0*M_PI));});
}