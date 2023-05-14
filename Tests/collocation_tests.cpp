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

	// collocation_doubling_test(2);
	// std::cout << std::endl;
	// collocation_doubling_test(3);

	Evaluator eval(tcont);
	eval.test(tcont);

	return;
}