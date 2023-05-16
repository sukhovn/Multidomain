#include <cmath>
#include <limits>


class Test_Equations_3D{
public:
	double solution_check(int fn, double *x, double *phi, double **dphi, double **ddphi) const{
		switch(fn){
			case 0:
				return phi[0] - std::cos(M_PI*x[0])*std::sin(M_PI*x[1])*std::cos(M_PI*x[2]);
				break;
			default:
				std::cout << "Wrong equation number" << std::endl;
				throw("Wrong equation number");
		}
	}

	double operator()(double *x, int fn) const{
		double phi = 0.0;
		return -solution_check(fn, x, &phi, NULL, NULL);
	}

	double equations(int fn, double *x, double *phi, double **dphi, double **ddphi) const{
		switch(fn){
			case 0:
				return ddphi[0][0] + ddphi[1][0] + ddphi[2][0] - phi[0] + (1.0 + 3.0 * M_PI * M_PI)*std::cos(M_PI*x[0])*std::sin(M_PI*x[1])*std::cos(M_PI*x[2]);
				break;
			default:
				std::cout << "Wrong equation number" << std::endl;
				throw("Wrong equation number");
		}
	}

	double boundary_condition_left(int dm, int fn, double *phi, double *dphi) const{
		switch(fn){
			case 0:
				if(dm == 0) return dphi[0];
				else if(dm == 1) return phi[0];
				else if(dm == 2) return dphi[0];
				break;
		}
		std::cout << "Wrong equation number" << std::endl;
		throw("Wrong equation number");
	}

	double boundary_condition_right(int dm, int fn, double *phi, double *dphi) const{
		switch(fn){
			case 0:
				if(dm == 0) return dphi[0];
				else if(dm == 1) return phi[0];
				else if(dm == 2) return dphi[0];
				break;
		}
		std::cout << "Wrong equation number" << std::endl;
		throw("Wrong equation number");
	}
};