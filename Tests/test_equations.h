#include <cmath>
#include <limits>

class No_Compactification{
public:
	double xc(int dm, double *x) const {return x[dm];}
	double dxc(int dnum, int dm, double *x) const {return 1.0;}
	double ddxc(int dnum, int dm, double *x)const {return 0.0;}

	double xc_inv(int dm, double *x) const{return x[dm];}
};


class Test_Equations{
public:
	double solution_check(int fn, double *x, double *phi, double **dphi, double **ddphi) const{
		switch(fn){
			case 0:
				return phi[0] - std::sin(M_PI*x[0])*std::sin(M_PI*x[1]);
				break;
			case 1:
				return phi[1] - std::cos(2.0*M_PI*x[0])*std::cos(2.0*M_PI*x[1]);
				break;
			default:
				std::cout << "Wrong equation number" << std::endl;
				throw("Wrong equation number");
		}
	}

	double equations(int fn, double *x, double *phi, double **dphi, double **ddphi) const{
		switch(fn){
			case 0:
				return ddphi[0][0] + ddphi[1][0] - phi[0] + (1.0 + 2.0 * M_PI * M_PI)*std::sin(M_PI*x[0])*std::sin(M_PI*x[1]);
				break;
			case 1:
				return ddphi[0][1] + ddphi[1][1] - phi[1] + (1.0 + 8.0 * M_PI * M_PI)*std::cos(2.0*M_PI*x[0])*std::cos(2.0*M_PI*x[1]);
				break;
			default:
				std::cout << "Wrong equation number" << std::endl;
				throw("Wrong equation number");
		}
	}

	double boundary_condition_left(int dm, int fn, double *phi, double *dphi) const{
		switch(fn){
			case 0:
				if(dm == 0) return phi[0];
				else if(dm == 1) return phi[0];
				break;
			case 1:
				if(dm == 0) return dphi[1];
				else if(dm == 1) return dphi[1];
				break;
		}
		std::cout << "Wrong equation number" << std::endl;
		throw("Wrong equation number");
	}

	double boundary_condition_right(int dm, int fn, double *phi, double *dphi) const{
		switch(fn){
			case 0:
				if(dm == 0) return phi[0];
				else if(dm == 1) return phi[0];
				break;
			case 1:
				if(dm == 0) return dphi[1];
				else if(dm == 1) return dphi[1];
				break;
		}
		std::cout << "Wrong equation number" << std::endl;
		throw("Wrong equation number");
	}
};