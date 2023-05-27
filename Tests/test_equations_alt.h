#include <cmath>
#include <limits>
#include <iostream>
#include <vector>
#include <container.h>

class Single_Compactification{
private:
	int dim;
	double xst;
	double cc;

	std::vector<bool> domains_compactified;
public:
	Single_Compactification(const Container &cont, int _dim = 0, bool print_selection = false): dim(_dim), xst(cont.partition.partition_start(_dim)){
		int ndomains = 0;
		double xend = cont.partition.partition_end(dim);
		for(const Bulk_Domain *dom = cont.partition.first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			if(dom->domain_start(dim) > xst) xst = dom->domain_start(dim);
			++ndomains;
			domains_compactified.push_back(false);
		}
		if(print_selection){
			std::cout << "Selected compactification region:" << std::endl;
			std::cout << "Dimension " << dim << ": [" << xst << ',' << xend << ']' << std::endl;
		}
		cc = 1.0/(xend - xst);
		for(const Bulk_Domain *dom = cont.partition.first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			if(dom->domain_end(dim) > xst + COMP_TOL)
				domains_compactified[dom->bulk_domain_number] = true;
		}
	}

	double xc(int dm, double *x) const{
		if(dm == dim){
			double xshift = x[dim] - xst;
			if(xshift <= 0.0) return x[dim];
			if(xshift*cc >= 1.0) return std::numeric_limits<double>::infinity();
			return xst + xshift/(1 - xshift*cc);
		}
		else
			return x[dm];
	}

	double dxc(int dnum, int dm, double *x) const{
		if(!domains_compactified[dnum])
			return 1.0;
		if(dm == dim){
			double xshift = x[dim] - xst;
			if(xshift*cc >= 1.0) return std::numeric_limits<double>::infinity();
			return 1.0/std::pow(1 - xshift*cc, 2.0);
		}
		else
			return 1.0;
	}

	double ddxc(int dnum, int dm, double *x) const{
		if(!domains_compactified[dnum])
			return 0.0;
		if(dm == dim){
			double xshift = x[dim] - xst;
			if(xshift*cc == 1.0) return std::numeric_limits<double>::infinity();
			return 2.0 * cc/std::pow(1 - xshift*cc, 3.0);
		}
		else
			return 0.0;
	}

	double xc_inv(int dm, double *x) const{
		if(dm == dim){
			double xshift = x[dim] - xst;
			if(xshift <= 0.0) return x[dm];
			return xst + xshift/(1 + xshift*cc);
		}
		else
			return x[dm];
	}
};


class Full_Compactification{
private:
	std::vector<double> xst;
	std::vector<double> cc;

public:
	void adjust_compactification(const Domain_Partitioning &part, bool print_selection = false){
		if(part.n_dim != cc.size()){
			std::cout << "Changing number of dimensions while adjusting partitioning is not permitted" << std::endl;
			throw("Changing dimensions prohibited");
		}
		xst = part.partition_start();
		std::vector<double> xend(part.partition_end());
		for(const Bulk_Domain *dom = part.first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			for(int dim = 0; dim < part.n_dim; ++dim)
				if(dom->domain_start(dim) > xst[dim]) xst[dim] = dom->domain_start(dim);
		}
		if(print_selection){
			std::cout << "Selected compactification region:" << std::endl;
			for(int i = 0; i < xst.size(); ++i){
				std::cout << "Dimension " << i << ": [" << xst[i] << ',' << xend[i] << ']' << std::endl;
			}
		}
		for(int i = 0; i < xst.size(); ++i)
			cc[i] = 1.0/(xend[i] - xst[i]);
	}

	// Automatically selects compactification region to match partition's outmost domains
	Full_Compactification(const Domain_Partitioning &part, bool print_selection = false): cc(part.n_dim){
		adjust_compactification(part, print_selection);
	}
	Full_Compactification(const Full_Compactification &compact): xst(compact.xst), cc(compact.cc) {};

	double xc(int dm, double *x) const{
		double xshift = x[dm] - xst[dm];
		if(xshift <= 0.0) return x[dm];
		if(xshift*cc[dm] >= 1.0) return std::numeric_limits<double>::infinity();
		return xst[dm] + xshift/(1 - xshift*cc[dm]);
	}
	double dxc(int dnum, int dm, double *x) const{
		double xshift = x[dm] - xst[dm];
		if(xshift <= 0.0) return 1.0;
		if(xshift*cc[dm] >= 1.0) return std::numeric_limits<double>::infinity();
		return 1.0/std::pow(1 - xshift*cc[dm], 2.0);
	}
	double ddxc(int dnum, int dm, double *x) const{
		double xshift = x[dm] - xst[dm];
		if(xshift <= 0.0) return 0.0;
		if(xshift*cc[dm] >= 1.0) return std::numeric_limits<double>::infinity();
		return 2.0 * cc[dm]/std::pow(1 - xshift*cc[dm], 3.0);
	}
	double xc_inv(int dm, double *x) const{
		double xshift = x[dm] - xst[dm];
		if(xshift <= 0.0) return x[dm];
		return xst[dm] + xshift/(1 + xshift*cc[dm]);
	}
};


class Cylindric_Test{
	double a;

	double corr_func(double x, double y, int fn) const{
		switch(fn){
			case 0:
				return y*y*std::exp(- 1.5*a*(x*x + y*y)) + std::exp(- 0.5*a*(x*x + y*y))*(4.0 - 7.0*a*y*y + a*a*y*y*(x*x + y*y));
				break;
			case 1:
				return std::exp(- a*(x*x + y*y))*(-6.0*a + y*y*y*y + 4.0*a*a*(x*x + y*y));
				break;
			default:
				std::cout << "Wrong equation number" << std::endl;
				throw("Wrong equation number");
		}
	}

public:
	Cylindric_Test(double ai): a(ai) {};

	double solution_check(int fn, double *x, double *phi, double **dphi, double **ddphi) const{
		switch(fn){
			case 0:
				return phi[0] - x[1]*x[1]*std::exp(- 0.5*a*(x[0]*x[0] + x[1]*x[1]));
				break;
			case 1:
				return phi[1] - std::exp(- a*(x[0]*x[0] + x[1]*x[1]));
				break;
			default:
				std::cout << "Wrong equation number" << std::endl;
				throw("Wrong equation number");
		}
	}

	double equations(int fn, double *x, double *phi, double **dphi, double **ddphi) const{
		switch(fn){
			case 0:
				return ddphi[1][0] + dphi[1][0]/x[1] + ddphi[0][0] + phi[1]*phi[0] - corr_func(x[0], x[1], 0);
				break;
			case 1:
				return ddphi[1][1] + dphi[1][1]/x[1] + ddphi[0][1] + phi[0]*phi[0] - corr_func(x[0], x[1], 1);
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
				return phi[0];
				break;
			case 1:
				return phi[1];
				break;
			default:
				std::cout << "Wrong equation number" << std::endl;
				throw("Wrong equation number");
		}
	}
};
