#include <fftw3.h>
#include <container.h>
#include <point_class.h>

#ifndef COLLOCATION_SOLUTION_H
#define COLLOCATION_SOLUTION_H

template<class Compactification>
class Collocation_Solution: public Domain_Partitioning{
protected:
	Point_Compactified<Compactification> pnt;
	Container *cont;

public:
	void fill_boundary_domains(Domain_Partitioning &part){
		std::unordered_set<int> exclude_left;
		std::unordered_set<int> exclude_right;
		for(int dim = 0; dim < cont->n_dim; ++dim){
			if(!cont->parity_st[dim].empty()) exclude_left.insert(dim);
			if(!cont->parity_end[dim].empty()) exclude_right.insert(dim);
		}
		part.insert_boundary_domains(exclude_left, exclude_right);

		auto boundary_function = [this](int dim, double xst, double xend, int n, double x) -> int{
			if(!cont->parity_st[dim].empty() && std::fabs(xst - cont->partition_start(dim)) < COMP_TOL){
				return ch_x_inv_low(2*xst-xend, xend, 2*n-1, x)-n;
			}
			else if(!cont->parity_end[dim].empty() && std::fabs(xend -  cont->partition_end(dim)) < COMP_TOL){
				return ch_x_inv_low(xst, 2*xend-xst, 2*n-1, x);
			}
			else{
				return ch_x_inv_low(xst, xend, n-1, x);
			}
		};
		part.fill_boundary_point_range(boundary_function, exclude_left, exclude_right);
	}

public:
	const Domain_Partitioning &partition(){return *this;}

	Collocation_Solution(const Compactification &compact, Container *conti):
		Domain_Partitioning(*conti), cont(conti), pnt(compact, conti) {};

	bool verify_container_solution_correspondence(void){
		if(compare_bulk_domains(*cont)) return false;
		Domain_Partitioning part(*cont);
		if(boundary_domains_present()){
			fill_boundary_domains(*this);
		}
		this->Domain_Partitioning::operator=(part);
		return true;
	}

	void change_container(Container &conti){
		cont = &conti;
		pnt.change_container(cont);
		verify_container_solution_correspondence();
	}

	// template<class Callable>
	// void fill_seed(Callable func, double value_at_infinity = 0.0){
	// 	verify_container_solution_correspondence();
	// 	unsigned int itmp, itr = 0;
	// 	std::vector<unsigned int> ind(n_dim);
	// 	double *x_c = new double[n_dim];
	// 	double *x_tmp = new double[n_dim];

	// 	for(const Bulk_Domain *bulk = first_bulk_domain; bulk != NULL; bulk = bulk->next_bulk_domain()){
	// 		int dm = bulk->bulk_domain_number;

	// 		for(int indx = 0; indx < bulk->size; indx++){
	// 			for(int i = n_dim-1, tmpi = indx; i >= 0; i--){
	// 				itmp = tmpi % bulk->dim_num(i);
	// 				x_c[i] = ch_x(bulk->domain_start(i), bulk->domain_end(i), bulk->dim_num(i)-1, itmp);
	// 				if(itmp) break;
	// 				tmpi /= bulk->dim_num(i);
	// 			}

	// 			bool is_inf = false;
	// 			for(int i = 0; i < n_dim; i++){
	// 				x_tmp[i] = pnt.compact.xc(i, x_c);
	// 				if(std::isinf(x_tmp[i])){
	// 					is_inf = true;
	// 					break;
	// 				}
	// 			}

	// 			if(is_inf){
	// 				for(int fn = 0; fn < cont->nfunc; fn++) cont->data[itr++] = value_at_infinity;
	// 				continue;
	// 			}
	// 			for(int fn = 0; fn < cont->nfunc; fn++)
	// 				cont->data[itr++] = func(fn, x_tmp);
	// 		}
	// 	}

	// 	delete [] x_c;
	// 	delete [] x_tmp;
	// }
};

#endif