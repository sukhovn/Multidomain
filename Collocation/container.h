#include <vector>
#include <boundaries.h>
#include <spectral.h>

#define EIGEN_SUPERLU_SUPPORT
#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE

#define LAPACK_COMPLEX_CUSTOM
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

#include <eigen3/Eigen/Dense>

#ifndef CONTAINER_H
#define CONTAINER_H

int count_func_files(const std::string &folder_name);
Domain_Partitioning get_partitioning_from_data(const std::string &folder_name);


class Indexer{
protected:
	int ndim;
	int nfunc;
	int ndomains;
	std::vector<int> index_st;
	std::vector<std::vector<int>> dim_nums;
	int npoints;

	std::vector<std::vector<double>> doms_start;
	std::vector<std::vector<double>> doms_end;

public:
	Indexer(const Domain_Partitioning &part, int nfunci): ndim(part.n_dim), nfunc(nfunci), npoints(0), ndomains(0){
		for(const Bulk_Domain *dom = part.first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			dim_nums.push_back(dom->dim_num());
			index_st.push_back(npoints);
			npoints += nfunc*dom->size;
			doms_start.push_back(dom->domain_start());
			doms_end.push_back(dom->domain_end());
			if(dom->bulk_domain_number != ndomains++)
				std::cout << "Warning! Indexer numbering and partition numbering deviate" << std::endl;
		}
	}

	int number_of_functions() const {return nfunc;}
	int size() const {return npoints;};

	int index(unsigned int dm, const std::vector<unsigned int> &ind, unsigned int fn) const{
		if(dm >= ndomains){
			std::cout << "Wrong domain number" << std::endl;
			throw("Wrong domain number");
		}
		for(int i = 0; i < ndim; i++){
			if(ind[i] >= dim_nums[dm][i]){
				std::cout << "Wrong point number" << std::endl;
				throw("Wrong point number");
			}
		}
		if(fn >= nfunc){
			std::cout << "Wrong function number" << std::endl;
			throw("Wrong function number");
		}

		int posit = ind[0];
		for(int i = 1; i < ndim; i++){
			posit = dim_nums[dm][i]*posit + ind[i];
		}

		return index_st[dm] + nfunc*posit + fn;
	}

	void indexer_regrid(const Domain_Partitioning &part){
		if(part.n_dim != ndim){
			std::cout << "Regrid cannot change dimensions" << std::endl;
			throw("Regrid dimension change");
		}

		npoints = 0;
		ndomains = 0;
		index_st.resize(0);
		dim_nums.resize(0);
		doms_start.resize(0);
		doms_end.resize(0);
		for(const Bulk_Domain *dom = part.first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			dim_nums.push_back(dom->dim_num());
			index_st.push_back(npoints);
			npoints += nfunc*dom->size;
			doms_start.push_back(dom->domain_start());
			doms_end.push_back(dom->domain_end());
			if(dom->bulk_domain_number != ndomains++)
				std::cout << "Warning! Indexer numbering and partition numbering deviate" << std::endl;
		}
	}
};


class Collocation_Matrix: public Indexer{
public:
	Eigen::MatrixXd matr;
	Eigen::VectorXd rhs;
	Eigen::VectorXd delta;

public:
	Collocation_Matrix(const Domain_Partitioning &part, int nfunci): Indexer(part, nfunci)
	{
		matr.resize(npoints, npoints);
		rhs.resize(npoints);
		delta.resize(npoints);
		matr.setZero();
		rhs.setZero();
		delta.setZero();
	}

	//This functions zeroes the matrix
	void setZero_jacobian(void){
		matr.setZero();
	}

	void setZero_rhs(void){
		rhs.setZero();
	}

	//This function zeroes the matrix row
	void setZero_row(unsigned int dm, const std::vector<unsigned int> &ind, unsigned int fn){
		matr.row(index(dm, ind, fn)).setZero();
	}

	//This function acceses right hand side array
	double &right_hand_side(unsigned int dm, const std::vector<unsigned int> &ind, unsigned int fn){
		return rhs(index(dm, ind, fn));
	}

	//This function adds a differention operator to the (row_dm, row_ind) row, multiplying it by coef
	//der1 and der2 indicate derivative dimensions, der1 = -1 and/or der2 = -1 are used to indicate no differention
	//This is the main function used to construct the Jacobian
	void row_add_Dn_factor(unsigned int row_dm, const std::vector<unsigned int> &row_ind, unsigned int row_fn,
		int der1, int der2, unsigned int dm, const std::vector<unsigned int> &ind, unsigned int fn, double coef);


	//Computes the next step given the Jacobian and the equation difference
	void compute_delta(void){
		delta = matr.partialPivLu().solve(rhs);
	}

	// Changes domains configuration
	void regrid(const Domain_Partitioning &part){
		indexer_regrid(part);
		matr.resize(npoints, npoints);
		rhs.resize(npoints);
		delta.resize(npoints);
		matr.setZero();
		rhs.setZero();
		delta.setZero();
	}
};


// Calculates coordinate based on index in simple ordering, only works when started with ind = 0
inline void update_coordinate(const Bulk_Domain *dom, int ind, double *coords){
	int tmp;
	for(int dim_ch = dom->n_dim-1; dim_ch >= 0; dim_ch--){
		if(tmp = (ind % dom->dim_num(dim_ch))){
			coords[dim_ch] = ch_x(dom->domain_start(dim_ch), dom->domain_end(dim_ch), dom->dim_num(dim_ch)-1, tmp);
			break;
		}
		coords[dim_ch] = ch_x(dom->domain_start(dim_ch), dom->domain_end(dim_ch), dom->dim_num(dim_ch)-1, 0);
		ind /= dom->dim_num(dim_ch);
	}
}


class Container: protected Domain_Partitioning{
private:
	int nfunc;

	int npoints;
	double *data;
	std::vector<int> index_st;

	std::vector<std::vector<bool>> parity_st; //false - even, true - odd, empty - no doubling
	std::vector<std::vector<bool>> parity_end;

	void reallocate_partitioning(const Domain_Partitioning &part){
		this->Domain_Partitioning::operator=(part);
		int ndomains = 0, npoints_old = npoints;
		npoints = 0;
		index_st.resize(0);
		for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			index_st.push_back(npoints);
			npoints += nfunc*dom->size;
			if(dom->bulk_domain_number != ndomains++)
				std::cout << "Warning! Container numbering and partition numbering deviate" << std::endl;
		}
		if(npoints_old != npoints){
			delete [] data;
			data = new double[npoints];
		}
	}

public:
	const Domain_Partitioning &partition;

	Container(int nfunci, const Domain_Partitioning &parti):
		Domain_Partitioning(parti), partition(*this), nfunc(nfunci), npoints(0),
		parity_st(parti.n_dim), parity_end(parti.n_dim)
	{
		int ndomains = 0;
		for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			index_st.push_back(npoints);
			npoints += nfunc*dom->size;
			if(dom->bulk_domain_number != ndomains++)
				std::cout << "Warning! Container numbering and partition numbering deviate" << std::endl;
		}
		data = new double[npoints];
		for(int i = 0; i < npoints; i++) data[i] = 0.0;
	};

	Container(const std::string &folder_name):
		Domain_Partitioning(get_partitioning_from_data(folder_name)), partition(*this), npoints(0),
		nfunc(count_func_files(folder_name))
	{
		parity_st.resize(n_dim);
		parity_end.resize(n_dim);
		int ndomains = 0;
		for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			index_st.push_back(npoints);
			npoints += nfunc*dom->size;
			if(dom->bulk_domain_number != ndomains++)
				std::cout << "Warning! Container numbering and partition numbering deviate" << std::endl;
		}
		data = new double[npoints];
		load(folder_name);
	}

	Container(const Container &cont):
		Domain_Partitioning(cont), partition(*this), nfunc(cont.nfunc), npoints(cont.npoints), index_st(cont.index_st),
		parity_st(cont.parity_st), parity_end(cont.parity_end)
	{
		data = new double[npoints];
		for(int i = 0; i < npoints; i++) data[i] = cont.data[i];
	};

	void set_doubling_begin(unsigned int dim, const std::vector<bool> &parity){
		if(dim >= n_dim){
			std::cout << "Wrong dimension number" << std::endl;
			throw("Wrong domain number");
		}
		if(parity.size() != nfunc){
			std::cout << "Wrong number of functions" << std::endl;
			throw("Wrong number of functions");
		}
		if(!parity_end[dim].empty()){
			for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
				if(std::fabs(dom->domain_start(dim) - partition_start(dim)) < COMP_TOL && std::fabs(dom->domain_end(dim) - partition_end(dim)) < COMP_TOL){
					std::cout << "Cannot set doubling on both sized of the same domain" << std::endl;
					throw("Cannot set doubling on both sized of the same domain");
				}
			}
		}
		parity_st[dim] = parity;
	}

	void set_doubling_end(unsigned int dim, const std::vector<bool> &parity){
		if(dim >= n_dim){
			std::cout << "Wrong dimension number" << std::endl;
			throw("Wrong domain number");
		}
		if(parity.size() != nfunc){
			std::cout << "Wrong number of functions" << std::endl;
			throw("Wrong number of functions");
		}
		if(!parity_st[dim].empty()){
			for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
				if(std::fabs(dom->domain_start(dim) - partition_start(dim)) < COMP_TOL && std::fabs(dom->domain_end(dim) - partition_end(dim)) < COMP_TOL){
					std::cout << "Cannot set doubling on both sized of the same domain" << std::endl;
					throw("Cannot set doubling on both sized of the same domain");
				}
			}
		}
		parity_end[dim] = parity;
	}

	Container &operator=(const Container &cont){
		this->Domain_Partitioning::operator=(cont);
		nfunc = cont.nfunc;
		if(npoints != cont.npoints){
			delete [] data;
			npoints = cont.npoints;
			data = new double[npoints];
		}
		index_st = cont.index_st;
		parity_st.assign(cont.parity_st.begin(), cont.parity_st.end());
		parity_end.assign(cont.parity_end.begin(), cont.parity_end.end());
		for(int i = 0; i < npoints; i++) data[i] = cont.data[i];
		return *this;
	}

	~Container(){
		delete [] data;
	}

	int number_of_functions() const {return nfunc;}
	int size() const {return npoints;}

	double norm() const{
		double tmp = 0.0;
		for(int i = 0; i < npoints; i++){
			tmp += data[i]*data[i];
		}
		return std::sqrt(tmp);
	}

	const double &operator[](int i) const {return data[i];}
	double &operator[](int i){return data[i];}

	template<class Callable>
	void fill_function(const Callable &func){
		int indx = 0;
		double *coords = new double[n_dim];
		for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			for(int ind = 0; ind < dom->size; ind++){
				update_coordinate(dom, ind, coords);
				for(int fn = 0; fn < nfunc; fn++)
					data[indx++] = func(fn, coords);
			}
		}
		delete [] coords;
	}

	template<class Callable>
	double compare_function(const Callable &func){
		int indx = 0;
		double *coords = new double[n_dim];
		double tmp, max_diff = 0.0;
		for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			for(int ind = 0; ind < dom->size; ind++){
				update_coordinate(dom, ind, coords);
				for(int fn = 0; fn < nfunc; fn++){
					tmp = std::fabs(data[indx++] - func(fn, coords));
					if(tmp > max_diff) max_diff = tmp;
				}
			}
		}
		delete [] coords;
		return max_diff;
	}

	void subtract_delta(const Collocation_Matrix &colm, double coeff = 1.0){
		if(npoints != colm.size()){
			std::cout << npoints << "," << colm.size() << std::endl;
			std::cout << "Error in Container delta subtraction, sizes don't match\n";
			return;
		}

		int indx = 0;
		std::vector<unsigned int> indv(n_dim, 0);
		for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			for(int ind = 0; ind < dom->size; ind++){
				for(int i = n_dim-1, tmp = ind; i >= 0; i--){
					if(indv[i] = (tmp % dom->dim_num(i))) break;
					tmp /= dom->dim_num(i);
				}
				for(int fn = 0; fn < nfunc; fn++){
					data[indx++] -= coeff*colm.delta(colm.index(dom->bulk_domain_number, indv, fn));
				}
			}
		}
		return;
	}

	void copy_to_vector(const Indexer &indexer, Eigen::VectorXd &vec) const{
		if(npoints != indexer.size()){
			vec.resize(indexer.size());
		}

		int indx = 0;
		std::vector<unsigned int> indv(n_dim, 0);
		for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			for(int ind = 0; ind < dom->size; ind++){
				for(int i = n_dim-1, tmp = ind; i >= 0; i--){
					if(indv[i] = (tmp % dom->dim_num(i))) break;
					tmp /= dom->dim_num(i);
				}
				for(int fn = 0; fn < nfunc; fn++){
					vec(indexer.index(dom->bulk_domain_number, indv, fn)) = data[indx++];
				}
			}
		}
	}

	//Changes partitioning while zeroing solution
	void adjust_partitioning(const Domain_Partitioning &part_init){
		if(this->compare_bulk_domains(part_init)) return;
		Domain_Partitioning part(part_init);
		if(part.boundary_domains_present()) part.clear_boundary_domains();
		reallocate_partitioning(part);
		for(int i = 0; i < npoints; i++) data[i] = 0.0;
	}

	void save(const std::string &folder_name) const;
	void load(const std::string &folder_name);

	//Interpolates solution to a new partitioning
	void regrid(const Domain_Partitioning &part);

	//Changes compactification
	template<class Compactification>
	void change_compactification(
			const Domain_Partitioning &part_new,
			const Compactification &comp_old,
			const Compactification &comp_new);

	friend class Evaluator;
	template<class>
	friend class Point;
	template<class>
	friend class Collocation_Solution;
	template<class, class>
	friend class Collocation_Solver_Full;
};


class Evaluator: public Domain_Partitioning{
private:
	int nfunc;

	int npoints;
	double *data;
	std::vector<int> index_st;

	void array_fft(void);

public:
	Evaluator(const Container &cont):
		Domain_Partitioning(cont), nfunc(cont.nfunc), npoints(cont.npoints), index_st(cont.index_st)
	{
		data = new double[npoints];
		for(int i = 0; i < npoints; i++) data[i] = cont.data[i];
		array_fft();
	}

	template<class Callable>
	Evaluator(const Domain_Partitioning &cont, Callable &func):
		Domain_Partitioning(cont), nfunc(1), npoints(0)
	{
		for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			index_st.push_back(npoints);
			npoints += dom->size;
		}
		data = new double[npoints];

		std::vector<unsigned int> indv(n_dim, 0);
		int indx = 0;
		for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			for(int ind = 0; ind < dom->size; ind++){
				for(int i = n_dim-1, tmp = ind; i >= 0; i--){
					if(indv[i] = (tmp % dom->dim_num(i))) break;
					tmp /= dom->dim_num(i);
				}
				data[indx++] = func(dom->bulk_domain_number, indv);
			}
		}

		array_fft();
	}

	Evaluator(Evaluator &&evali): Domain_Partitioning(std::move(evali)),
		nfunc(evali.nfunc), npoints(evali.npoints), data(evali.data), index_st(evali.index_st)
	{
		evali.data = NULL;
	}

	~Evaluator(){
		if(data != NULL) delete [] data;
	}

	double operator()(const std::vector<double> &coord, int fn) const;

	void test(const Container &cont){
		int indx = 0, i;
		double max = 0.0, tmp;

		std::vector<double> coord(n_dim);
		for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			for(int ind = 0; ind < dom->size; ind++){
				for(int dim = n_dim-1, tmpi = ind; dim >= 0; dim--){
					i = tmpi % dom->dim_num(dim);
					coord[dim] = ch_x(dom->domain_start(dim), dom->domain_end(dim), dom->dim_num(dim)-1, i);
					if(i) break;
					tmpi /= dom->dim_num(dim);
				}
				for(int fn = 0; fn < nfunc; fn++){
					tmp = std::fabs(cont.data[indx++] - (*this)(coord, fn));
					if(tmp > max) max = tmp;
				}
			}
		}

		std::cout << "Maximal discrepancy is " << max << std::endl;
	}

	friend class Container;
};

template<class Compactification>
void Container::change_compactification(
		const Domain_Partitioning &part_new,
		const Compactification &comp_old,
		const Compactification &comp_new)
{
	Evaluator eval_old(*this);
	adjust_partitioning(part_new);
	//Reinterpolating function

	unsigned int itmp, itr = 0;
	std::vector<unsigned int> ind(n_dim);
	double *xc_new = new double[n_dim];
	double *coord = new double[n_dim];
	std::vector<double> xc_old(n_dim);

	for(const Bulk_Domain *bulk = first_bulk_domain; bulk != NULL; bulk = bulk->next_bulk_domain()){
		int dm = bulk->bulk_domain_number;

		for(int indx = 0; indx < bulk->size; indx++){
			for(int i = n_dim-1, tmpi = indx; i >= 0; i--){
				itmp = tmpi % bulk->dim_num(i);
				xc_new[i] = ch_x(bulk->domain_start(i), bulk->domain_end(i), bulk->dim_num(i)-1, itmp);
				if(itmp) break;
				tmpi /= bulk->dim_num(i);
			}

			bool is_inf = false;
			for(int i = 0; i < n_dim; i++){
				coord[i] = comp_new.xc(i, xc_new);
				if(std::isinf(coord[i])){
					is_inf = true;
					break;
				}
				else xc_old[i] = comp_old.xc_inv(i, coord);
			}

			if(is_inf){
				itr += nfunc;
				continue;
			}
			for(int fn = 0; fn < nfunc; fn++)
				data[itr++] = eval_old(xc_old, fn);
		}
	}

	delete [] xc_new;
	delete [] coord;
}

#endif