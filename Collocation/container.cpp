#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

#include <fftw3.h>
#include <container.h>

#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

void Container::save(const std::string &folder_name) const{
	struct stat sttus = {0};
	if (stat(folder_name.c_str(), &sttus) == -1) {
		mkdir(folder_name.c_str(), 0700);
	}

	std::ofstream *out_f = new std::ofstream[nfunc];
	for(int fn = 0; fn < nfunc; fn++){
		out_f[fn].open(folder_name + "/f" + std::to_string(fn) + ".dat");
		//Recording parity
		for(int dm = 0; dm < n_dim; ++dm){
			if(dm) out_f[fn] << ',';
			if(parity_st[dm].empty())
				out_f[fn] << 'n' << ',';
			else if(parity_st[dm][fn])
				out_f[fn] << 'o' << ',';
			else
				out_f[fn] << 'e' << ',';
			if(parity_end[dm].empty())
				out_f[fn] << 'n';
			else if(parity_end[dm][fn])
				out_f[fn] << 'o';
			else
				out_f[fn] << 'e';
		}
		out_f[fn] << std::setprecision(15) << std::fixed << std::endl;
	}

	int indx = 0;
	const Bulk_Domain *dom = first_bulk_domain;
	while(dom != NULL){
		for(int fn = 0; fn < nfunc; fn++){
			out_f[fn] << dom->dim_num(0) << "," << dom->domain_start(0) << "," << dom->domain_end(0);
			for(int dm = 1; dm < n_dim; dm++)
				out_f[fn] << "," << dom->dim_num(dm) << "," << dom->domain_start(dm) << "," << dom->domain_end(dm);
			out_f[fn] << std::endl;
		}
		for(int i = 0; i < dom->size-1; i++){
			for(int fn = 0; fn < nfunc; fn++)
				out_f[fn] << data[indx++] << ",";
		}
		for(int fn = 0; fn < nfunc; fn++)
			out_f[fn] << data[indx++] << std::endl;

		dom = dom->next_bulk_domain();
	}

	for(int fn = 0; fn < nfunc; fn++)
		out_f[fn].close();
	delete [] out_f;
}

int count_func_files(const std::string &folder_name){
	std::string path;
	int nfunc_load = 0;
	struct stat sttus = {0};
	while(1){
		path = folder_name + "/f" + std::to_string(nfunc_load) + ".dat";
		if(stat(path.c_str(), &sttus) == 0) nfunc_load++;
		else break;
	}
	if(nfunc_load == 0){
		std::cout << "No function files found in folder: " << folder_name << std::endl;
		throw("reading error");
	}
	return nfunc_load;
}

Domain_Partitioning get_partitioning_from_data(const std::string &folder_name){
	std::ifstream in(folder_name + "/f0.dat");
	if(in.is_open() == false){
		std::cout << "No data file 0 found" << std::endl;
		throw("No data file found");
	}

	std::string::size_type sztp;
	std::vector<std::string> domain_data;
	int nlines = 0;
	std::string line, token;
	while(1){
		if(!std::getline(in, line)) break;
		if(line.empty()) break;

		//Skipping parity description line
		std::istringstream iss(line);
		std::getline(iss, token, ',');
		try{
			std::stod(token, &sztp);
		}
		catch(std::invalid_argument const& ex){
			continue;
		}

		if((nlines++)%2 == 0) domain_data.push_back(line);
	}

	return Domain_Partitioning(domain_data);
}

void Collocation_Matrix::row_add_Dn_factor(unsigned int row_dm, const std::vector<unsigned int> &row_ind, unsigned int row_fn,
		int der1, int der2, unsigned int dm, const std::vector<unsigned int> &ind, unsigned int fn, double coef)
{
	int row_index = index(row_dm, row_ind, row_fn);

	if(der1 < 0 && der2 < 0){
		matr(row_index, index(dm, ind, fn)) += coef;
		return;
	}

	int ord_x = 0, ord_y = 0;
	std::vector<unsigned int> ind_tmp(ind);
	if(der2 < der1){int tmp = der1; der1 = der2; der2 = tmp;}


	double factor;
	int dim2_n = dim_nums[dm][der2];
	if(der1 < 0 || der1 == der2){
		int ord = (der1 == der2 ? 2 : 1);
		std::pair<int, bool> parity = check_parity(dm, der2, fn);
		if(parity.first){
			factor = coef*std::pow(1.0/(doms_end[dm][der2] - doms_start[dm][der2]), ord);
			for(ind_tmp[der2] = 0; ind_tmp[der2] < dim2_n; ind_tmp[der2]++){
				matr(row_index, index(dm, ind_tmp, fn)) += factor*cheb_Dn_double(parity.first == 2, parity.second, dim2_n, ord, ind[der2], ind_tmp[der2]);
			}
		}
		else{
			factor = coef * std::pow(2.0/(doms_end[dm][der2] - doms_start[dm][der2]), ord);
			for(ind_tmp[der2] = 0; ind_tmp[der2] < dim2_n; ind_tmp[der2]++){
				matr(row_index, index(dm, ind_tmp, fn)) += factor*cheb_Dn_el(dim2_n-1, ord, dim2_n-ind[der2]-1, dim2_n-ind_tmp[der2]-1);
			}
		}
		return;
	}
	else{
		int dim1_n = dim_nums[dm][der1];
		std::pair<int, bool> parity1 = check_parity(dm, der1, fn);
		std::pair<int, bool> parity2 = check_parity(dm, der2, fn);
		double tmp;
		factor = coef*4.0/(doms_end[dm][der1] - doms_start[dm][der1])/(doms_end[dm][der2] - doms_start[dm][der2]);
		for(ind_tmp[der1] = 0; ind_tmp[der1] < dim1_n; ind_tmp[der1]++){
			if(parity1.first)
				tmp = 0.5*factor*cheb_Dn_double(parity1.first == 2, parity1.second, dim1_n, 1, ind[der1], ind_tmp[der1]);
			else
				tmp = factor*cheb_Dn1_el(dim1_n-1, dim1_n-ind[der1]-1, dim1_n-ind_tmp[der1]-1);
			if(parity2.first){
				for(ind_tmp[der2] = 0; ind_tmp[der2] < dim2_n; ind_tmp[der2]++){
					matr(row_index, index(dm, ind_tmp, fn)) += 0.5*tmp*cheb_Dn_double(parity2.first == 2, parity2.second, dim2_n, 1, ind[der2], ind_tmp[der2]);
				}
			}
			else{
				for(ind_tmp[der2] = 0; ind_tmp[der2] < dim2_n; ind_tmp[der2]++){
					matr(row_index, index(dm, ind_tmp, fn)) += tmp*cheb_Dn1_el(dim2_n-1, dim2_n-ind[der2]-1, dim2_n-ind_tmp[der2]-1);
				}
			}
		}
		return;
	}
}

void Container::load(const std::string &folder_name){
	std::ifstream *in_f = new std::ifstream[nfunc];
	for(int fn = 0; fn < nfunc; fn++){
		in_f[fn].open(folder_name + + "/f" + std::to_string(fn) + ".dat");
		if(in_f[fn].is_open() == false){
			delete [] in_f;
			std::cout << "Not enough function files to load" << std::endl;
			throw("Not enough function files to load");
		}
	}

	std::string line, piece;
	std::string::size_type sztp;
	//Reading parity description
	int ndim_est = 0;
	parity_st.clear(); parity_end.clear();
	for(int fn = 0; fn < nfunc; fn++){
		std::getline(in_f[fn], line);
		std::istringstream iss(line);
		std::getline(iss, piece, ',');
		try{
			std::stod(piece, &sztp);
		}
		catch(std::invalid_argument const& ex){
			iss.seekg(0, iss.beg);
			int cnt = 0;
			while(std::getline(iss, piece, ',')){
				if(piece.compare("n") == 0){
					if(cnt%2 && cnt/2 >= parity_end.size()) parity_end.push_back({});
					else if(cnt/2 >= parity_end.size()) parity_st.push_back({});
				}
				else if(piece.compare("e") == 0){
					if(cnt%2){
						if(cnt/2 >= parity_end.size()) parity_end.push_back({false});
						else parity_end[cnt/2].push_back(false);
					}
					else{
						if(cnt/2 >= parity_st.size()) parity_st.push_back({false});
						else parity_st[cnt/2].push_back(false);
					}
				}
				else if(piece.compare("o") == 0){
					if(cnt%2){
						if(cnt/2 >= parity_end.size()) parity_end.push_back({true});
						else parity_end[cnt/2].push_back(true);
					}
					else{
						if(cnt/2 >= parity_st.size()) parity_st.push_back({true});
						else parity_st[cnt/2].push_back(true);
					}
				}
				else{
					std::cout << "Wrong parity flag read: " << piece << std::endl;
					throw("Wrong parity flag read");
				}
				++cnt;
			}
			if(ndim_est == 0 && cnt) ndim_est = cnt;
			if(ndim_est != cnt){
				std::cout << "Parity data corrupted" << std::endl;
				throw("Parity data corrupted");
			}
			continue;
		}
		in_f[fn].seekg(0, in_f[fn].beg);
	}
	for(int dm = 0; dm < parity_st.size(); ++dm){
		if((parity_st[dm].size() && parity_st[dm].size() != nfunc) || (parity_end[dm].size() && parity_end[dm].size() != nfunc)){
			std::cout << "Parity data corrupted" << std::endl;
			throw("Parity data corrupted");
		}
	}

	for(int fn = 0; fn < nfunc; fn++){
		std::vector<std::string> domain_data;
		std::vector<std::string> function_data;

		int nlines = 0;
		while(1){
			if(!std::getline(in_f[fn], line)) break;
			if(line.empty()) break;
			if((nlines++)%2) function_data.push_back(line);
			else domain_data.push_back(line);
		}

		if(domain_data.size() != function_data.size()){
			std::cout << "Data file is " << fn << " broken" << std::endl;
			throw("Data file is broken");
		}

		if(fn == 0){
			Domain_Partitioning part(domain_data);
			if(!(*this == part)) reallocate_partitioning(part);
		}
		else if(!(*this == Domain_Partitioning(domain_data))){
			std::cout << "Data file " << fn << " domain structure does not match container structure" << std::endl;
			throw("Data file domain structure does not match container structure");
		}
		else if(parity_st.size() != n_dim){
			if(parity_st.size() == 0){
				parity_st.resize(n_dim);
				parity_end.resize(n_dim);
				std::fill<std::vector<std::vector<bool>>::iterator,std::vector<bool>>(parity_st.begin(), parity_st.end(), {});
				std::fill<std::vector<std::vector<bool>>::iterator,std::vector<bool>>(parity_end.begin(), parity_end.end(), {});
			}
			else{
				std::cout << "Data file " << fn << " parity structure does not match number of dimensions" << std::endl;
				throw("Data file parity structure does not match number of dimensions");
			}
		}

		const Bulk_Domain *dom = first_bulk_domain;
		while(dom != NULL){
			int dm = dom->bulk_domain_number;
			std::istringstream iss(function_data[dm]);
			int index_end = index_st[dm]+nfunc*dom->size;
			for(int i = index_st[dm]+fn; i < index_end; i+= nfunc){
				std::getline(iss, piece, ',');
				data[i] = std::stod(piece, &sztp);
			}
			dom = dom->next_bulk_domain();
		}
	}

	for(int fn = 0; fn < nfunc; fn++)
		in_f[fn].close();
	delete [] in_f;
}

void Container::subtract_delta(const Collocation_Matrix &colm, double coeff){
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

void Container::copy_to_vector(const Indexer &indexer, Eigen::VectorXd &vec) const{
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

void Container::regrid(const Domain_Partitioning &part_init){
	if(part_init.n_dim != n_dim){
		std::cout << "Regrid cannot change dimensions" << std::endl;
		throw("Regrid dimension change");
	}
	if(partition_start() != part_init.partition_start() || partition_end() != part_init.partition_end()){
		std::cout << "Regrid cannot change outside boundaries" << std::endl;
		throw("Regrid boundary change");	
	}

	if(this->compare_bulk_domains(part_init)){
		std::cout << "No container regridding is needed" << std::endl;
		return;
	}
	Domain_Partitioning part(part_init);
	if(part.boundary_domains_present()) part.clear_boundary_domains();
	Evaluator eval(*this);
	
	reallocate_partitioning(part);

	int tmpi, indx = 0;
	std::vector<double> coord(n_dim);
	for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
		for(int ind = 0; ind < dom->size; ind++){
			for(int i = n_dim-1, tmp = ind; i >= 0; i--){
				tmpi = tmp % dom->dim_num(i);
				coord[i] = ch_x(dom->domain_start(i), dom->domain_end(i), dom->dim_num(i)-1, tmpi);
				if(tmpi) break;
				tmp /= dom->dim_num(i);
			}
			for(int fn = 0; fn < nfunc; fn++)
				data[indx++] = eval(coord, fn);
		}
	}
}

void Evaluator::array_fft(void){
	double *ptr;
	
	double fft_coef;
	fftw_plan planCOS;
	int *n_array = new int[n_dim];
	fftw_r2r_kind *kind = new fftw_r2r_kind[n_dim];
	for(int dim = 0; dim < n_dim; dim++) kind[dim] = FFTW_REDFT00;

	for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
		ptr = data + index_st[dom->bulk_domain_number];
		
		fft_coef = 1.0;
		for(int dim = 0; dim < n_dim; dim++){
			fft_coef *= dom->dim_num(dim)-1;
			n_array[dim] = dom->dim_num(dim);
		}
		planCOS = fftw_plan_many_r2r(n_dim, n_array, nfunc, ptr, n_array, nfunc, 1, ptr, n_array, nfunc, 1, kind, FFTW_ESTIMATE);

		for(int i = 0; i < nfunc*dom->size; i++) ptr[i] /= fft_coef;

		fftw_execute_r2r(planCOS, ptr, ptr);

		fftw_destroy_plan(planCOS);
	}

	delete [] kind;
	delete [] n_array;
}

double Evaluator::operator()(const std::vector<double> &coord, int fn) const{
	const Bulk_Domain *dom = this->Domain_Partitioning::operator()(coord);
	int indx = index_st[dom->bulk_domain_number] + (dom->size - 1)*nfunc + fn;

	double sv, reslt;
	std::vector<double> y2(n_dim), d(n_dim, 0.0), dd(n_dim, 0.0);
	for(int dim = 0; dim < n_dim; dim++) y2[dim] = 2.0 * (dom->domain_start(dim) + dom->domain_end(dim) - 2.0 * coord[dim])
													/ (dom->domain_end(dim) - dom->domain_start(dim));

	for(int ind = dom->size-1; ind >= 0; ind--){
		if(ind % dom->dim_num(n_dim-1)){
			sv = d[n_dim-1];
			d[n_dim-1] = y2[n_dim-1] * d[n_dim-1] - dd[n_dim-1] + data[indx];
			dd[n_dim-1] = sv;
		}
		else{
			reslt = 0.5 * y2[n_dim-1] * d[n_dim-1] - dd[n_dim-1] + 0.5 * data[indx];

			for(int dim = n_dim-2, tmpi = ind; dim >= 0; dim--){
				tmpi /= dom->dim_num(dim+1);
				d[dim+1] = 0.0;
				dd[dim+1] = 0.0;

				if(tmpi % dom->dim_num(dim)){
					sv = d[dim];
					d[dim] = y2[dim] * d[dim] - dd[dim] + reslt;
					dd[dim] = sv;
					break;
				}
				else{
					reslt = 0.5 * y2[dim] * d[dim] - dd[dim] + 0.5 * reslt;
				}
			}
		}
		indx -= nfunc;
	}

	return reslt; 
}