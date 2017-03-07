// PageRank.cpp : Defines the entry point for the console application.
//

//For AVX2 instructions
#define EIGEN_MAX_ALIGN_BYTES 32

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <fstream>
#include <iostream>
#include <set>
#include <string>

namespace Eigen
{
	using Matrix2Xi64 = Matrix<int64_t, 2, -1>;
	using VectorXi64 = Matrix<int64_t, -1, 1>;
}

Eigen::Matrix2Xi64 parse_association_list(const std::string& al_path)
{
	int n_edges = 0;

	{
		std::ifstream file{ al_path };

		//Determine length of file -- each line is an edge, so this can be used to size matrix accordingly 
		char buffer[1000];
		while (file.getline(buffer, 1000))
		{
			n_edges++;
		}
	}

	//Start from beginning of file again
	std::ifstream file{ al_path };

	Eigen::Matrix2Xi64 al{ 2, n_edges };
	//Add edges into AL matrix

	Eigen::Index index = 0;
	std::string v1;
	std::string v2;

	while (file >> v1 >> v2)
	{
		const auto paper1 = std::stoll(v1);
		const auto paper2 = std::stoll(v2);
		al(0, index) = paper1;
		al(1, index) = paper2;
		index++;
	}

	return al;
}

//Generate a vector of unique IDs for each paper
//The position of the paper in this vector is its unique ID
//This is used to map rows and columns of the M matrix, since the paper IDs are not contiguous and do not start from 0
//returned vector is sorted in ascending order
Eigen::VectorXi64 unique_ids(const Eigen::Matrix2Xi64& al)
{
	std::set<int64_t> uids;

	for (Eigen::Index i = 0; i < al.cols(); i++)
	{
		uids.emplace(al(0, i));
		uids.emplace(al(1, i));
	}

	Eigen::VectorXi64 uid_vector{ static_cast<Eigen::Index>(uids.size()) };

	Eigen::Index index = 0;
	for (const auto i : uids)
	{
		uid_vector(index) = i;
		index++;
	}

	return uid_vector;
}

Eigen::Index paperid_to_uid(const Eigen::VectorXi64& uids, int64_t paper)
{
	//Take advantage of sorted vector -- do binary search
	const auto n = uids.size();
	Eigen::Index lowerBound = 0;
	Eigen::Index upperBound = n - 1;

	Eigen::Index midpoint = 0;

	do {
		midpoint = (lowerBound + upperBound) / 2;
		if (uids(midpoint) < paper)
		{
			//Take second half
			lowerBound = midpoint + 1;
		}
		else if (paper < uids(midpoint))
		{
			//Take first half
			upperBound = midpoint - 1;
		}

	} while (uids(midpoint) != paper);

	return midpoint;
}


namespace dense
{

	//Generate adjacency matrix from association list and unique IDs
	//This is a dense matrix.  If many papers are used, this will *not* scale!
	Eigen::MatrixXd adjacency_matrix(const Eigen::Matrix2Xi64& al, const Eigen::VectorXi64& uids)
	{
		//matrix must only be N x N in size, where N is the number of papers
		Eigen::MatrixXd A = Eigen::MatrixXd::Zero( uids.size(), uids.size() );

		for (Eigen::Index i = 0; i < al.cols(); i++)
		{
			const auto paper1 = al(0, i);
			const auto paper2 = al(1, i);
			const auto uid1 = paperid_to_uid(uids, paper1);
			const auto uid2 = paperid_to_uid(uids, paper2);
			A(uid1, uid2) = 1;
		}

		return A;
	}

	//Divide the contents of each row by the number of outlinks (or the L1 norm)
	//Fix spider trap in here as well
	//Second argument is the "teleport" chance for random surfer
	void convert_to_link_matrix(Eigen::MatrixXd& A, const double t)
	{
		const auto C = A.rowwise().sum();

		const auto uniform = 1/static_cast<double>(A.rows());
		const auto uniformT = t*uniform;
		for (Eigen::Index i = 0; i < A.rows(); i++)
		{
			if (C(i) != 0)
			{
				A.row(i).array() /= C(i);
				A.row(i).array() *= (1 - t);
				A.row(i).array() += uniformT;
				//std::cout << "i: " << A.row(i) << std::endl;
			}
			else
			{
				//Fix spider trap
				A.row(i).setConstant(uniform);
			}
		}
	}

	Eigen::VectorXd principal_eigenvector(const Eigen::MatrixXd& A)
	{

		//Initialize surfer probability starting at paper 1
		Eigen::RowVectorXd P = Eigen::VectorXd::Zero(A.rows());
		P(0) = 1;
		double diff = 1;

		int nIterations = 0;
		while (diff != 0)
		{
			Eigen::RowVectorXd PN = P * A;
			diff = (PN - P).squaredNorm();
			P = PN;
			nIterations++;
		}

		std::cout << "Convergence after " << nIterations << " iterations" << std::endl;

		return P;
	}

	std::map<double, int64_t> pagerank(const Eigen::Matrix2Xi64& al, const double t)
	{
		std::cout << "Dense matrix pagerank selected" << std::endl;
		const auto uids = unique_ids(al);
		std::cout << "Number of edges: " << al.cols() << std::endl;
		std::cout << "Number of papers: " << uids.size() << std::endl;
		auto A = adjacency_matrix(al, uids);
		std::cout << "Adjacency matrix generated" << std::endl; //\n" << A << std::endl;
		convert_to_link_matrix(A, t);
		std::cout << "Link matrix generated" << std::endl; // \n" << std::endl;
		const auto P = principal_eigenvector(A);
		std::map<double, int64_t> pageRanks;
		for (Eigen::Index i = 0; i < P.size(); i++)
		{
			pageRanks.emplace(P(i), uids(i));
		}

		//M is now of proper form -- iterate to find principal eigenvector

		return pageRanks;


	}

}

namespace sparse
{

	//Row major matrices used for efficiency: left eigenvalue 
	Eigen::SparseMatrix<double, Eigen::RowMajor> adjacency_matrix(const Eigen::Matrix2Xi64& al, const Eigen::VectorXi64& uids)
	{
		Eigen::SparseMatrix<double, Eigen::RowMajor> A{ uids.size(), uids.size() };
		Eigen::SparseVector<double> C{ uids.size() };

		std::vector<Eigen::Triplet<double>> triplets;
		for (Eigen::Index i = 0; i < al.cols(); i++)
		{
			const auto paper1 = al(0, i);
			const auto paper2 = al(1, i);
			triplets.emplace_back(paperid_to_uid(uids, paper1), paperid_to_uid(uids, paper2), 1);
		}


		A.setFromTriplets(triplets.begin(), triplets.end());

		return A;
	}

	//Switch to column major matrices for computing left eigenvector (test)
	std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd> adj_to_link_matrix(Eigen::SparseMatrix<double, Eigen::RowMajor>&& A, const double t)
	{
		//Also return a dense vector d that element i is 1 if row i is all zeros and 0 otherwise
		Eigen::VectorXd d = Eigen::VectorXd::Zero(A.rows());
		
		for (Eigen::Index i = 0; i < A.rows(); i++)
		{
			const auto Ci = A.row(i).nonZeros();
			if (Ci > 0)
			{
				A.row(i) /= Ci;
			}
			else {
				d(i) = 1;
			}
		}

		return { Eigen::SparseMatrix<double> { A*(1-t) }, d };
	}

	Eigen::RowVectorXd principal_eigenvector(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& d, const double t)
	{
		const auto n = A.rows();
		const auto teT_n = Eigen::RowVectorXd::Constant(n, (1-t) / n);

		//const auto tE = Eigen::MatrixXf::Constant(n, n, t / n);
		const auto t_n = t / n;

		//M = tE + (1 - t)D + (1 - t)A
		//A is prescaled

		//E = uniform(1/n)
		//D = deT/n
		//n = A.rows()

		Eigen::RowVectorXd p = Eigen::RowVectorXd::Zero(n);
		p(0) = 1;

		double dist = 1;
		int nIterations = 0;
		while (dist != 0)
		{
			Eigen::RowVectorXd p2 = p * A + ((p*d)*(teT_n).array() + p.sum()*(t_n)).matrix();
			dist = (p2 - p).squaredNorm();
			nIterations++;
			p = p2;
		}

		std::cout << "Convergence after " << nIterations << " iterations" << std::endl;

		return p;
	}

	std::map<double, int64_t> pagerank(const Eigen::Matrix2Xi64& al, const double t)
	{
		std::cout << "Sparse matrix pagerank selected" << std::endl;
		const auto uids = unique_ids(al);
		std::cout << "Number of edges: " << al.cols() << std::endl;
		std::cout << "Number of papers: " << uids.size() << std::endl;

		//double t = 0.15f;
		const auto pair = adj_to_link_matrix(adjacency_matrix(al, uids), t);
		const auto& A = pair.first;
		const auto& d = pair.second;

		const auto P = principal_eigenvector(A, d, t);

		std::map<double, int64_t> pageRanks;
		for (Eigen::Index i = 0; i < P.size(); i++)
		{
			pageRanks.emplace(P(i), uids(i));
		}
		//std::cout << "A:\n" << A << std::endl;

		return pageRanks;
	}
}

enum class ComputationType
{
	Sparse, Dense
};

ComputationType parseComputationArg(const std::string& arg)
{
	if (arg == "--dense")
	{
		return ComputationType::Dense;
	}
	else if (arg == "--sparse")
	{
		return ComputationType::Sparse;
	}
	else
	{
		std::cerr << "Error: unrecognized argument " << arg << std::endl;
		throw std::runtime_error{ "unrecognized argument" };
	}
}

int __cdecl main(int argc, char*argv[])
{

	if (argc < 3)
	{
		std::cout << "Usage: " << argv[0] << " <association list file> (--dense | --sparse) (-t <n>)" << std::endl;
		return 0;
	}

	const std::string path{ argv[1] };

	const ComputationType ctype = parseComputationArg(argv[2]);

	double t = 0.15;

	if (argc == 5)
	{
		t = std::stod(argv[4]);
	}

	const auto al = parse_association_list(path);
	std::map<double, int64_t> pageRanks;

	if (ctype == ComputationType::Dense)
	{
		pageRanks = dense::pagerank(al, t);
	}
	else if (ctype == ComputationType::Sparse)
	{
		pageRanks = sparse::pagerank(al, t);
	}

	std::cout << "PageRanks:\n" << std::endl;
	for (auto it = pageRanks.crbegin(); it != pageRanks.crend(); it++)
	{
		const auto& rank = *it;
		std::cout << rank.first << " " << rank.second << std::endl;
	}



    return 0;
}

//Dense PageRank
