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

namespace dense
{
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
		bool found = false;
		Eigen::Index lowerBound = 0;
		Eigen::Index upperBound = n-1;
		
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

	//Generate adjacency matrix from association list and unique IDs
	//This is a dense matrix.  If many papers are used, this will *not* scale!
	Eigen::MatrixXf adjacency_matrix(const Eigen::Matrix2Xi64& al, const Eigen::VectorXi64& uids)
	{
		//matrix must only be N x N in size, where N is the number of papers
		Eigen::MatrixXf A = Eigen::MatrixXf::Zero( uids.size(), uids.size() );

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
	void convert_to_link_matrix(Eigen::MatrixXf& A, const float t)
	{
		const auto C = A.rowwise().sum();

		const auto uniform = 1/static_cast<float>(A.rows());
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

	Eigen::VectorXf principal_eigenvector(const Eigen::MatrixXf& A)
	{

		//Initialize surfer probability starting at paper 1
		Eigen::RowVectorXf P = Eigen::VectorXf::Zero(A.rows());
		P(0) = 1;
		float diff = 1;

		int nIterations = 0;
		while (diff != 0)
		{
			Eigen::RowVectorXf PN = P * A;
			diff = (PN - P).squaredNorm();
			P = PN;
			nIterations++;
		}

		std::cout << "Convergence after " << nIterations << " iterations" << std::endl;

		return P;
	}

	Eigen::VectorXf pagerank(const Eigen::Matrix2Xi64& al)
	{
		std::cout << "Dense matrix pagerank selected" << std::endl;
		const auto uids = unique_ids(al);
		std::cout << "Number of edges: " << al.cols() << std::endl;
		std::cout << "Number of papers: " << uids.size() << std::endl;
		auto A = adjacency_matrix(al, uids);
		std::cout << "Adjacency matrix generated" << std::endl; //\n" << A << std::endl;
		convert_to_link_matrix(A, 0.15);
		std::cout << "Link matrix generated" << std::endl; // \n" << std::endl;
		const auto P = principal_eigenvector(A);
		std::map<float, int64_t> pageRanks;
		for (Eigen::Index i = 0; i < P.size(); i++)
		{
			pageRanks.emplace(P(i), uids(i));
		}
		std::cout << "PageRanks:\n" << std::endl;
		for (auto it = pageRanks.crbegin(); it != pageRanks.crend(); it++)
		{
			const auto& rank = *it;
			std::cout << rank.first << " " << rank.second << std::endl;
		}

		//M is now of proper form -- iterate to find principal eigenvector

		return {};


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
		std::cout << "Usage: " << argv[0] << " <association list file> (--dense | --sparse)" << std::endl;
		return 0;
	}

	const std::string path{ argv[1] };

	const ComputationType ctype = parseComputationArg(argv[2]);

	const auto al = parse_association_list(path);

	if (ctype == ComputationType::Dense)
	{
		dense::pagerank(al);
	}
	else if (ctype == ComputationType::Sparse)
	{
		throw std::runtime_error{ "sparce matrices not supported yet" };
	}



    return 0;
}

//Dense PageRank
