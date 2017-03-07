// PageRank.cpp : Defines the entry point for the console application.
//

//Important: Want indices to be 32 bits signed for speed -- won't have matrices with more than 2 billion elements
#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int32_t

//For AVX2 instructions
#define EIGEN_MAX_ALIGN_BYTES 32

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <fstream>
#include <iostream>
#include <string>

Eigen::Matrix2Xi parse_association_list(const std::string& al_path)
{
	std::ifstream file{ al_path };
	file.exceptions(std::ifstream::failbit | std::ifstream::badbit);

	//Determine length of file -- each line is an edge, so this can be used to size matrix accordingly 
	int n_edges = 0;
	char buffer[1000];
	while (file.getline(buffer, 1000))
	{
		n_edges++;
	}

	//Start from beginning of file again
	file.seekg(0);

	Eigen::Matrix2Xi al{ 2, n_edges };
	//Add edges into AL matrix

	int32_t index = 0;
	std::string v1;
	std::string v2;

	while (file >> v1 >> v2)
	{
		int paper1 = std::stoi(v1);
		int paper2 = std::stoi(v2);
		al(0, index) = paper1;
		al(1, index) = paper2;
		index++;
	}

	return al;
}

namespace dense
{

	Eigen::MatrixXf generate_link_matrix(const Eigen::Matrix2Xi& al)
	{
		
	}

	void fixup_spider_trap(Eigen::MatrixXf& link_matrix)
	{

	}

	void divide_by_outgoing_edges(Eigen::MatrixXf& matrix)
	{

	}

	Eigen::VectorXf pagerank(const Eigen::Matrix2Xi& al)
	{
		auto M = generate_link_matrix(al);
		fixup_spider_trap(M);
		divide_by_outgoing_edges(M);

		//M is now of proper form -- iterate to find principal eigenvector


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
