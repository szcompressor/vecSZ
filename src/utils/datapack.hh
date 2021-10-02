#ifndef UTILS_DATAPACK_HH
#define UTILS_DATAPACK_HH

#include "io.hh"
#include <cstring>

namespace datapack {

// generate metadata required for decompression
template <typename C, typename O>
void generate_metadata()
{

}

// pack the data array to prep for write to file
template <typename C, typename O>
unsigned char* pack(size_t dataSize, int* dimsMeta, int ndims, unsigned int* lensMeta, int nlens, double eb, O* padValues, C* code, O* outlier)
{
	auto dataOut = new unsigned char[dataSize];

	memcpy(dataOut, dimsMeta, sizeof(int) * (ndims + 1));
	memcpy(dataOut + sizeof(int) * (ndims + 1), &eb, sizeof(double));
	memcpy(dataOut + sizeof(int) * (ndims + 1) + sizeof(double), lensMeta, sizeof(unsigned int) * nlens);
	memcpy(dataOut + sizeof(int) * (ndims + 1) + sizeof(double) + sizeof(unsigned int) * nlens, padValues, sizeof(O) * lensMeta[0]);
	memcpy(dataOut + sizeof(int) * (ndims + 1) + sizeof(double) + sizeof(unsigned int) * nlens + sizeof(O) * lensMeta[0], code, sizeof(C) * lensMeta[1]);
	memcpy(dataOut + sizeof(int) * (ndims + 1) + sizeof(double) + sizeof(unsigned int) * nlens + sizeof(O) * lensMeta[0] + sizeof(C) * lensMeta[1], outlier, sizeof(O) * lensMeta[2]);

	return dataOut;
}

// unpack the data array read from file
template <typename C, typename O>
void unpack(unsigned char* dataIn, int* ndims, int** dimsMeta, unsigned int **lensMeta, double* eb, O** padValues, unsigned int* padLength, C** code, O** outlier)
{
	int nlens = 4; // store len of code(1) and outlier(2)

	memcpy(ndims, dataIn, sizeof(int));

	*dimsMeta = new int[*ndims];
	memcpy(*dimsMeta, dataIn + sizeof(int), sizeof(int) * (*ndims));

	memcpy(eb, dataIn + sizeof(int) * (*ndims + 1), sizeof(double));

	*lensMeta = new unsigned int[nlens];
	memcpy(*lensMeta, dataIn + sizeof(int) * (*ndims + 1) + sizeof(double), sizeof(unsigned int) * nlens);

    *padLength = (*lensMeta)[0];
	*padValues = new O[(*lensMeta)[0]];
	memcpy(*padValues, dataIn + sizeof(int) * (*ndims + 1) + sizeof(double) + sizeof(unsigned int) * nlens, sizeof(O) * (*lensMeta)[0]);

	*code = new C[(*lensMeta)[1]];
	memcpy(*code, dataIn + sizeof(int) * (*ndims + 1) + sizeof(double) + sizeof(unsigned int) * nlens + sizeof(O) * (*lensMeta)[0], sizeof(C) * (*lensMeta)[1]);

	*outlier = new O[(*lensMeta)[2]];
	memcpy(*outlier, dataIn + sizeof(int) * (*ndims + 1) + sizeof(double) + sizeof(unsigned int) * nlens + sizeof(O) * (*lensMeta)[0] + sizeof(C) * (*lensMeta)[1], sizeof(O) * (*lensMeta)[2]);
}

} // namespace datapack

#endif //UTILS_DATAPACK_HH
