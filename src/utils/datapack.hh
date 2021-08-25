#ifndef UTILS_DATAPACK_HH
#define UTILS_DATAPACK_HH

#include "io.hh"

namespace datapack {

// generate metadata required for decompression
template <typename C, typename O>
void generate_metadata()
{

}

// pack the data array to prep for write to file
template <typename C, typename O>
unsigned char* pack(size_t dataSize, int* dimsMeta, int ndims, size_t* lensMeta, int nlens, double eb, O* padValues, C* code, O* outlier)
{
	auto dataOut = new unsigned char[dataSize];

	memcpy(dataOut, dimsMeta, sizeof(int) * (ndims + 1));
	memcpy(dataOut + sizeof(int) * (ndims + 1), &eb, sizeof(double));
	memcpy(dataOut + sizeof(int) * (ndims + 1) + sizeof(double), lensMeta, sizeof(size_t) * nlens);
	memcpy(dataOut + sizeof(int) * (ndims + 1) + sizeof(double) + sizeof(size_t) * nlens, padValues, sizeof(O) * lensMeta[0]);
	memcpy(dataOut + sizeof(int) * (ndims + 1) + sizeof(double) + sizeof(size_t) * nlens + sizeof(O) * lensMeta[0], code, sizeof(C) * lensMeta[1]);
	memcpy(dataOut + sizeof(int) * (ndims + 1) + sizeof(double) + sizeof(size_t) * nlens + sizeof(O) * lensMeta[0] + sizeof(C) * lensMeta[1], outlier, sizeof(O) * lensMeta[2]);

	return dataOut;
}

// unpack the data array read from file
template <typename C, typename O>
void unpack(unsigned char* dataIn, int* ndims, int** dimsMeta, double* eb, O** padValues, size_t* padLength, C** code, O** outlier)
{
	int nlens = 3; // store len of code(1) and outlier(2)

	memcpy(ndims, dataIn, sizeof(int));

	*dimsMeta = new int[*ndims];
	memcpy(*dimsMeta, dataIn + sizeof(int), sizeof(int) * (*ndims));

	memcpy(eb, dataIn + sizeof(int) * (*ndims + 1), sizeof(double));

	auto lensMeta = new size_t[nlens];
	memcpy(lensMeta, dataIn + sizeof(int) * (*ndims + 1) + sizeof(double), sizeof(size_t) * nlens);

    *padLength = lensMeta[0];
	*padValues = new O[lensMeta[0]];
	memcpy(*padValues, dataIn + sizeof(int) * (*ndims + 1) + sizeof(double) + sizeof(size_t) * nlens, sizeof(O) * lensMeta[0]);

	*code = new C[lensMeta[1]];
	memcpy(*code, dataIn + sizeof(int) * (*ndims + 1) + sizeof(double) + sizeof(size_t) * nlens + sizeof(O) * lensMeta[0], sizeof(C) * lensMeta[1]);

	*outlier = new O[lensMeta[2]];
	memcpy(*outlier, dataIn + sizeof(int) * (*ndims + 1) + sizeof(double) + sizeof(size_t) * nlens + sizeof(O) * lensMeta[0] + sizeof(C) * lensMeta[1], sizeof(O) * lensMeta[2]);
}

} // namespace datapack

#endif //UTILS_DATAPACK_HH
