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
unsigned char* pack(size_t dataSize, int* dimsMeta, int ndims, size_t* lensMeta, int nlens, double eb, C* code, O* outlier)
{
	FILE* fp = fopen("datapack_out.txt", "wb");
	fprintf(fp, "dimsMeta: %d ", ndims);
	for (int i = 1; i < ndims + 1; i++) fprintf(fp, "%d ",dimsMeta[i]);
	fprintf(fp, "\n");
	fprintf(fp, "EB: %f\n", eb);
	fprintf(fp, "lensMeta: %ld %ld\n", lensMeta[0], lensMeta[1]);
	fprintf(fp, "code: \n", ndims);
	for (int i = 0; i < lensMeta[0]; i++) {fprintf(fp, "%d ",code[i]); if(i > 0 && (i%10==0)) fprintf(fp,"\n");}
	fprintf(fp, "\n");
	fprintf(fp, "outlier: \n", ndims);
	for (int i = 0; i < lensMeta[1]; i++) {fprintf(fp, "%f ",outlier[i]); if(i > 0 && (i%10==0)) fprintf(fp,"\n");}
	fprintf(fp, "\n");
	fclose(fp);

	fp = fopen("dataout_struct.txt", "wb");

	auto dataOut = new unsigned char[dataSize];

	memcpy(dataOut, dimsMeta, sizeof(int) * (ndims + 1));
	fprintf(fp, "dimsMeta: %d ", *dataOut);
	for (int i = 1; i < ndims + 1; i++) fprintf(fp, "%d ", (int)*(dataOut + i * sizeof(int)));
	fprintf(fp, "\n");
	memcpy(dataOut + sizeof(int) + (ndims + 1), &eb, sizeof(double));
	fprintf(fp, "EB: %f\n", (double)*(dataOut + (ndims + 1) * sizeof(int)));
	memcpy(dataOut + sizeof(int) + (ndims + 1) + sizeof(double), lensMeta, sizeof(size_t) * nlens);
	fprintf(fp, "lensMeta: %ld %ld\n", (size_t)*(dataOut + sizeof(int) + (ndims + 1) + sizeof(double)) , (size_t)*(dataOut + sizeof(int) + (ndims + 1) + sizeof(double) + sizeof(size_t)));
	memcpy(dataOut + sizeof(int) + (ndims + 1) + sizeof(double) + sizeof(size_t) * nlens, code, sizeof(C) * lensMeta[0]);
	fprintf(fp, "code: \n", ndims);
	for (int i = 0; i < lensMeta[0]; i++) {fprintf(fp, "%d ", (C)*(dataOut + sizeof(int) + (ndims + 1) + sizeof(double) + sizeof(size_t) * nlens + sizeof(C) * i)); if(i > 0 && (i%10==0)) fprintf(fp,"\n");}
	fprintf(fp, "\n");
	memcpy(dataOut + sizeof(int) + (ndims + 1) + sizeof(double) + sizeof(size_t) * nlens + sizeof(C) * lensMeta[0], outlier, sizeof(O) * lensMeta[1]);
	fprintf(fp, "outlier: \n", ndims);
	for (int i = 0; i < lensMeta[1]; i++) {fprintf(fp, "%f ", (O)*(dataOut + sizeof(int) + (ndims + 1) + sizeof(double) + sizeof(size_t) * nlens + sizeof(C) * lensMeta[0] + sizeof(O) * i)); if(i > 0 && (i%10==0)) fprintf(fp,"\n");}
	fprintf(fp, "\n");
	fclose(fp);

	return dataOut;
}

// unpack the data array read from file
template <typename C, typename O>
void unpack(unsigned char* dataIn, int* ndims, int** dimsMeta, int* nlens, size_t** lensMeta, double* eb, C** code, O** outlier)
{
	*nlens = 2; // store len of code(1) and outlier(2)

	memcpy(ndims, dataIn, sizeof(int));

	*dimsMeta = new int[*ndims];
	memcpy(*dimsMeta, dataIn + sizeof(int), sizeof(int) * (*ndims));

	memcpy(eb, dataIn + sizeof(int) * (*ndims + 1), sizeof(double));

	*lensMeta = new size_t[*nlens];
	memcpy(*lensMeta, dataIn + sizeof(int) * (*ndims + 1) + sizeof(double), sizeof(size_t) * (*nlens));

	*code = new C[(*lensMeta)[0]];
	memcpy(*code, dataIn + sizeof(int) * (*ndims + 1) + sizeof(double) + sizeof(size_t) * 2, sizeof(C) * (*lensMeta)[0]);

	*outlier = new O[(*lensMeta)[1]];
	memcpy(*outlier, dataIn + sizeof(int) * (*ndims + 1) + sizeof(double) + sizeof(size_t) * 2 + sizeof(C) * (*lensMeta)[0], sizeof(O) * (*lensMeta)[1]);
}

} // namespace datapack

#endif //UTILS_DATAPACK_HH
