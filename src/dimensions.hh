#ifndef DIMENSIONS_HH
#define DIMENSIONS_HH

#include <cstring>
#include <string>
#include <vector>

#include "argument_parser/argparse.hh"
#include "constants.hh"

size_t* InitializeDims(ArgParse* ap)
{
	std::unordered_map<string, std::vector<int>> dataset_entries = {
        {std::string("hacc"), {280953867, 1, 1, 1, 1}},    {std::string("hacc1b"), {1073726487, 1, 1, 1, 1}},
        {std::string("cesm"), {3600, 1800, 1, 1, 2}},      {std::string("hurricane"), {500, 500, 100, 1, 3}},
        {std::string("nyx"), {512, 512, 512, 1, 3}},       {std::string("nyx-m"), {1024, 1024, 1024, 1, 3}},
        {std::string("qmc"), {288, 69, 7935, 1, 3}},       {std::string("qmcpre"), {69, 69, 33120, 1, 3}},
        {std::string("exafel"), {388, 59200, 1, 1, 2}},    {std::string("aramco"), {235, 849, 849, 1, 3}},
        {std::string("parihaka"), {1168, 1126, 922, 1, 3}}};

	if (not ap->demo_dataset.empty())
	{
		auto demo_dims = dataset_entries.at(ap->demo_dataset);

		ap->dim4._0 = (int)demo_dims[0];
		ap->dim4._1 = (int)demo_dims[1];
		ap->dim4._2 = (int)demo_dims[2];
		ap->dim4._3 = (int)demo_dims[3];
		ap->ndim    = (int)demo_dims[4];
	}

	auto get_nblk = [&](int d) { return (d + ap->block_size + 1) / ap->block_size; };

	ap->nblk4._0 = get_nblk(ap->dim4._0);
	ap->nblk4._1 = get_nblk(ap->dim4._1);
	ap->nblk4._2 = get_nblk(ap->dim4._2);
	ap->nblk4._3 = get_nblk(ap->dim4._3);

	ap->len = ap->dim4._0 * ap->dim4._1 * ap->dim4._2 * ap->dim4._3;

	ap->stride4 = {
		1,
		ap->dim4._0,
		ap->dim4._0 * ap->dim4._1,
		ap->dim4._0 * ap->dim4._1 * ap->dim4._2};


    auto dims_L16 = new size_t[16]();

	dims_L16[nDIM]   = ap->ndim;
	dims_L16[DIM0]   = ap->dim4._0;
	dims_L16[DIM1]   = ap->dim4._1;
	dims_L16[DIM2]   = ap->dim4._2;
	dims_L16[DIM3]   = ap->dim4._3;
    dims_L16[nBLK0]  = ap->nblk4._0;
    dims_L16[nBLK1]  = ap->nblk4._1;
    dims_L16[nBLK2]  = ap->nblk4._2;
    dims_L16[nBLK3]  = ap->nblk4._3;
    dims_L16[LEN]    = dims_L16[DIM0] * dims_L16[DIM1] * dims_L16[DIM2] * dims_L16[DIM3];
    dims_L16[CAP]    = ap->dict_size;
    dims_L16[RADIUS] = ap->radius;

	return dims_L16;
}

#endif // DIMENSIONS_HH