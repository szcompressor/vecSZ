#ifndef ARGPARSE_HH
#define ARGPARSE_HH

#ifndef MAX_VECTOR_LENGTH
#define MAX_VECTOR_LENGTH -1
#endif

#include <cstdlib>
#include <iostream>
#include <regex>
#include <string>

#include "../types.hh"
#include "utils/format.hh"

using std::string;

struct SZWorkflow {
	bool use_demo{false};
	bool lossy_compress{false};
	bool lossy_decompress{false};
	bool lossy_dryrun{false};
	bool autotune{false};
};

struct fileNames {
	string input_file;
	string output_file;
};

class ArgParse {
	public:

		struct fileNames files;
		struct SZWorkflow szwf;

		stat_t stat;
		int read_args_status{0};

		string mode;
		string demo_dataset;
		string opath;
		string dtype;

		int ndim{-1};

		double eb{0.0001};
		double sample_percentage{-1};

		size_t   len{1};
		Integer4 dim4{1, 1, 1, 1};
		Integer4 nblk4{1, 1, 1, 1};
		Integer4 stride4{1, 1, 1, 1};
		int      dict_size{4096};
		int      radius{2048};
		int      vector_length{MAX_VECTOR_LENGTH};
		int      block_size{8};
		int      num_iterations{-1};

		bool verbose{false};

		static void vecszDoc();

		static void vecszFullDoc();

		static string format(const string& s);

		int trap(int _status);

		ArgParse() = default;

    		static int SelfMultiple4(Integer4 i) { return i._0 * i._1 * i._2 * i._3; }

		void ParseVecszArgs(int argc, char** argv);

		void CheckArgs();

		void PrintArgs();
};

typedef ArgParse argparse;

#endif //ARGPARSE_HH
