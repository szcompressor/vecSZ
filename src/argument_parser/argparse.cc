/*
 * @file argparse.cc
 * @author Griffin Dube
 * @brief Command Line argument parser
 * @version 1.0
 * @date 2021-06-01
 * Created on: 21-06-01
 *
 * @copyright (C) 2021 by Clemson University, Washington State University, Argonne National Laboratory
 * See LICENSE in top-level directory
 *
 */

#include "argparse.hh"
#include "document.hh"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <regex>
#include <string>
#include "../utils/format.hh"

using std::cerr;
using std::cout;
using std::endl;
using std::string;

const char* version_text  = "version: pre-alpha, build: 2021-06-01";
const int   version       = 210601;
const int   compatibility = 0;

bool isNumeric(string str) {
   for (int i = 0; i < str.length(); i++)
      if (isdigit(str[i]) == false)
         return false;
      return true;
}

void ArgParse::vecszDoc()
{
	cout << vecsz_short_doc << endl;
}

void ArgParse::vecszFullDoc()
{
	cout << format(vecsz_full_doc) << endl;
}

void ArgParse::PrintArgs()
{
	cout << log_cfg << "Input File:          " << files.input_file                                             << endl;
	cout << log_cfg << "Output File:         " << files.output_file                                            << endl;
	cout << log_cfg << "Compress?            " << string((szwf.lossy_compress) ? "Yes" : "No")                 << endl;
	cout << log_cfg << "Decompress?          " << string((szwf.lossy_decompress) ? "Yes" : "No")               << endl;
	cout << log_cfg << "Dry-Run?             " << string((szwf.lossy_dryrun) ? "Yes" : "No")                   << endl;
	cout << log_cfg << "Demo Dataset?        " << string((demo_dataset.empty()) ? "No" : demo_dataset)         << endl;
	cout << log_cfg << "Error Mode:          " << mode                                                         << endl;
	cout << log_cfg << "Data Type:           " << dtype                                                        << endl;
	cout << log_cfg << "Num. Dimensions:     " << ndim                                                         << endl;
	cout << log_cfg << "Error Bound:         " << eb                                                           << endl;
	cout << log_cfg << "Dimensions (x,y,z):  " << "(" << dim4._0 << "," << dim4._1 << "," << dim4._2 << ")"    << endl;
	cout << log_cfg << "Num. Blocks (x,y,z): " << "(" << nblk4._0 << "," << nblk4._1 << "," << nblk4._2 << ")" << endl;
	cout << log_cfg << "Dict. Size:          " << dict_size                                                    << endl;
	cout << log_cfg << "Radius:              " << radius                                                       << endl;
	cout << log_cfg << "Vector Length:       " << vector_length                                                << endl;
	cout << log_cfg << "Block Size:          " << block_size                                                   << endl;
	cout << log_cfg << "Verbose?             " << string((verbose) ? "Yes" : "No")                             << endl;
	cout << log_cfg << "Autotune?            " << string((szwf.autotune) ? "Yes" : "No")                       << endl;
	cout << log_cfg << "Num. Iterations:     " << num_iterations                                               << endl;
	cout << log_cfg << "Sample Percentage:   " << sample_percentage                                            << endl;
}

void ArgParse::CheckArgs()
{
	bool abort = false;

	// check input file
	if (files.input_file.empty())
	{
		cerr << log_err << "Not specifying input file!" << endl;
		abort = true;
	}

	// check dimensions
	if (SelfMultiple4(dim4) == 1 and not szwf.use_demo)
	{
		if (szwf.lossy_compress or szwf.lossy_dryrun)
		{
			cerr << log_err << "Incorrect input size(s)!" << endl;
			abort = true;
		}
	}

	// check compression vs. decompression vs. dry-run
	if (not szwf.lossy_compress and not szwf.lossy_decompress and not szwf.lossy_dryrun)
	{
		cerr << log_err << "Select compress (-z), decompress (-x), or dry-run (-r)!" << endl;
		abort = true;
	}

	// check datatype
	if (dtype != "f32" and dtype != "f64")
	{
		if (szwf.lossy_compress or szwf.lossy_dryrun)
		{
			cout << dtype << endl;
			cerr << log_err << "Not specifying data type!" << endl;
			abort = true;
		}
	}

	// check for valid vector length in bits
	if ((vector_length != -1   and
	     vector_length != 128  and
	     vector_length != 256  and
	     vector_length != 512) or
	     vector_length > MAX_VECTOR_LENGTH)
	{
		cerr << log_err << "Invalid vectorization parameter (-s)! Max Vector Length is " << MAX_VECTOR_LENGTH << "." << endl;
		abort = true;
	}

	// ensure reasonable dictionary size
	assert(dict_size <= 65536);

	// make sure dry-run and compress/decompress do not occur at same time
	if (szwf.lossy_dryrun and szwf.lossy_compress and szwf.lossy_decompress)
	{
		cerr << log_warn << "No need to dry-run, compress and decompress at the same time!" << endl;
		cerr << log_warn << "Will dry-run only." << endl << endl;
		szwf.lossy_compress   = false;
		szwf.lossy_decompress = false;
	}

	if (szwf.lossy_dryrun and szwf.lossy_compress)
	{
		cerr << log_warn << "No need to dry-run and compress at the same time!" << endl;
		cerr << log_warn << "Will dry-run only." << endl << endl;
		szwf.lossy_compress   = false;
	}

	if (szwf.lossy_dryrun and szwf.lossy_decompress)
	{
		cerr << log_warn << "No need to dry-run and decompress at the same time!" << endl;
		cerr << log_warn << "Will dry-run only." << endl << endl;
		szwf.lossy_decompress   = false;
	}

	if (not szwf.autotune and ((num_iterations >= 0) or (sample_percentage >= 0)))
	{
		cerr << log_err << "Autotune option not specified! Specify `-a` or `--autotune` to use `--sample` or `--num-iter` options." << endl;
		abort = true;
	}

	if (szwf.autotune)
	{
		if ((num_iterations <= 0) and (sample_percentage <= 0))
		{
			cerr << log_warn << "Invalid autotune configuration! Using default value." << endl;
			num_iterations = 10;
			sample_percentage = 1.0;
		}
		else if (num_iterations <= 0)
		{
			cerr << log_warn << "Invalid autotune configuration! Using default value." << endl;
			num_iterations = 10;
		}
		else if (sample_percentage <= 0)
		{
			cerr << log_warn << "Invalid autotune configuration! Using default value." << endl;
			sample_percentage = 1.0;
		}
	}

	if (abort)
	{
		vecszDoc();
		exit(-1);
	}

}

string  ArgParse::format(const string& s)
{
    std::regex  bful("@(.*?)@");
    std::string bful_text("\e[1m\e[4m$1\e[0m");
    std::regex  bf("\\*(.*?)\\*");
    std::string bf_text("\e[1m$1\e[0m");
    std::regex  ul(R"(_((\w|-|\d|\.)+?)_)");
    std::string ul_text("\e[4m$1\e[0m");
    std::regex  red(R"(\^\^(.*?)\^\^)");
    std::string red_text("\e[31m$1\e[0m");
    auto        a = std::regex_replace(s, bful, bful_text);
    auto        b = std::regex_replace(a, bf, bf_text);
    auto        c = std::regex_replace(b, ul, ul_text);
    auto        d = std::regex_replace(c, red, red_text);
    return d;
}

int  ArgParse::trap(int _status)
{
    this->read_args_status = _status;
    return read_args_status;
}

void ArgParse::ParseVecszArgs(int argc, char** argv)
{
	if (argc == 1) {
		vecszDoc();
		exit(0);
	}

	opath = "";

        auto str2int = [&](const char* s) {
            char* end;
            auto  res = std::strtol(s, &end, 10);
            if (*end) {
                const char* notif = "invalid option value, non-convertible part: ";
                cerr << log_err << notif << "\e[1m" << end << "\e[0m" << endl;
                cerr << string(log_null.length() + strlen(notif), ' ') << "\e[1m"  //
                     << string(strlen(end), '~')                                   //
                     << "\e[0m" << endl;
                trap(-1);
                return 0;  // just a placeholder
            }
            return (int)res;
        };

        auto str2fp = [&](const char* s) {
            char* end;
            auto  res = std::strtod(s, &end);
            if (*end) {
                const char* notif = "invalid option value, non-convertible part: ";
                cerr << log_err << notif << "\e[1m" << end << "\e[0m" << endl;
                cerr << string(log_null.length() + strlen(notif), ' ') << "\e[1m"  //
                     << string(strlen(end), '~')                                   //
                     << "\e[0m" << endl;
                trap(-1);
                return 0;  // just a placeholder
            }
            return (int)res;
        };

	int i = 1;
	while (i < argc)
	{
		if (argv[i][0] == '-')
		{
			auto long_opt = string(argv[i]);
			switch (argv[i][1])
			{
		                // ----------------------------------------------------------------
                                case '-':
                                    if (long_opt == "--help") goto tag_help;              // DOCUMENT
                                    if (long_opt == "--version") goto tag_version;        //
                                    if (long_opt == "--verbose") goto tag_verbose;        //
                                    if (long_opt == "--mode") goto tag_mode;              // COMPRESSION CONFIG
                                    if (long_opt == "--eb") goto tag_error_bound;         //
                                    if (long_opt == "--dict-size") goto tag_dict;         //
                                    if (long_opt == "--block-size") goto tag_block;       //
                                    if (long_opt == "--dtype") goto tag_type;             //
                                    if (long_opt == "--input") goto tag_input;            // INPUT
                                    if (long_opt == "--demo") goto tag_demo;              //
                                    if (long_opt == "--len") goto tag_len;                //
                                    if (long_opt == "--compress") goto tag_compress;      // WORKFLOW
                                    if (long_opt == "--zip") goto tag_compress;           //
                                    if (long_opt == "--decompress") goto tag_decompress;  //
                                    if (long_opt == "--unzip") goto tag_decompress;       //
                                    if (long_opt == "--dry-run") goto tag_dryrun;         //
                                    if (long_opt == "--output") goto tag_x_out;           //
				    if (long_opt == "--vector") goto tag_vector;          // VECTOR LENGTH
				    if (long_opt == "--autotune") goto tag_autotune;      //
				    if (long_opt == "--num-iter") {
					    if (i + 1 <= argc) {
                                            	num_iterations = str2int(argv[++i]);
					    }
				    }
				    if (long_opt == "--sample") {
					    if (i + 1 <= argc) {
                                        	char* end;
                                        	this->sample_percentage = std::strtod(argv[++i], &end);
					    }
				    }
                                case 'z':
                                tag_compress:
                                    szwf.lossy_compress = true;
                                    break;
                                case 'x':
                                tag_decompress:
                                    szwf.lossy_decompress = true;
                                    break;
                                case 'r':
                                tag_dryrun:
                                    szwf.lossy_dryrun = true;
                                    break;
                                // COMPRESSION CONFIG
                                case 'm':  // mode
                                tag_mode:
                                    if (i + 1 <= argc) mode = string(argv[++i]);
                                    break;
                                // INPUT
                                case 'l':
                                tag_len:
                                    if (i + 1 <= argc) {
                                        std::stringstream   datalen(argv[++i]);
                                        std::vector<string> dims;
                                        while (datalen.good()) {
                                            string substr;
                                            getline(datalen, substr, ',');
                                            dims.push_back(substr);
                                        }
                                        ndim = dims.size();
                                        if (ndim == 1) {  //
                                            auto d0 = str2int(dims[0].c_str());
                                            dim4    = {d0, 1, 1, 1};
                                        }
                                        if (ndim == 2) {  //
                                            auto d0 = str2int(dims[0].c_str()), d1 = str2int(dims[1].c_str());
                                            dim4 = {d0, d1, 1, 1};
                                        }
                                        if (ndim == 3) {
                                            auto d0 = str2int(dims[0].c_str()), d1 = str2int(dims[1].c_str());
                                            auto d2 = str2int(dims[2].c_str());
                                            dim4    = {d0, d1, d2, 1};
                                        }
                                        if (ndim == 4) {
                                            auto d0 = str2int(dims[0].c_str()), d1 = str2int(dims[1].c_str());
                                            auto d2 = str2int(dims[2].c_str()), d3 = str2int(dims[3].c_str());
                                            dim4 = {d0, d1, d2, d3};
                                        }
                                    }
                                    break;
                                case 'i':
                                tag_input:
                                    if (i + 1 <= argc)
				    {
					    files.input_file = string(argv[++i]);
					    if (files.output_file == "") files.output_file = files.input_file + ".sz";
				    }
                                    break;
                                    // alternative output
                                case 'o':
                                tag_x_out:
                                    if (i + 1 <= argc) files.output_file = string(argv[++i]);
                                    break;
                                // demo datasets
                                case 'D':
                                tag_demo:
                                    if (i + 1 <= argc) {
                                        szwf.use_demo = true;
                                        demo_dataset        = string(argv[++i]);
                                    }
                                    break;
                                // DOCUMENT
                                case 'h':
                                tag_help:
                                    vecszFullDoc();
                                    exit(0);
                                    break;
                                case 'v':
                                tag_version:
                                    // TODO
                                    cout << log_info << version_text << endl;
				    exit(0);
                                    break;
                                // COMPRESSION CONFIG
                                case 't':
                                tag_type:
                                    if (i + 1 <= argc) {
                                        string s = string(string(argv[++i]));
                                        if (s == "f32" or s == "fp4")
                                            dtype = "f32";
                                        else if (s == "f64" or s == "fp8")
                                            dtype = "f64";
                                    }
                                    break;
                                case 'e':
                                tag_error_bound:
                                    if (i + 1 <= argc) {
                                        char* end;
                                        this->eb = std::strtod(argv[++i], &end);
                                    }
                                    break;
                                case 'V':
                                tag_verbose:
                                    verbose = true;
                                    break;
                                case 'd':
                                tag_dict:
                                    if (i + 1 <= argc) {
                                        dict_size = str2int(argv[++i]);
                                        radius    = dict_size / 2;
                                    }
                                    break;
				case 'b':
				tag_block:
				    if (i + 1 <= argc) {
					block_size = str2int(argv[++i]);
				    }
				    break;
				case 's':
				tag_vector:
				    if (i + 1 <= argc) {
					string s = string(string(argv[++i]));
					if (isNumeric(s)) vector_length = std::stoi(s);
					else if (s == "None"  ) vector_length = -1;
					else if (s == "AVX"   ) vector_length = 128;
					else if (s == "AVX2"  ) vector_length = 256;
					else if (s == "AVX512") vector_length = 512;
					else vector_length = 0;
				    }
				    break;
				case 'a':
				tag_autotune:
				    szwf.autotune = true;
				    num_iterations = 10;
				    sample_percentage = 1.0;
				    break;
                                default:
                                    const char* notif_prefix = "invalid option value at position ";
                                    char*       notif;
                                    int         size = asprintf(&notif, "%d: %s", i, argv[i]);
                                    cerr << log_err << notif_prefix << "\e[1m" << notif << "\e[0m"
                                         << "\n";
                                    cerr << string(log_null.length() + strlen(notif_prefix), ' ');
                                    cerr << "\e[1m";
                                    cerr << string(strlen(notif), '~');
                                    cerr << "\e[0m\n";
                                    trap(-1);
                            }
        	}
        else {
            const char* notif_prefix = "invalid option at position ";
            char*       notif;
            int         size = asprintf(&notif, "%d: %s", i, argv[i]);
            cerr << log_err << notif_prefix << "\e[1m" << notif
                 << "\e[0m"
                    "\n"
                 << string(log_null.length() + strlen(notif_prefix), ' ')  //
                 << "\e[1m"                                                //
                 << string(strlen(notif), '~')                             //
                 << "\e[0m\n";
            trap(-1);
        }
        i++;
    }

    // phase 1: check grammar
    if (read_args_status != 0) {
        cout << log_info << "Exiting..." << endl;
        // after printing ALL argument errors
        exit(-1);
    }

    // phase 2: check if meaningful
    CheckArgs();

    LogAll(log_info,"parsing command-line args");
}
