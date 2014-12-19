/*
 *	Author : Shujia Huang
 *	Date   : 2012-03-22
 *
 *	Modify : 2012-07-31 As fosmid reads have been aliged to the reference. I don't need to analyse the pooling ID any more, 
			 which will been delete in this version, and just finding the continue region as fosmid by bamfomat.
			 Add mapping quality filter parameter.
			 Add the function to read BAM file
 *	Modify : 2012-03-28 Upgrade the function, and add one more query_length_list_file as input!
 *	Modify : 2012-03-26	Anchor fosmid not pooling id!!!
 *	Just need one input file.
 *	Input file format: BAM
 */
#include <iostream>
#include <assert.h>
#include <iomanip>
#include <fstream>
#include <getopt.h>
#include <vector>
#include <map>
#include <set>
#include "anchor_fosmid.h"
#include "api/BamReader.h"

using namespace std;
using namespace BamTools;

void usage( const char* prog ) {
	cerr << "Version 1.0.1\nUsage: " << prog << " [Options] -i [in_bam_file]" << endl;
	cerr << "   Options:                                                          \n"
         << "       -i\tInput Sorted BAM file.                                    \n"
		 << "       -o\tOutput file.											  \n"
         << "       -p\tPooling ID.                                               \n"
		 << "       -q\tMappinp quality threshold. [20]							  \n"
		 << "       -r\tCoverage Ratio. [0.5]									  \n"
         << "       -l\tContinue region length threshold. [10000]                 \n"
         << "       -d\tThe average depth threshold. [5]                          \n"
         << "       -g\tThe gap length between two nearby continue region. [20000]\n"
		 << "       -t\tAdd Header label in output file. [false]				  \n"
		 << "       -h\tOutput this help. 										  \n"
         << endl;
	exit(1);
}

void decid_fosmid ( vector< SamExt >& sam_array, int& fosmid_num, const double ratio, double depth, ofstream& O1, ofstream& O2 );
void out_filter   ( SamExt & sam );

int main( int argc, char* argv[] ) {

	bool no_title( true );
	char c;
	double ratio(0.5);
	int length_threshold(10000), mapping_quality_threshold( 20 ), gap_threshold(20000), depth_threshold(5);
	string infile, qlength_file, pooling_id, outfile;
	while ( ( c = getopt( argc, argv, "i:l:r:p:q:d:g:o:th" ) ) != -1 ) {
		switch ( c ) {
			case 'i' : infile       = optarg; break;
			case 'o' : outfile      = optarg; break;
			case 'p' : pooling_id   = optarg; break;
			case 't' : no_title     = false;  break;
			case 'h' : usage( argv[0] );      break;
			case 'q' : mapping_quality_threshold = atoi(optarg); break;
			case 'l' : length_threshold          = atoi(optarg); break;
			case 'd' : depth_threshold			 = atoi(optarg); break;
			case 'g' : gap_threshold			 = atoi(optarg); break;
			case 'r' : ratio                     = atof(optarg); break;
			default  : 
				usage ( argv[0] );
		}
	}
	if ( infile.empty() || outfile.empty() || pooling_id.empty() ) usage( argv[0] );
	cerr << "Important Parameter: -L " << length_threshold << " -r "<< ratio << " -q " << mapping_quality_threshold
		 << " -p " << pooling_id << " -g " << gap_threshold << " -d "<< depth_threshold << endl;
	string outfile2 = outfile + ".region";
	
	vector< SamExt > sam_array;
	SamExt  sam;					// SamExt class
	SamLine samline;				// SamLine class
	string tmp;
	int fosmid_num(0), max_end(0);

	ofstream O1( outfile.c_str() ), O2( outfile2.c_str() );
	if ( !O1 ) { cerr << "Cannot write to file : " << outfile  << endl; exit(1); }
	if ( !O2 ) { cerr << "Cannot write to file : " << outfile2 << endl; exit(1); }
	if ( !no_title ) {
		O1 << "#FosmidID\tRefID\tRegStart\tRegEnd\tReadID\tReadStart\tReadEnd\tMappingInf" << endl;
		O2 << "#FosmidID\tRefID\tRegStart\tRegEnd\tRegLength\tPoolingID\tCoverLength\tAverageDepth\tCoverReadNum" << endl;
	}

	BamReader I; // bam input file handle
	BamAlignment al;
	if ( !I.Open( infile ) ) cerr << "[ERROR]: " << I.GetErrorString() << endl;
	while ( I.GetNextAlignment( al ) ) {
	//HUMuvjD320poolingDFAAPEI-10_Index8_Index2_C958  16   scaffold27   24677427   37  116M1D67M    *    0 0   ...	

		if ( !al.IsMapped() || al.MapQuality < mapping_quality_threshold ) continue;
		samline._RID      = al.Name; samline._Flag = al.AlignmentFlag;
        samline._ref_id   = I.GetReferenceData()[al.RefID].RefName;
        samline._position = al.Position + 1;     // Position (0-base starts in BamTools), but I need 1-base starts
        samline._mapQ     = al.MapQuality;
        // MateRefID == -1 means mate read is unmapping
        samline._XorD     = ( al.MateRefID > -1 ) ? I.GetReferenceData()[al.MateRefID].RefName : "*";
        samline._coor     = al.MatePosition + 1; // Position (0-base starts in BamTools), but I need 1-base starts
        samline._seq      = al.QueryBases;
        samline._insert_size = abs (al.InsertSize);
		if ( samline._ref_id == "BIG_ID_CAT" ) continue; // Ignore "BIG_ID_CAT"

		// get cigar;
        samline._cigar = itoa(al.CigarData[0].Length); samline._cigar.append( 1, al.CigarData[0].Type );
        for ( size_t i(1); i < al.CigarData.size(); ++i ) {
            samline._cigar += itoa(al.CigarData[i].Length);
            samline._cigar.append( 1, al.CigarData[i].Type );
        }

		sam.assign( &samline, pooling_id );

		if ( !sam_array.empty() ) {

			if ( sam._ref_id ==  sam_array.back()._ref_id ) {

				if ( sam._position < sam_array.back()._position ) {
					cerr << "[ERROR]: Your file " << infile << " should be sorted!!" << endl;
					exit(1);
				}

				max_end = ( max_end < sam_array.back().ref_end() ) ? sam_array.back().ref_end() : max_end;
				if ( sam._position - max_end > gap_threshold ) {

					if ( max_end - sam_array.front().ref_start() + 1 >= length_threshold ) {

						while ( !sam_array.empty() && ( max_end - sam_array.front().ref_start() + 1 >= length_threshold) ) {
							decid_fosmid(sam_array, fosmid_num, ratio, depth_threshold, O1, O2);
							for ( size_t i(0); i < sam_array.size(); ++i ) {
								if ( i == 0 ) { max_end = sam_array[i].ref_end(); continue; }
								if ( max_end < sam_array[i].ref_end() ) max_end = sam_array[i].ref_end();
							}
						}
						if ( !sam_array.empty() && ( max_end - sam_array.front().ref_start() + 1 < length_threshold) ) {
							for ( size_t j(0); j < sam_array.size(); ++j ) out_filter ( sam_array[j] );
							sam_array.clear();
						} 
					} else {
						for ( size_t j(0); j < sam_array.size(); ++j ) out_filter ( sam_array[j] );
						sam_array.clear();
					}
					assert( sam_array.empty() );
				}
			} else {

				if ( max_end - sam_array.front().ref_start() + 1 >= length_threshold ) { 

					while ( !sam_array.empty() && ( max_end - sam_array.front().ref_start() + 1 >= length_threshold) ) {
						decid_fosmid(sam_array, fosmid_num, ratio, depth_threshold, O1, O2);
						for ( size_t i(0); i < sam_array.size(); ++i ) {
							if ( i == 0 ) { max_end = sam_array[i].ref_end(); continue; }
							if ( max_end < sam_array[i].ref_end() ) max_end = sam_array[i].ref_end();
						}
					}
					if ( !sam_array.empty() && ( max_end - sam_array.front().ref_start() + 1 < length_threshold) ) {
						for ( size_t j(0); j < sam_array.size(); ++j ) out_filter ( sam_array[j] );
						sam_array.clear();
					}
				} else {
					for ( size_t j(0); j < sam_array.size(); ++j ) out_filter ( sam_array[j] );
					sam_array.clear();
				}
				assert( sam_array.empty() );
				max_end = sam.ref_end();
			}
			sam_array.push_back( sam );
		} else {
			max_end = sam.ref_end();
			sam_array.push_back( sam );
		}
	}
	I.Close();

	if ( !sam_array.empty() ) { 
		if ( max_end - sam_array.front().ref_start() + 1 >= length_threshold ) { 

			while ( !sam_array.empty() && ( max_end - sam_array.front().ref_start() + 1 >= length_threshold) ) {
				decid_fosmid(sam_array, fosmid_num, ratio, depth_threshold, O1, O2);
				for ( size_t i(0); i < sam_array.size(); ++i ) {
					if ( i == 0 ) { max_end = sam_array[i].ref_end(); continue; }
					if ( max_end < sam_array[i].ref_end() ) max_end = sam_array[i].ref_end();
				}
			}
			if ( !sam_array.empty() && ( max_end - sam_array.front().ref_start() + 1 < length_threshold) ) {
				for ( size_t j(0); j < sam_array.size(); ++j ) out_filter ( sam_array[j] );
				sam_array.clear();
			}
		} else {
			for ( size_t j(0); j < sam_array.size(); ++j ) out_filter ( sam_array[j] );
			sam_array.clear();
		}
		assert( sam_array.empty() );
	}

	O1.close();
	O2.close();

	cerr << "\nDone.\nOutput file list:\n" << outfile << "\n" << outfile2 << endl;

	return 0;
}

void decid_fosmid( vector< SamExt >& array, int& fosmid_num, const double ratio, double depth_threshold, ofstream& O1, ofstream& O2 ) {

	assert( !array.empty() );

	string pool_id ( array[0].pool_id()    ), ref_id( array[0]._ref_id );
	int cov_length ( array[0].map_length() );
	int region_start( array[0].ref_start() ); // 2012-4-16
    int region_end  ( array[0].ref_end()   ); // 2012-4-16

	int query2endp ( array[0].ref_end() ); 

	double accum_depth(0.0);
	for ( size_t i(0); i < array.size(); ++i ) {

		if ( array[i].ref_start() < region_start ) region_start = array[i].ref_start(); // 2012-4-16
		if ( array[i].ref_end()   > region_end   ) region_end   = array[i].ref_end();   // 2012-4-16

		// Skip the overlap part!
		cov_length += array[i].cov_length_after( query2endp );
		if ( array[i].ref_end() > query2endp ) query2endp = array[i].ref_end();

		// Calculate the coverage depth
		for ( size_t j(0); j < array[i].cigar_seq().length(); ++j ) {
			if ( array[i].cigar_seq()[j] == 'M' ) ++accum_depth;
		}
	}

	int region_length = region_end - region_start + 1; // 2012-4-16
	double covg       = double( cov_length ) / region_length;
	double depth      = accum_depth / region_length;
	/* Debug
	if ( covg > 1 ) 
		cerr << "[WARNING] Coverage length is bigger than the region!! " 
             << cov_length << " > " << region_length << endl;
	*/
	if ( covg >= ratio && depth >= depth_threshold ) {

		++fosmid_num;

		O2.setf(ios::fixed);
		O2 << fosmid_num << "\t" << ref_id  << "\t" << region_start         << "\t"       << region_end << "\t"
		   << region_end - region_start + 1 << "\t" << pool_id      << "\t" << cov_length << "\t"
		   << setprecision(2)    << depth   << "\t" << array.size() << endl;
		for ( size_t j(0); j < array.size(); ++j ) {

			O1  << fosmid_num << "\t" << array[j]._ref_id			     << "\t" 
				<< array[j].ref_start() << "\t" << array[j].ref_end()    << "\t" 
				<< array[j].read_id()   << "\t" << array[j].read_start() << "\t" 
				<< array[j].read_end()  << "\t"
                << array[j]._RID        << ":" << array[j]._ref_id << ":" << array[j]._position << ":"
				<< array[j]._Flag       << ":" << array[j]._mapQ   << ":" << array[j]._cigar    << ":"
                << array[j]._XorD       << ":" << array[j]._coor   << ":"
                << array[j]._insert_size
                << endl;
		}
		array.clear();
	} else { 
		// Skip the first line
		out_filter ( array.front() );
		array.erase( array.begin() );
	}

	return;
}

void out_filter ( SamExt & sam ) {

	cout << sam._ref_id << "\t"
         << sam.ref_start() << "\t" << sam.ref_end()    << "\t"
         << sam.read_id()   << "\t" << sam.read_start() << "\t"
         << sam.read_end()  << "\t"
         << sam._RID        << ":" << sam._ref_id << ":" << sam._position << ":"
         << sam._Flag       << ":" << sam._mapQ   << ":" << sam._cigar    << ":"
         << sam._XorD       << ":" << sam._coor   << ":"
         << sam._insert_size
         << endl;

	return;
}












