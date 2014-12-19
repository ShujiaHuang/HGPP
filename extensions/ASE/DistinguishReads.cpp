/*
 *	Author : Shujia Huang
 *	Date   : 2012-06-26
 *
 *	
 *	Distinguish the mapping reads into two parts by the haplotype phasing result.
 *	Need two input file for this program.	
 *	1) Input  Format: bam 
 *  2) Phasing haplotype: 
		scaffold99_9    2980260 2980453 2/2     2980304,2980438 G,A     R,R     1       @A81EE0ABXX_5_22
 *	Just one output file.
 *	Output Format: bam
 *
 *	Modify: 2013-8-16 16:19:45	Add to output the mate pair reads, which had been missed for a long time!!
 *	Modify: 2013-8-15 18:35:31	Add homo read output into a independent file, names: "perfix.homo.bam"
 *	Modify: version0.1(2012-06-29) Add an parameter(-B) for BS reads, for this matter we have to add the other parameter for reference fa(-r)
 *	Modify: version0.1(2012-06-26) Create the programming functions. Use new program type since today(Camel and Pascal)!!
 */
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <map>
#include <assert.h>
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "anchor_variation.h"

using namespace std;
using namespace BamTools;

void Usage( const char* prog ) {
	cerr << "Usage: " << prog << " [-i input_bamfile] [-f input_phasing_hap] [-o output_bamfile_prefix] <-B For BS read.> <-r reference. Use with -B parameter.>" << endl;
	exit(1);
}
void readRef( const char* file, map< string, string >& refSeq );
void LoadHapSNP ( const char* phasingFile, map< string, bool > & hapSNP );
void Distinguish( string mappingFile, map< string, bool > & hapSNP, map< string, string > & refSeq, string outfilePrefix );
void CreateMap  ( const string refID, const string position, const string allele, const string genotype, map< string, bool >& hapSNP );
void modifyBSreadBases ( string refID, int refPos, int readPos, string cigarStr, string refStr, string & queryBases, map< string, bool > & hapSNP, bool isC2T ); 
int Decide( SamExt & sam, map< string, bool > & hapSNP );

int main ( int argc, char* argv[] ) {

	char c;
	bool isBS( false );
	string mappingFile, phasingFile, outfilePrefix, referenceFile;
	while ( ( c = getopt( argc, argv, "i:f:o:r:Bh" ) ) != -1 ) {
		switch ( c ) {
			case 'i' : mappingFile   = optarg; break;
			case 'f' : phasingFile   = optarg; break;
			case 'o' : outfilePrefix = optarg; break;
			case 'r' : referenceFile = optarg; break;
			case 'B' : isBS          = true;   break;
			case 'h' : Usage( argv[0] );
			default  : cerr << "Parameter error! -" << c << endl; Usage( argv[0] );
		}
	}
	if ( argc == 1 || mappingFile.empty() || phasingFile.empty() || outfilePrefix.empty() ) Usage( argv[0] );
	if ( (isBS && referenceFile.empty()) || ( !isBS && !referenceFile.empty() ) ) Usage( argv[0] );

	map< string, string > refSeq;
	if ( isBS ) readRef( referenceFile.c_str(), refSeq);

	map< string, bool > hapSNP; // Define: 'true' for haplotype1 and 'false' for haplotype2
	LoadHapSNP ( phasingFile.c_str(), hapSNP );
	Distinguish( mappingFile, hapSNP, refSeq, outfilePrefix );
	return 0;
}

// Input File Format: fa format
void readRef( const char* file, map< string, string >& refSeq ) {
// fa file 
    cerr << "Loading reference ... " << endl;
    ifstream I ( file );
    if ( !I ) {
        cerr << "Cannot open the file : " << file << endl;
        exit(1);
    }

    string tmp, key, seq;
    int refNum (0);
    while ( 1 ) {
        I >> tmp;
        if ( I.eof() ) break;
        if ( tmp[0] != '>' ) {

            seq += tmp;
        } else {

            ++refNum;
            if ( !seq.empty() ) refSeq.insert( make_pair( key, seq) );
            seq.clear();
            key.assign( tmp, 1, string::npos );
        }
        getline( I, tmp, '\n' );
    }
    I.close();

    if ( !seq.empty() ) refSeq.insert( make_pair( key, seq) );
    cerr << "Reference loaded done. And there's  " << refNum << " reference fa sequeces. " << refSeq.size() << endl;
}

/* 
Input file format: 
scaffold99_13   3110605 3152184 4/5     3124348,3132428,3133369,3133778,3147227 -,T,T,G,+T       Y,Y,K,R,+T|-C       1       HUMuv
*/
void LoadHapSNP( const char* phasingFile, map< string, bool >& hapSNP ) {
	
	ifstream I ( phasingFile );
	assert ( I.is_open() );      // Assert phasingFile is open, if open fail program will exit.

	string tmp, refID;
	vector< string > genotype, allele, position;
	while ( 1 ) {
		
		I >> tmp;
		if ( I.eof() ) break;
		refID.assign( tmp, 0, tmp.find_first_of( "_" ) );
		I >> tmp >> tmp >> tmp >> tmp;  split( ",", tmp, position );//Position: 3124348,3132428
		I >> tmp;						split( ",", tmp, allele   );//Allele  : T,T
		I >> tmp;						split( ",", tmp, genotype );//genotype: Y,Y
		getline( I, tmp, '\n');
		assert( (position.size() == allele.size()) && (genotype.size() == allele.size()) );//Size should be the same or input file error

		for ( size_t i(0); i < genotype.size(); ++i ) {
			if ( allele[i][0] == '-' || genotype[i][0] == '+' || genotype[i][0] == '-' ) continue;
			CreateMap( refID, position[i], toupper(allele[i]), toupper(genotype[i]), hapSNP );
		}
		/* To be continue here */
	}

	I.close();
}

void CreateMap( const string refID, const string position, const string allele, const string genotype, map< string, bool >& hapSNP ) {

	map< char, map<char, char> > alleleTable;
    alleleTable['S']['C'] = 'G'; alleleTable['S']['G'] = 'C';
    alleleTable['M']['A'] = 'C'; alleleTable['M']['C'] = 'A';
    alleleTable['K']['G'] = 'T'; alleleTable['K']['T'] = 'G';
    alleleTable['R']['A'] = 'G'; alleleTable['R']['G'] = 'A';
    alleleTable['W']['A'] = 'T'; alleleTable['W']['T'] = 'A';
    alleleTable['Y']['C'] = 'T'; alleleTable['Y']['T'] = 'C';

	assert (alleleTable.count( genotype[0] ) && alleleTable[genotype[0]].count( allele[0] ) );//Assert the allele match the alleleTable
	string key;
	key = refID + ":" + position + ":" + allele; // Key=scaffold20:210981:A

	if ( hapSNP.count( key ) ) return;

	hapSNP[key] = true;  // for Hap1
	key = refID + ":" + position + ":" + char2str( alleleTable[genotype[0]][allele[0]] );
	hapSNP[key] = false; // for Hap2
	
	return;
}

/*
Input file format: BAM file
*/
void Distinguish( string mappingFile, map< string, bool > & hapSNP, map< string, string > & refSeq, string outfilePrefix ) {
// [WARNING] If refSeq is not empty than it's BS data, otherwirse is resequencing data.
	string outfile1 = outfilePrefix + ".1.bam";
	string outfile2 = outfilePrefix + ".2.bam";
	string outHomoF = outfilePrefix + ".homo.bam";
	string outUdeci = outfilePrefix + ".ambiguity.bam"; // The ambiguity reads or the reads which has't any SNP.

	BamReader h_I; // bam input file handle
	if ( !h_I.Open( mappingFile ) ) cerr << "[ERROR]: " << h_I.GetErrorString() << endl;

	// "header" and "references" from BAM files, these are required by BamWriter
	const SamHeader header     = h_I.GetHeader();
	const RefVector references = h_I.GetReferenceData();
	
	BamWriter h_O1, h_O2, h_U, h_H;
	if ( !h_O1.Open( outfile1, header, references ) ) { cerr << "Cannot open output BAM file: " << outfile1 << endl; exit(1); }
	if ( !h_O2.Open( outfile2, header, references ) ) { cerr << "Cannot open output BAM file: " << outfile2 << endl; exit(1); }
	if ( !h_U.Open ( outUdeci, header, references ) ) { cerr << "Cannot open output BAM file: " << outUdeci << endl; exit(1); }
	if ( !h_H.Open ( outHomoF, header, references ) ) { cerr << "Cannot open output BAM file: " << outHomoF << endl; exit(1); }

	int readsNumberRecord(0);

	SamLine samline;                                       // Samline class
	SamExt      sam;
	BamAlignment al;
	map< string, pair<BamAlignment, SamExt> > firstMateAl; // record the first mate reads alignment, HIstory problem to be like this struct!!

	string refstr; 					                       // Just For BS data.
	bool isC2T ( false );			                       // Just For BS data.
	while ( h_I.GetNextAlignment( al ) ) {

		++readsNumberRecord;
		if ( readsNumberRecord % 1000000 == 0 ) 
			cerr << "Have been dealed " << readsNumberRecord << " lines. " << local_time ();

		if ( !al.IsMapped() ) continue;
        //if ( al.InsertSize == 0 || al.RefID != al.MateRefID ) continue;

		samline._RID      = al.Name; samline._Flag= al.AlignmentFlag;
        samline._ref_id   = h_I.GetReferenceData()[al.RefID].RefName;
        samline._position = al.Position + 1;     // Position (0-base starts in BamTools), but I need 1-base starts
        samline._mapQ     = al.MapQuality;
		// MateRefID == -1 means mate read is unmapping
        samline._XorD     = ( al.MateRefID > -1 ) ? h_I.GetReferenceData()[al.MateRefID].RefName : "*"; 
        samline._coor     = al.MatePosition + 1; // Position (0-base starts in BamTools), but I need 1-base starts
        samline._seq      = al.QueryBases;
        samline._insert_size = abs (al.InsertSize);
		if ( samline._ref_id.compare( "BIG_ID_CAT" ) == 0 ) continue; // Ignore "BIG_ID_CAT"

		// get cigar;
        samline._cigar = itoa(al.CigarData[0].Length); samline._cigar.append( 1, al.CigarData[0].Type );
		for ( size_t i(1); i < al.CigarData.size(); ++i ) {
            samline._cigar += itoa(al.CigarData[i].Length);
            samline._cigar.append( 1, al.CigarData[i].Type );
        }

        sam.assign( &samline );
		/*********************************** For BS Data *********************************************/
		if ( !refSeq.empty() ) { // If the data is BS data, we should modify the QueryBases.
			if ( !refSeq.count( samline._ref_id ) ) { 
				cerr << "[ERROR]There's no such reference in the reference file. " << samline._ref_id << endl;
				exit(1);
			}
			if ( al.IsFirstMate() && !al.IsReverseStrand() ) { 
				isC2T = true; 
			} 
			else if ( al.IsFirstMate() && al.IsReverseStrand() ) { 
				isC2T = false; 
			} 
			else if ( al.IsSecondMate() && !al.IsReverseStrand() ) {
				isC2T = false;
			}
			else if ( al.IsSecondMate() && al.IsReverseStrand() ) {
				isC2T = true;
			} else {
				cerr << "[ERROR MATCH] " << endl; exit(1);
			}
			refstr.assign( refSeq[samline._ref_id], sam.ref_start() - 1, sam.ref_end() - sam.ref_start() + 1 );
			modifyBSreadBases( samline._ref_id, sam.ref_start (), sam.read_start(), sam.cigar_seq(), refstr, sam._seq, hapSNP, isC2T );
		}
		/********************************** End For BS Data *******************************************/

		// Consider the mate pair reads
		if ( !firstMateAl.count(al.Name) && (al.MateRefID > -1) ) { 

			firstMateAl[al.Name] = std::make_pair( al, sam );
		} else { // Consider the mate pair reads

			if ( !firstMateAl.count(al.Name) ) {

				switch ( Decide( sam, hapSNP ) ) { 
					case 1 : h_O1.SaveAlignment( al ); break; // Hap1
					case 2 : h_O2.SaveAlignment( al ); break; // Hap2
					case 0 : h_U.SaveAlignment ( al ); break; // Ambiguity
					default: // This alignment didn't contain any hete SNP.
							 h_H.SaveAlignment ( al );        // Homozygous reads
				}
			} else {
				
				int mark1 = Decide( firstMateAl[al.Name].second, hapSNP );
				int mark2 = Decide( sam, hapSNP );

				if ( mark1 == 1 && mark2 == 1 ) {
					
					h_O1.SaveAlignment( firstMateAl[al.Name].first );
					h_O1.SaveAlignment( al );
				} else if ( (mark1 == 1 && mark2 ==  0) || (mark1 == 0  && mark2 == 1) ) {
					
					h_O1.SaveAlignment( firstMateAl[al.Name].first );
					h_O1.SaveAlignment( al );
				} else if ( (mark1 == 1 && mark2 == -1) || (mark1 == -1 && mark2 == 1) ) {

					h_O1.SaveAlignment( firstMateAl[al.Name].first );
					h_O1.SaveAlignment( al );
				} else if ( mark1 == 2 && mark2 == 2 ) {

					h_O2.SaveAlignment( firstMateAl[al.Name].first );
					h_O2.SaveAlignment( al );
				} else if ( (mark1 == 2 && mark2 == 0 ) || (mark1 == 0  && mark2 == 2) ) {
					
					h_O2.SaveAlignment( firstMateAl[al.Name].first );
					h_O2.SaveAlignment( al );
				} else if ( (mark1 == 2 && mark2 == -1) || (mark1 == -1 && mark2 == 2) ) {

					h_O2.SaveAlignment( firstMateAl[al.Name].first );
					h_O2.SaveAlignment( al );
				} else if ( mark1 == -1 && mark2 == -1 ) {

					h_H.SaveAlignment ( firstMateAl[al.Name].first );
					h_H.SaveAlignment ( al );
				} else {
					
					h_U.SaveAlignment ( firstMateAl[al.Name].first );
					h_U.SaveAlignment ( al );
				}

				firstMateAl.erase( al.Name );
			}
		}
	}

	cerr << "------------  Remaind size: " << firstMateAl.size() << endl;
	for (map< string, pair<BamAlignment, SamExt> >::iterator it( firstMateAl.begin() ); it != firstMateAl.end(); ++it ) {

		switch ( Decide( it->second.second, hapSNP ) ) {
			case 1 : h_O1.SaveAlignment( it->second.first ); break; // Hap1
			case 2 : h_O2.SaveAlignment( it->second.first ); break; // Hap2
			case 0 : h_U.SaveAlignment ( it->second.first ); break; // Ambiguity
			default: // This alignment didn't contain any hete SNP.
					 h_H.SaveAlignment ( it->second.first );        // Homozygous reads
		}
	}

	h_I.Close();
	h_U.Close();
	h_H.Close();
	h_O1.Close();
	h_O2.Close();

	cerr << ">>>>>>>>>>>>> All Done <<<<<<<<<<<<<<" << endl;
	cerr << "Write to output file: " << outfile1    << endl;
	cerr << "Write to output file: " << outfile2    << endl;
	cerr << "Write to output file: " << outHomoF    << endl;
	cerr << "Write to output file: " << outUdeci    << endl;

	return;
}

void modifyBSreadBases ( string refID, int refPos, int readPos, string cigarStr, string refBases, string & queryBases, map< string, bool > & hapSNP, bool isC2T ) 
{
	// "hapSNP" uses for checking the base whether is SNP, if the base is a SNP, then we should not modify it. 
	int querySeqIndex(0), refSeqIndex(0); 
	string key;
	for ( size_t i(0); i < cigarStr.length(); ++i ) {
		switch ( cigarStr[i] ){
			case 'M' :

				if ( querySeqIndex >= queryBases.length() ) { cerr << "[BUG] querySeqIndex Overflow. " << endl; exit(1); }
				if ( refSeqIndex   >= refBases.length()   ) { cerr << "[BUG] refSeqIndex   Overflow. " << endl; exit(1); }

				queryBases[querySeqIndex] = toupper(queryBases[querySeqIndex]);
				refBases[refSeqIndex]     = toupper(refBases[refSeqIndex]    );
				key = refID + ":" + itoa( refPos ) + ":" + char2str(queryBases[querySeqIndex]); // scaffold20:201290:A
				if ( hapSNP.count( key ) ) continue; // SNP base should not modify it.
				/* To Continue Here */
				if ( isC2T ) { // For C->T reads
					if ( queryBases[querySeqIndex] == 'T' && refBases[refSeqIndex] == 'C' ) { queryBases[querySeqIndex] = 'C'; }
				} else {       // For G->A reads
					if ( queryBases[querySeqIndex] == 'A' && refBases[refSeqIndex] == 'G' ) { queryBases[querySeqIndex] = 'G'; }
				}
				++readPos; ++refPos; ++querySeqIndex; ++refSeqIndex;
				break;
			case 'I' : ++querySeqIndex; ++readPos; break;
			case 'S' : ++querySeqIndex; ++readPos; break;
			case 'D' : ++refSeqIndex  ; ++refPos ; break;
			case 'N' : ++refSeqIndex  ; ++refPos ; break;
			default  :
                cerr << "[BUG] Undefined to deal with this mark : " << cigarStr[i] << "\t" << cigarStr[i] << endl;
                exit(1);
		}
	}
}

int Decide( SamExt & sam, map< string, bool > & hapSNP ) {

	int readPos( sam.read_start() );
    int refPos ( sam.ref_start () );
    int index(0), one(0), two(0);
    string key;
	for ( size_t i(0); i < sam.cigar_seq().length(); ++i ) {
		
		switch( sam.cigar_seq()[i] ) {
			case 'M' : 
				key = sam._ref_id + ":" + itoa( refPos ) + ":" + char2str(toupper(sam._seq[index])); // scaffold20:201290:A
				if ( hapSNP.count( key ) ) {
					if ( hapSNP[key] ) {
						++one;
					} else {
						++two;
					}
				}
				++readPos; ++refPos; ++index;
				break;
			case 'I' : ++index; ++readPos; break;
			case 'S' : ++index; ++readPos; break;
			case 'D' : ++refPos;           break;
			case 'N' : ++refPos;           break;
			default  :
                cerr << "[BUG] Undefined to deal with this mark : " << sam.cigar_seq()[i] << "\t" << sam._cigar << endl;
                exit(1);
		}
	}

	int mark;
	if ( one > two ) {
		mark = 1;
	} 
	else if ( one < two ) {
		mark = 2;
	}
	else if ( one > 0 && one == two ) {
		mark = 0;
	} else { // one == two = 0
		if ( one > 0 || two > 0 ) { cerr << "[BUG] one and two are not equare to 0. " << one << " " << two << endl; exit(1); }
		mark = -1;
	}

	return mark;
}












