/*
 *	Author : Shujia Huang
 *	Date   : 2012-08-24
 *
 *	Modify : 2012-08-31	Add the error correction to correcte the fosmid alelles which contain error.
 *	This program uses for evaluating the consistency of SNP alleles of fosmid with others.
 *	Just need one input file and whill create one output file.
 *	Input file format (The file should be sorted by the fisrt SNP position ) : (We just use column 1~4) 
 *	e.g : chr6    28527374,28537251,28538474      W,R,Y   A:1,A:1,C:1     chr6-28519579-28538577:28519579-28538577
 	or  : chr6    28527374,28537251,28538474      W,R,Y   A,A,C           chr6-28519579-28538577:28519579-28538577 

 *	Output file format:
 *
 */
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <assert.h>
#include <vector>
#include <list>
#include <map>

#include "utility.h"

using namespace std;

typedef struct {
	int consistency;
	int inconsistency;
} OvlpSite;

typedef struct {

    string refId;
    string fragmentId;
    list < int     > position;
    list < string  > genotype;
    list < string  > allele;
	map  < int, int> posOvlpDep;
	map  < int, double > aleErrWgt; 	   // Allele error weight
	map  < int, vector<OvlpSite> > overlap;// overlap[2].push_back( OvlpSite )
										   // overlap[2][0].consistency;
										   // overlap[2][0].inconsistency;
} Fragment;
void Usage ( const char *prog ) {
	cerr << "Usage: " << prog << " -i [fosmid region file] -o [outputFile] -s <Reference ID, [All]> -n <Overlap number, default[1]> <-r [this parameter means the program will correct the error alleles ]> " << endl;
	exit(1);
}
void output ( Fragment & frag ) ;
void output ( Fragment & frag, ofstream & O );
void Calculate ( list< Fragment > & fragments, int ovlpNumThreshold, bool recorrect, /*bool mask, */ ofstream & O );
void Statistic ( vector< Fragment > & fragmentArray, bool recorrect, /* bool mask, */ ofstream & O );
OvlpSite RecordAgree ( Fragment & fragment1, Fragment & fragment2 );
bool IsOverlap ( list< int > & position, vector< Fragment > & fragments, int ovlpNumThreshold );
bool CheckOverlap ( list < int > &position1, list < int > &position2, int ovlpNumThreshold );

void GetAleRelate ( vector< Fragment > & fragments, map< string, int > & aleCombi, map< string, int > & posCombiNum, int combiDistance );
void CorrErr ( Fragment & fragment, map< string, int > & aleCombi, map< string, int > & posCombiNum, int combiDistance );
void GetErr ( Fragment & fragment1, Fragment & fragment2, double weight, bool isSameHap );

map< char, map<char, string> > ALLELETABLE; // A Globle variant!

int main ( int argc, char* argv[] ) {

	// init the globle variant
    ALLELETABLE['S']['C'] = "G"; ALLELETABLE['S']['G'] = "C";
    ALLELETABLE['M']['A'] = "C"; ALLELETABLE['M']['C'] = "A";
    ALLELETABLE['K']['G'] = "T"; ALLELETABLE['K']['T'] = "G";
    ALLELETABLE['R']['A'] = "G"; ALLELETABLE['R']['G'] = "A";
    ALLELETABLE['W']['A'] = "T"; ALLELETABLE['W']['T'] = "A";
    ALLELETABLE['Y']['C'] = "T"; ALLELETABLE['Y']['T'] = "C";

	char c;
	int  ovlpNum(1);
	bool recorrect( false ); // Determine to recorrect or not
	//bool mask     ( false ); // The alleles will be mask, whether you choose error correction or not, if you choose to mask the incorrect alleles.
	string infile, outfile, refID( "All" );
	while ( ( c = getopt( argc, argv, "i:o:s:n:mrh" ) ) != -1 ) {

		switch ( c ) {
			case 'i' : infile    = optarg; 		   break;
			case 'o' : outfile   = optarg;         break;
			case 's' : refID     = optarg;         break;
			case 'r' : recorrect = true;		   break;
			//case 'm' : mask      = true;		   break;
			case 'n' : ovlpNum   = atoi( optarg ); break;
			case 'h' : Usage ( argv[0] );
			default  : Usage ( argv[0] );
		}
	}
	if ( infile.empty() || outfile.empty() ) Usage( argv[0] );
	cerr << "[Command Parameters]\n" << argv[0] << " -i " << infile << " -o " << outfile << " -n " << ovlpNum << "\trecorrect: " << recorrect /* << "\t" << mask */ << " -s " << refID     << endl;

	ofstream O;
    ifstream I ( infile.c_str() );
	if ( !I ) {
		cerr << "No such file or diretory: " << infile << endl;
		exit(1);
	}

	Fragment fragment;
	list  < Fragment > fragments;
    vector< string   > alleleTmp;

	string tmp;
    int regionEnd;

	while ( 1 ) {
		//scaffold10036   18950,18988,19121       R,M,Y   G,A,C           @A81EVJABXX_4_1105_17480_162292:1-100#scaff
        I >> fragment.refId; if ( I.eof() ) break;
        I >> tmp;                     split ( ",", tmp, fragment.position ); // "8209527,8209528,8209531"
        I >> tmp; tmp = toupper(tmp); split ( ",", tmp, fragment.genotype ); // "Y,S,W"
        I >> tmp; tmp = toupper(tmp); split ( ",", tmp, fragment.allele   ); // "C:3,C:1,A:1" or "G,A,C"
        I >> fragment.fragmentId;
        getline ( I, tmp, '\n' );
		if ( refID != "All" && fragment.refId != refID ) continue;
        if ( fragment.position.size() < 2 ) continue;

        if ( fragment.fragmentId == "-"   ) fragment.fragmentId.clear();
        if ( fragment.allele.front().find( ":" ) != string::npos ) // Ignore the weight!
            for ( list< string >::iterator it( fragment.allele.begin() ); it != fragment.allele.end(); ++it ) {
                split( ":", *it, alleleTmp ); *it = alleleTmp[0];
            }
		for ( list< int >::iterator it( fragment.position.begin() ); it != fragment.position.end(); ++it ) { // Just use for correct alleles or mask alleles
			fragment.aleErrWgt[*it]  = 0.0;
			fragment.posOvlpDep[*it] = 0;
		}

        if ( !fragments.empty() ) {
			if ( !O.is_open() ) { O.open( outfile.c_str()); }
            if ( fragment.refId == fragments.back().refId && regionEnd >= fragment.position.front() ){

                if ( fragment.position.front() < fragments.back().position.front() ) {
                    cerr << "[ERROR] Your file should been sorted." << endl; 
					cerr << fragment.refId << "\t" << join (",", fragment.position ) << "\t" 
						 << join (",", fragment.genotype ) << "\t"
						 << join (",", fragment.allele   ) << "\t" << fragment.fragmentId << endl;
					exit(1);
                }
                if ( regionEnd < fragment.position.back() ) regionEnd = fragment.position.back();
            } else {

                Calculate( fragments, ovlpNum, recorrect, /*mask,*/ O );
                assert( fragments.empty() );
                regionEnd = fragment.position.back();
            }
        } else {
            regionEnd = fragment.position.back();
        }
        fragments.push_back( fragment );
	}
	I.close();

	if ( !fragments.empty() ) { 
		if ( !O.is_open() ) { O.open ( outfile.c_str()); } 
		Calculate( fragments, ovlpNum, recorrect, /*mask*/ O ); 
	}

	if ( O.is_open() ) O.close();

	return 0;
}

void Calculate ( list< Fragment > & fragments, int ovlpNumThreshold, bool recorrect, /*bool mask, */  ofstream & O ) {

	vector< Fragment > fragmentArray;
	while ( !fragments.empty() ) {

		fragmentArray.clear();
		fragmentArray.push_back( fragments.front() ); fragments.pop_front();

		for(list< Fragment >::iterator itFrag( fragments.begin() ); itFrag != fragments.end(); ) {

			if ( IsOverlap( itFrag->position, fragmentArray, ovlpNumThreshold ) ) {
				fragmentArray.push_back( *itFrag );
				itFrag = fragments.erase( itFrag );
			} else {
				++itFrag;
			}
		}
		Statistic ( fragmentArray, recorrect, /* mask, */ O ); // fragmentArray is one of the overlap block
	}
}

void Statistic ( vector< Fragment > & fragmentArray, bool recorrect, /* bool mask, */ ofstream & O ) {

	int THRESHOLD = 3;
	if ( fragmentArray.empty() ) return;

	map< string, int > aleCombi;
	map< string, int > posCombiNum;

	GetAleRelate ( fragmentArray, aleCombi, posCombiNum, THRESHOLD );// For alleles correcting or alleles masking

	OvlpSite ovlp;
	int ovlpNum;
	double weight;
	bool isSameHap;
	for ( size_t i(0); i < fragmentArray.size(); ++i ) {

		for ( size_t j(i+1); j < fragmentArray.size(); ++j ) {

			ovlp    = RecordAgree ( fragmentArray[i], fragmentArray[j] );
			ovlpNum = ovlp.consistency + ovlp.inconsistency;
			if ( ovlpNum == 0 ) continue;

			fragmentArray[i].overlap[ovlpNum].push_back( ovlp );
			//fragmentArray[j].overlap[ovlpNum].push_back( ovlp );
			fragmentArray[j].overlap[ovlpNum].push_back( ovlp );
			isSameHap = ( ovlp.consistency > ovlp.inconsistency ) ? true : false;
			weight    = ( isSameHap ) ? double (ovlp.consistency) / ovlpNum : double ( ovlp.inconsistency ) / ovlpNum;
			GetErr ( fragmentArray[i], fragmentArray[j], weight, isSameHap );
		}

		if ( /* !mask && */ recorrect ) { /* cerr << "Re-CorrERR Aleles" << endl; */ CorrErr ( fragmentArray[i], aleCombi, posCombiNum, THRESHOLD ); }
		output ( fragmentArray[i], O );
	}
	fragmentArray.clear();
}

void GetAleRelate ( vector< Fragment > & fragment, map< string, int > & aleCombi, map< string, int > & posCombiNum, int combiDistance ) {

	string posCombiTyp, aleCombiTyp;
	for ( size_t i(0); i < fragment.size(); ++i ) {

		list< int    >::iterator itPos ( fragment[i].position.begin() );
		list< string >::iterator itGty ( fragment[i].genotype.begin() );
		list< string >::iterator itAle ( fragment[i].allele.begin()   );
		for ( ; itPos != fragment[i].position.end(); ++itPos, ++itGty, ++itAle ) {

			if ( (*itGty)[0] == '+' || (*itGty)[0] == '-'  ) continue; // Just for SNP

			list< int    >::iterator itPosNext ( itPos ); ++itPosNext;
			list< string >::iterator itGtyNext ( itGty ); ++itGtyNext;
			list< string >::iterator itAleNext ( itAle ); ++itAleNext;
			for ( int n (0); itPosNext != fragment[i].position.end() && n < combiDistance; ++itPosNext, ++itGtyNext, ++itAleNext, ++n ) {
				if ( (*itGtyNext)[0] == '+' || (*itGtyNext)[0] == '-'  ) continue; // Just for SNP
				posCombiTyp = fragment[i].refId + ":" + itoa( *itPos ) + ":" + itoa( *itPosNext );
				aleCombiTyp = *itAle + ":" + *itAleNext; 

				//++aleCombi[ posCombiTyp ][ aleCombiTyp ];
				++aleCombi[ posCombiTyp + "|" + aleCombiTyp ];
				++posCombiNum[ posCombiTyp ];

				aleCombiTyp = ALLELETABLE[ (*itGty)[0] ][ (*itAle)[0] ] + ":" + ALLELETABLE[ (*itGtyNext)[0] ][ (*itAleNext)[0] ];
				//++aleCombi[ posCombiTyp ][ aleCombiTyp ];
				++aleCombi[ posCombiTyp + "|" + aleCombiTyp ];
			}
		}
	}

	return;
}

void CorrErr ( Fragment & fragment, map< string, int > & aleCombi, map< string, int > & posCombiNum, int combiDistance ) {

	const double RATIO = 0.5;
	cerr << "Start correcting. " << local_time();

	list< int    >::iterator itPos ( fragment.position.begin() );
	list< string >::iterator itGty ( fragment.genotype.begin() );
	list< string >::iterator itAle ( fragment.allele.begin()   );

	string posCombiTyp, aleCombiTyp, oldAleCombiTyp;
	double ratio, r;
	int un, total;
	for ( ; itPos != fragment.position.end(); ++itPos, ++itGty, ++itAle ) {

		if ( (*itGty)[0] == '+' || (*itGty)[0] == '-' ) continue; // Just for SNPs
		list< int    >::iterator itPosNext ( itPos ); ++itPosNext;
		list< string >::iterator itGtyNext ( itGty ); ++itGtyNext;
		list< string >::iterator itAleNext ( itAle ); ++itAleNext;

		un = 0; total = 0;
		for ( int n(0); itPosNext != fragment.position.end() && n < combiDistance; ++itPosNext, ++itGtyNext, ++itAleNext, ++n ) {

			if ( (*itGtyNext)[0] == '+' || (*itGtyNext)[0] == '-' ) continue; // Just for SNPs
			posCombiTyp = fragment.refId + ":" + itoa( *itPos ) + ":" + itoa( *itPosNext );
            aleCombiTyp = *itAle + ":" + *itAleNext;

			if ( !posCombiNum.count( posCombiTyp ) ) continue;
			assert ( aleCombi.count( posCombiTyp + "|" + aleCombiTyp ) );
			++total;

			ratio = double ( aleCombi[posCombiTyp + "|" + aleCombiTyp] ) / double( posCombiNum[posCombiTyp] );
			if ( ratio <= RATIO ) {
				++un;
			}
		}
		if ( total == 0 ) continue;
		if ( !ALLELETABLE.count( (*itGty)[0] ) || !ALLELETABLE[ (*itGty)[0] ].count( (*itAle)[0] ) ) {
			cerr << "[ERROR]ERROR Genotype. or contain Indel! " << *itPos << "\t" << *itGty << "\t" << *itAle << "\t" << fragment.fragmentId << endl;
			exit(1);
		}
		r = double ( un ) / total;

		if ( r > 0.5 ) {
			
			*itAle = ALLELETABLE[ (*itGty)[0] ][ (*itAle)[0] ];
			// Check if it's OK! 
			itPosNext = itPos; ++itPosNext;
            itGtyNext = itGty; ++itGtyNext;
            itAleNext = itAle; ++itAleNext;
			un = 0; total = 0;
			for ( int n(0); itPosNext != fragment.position.end() && n < combiDistance; ++itPosNext, ++itGtyNext, ++itAleNext, ++n ) {
				if ( (*itGtyNext)[0] == '+' || (*itGtyNext)[0] == '-' ) continue; // Just for SNPs
				posCombiTyp    = fragment.refId + ":" + itoa( *itPos ) + ":" + itoa( *itPosNext );
				aleCombiTyp    = *itAle + ":" + *itAleNext;
				oldAleCombiTyp = ALLELETABLE[(*itGty)[0]][(*itAle)[0]] + ":" + *itAleNext;

				if ( !posCombiNum.count( posCombiTyp ) ) continue;
				assert ( aleCombi.count( posCombiTyp + "|" + oldAleCombiTyp ) );

				//if ( !aleCombi[posCombiTyp].count(aleCombiTyp) ) { 

					--aleCombi[posCombiTyp + "|" + oldAleCombiTyp]; if ( aleCombi[posCombiTyp + "|" + oldAleCombiTyp] < 0 ) aleCombi[posCombiTyp + "|" + oldAleCombiTyp] = 0;
					--posCombiNum[ posCombiTyp ];

				//	++aleCombi[ posCombiTyp ][ aleCombiTyp ];
				//	++posCombiNum[ posCombiTyp ];

				//	++aleCombi[ posCombiTyp ][ ALLELETABLE[(*itGty)[0]][(*itAle)[0]]+":"+ALLELETABLE[(*itGtyNext)[0]][(*itAleNext)[0]] ];
				//	break; // End the loop
				//}
				++aleCombi[ posCombiTyp + "|" + aleCombiTyp ];
				++posCombiNum[ posCombiTyp ];
				++aleCombi[ posCombiTyp + "|" + ALLELETABLE[(*itGty)[0]][(*itAle)[0]]+":"+ALLELETABLE[(*itGtyNext)[0]][(*itAleNext)[0]] ];
			}

			itPosNext = itPos; ++itPosNext;
            itGtyNext = itGty; ++itGtyNext;
            itAleNext = itAle; ++itAleNext;
			for ( int n(0); itPosNext != fragment.position.end() && n < combiDistance; ++itPosNext, ++itGtyNext, ++itAleNext, ++n ) {
				if ( (*itGtyNext)[0] == '+' || (*itGtyNext)[0] == '-' ) continue; // Just for SNPs
				posCombiTyp    = fragment.refId + ":" + itoa( *itPos ) + ":" + itoa( *itPosNext );
                aleCombiTyp    = *itAle + ":" + *itAleNext;
                oldAleCombiTyp = ALLELETABLE[(*itGty)[0]][(*itAle)[0]] + ":" + *itAleNext;

                if ( !posCombiNum.count( posCombiTyp ) ) continue;
				assert ( aleCombi.count( posCombiTyp + "|" + aleCombiTyp ) );

				ratio = double ( aleCombi[posCombiTyp + "|" + aleCombiTyp] ) / double( posCombiNum[posCombiTyp] );
				//cerr << "[RATIO] " << ratio << endl;
				if ( ratio <= RATIO ) {
					++un;
				}
			}
			if ( un > 0 ) {
				cerr << "[Mask] " << fragment.refId << "\t" << *itPos << "\t" << *itGty << "\t" 
       				 << *itAle << "\t" << fragment.fragmentId << endl;
				itPos = fragment.position.erase ( itPos ); --itPos; // Rock back
       			itGty = fragment.genotype.erase ( itGty ); --itGty; // Rock back
            	itAle = fragment.allele.erase ( itAle )  ; --itAle; // Rock back
				fragment.posOvlpDep.erase ( *itPos );
            	fragment.aleErrWgt.erase  ( *itPos );
			} else {
				
				cerr << "[ErrCorr] " << fragment.refId << "\t" << *itPos << "\t" << *itGty << "\t"
					 << ALLELETABLE[ (*itGty)[0] ][ (*itAle)[0] ] << "=>" << *itAle << "\t" << fragment.fragmentId << endl;
			}
		} else if ( r == 0.5 ) { // Candidate allele , which need to be correcte. But I do nothing now!
			cerr << "Ratio: " << r << " , but I'm doing nothing here." << posCombiTyp << "\t" << aleCombiTyp << "\t" << posCombiNum[posCombiTyp] << "\t" << aleCombi[posCombiTyp + "|" + aleCombiTyp] << endl;
		}
	}

	cerr << "Done correcting. " << local_time();
}

void GetErr ( Fragment & fragment1, Fragment & fragment2, double weight, bool isSameHap ) {
// fragment1 need to be correcte! fragment2 is just a control.

	list< int >::iterator    itPos1( fragment1.position.begin() ), itPos2( fragment2.position.begin() );
	list< string >::iterator itAle1( fragment1.allele.begin()   ), itAle2( fragment2.allele.begin()   );
	for ( ; itPos1 != fragment1.position.end(); ++itPos1, ++itAle1 ) {
		assert( itAle1 != fragment1.allele.end() );
		while( (itPos2 != fragment2.position.end()) && ( *itPos1 > *itPos2 ) ) { ++itPos2; ++itAle2; }
		if ( itPos2 == fragment2.position.end() ) break;

		if ( *itPos1 == *itPos2 ) {

			++fragment1.posOvlpDep[*itPos1];
			++fragment2.posOvlpDep[*itPos2];

			if ( isSameHap ) {
				if ( *itAle1 != *itAle2 ) {
					fragment1.aleErrWgt[*itPos1] += weight;
				}
			} else {
				if ( *itAle1 == *itAle2 ) {
					fragment1.aleErrWgt[*itPos1] += weight;
				}
			}
		}
	}
}

OvlpSite RecordAgree ( Fragment & fragment1, Fragment & fragment2 ) {

	OvlpSite agree; agree.consistency = 0; agree.inconsistency = 0;
	list< int    >::iterator itPos1( fragment1.position.begin() ), itPos2( fragment2.position.begin() );
	list< string >::iterator itAle1( fragment1.allele.begin()   ), itAle2( fragment2.allele.begin()   );
	for (; itPos1 != fragment1.position.end(); ++itPos1, ++itAle1 ) {

		assert( itAle1 != fragment1.allele.end() );
		while( (itPos2 != fragment2.position.end()) && ( *itPos1 > *itPos2 ) ) { ++itPos2; ++itAle2; }
        if ( itPos2 == fragment2.position.end() ) break;

        if ( *itPos1 == *itPos2 ) {
            if ( *itAle1 == *itAle2 ) {
				++agree.consistency;
			} else {
				++agree.inconsistency;
			}
        }
    }

	if ( agree.inconsistency != 0 && agree.consistency != 0 ) {

		double ratio = double ( agree.consistency ) / ( agree.consistency + agree.inconsistency );
		if ( ratio < 0.5 ) ratio = 1 - ratio;
/*
		cerr << agree.consistency << " / (" << agree.consistency << "+" << agree.inconsistency << ") " << ratio << endl;
		cerr << fragment1.refId << "\t" << join (",", fragment1.position ) << "\t" << join (",", fragment1.genotype ) << "\t" 
			 << join (",", fragment1.allele ) << "\t" << fragment1.fragmentId << endl;
		cerr << fragment2.refId << "\t" << join (",", fragment2.position ) << "\t" << join (",", fragment2.genotype ) << "\t" 
			 << join (",", fragment2.allele ) << "\t" << fragment2.fragmentId << endl;
		cerr << "********************************************************" << endl; 
*/
	}

    return agree;
}

bool IsOverlap ( list< int > & position, vector< Fragment > & fragments, int ovlpNumThreshold ) {

    bool isOverlap( false );
    for ( vector< Fragment >::iterator itFrag( fragments.begin() ); itFrag != fragments.end(); ++itFrag ) {
        isOverlap = CheckOverlap( itFrag->position, position, ovlpNumThreshold );
        if ( isOverlap ) break;
    }
    return isOverlap;
}

bool CheckOverlap ( list < int > &position1, list < int > &position2, int ovlpNumThreshold ) {

    if ( (position1.size() < ovlpNumThreshold) || (position2.size() < ovlpNumThreshold) ) return false;
    bool overlap ( false );
    int overlapNum (0);
    list< int >::iterator itPos1( position1.begin() ), itPos2( position2.begin() );
    for (; itPos1 != position1.end(); ++itPos1 ) {

        while ( (itPos2 != position2.end()) && (*itPos1 > *itPos2) ) {
            ++itPos2;
        }
        if ( itPos2 == position2.end() ) break;

        if ( *itPos1 == *itPos2 ) {
            ++itPos2;
            ++overlapNum;
            if ( overlapNum >= ovlpNumThreshold ) { overlap = true; break; }
        }
    }
    return overlap;
}

void output ( Fragment & frag ) {

    cout << "--: "<< frag.refId << "\t" << join ( ",", frag.position ) << "\t" << join ( ",", frag.genotype ) << "\t"
      << join ( ",", frag.allele ) << "\t" << frag.fragmentId << endl;
}


void output ( Fragment & frag, ofstream & O ) {

	if ( frag.overlap.empty() ) return;
	O << frag.refId << "\t" << join ( ",", frag.position ) << "\t" << join ( ",", frag.genotype ) << "\t"
	  << join ( ",", frag.allele ) << "\t" << frag.fragmentId << endl;
}























