/*
 *	Author : Shujia Huang
 *	Date   : 2012-07-30
 *
 *	Modify : 2012-08-21 Fix bugs, but no bugs! So fix nothing.
 *	Modify : 2012-08-15	Improve
 *
 *	Using a direct connecting method to build the inital haplotype. 
 *	Just need one input file and will create two output file.
 *	Input  file format: (We just use column 1~4)
 *	e.g. : scaffold10036   18898,18950     R,R     G,A @A81EVJABXX_4_1106_9636_61942:1-100#scaffold10036-18882-18981
	or   : scaffold10036   18950,18988,19121       R,M,Y   G:1,A:1,C:1     @A81X_4_05_10_162292:1-100#scaffold10036-18882-18981

 *	Output file format:
 *	1) "inital haplotype": scaffold10036   18898,18950  R,R  G,A @A81_4_1106_9636_642:1-100#scaffold10036-18882-18981
 *	2) "Ambiguity region": scaffold10036   18950	19121
 *	
 */
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <assert.h>
#include <vector>
#include <list>
#include <map>
#include <set>
#include "utility.h"

using namespace std;

typedef struct {

    string refId;
    string fragmentId;
    list< int     > position;
    list< string  > genotype;
    list< string  > allele;

} Fragment;

void Usage( const char *prog ) {

    cerr << "Usage: " << prog << " -i [variation file] -o [outfile] -n <Overlap number, default[1]> > ambiguityRegion " << endl;
    exit(1);
}
void Insert ( Fragment & haplotype, Fragment & fragment );
void output ( Fragment & hapFrag, ofstream & O );
bool IsOverlap    ( list < int > & position,  vector< Fragment > & fragments, int ovlpNumThreshold ); 
bool CheckOverlap ( list < int > & position1, list < int > &position2,        int ovlpNumThreshold );
void ClusterAndOutput  ( list  < Fragment > & fragments, int ovlpNumThreshold, ofstream & O );
void CreateHaplotype   ( vector< Fragment > & fragmentArray, ofstream & O );
double CalcuAgreeRatio ( Fragment & fragment1, Fragment & fragment2 );
Fragment CreateOppositeHaplotype( Fragment & haplotype );

void output ( Fragment & fragment );

int main ( int argc, char* argv[] ) {

	char c;
	int ovlpNum(1);
	string infile, outfile;
	while ( ( c = getopt( argc, argv, "i:n:o:h" ) ) != -1 ) {
		
		switch ( c ) {
			case 'i' : infile  = optarg;         break;
			case 'o' : outfile = optarg;         break;
			case 'n' : ovlpNum = atoi( optarg ); break;
			case 'h' : Usage ( argv[0] );
			default  : Usage ( argv[0] );
		}
	}
	if ( infile.empty() || outfile.empty() ) Usage( argv[0] );
	cerr << "[Command Parameters]\n" << argv[0] << " -i " << infile << " -o " << outfile << " -n " << ovlpNum << endl;

	ofstream O ( outfile.c_str());
	ifstream I ( infile.c_str() );
	if ( ! I ) {
		cerr << "[ERROR] Cannot open " << infile << " : No such file or directory" << endl;
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
		I >> tmp;	                  split ( ",", tmp, fragment.position ); // "8209527,8209528,8209531"
		I >> tmp; tmp = toupper(tmp); split ( ",", tmp, fragment.genotype ); // "Y,S,W"
        I >> tmp; tmp = toupper(tmp); split ( ",", tmp, fragment.allele   ); // "C:3,C:1,A:1" or "G,A,C"
        I >> fragment.fragmentId;
		getline ( I, tmp, '\n' );
		if ( fragment.position.size() < 2 ) continue;
		if ( fragment.fragmentId == "-"   ) fragment.fragmentId.clear();
		if ( fragment.allele.front().find( ":" ) != string::npos ) // Ignore the weight! 
			for ( list< string >::iterator it( fragment.allele.begin() ); it != fragment.allele.end(); ++it ) {
				split( ":", *it, alleleTmp ); *it = alleleTmp[0];
			}

		if ( !fragments.empty() ) {
			if ( fragment.refId == fragments.back().refId && regionEnd >= fragment.position.front() ){

				if ( fragment.position.front() < fragments.back().position.front() ) {
					cerr << "[ERROR] Your file should been sorted.\n" << fragment.fragmentId << "\n" 
						 << fragments.back().fragmentId << endl; exit(1);
				}
				if ( regionEnd < fragment.position.back() ) regionEnd = fragment.position.back();
			} else {
				
				ClusterAndOutput( fragments, ovlpNum, O );
				assert( fragments.empty() );
				regionEnd = fragment.position.back();
			}
		} else {
			regionEnd = fragment.position.back();
		}
		fragments.push_back( fragment );
	}
	I.close();

	if ( !fragments.empty() ) ClusterAndOutput( fragments, ovlpNum, O );

	O.close();
}

void ClusterAndOutput ( list< Fragment > & fragments, int ovlpNumThreshold, ofstream & O ) {

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

//[for debug]for ( size_t i(0); i < fragmentArray.size(); ++i ) { output( fragmentArray[i] );} exit(1);
		// Cluster succeeful, Now build haplotype
		CreateHaplotype( fragmentArray, O );
	}
}

void CreateHaplotype ( vector< Fragment > & fragmentArray, ofstream & O ) {

	const double AGREETHRESHOLD = 0.79;

	/* inital haplotype1 and haplotype2 */
	Fragment haplotype1, haplotype2;
	haplotype1 = fragmentArray[0];
	haplotype2 = CreateOppositeHaplotype( haplotype1 );

	// Attention: The position of haplotype1 and haplotype2 will aways be the same.
	double agreeRatio1, agreeRatio2;
	int flag;
	Fragment h1, h2;
	for ( size_t i(1); i < fragmentArray.size(); ++i ) {

		agreeRatio1 = CalcuAgreeRatio( haplotype1, fragmentArray[i] );
		agreeRatio2 = CalcuAgreeRatio( haplotype2, fragmentArray[i] );

		if ( agreeRatio1 > AGREETHRESHOLD ) { 

			flag = 2;
			h2   = CreateOppositeHaplotype( fragmentArray[i] );
			Insert( haplotype1, fragmentArray[i]  );
		} else if ( agreeRatio2 > AGREETHRESHOLD ) {

			flag = 1;
			h1   = CreateOppositeHaplotype( fragmentArray[i] );
			Insert( haplotype2, fragmentArray[i] );
		} else { // Ambiguity region

			flag = 0;
		}

		switch ( flag ) {
			case 1 : Insert( haplotype1, h1 ); break;
			case 2 : Insert( haplotype2, h2 ); break;
			default: // Output ambiguity region
			/*	cout << fragmentArray[i].refId << "\t"   << fragmentArray[i].position.front() << "\t" 
					 << fragmentArray[i].position.back() << endl;
			*/
			// output( fragmentArray[i] )
				output( fragmentArray[i], O );
		}
	}
	output ( haplotype1, O );
	output ( haplotype2, O );

	return;
}

void Insert( Fragment & haplotype, Fragment & fragment ) {

	assert( haplotype.refId == fragment.refId );
	list< int    >::iterator itPosHap( haplotype.position.begin() ), itPosFrag( fragment.position.begin() );
	list< string >::iterator itGtyHap( haplotype.genotype.begin() ), itGtyFrag( fragment.genotype.begin() );
	list< string >::iterator itAleHap( haplotype.allele.begin()   ), itAleFrag( fragment.allele.begin()   );
	bool isEnd(false), overlap( false );
	for (; itPosFrag != fragment.position.end(); ++itPosFrag, ++itGtyFrag, ++itAleFrag ) {
		
		while ( itPosHap != haplotype.position.end() && *itPosHap < *itPosFrag ) { ++itPosHap; ++itGtyHap; ++itAleHap; }
		if ( itPosHap == haplotype.position.end() ) isEnd = true;
		if ( isEnd ) break;

		if ( *itPosHap > *itPosFrag ) {
			haplotype.position.insert( itPosHap, *itPosFrag );
			haplotype.genotype.insert( itGtyHap, *itGtyFrag );
			haplotype.allele.insert  ( itAleHap, *itAleFrag );
		} else { // *itPosHap == *itPosFrag 
			overlap = true;
			assert ( *itGtyHap == *itGtyFrag );
			if ( *itAleHap == "-" ) *itAleHap = *itAleFrag; // Uses the better one or Just the first one.
			++itPosHap; ++itGtyHap; ++itAleHap;
		}
	}

	if ( isEnd ) {

		assert( overlap );
		for (; itPosFrag != fragment.position.end(); ++itPosFrag, ++itGtyFrag, ++itAleFrag ) {

			assert( itGtyFrag != fragment.genotype.end() && itAleFrag != fragment.allele.end() );
			haplotype.position.push_back( *itPosFrag );
			haplotype.genotype.push_back( *itGtyFrag );
			haplotype.allele.push_back  ( *itAleFrag );
		}
	}

	if ( haplotype.fragmentId.empty() ) {
		haplotype.fragmentId = fragment.fragmentId;
	} else { 
		if ( !fragment.fragmentId.empty() ) haplotype.fragmentId += ( "," + fragment.fragmentId );
	}
}

Fragment CreateOppositeHaplotype( Fragment & haplotype ) {

	map< char, map<char, char> > alleleTable;
    alleleTable['S']['C'] = 'G'; alleleTable['S']['G'] = 'C';
    alleleTable['M']['A'] = 'C'; alleleTable['M']['C'] = 'A';
    alleleTable['K']['G'] = 'T'; alleleTable['K']['T'] = 'G';
    alleleTable['R']['A'] = 'G'; alleleTable['R']['G'] = 'A';
    alleleTable['W']['A'] = 'T'; alleleTable['W']['T'] = 'A';
    alleleTable['Y']['C'] = 'T'; alleleTable['Y']['T'] = 'C';

	Fragment oppHaplotype;
	oppHaplotype.refId      = haplotype.refId;
    oppHaplotype.position   = haplotype.position;
    oppHaplotype.genotype   = haplotype.genotype;
    vector< string > indel;
    list< string >::iterator itAle( haplotype.allele.begin() ), itGty( haplotype.genotype.begin() );
    for ( ; itAle != haplotype.allele.end(); ++itAle, ++itGty ) {
        assert( itGty != haplotype.genotype.end() );
        if ( (*itGty)[0] == '+' || (*itGty)[0] == '-' ) {
        // Indel
            split( "|", *itGty, indel );
            if ( indel.size() == 2 ) {
                if ( *itAle == indel[0] ) {
                    oppHaplotype.allele.push_back( indel[1] );
                } else if( *itAle == indel[1] ) {
                    oppHaplotype.allele.push_back( indel[0] );
                } else {
                    cerr << "[ERROR]Genotype error. " << *itGty << "\t" << *itAle << endl;
                }
            } else {
                if ( *itAle == "-" ) {
                    oppHaplotype.allele.push_back( indel[0] );
                } else if( *itAle == indel[0] ) {
                    oppHaplotype.allele.push_back( "-" );
                } else {
                    cerr << "[ERROR]Genotype error. " << *itGty << "\t" << *itAle << endl;
                }
            }
        } else {
        // SNP
			if ( !alleleTable.count( (*itGty)[0] ) || !alleleTable[(*itGty)[0]].count( (*itAle)[0] ) ) {
				cerr << "[ERROR]Your Genotype is unmatch!!"<< endl;
				cerr << (*itGty)[0] << "\n" << (*itAle)[0] << "\t" << haplotype.fragmentId << "\n" << endl;
			}
            assert (alleleTable.count( (*itGty)[0] ) && alleleTable[(*itGty)[0]].count( (*itAle)[0] ));
            oppHaplotype.allele.push_back( char2str(alleleTable[(*itGty)[0]][(*itAle)[0]]) );
        }
    }
	oppHaplotype.fragmentId.clear();

	return oppHaplotype;
}

double CalcuAgreeRatio ( Fragment & fragment1, Fragment & fragment2 ) {

	double agree(0), totalOvlp(0);

	list< int    >::iterator itPos1( fragment1.position.begin() ), itPos2( fragment2.position.begin() );
	list< string >::iterator itAle1( fragment1.allele.begin()   ), itAle2( fragment2.allele.begin()   );
	for (; itPos1 != fragment1.position.end(); ++itPos1, ++itAle1 ) {
		
		assert( itAle1 != fragment1.allele.end() );
		while( (itPos2 != fragment2.position.end()) && ( *itPos1 > *itPos2 ) ) { ++itPos2; ++itAle2; }
		if ( itPos2 == fragment2.position.end() ) break;

		if ( *itPos1 == *itPos2 ) {
			++totalOvlp;
			if ( *itAle1 == *itAle2 ) ++agree;
		}
	}
	double ratio = ( totalOvlp > 0 ) ? agree / totalOvlp : 0;
	return ratio;
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

void output ( Fragment & hapFrag, ofstream & O ) {

	if ( hapFrag.fragmentId.empty() ) return; // Ignore this type at today, but no forever. ( 2012-08-17 )
	//if ( hapFrag.fragmentId.empty() ) hapFrag.fragmentId = "-"; // This may be used next day.
	O << hapFrag.refId                  << "\t" << join ( ",", hapFrag.position ) << "\t"
	  << join ( ",", hapFrag.genotype ) << "\t" << join ( ",", hapFrag.allele   ) << "\t"
	  << hapFrag.fragmentId             << endl;
}

void output ( Fragment & fragment ) {

 cout << fragment.refId                  << "\t" << join ( ",", fragment.position ) << "\t"
      << join ( ",", fragment.genotype ) << "\t" << join ( ",", fragment.allele   ) << "\t"
      << fragment.fragmentId             << endl;

/*
	cerr << fragment.refId               << "\t" << join ( ",", fragment.position ) << "\t"
      << join ( ",", fragment.genotype ) << "\t" << join ( ",", fragment.allele   ) << endl;
*/
}













