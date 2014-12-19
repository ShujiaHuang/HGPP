/*
 *  Author : Shujia Huang
 *  Date   : 2012-04-06
 *  Anchor_fosmid header
 *
 */
#ifndef _ANCHOR_VARIATION_H
#define _ANCHOR_VARIATION_H

#include <vector>
#include "sam.h"
#include "utility.h"

typedef struct {
	char genotype;
    vector< char > allele;
    vector< int  > read_position;
	vector< string > reads;
} SNP_Point;

typedef struct {
	string genotype;
	vector< string > allele;
	vector< string > reads;
} Indel;

class SamExt : public Sam {

public:
	SamExt () : _read_start(0), _read_end(0) {};
	~SamExt() { _read_id.clear(); _pool_id.clear(); Sam::clear(); }

	int read_start() { return _read_start; }
	int read_end  () { return _read_end;   }
	string pool_id() { return _pool_id; }
	string read_id() { return _read_id; }

	void clear () { _read_id.clear(); _pool_id.clear(); Sam::clear(); }
	void assign( SamLine *B );

	//int cov_length_after ( int position );

private:

	int _read_start;
	int _read_end;
	string _pool_id;
	string _read_id;
};

void SamExt::assign( SamLine *B ) {

	Sam::assign( B );
	// query id format: HUMuvjD320poolingDFAAPEI-10_Index20_Index10_scaffold23_100
    // pool id:         HUMuvjD320poolingDFAAPEI-10_Index20_Index10
	// A81C7KABXX:3:63:11357:65868#AACTGC
	vector< string > tmp;
	_pool_id.assign( _RID, 0, _RID.find_first_of( '#' ) );
	split( ":", _pool_id, tmp );
	_pool_id = join( "_", tmp );
	_read_id = "@" + _pool_id;

	_read_start = 1;
	_read_end   = _read_start + _seq.length() - 1;

	return;
}

#endif

