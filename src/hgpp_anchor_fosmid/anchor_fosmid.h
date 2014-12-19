/*
 *  Author : Shujia Huang
 *  Date   : 2012-04-06
 *  Anchor_fosmid header
 *
 */
#ifndef _ANCHOR_FOSMID_H
#define _ANCHOR_FOSMID_H

#include <vector>
#include "sam.h"
#include "utility.h"

class SamExt : public Sam {

public:
	SamExt () : _read_start(0), _read_end(0) {};
	~SamExt() { _read_id.clear(); _pool_id.clear(); Sam::clear(); }

	int read_start() { return _read_start; }
	int read_end  () { return _read_end;   }
	string pool_id() { return _pool_id; }
	string read_id() { return _read_id; }

	void clear () { _read_id.clear(); _pool_id.clear(); Sam::clear(); }
	void assign( SamLine *B, string pooling_id );

	int cov_length_after ( int position );

private:

	int _read_start;
	int _read_end;
	string _pool_id;
	string _read_id;
};

void SamExt::assign( SamLine *B, string pooling_id ) {

	Sam::assign( B );
	_read_id.assign( _RID, 0, _RID.find_first_of( "#" ) );
	_pool_id = pooling_id;

	_read_start = 1;
	_read_end   = _read_start + _seq.length() - 1;

	return;
}

int SamExt::cov_length_after ( int position ) {

	int length(0);
	if ( position >= ref_end() ) { // have been contained! Don't need to calculate again!
		length = 0;
	} else if ( position <= ref_start() ){
		length = map_length();
	} else {

//cerr << "Start:     " << ref_start() << "\t" << ref_end() << endl;
//cerr << "Cigar:     " << cigar_seq2cigar() << endl;
//cerr << "Cigar_seq: " << cigar_seq() << endl;
//cerr << "Cov After: " << position << endl; 
		int pos ( ref_start() );
		for ( size_t i(0); i < cigar_seq().length(); ++i ) {
			if ( cigar_seq()[i] == 'M' && pos > position ) ++length;
			if ( cigar_seq()[i] == 'M' || cigar_seq()[i] == 'D' || cigar_seq()[i] == 'N' ) ++pos;
		}

//cerr << "Cov length: " << length << endl;
//cerr << "End POS   : " << pos    << endl;
//exit(1);
	}
	
	return length;
}

#endif

