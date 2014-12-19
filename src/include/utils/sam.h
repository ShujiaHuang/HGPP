/*
 *	Author : Shujia Huang
 *	Date   : 2012-04-06
 *	Sam/bam format header
 *
 */
#ifndef _SAM_H
#define _SAM_H

#include <iostream>
#include "utility.h"

using namespace std;

class SamLine {
// HUMuvjD320poolingDFAAPEI-10_Index8_Index2_C958  16   scaffold27   24677427   37  116M1D67M    *    0	0	...
public:
	string _RID;
	int    _Flag;
	string _ref_id;
	int    _position;
	int    _mapQ;
	string _cigar;
	string _XorD;
	int    _coor;
	int    _insert_size;
	string _seq;

	//string Qseq;
	SamLine (): _Flag(0), _position(0), _mapQ(0), _coor(0), _insert_size(0) {}
	~SamLine(){ clear(); }
	void clear() {
		_RID.clear();	_ref_id.clear();
		_cigar.clear();	_XorD.clear();
		_Flag = 0;		_position = 0;
		_mapQ = 0;		_coor	  = 0;
		_insert_size = 0; 
		_seq.clear();
	}
};

class Sam : public SamLine {

public:

	Sam (): _ref_start(0), _ref_end(0), _N_length(0), _map_length(0), _deletion_length(0) {}
	~Sam(){ clear();      }
	Sam   ( SamLine *B );
	void assign( SamLine *B );
	void clear() {
		SamLine::clear();
		_ref_start       = 0;
		_ref_end         = 0;
		_N_length		 = 0;
		_map_length 	 = 0;
		_deletion_length = 0;
		_cigar_seq.clear();
	}

	int n_length       () { return _N_length;        }
	int map_length     () { return _map_length;      }
	int deletion_length() { return _deletion_length; }
	int ref_start () { return _ref_start;  }
	int ref_end   () { return _ref_end;    }

	void para_cigar();
	string cigar_seq2cigar();
	string& cigar_seq() { return _cigar_seq; }
	void update_cigar( const string& cigar_seq  ) { _cigar_seq = cigar_seq; _cigar = cigar_seq2cigar(); para_cigar(); }
	// Will update cigar, cigar_seq and ref_end position at the same time
	inline void update_ref_start( const int pos ) { _ref_start = pos; _position = _ref_start; }
	inline void update_ref_end  ( const int pos ) { _ref_end   = pos; }

private:

	int _ref_start;
	int _ref_end;
	int _N_length; // Contain the sequence length, which skipped region from the reference, not just N
	int _map_length;
	int _deletion_length;

	string _cigar_seq;
};

#endif


