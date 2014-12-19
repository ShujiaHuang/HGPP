/*
 *  Author : Shujia Huang
 *  Date   : 2012-04-06
 *  Sam/bam format header
 *
 */

#include "sam.h"

Sam::Sam( SamLine *B ) {

	assign( B );
}

void Sam::assign( SamLine *B ) {

	_RID     = B->_RID ;    _Flag     = B->_Flag;
    _ref_id  = B->_ref_id;  _position = B->_position;
    _mapQ    = B->_mapQ;    _cigar    = B->_cigar;
    _XorD    = B->_XorD;    _coor     = B->_coor;
    _insert_size = B->_insert_size;
	_seq     = B->_seq;

    para_cigar();
}

void Sam::para_cigar() {

	_N_length        = 0;
	_map_length      = 0;
	_deletion_length = 0;
	_cigar_seq.clear();
	if ( _cigar.empty() ) return;

	int seq_length (0);
	string num; 
	for ( size_t i(0); i < _cigar.length(); ++i ) {

		switch ( _cigar[i] ) {
			case '0' : 
					// 0 must not be the header element of num
					if ( num.empty() ) {
						cerr << "[ERROR]The first Number is 0!! cigar = " << _cigar << endl;
						exit(1);
					}
			case '1' :
			case '2' :
			case '3' : 
			case '4' :
			case '5' :
			case '6' :
			case '7' :
			case '8' :
			case '9' : num.append( 1, _cigar[i] ); break;
			case 'H' : num.clear(); break;
			case 'M' :
				if ( !num.empty() ) {
					seq_length  += atoi( num.c_str() );
					_map_length += atoi( num.c_str() );
					_cigar_seq.append( atoi( num.c_str() ), 'M' );
				} else {
					cerr << "[ERROR]The Number is empty at M. !! cigar = " << _cigar << endl;
					exit(1);
				}
				num.clear();
				break;
			case 'D' :
				if ( !num.empty() ) {
					 _deletion_length += atoi( num.c_str() );
					 _cigar_seq.append( atoi( num.c_str() ), 'D' );
				} else {
					cerr << "[ERROR]The Number is empty at D. !! cigar = " << _cigar << endl;
                    exit(1);
				}
				num.clear();
				break;
			case 'I' : 
				if ( !num.empty() ) {
					seq_length  += atoi( num.c_str() );
                     _cigar_seq.append( atoi( num.c_str() ), 'I' );
                } else {
                    cerr << "[ERROR]The Number is empty at I. !! cigar = " << _cigar << endl;
                    exit(1);
                }
                num.clear();
				break;
			case 'S' : 
				if ( !num.empty() ) {
					seq_length  += atoi( num.c_str() );
                     _cigar_seq.append( atoi( num.c_str() ), 'S' );
                } else {
                    cerr << "[ERROR]The Number is empty at S. !! cigar = " << _cigar << endl;
                    exit(1);
                }
                num.clear();
				break;
			case 'N' :
				if ( !num.empty() ) {
					 _N_length += atoi( num.c_str() );
                     _cigar_seq.append( atoi( num.c_str() ), 'N' );
                } else {
                    cerr << "[ERROR]The Number is empty at N. !! cigar = " << _cigar << endl;
                    exit(1);
                }
                num.clear();
				break;

			default:
				cerr << "[BUG] Undefined to deal with this cigar mark: " << _cigar[i] << "\t" << num << "\nCigar: " << _cigar << endl;
				cerr << _RID      << "\t" << _Flag << "\t" << _ref_id      << "\t" 
					 << _position << "\t" << _mapQ << "\t" << _cigar	   << "\t" 
					 << _XorD     << "\t" << _coor << "\t" << _insert_size << "\t" << _seq << endl;
				num.clear();
		}
	}
	if ( _seq.length() != seq_length ) {
		cerr << "[ERROR] Sequence length Unmatch. Cigar_seq ERROR or input file error or forget to update the seq first, when programing.\n" 
			 << _RID << "\t" << _seq.length() << "\t" << seq_length << "\t" << _cigar << endl;
		exit(1);
	}
	_ref_start = _position;
	_ref_end   = _position + _map_length + _deletion_length + _N_length - 1;
//cerr << "*Merge: ref_end: " << this->ref_end() + 1 << " : " << _position << "\t" << _map_length << "\t" << _deletion_length << "\t" << _N_length<< endl;
	if ( _ref_end - _ref_start < 0 ) _ref_end = _ref_start;
}

string Sam::cigar_seq2cigar() {

	string cigar;
	if ( _cigar_seq.empty() ) return "*";

//cerr << "_cigar_seq : " << _cigar_seq << endl;
	int num(0);
	char type;
	for ( size_t i(0); i < _cigar_seq.length(); ++i ){
		
		if ( num > 0 ) {

			if ( type == _cigar_seq[i] ) {
				++num;
			} else {
				cigar += itoa( num ) + char2str( type );
//cerr << "*cigar : " << cigar << endl;
				type   = _cigar_seq[i];
				num    = 1;
			}
		} else {
			++num;
			type = _cigar_seq[i];
		}
	}
	cigar += itoa( num ) + char2str( type );
//cerr << "*cigar : " << cigar << endl;
//cerr << "Test cigar_seq2cigar: \ncigar_seq: " << _cigar_seq << "\ncigar    : " << cigar << endl; exit(1);
	return cigar;
}

