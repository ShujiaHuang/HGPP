// base.h
// Ver 1.0
// Author : Shujia Huang
// Date   : 2011/04/06

#ifndef _UTILITY_H
#define _UTILITY_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <deque>
#include <list>
#include <algorithm>

using namespace std;
// split function
void split( const char *delim, string strLine, vector<string>& tokens, bool is_add=0 );
void split( const char *delim, string strLine, vector<char>&   tokens, bool is_add=0 ); // New 2011/08/23
void split( const char *delim, string strLine, vector<int>&    tokens, bool is_add=0 );

void split( const char *delim, string strLine, list<string>& tokens, bool is_add=0 );
void split( const char *delim, string strLine, list<char>&   tokens, bool is_add=0 ); // New 2011/08/23
void split( const char *delim, string strLine, list<int>&    tokens, bool is_add=0 );

void join ( const char *delim, vector<string>& goven, string &strLine );
string join ( const char *delim, vector<string>& goven );
string join ( const char *delim, vector<char>& goven ); // New 2011/08/23
string join ( const char *delim, vector<int>& goven );
string join ( const char *delim, list<string> &goven );
string join ( const char *delim, list<char> &goven );
string join ( const char *delim, list<int>  &goven );
string join ( const char *delim, deque<string> &goven );
string join ( const char *delim, deque<char> &goven );
string join ( const char *delim, deque<int> &goven);
string join ( const char *delim, char* goven[], int size ); // New 2013-05-28 17:20:18

string itoa ( const int number );
string toupper ( string strLine );
string char2str( char c ); // New 2011/08/23

// 
char mytolower( const char c );
char mytoupper( const char c );

// open file 
void open_file ( ifstream &in, const string &file_name );
//Calculation locaktime
char* local_time ();

// Math function
template < class T >
double mean( vector<T> &values );
template < class T >
double mean( deque<T> &values );

template < class T >
double variance ( vector<T> &values );
template < class T >
double variance ( deque<T> &values  );

//Trim spaces among string	2011/09/28
string trim_space( string str );

#endif

