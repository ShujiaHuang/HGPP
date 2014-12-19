/*
 *    Filename:   cytosine.h
 *
 *    Copyright (c) 2010, Chang Yu (yuchang@genomics.org.cn)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef CYTOSINE_H_YC
#define CYTOSINE_H_YC

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

#define PER_ALLOC   0xFFFF
#define CHRN_LEN    256
#define PATT_LEN    4 
#define SEQ_LEN     3
#define REPEAT      1.5

typedef unsigned short depth_t;
//typedef void gzFile;

typedef struct {
	unsigned int pos;
	int ctype;
	/* effective detph for Tissue1 (methy, unmehty) Tissue2 (methy, unmethy)
	   for CpG those should be the total detph of both strains
	 */
	depth_t a, b, c, d;
	/* binormal distribution predict detpth for methylation if one cytosine would be called as mC  */
	int mCa, mCc;
	int flag;
	float mrate1, mrate2;
} cytosine_t;

typedef struct {
	char chr[CHRN_LEN];
	int n;
	cytosine_t *c_list;
} chr_cytosine_t;

typedef struct {
	int tid;
	unsigned int beg, end;
	depth_t cov1, cov2;
	float mrate1, mrate2;
	double p;
	float var[4];
} tdmr_t;

typedef struct {
	char chr[CHRN_LEN];
	int n;
	tdmr_t *tdmr_list;
} chr_tdmr_t;

chr_cytosine_t *load_cytosines(const char *, const char *, const char *, const int );

#endif    /*    End of cytonsine.h    */
