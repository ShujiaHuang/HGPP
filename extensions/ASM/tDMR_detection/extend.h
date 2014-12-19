/*
 *    Filename:   extend.h
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

#ifndef _EXTEND_H_YC
#define _EXTEND_H_YC

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include "cytosine.h"
#include "fisher.h"

#define PVALUE 0.05
#define DIFF 2
#define HYPO 0.1

#ifndef BOOTSTRAP
#define BOOTSTRAP 10000
#endif

#ifdef PTHREADS
#include <pthread.h>
#define N_PER_THREAD 0xF00
static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
#endif

typedef struct {
	int win, len, nearby, depth, mC_depth, bounder;
	int bdtest, btsp;
	int nthreads;
	float diff;
	double pvalue;
	char patt[PATT_LEN];
} ext_opt_t;


int stat_tdmr(cytosine_t *, cytosine_t *, float *, int *, int *, float *, float *, float *, float *);
void extend_detect(const int tid, const chr_cytosine_t *, const chr_tdmr_t *, const ext_opt_t *);
void extend_main(const char *, const char *, const char *, ext_opt_t * const);
int tdmr_extend(int , char *[]);

#endif
