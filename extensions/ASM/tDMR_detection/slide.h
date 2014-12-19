/*
 *    Filename:   slide.h
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

#ifndef SLIDE_H_YC
#define SLIDE_H_YC

#include <unistd.h>
#include <getopt.h>
#include "cytosine.h"
#include "fisher.h"


typedef struct {
	char patt[PATT_LEN];
	int win, slide;
	int len, nearby;
	int pthreads;
	int depth;
	double pvalue;
} slide_opt_t;

void sliding_core(chr_cytosine_t *, const slide_opt_t *);
void slide_main(const char *, const char *, const slide_opt_t *);
int tdmr_slide(int , char *[]);

#endif
