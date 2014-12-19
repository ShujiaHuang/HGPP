/*
 * =====================================================================================
 *
 *     Filename:  slide.c 
 *
 *     Copyright (c) 2010, Chang Yu (yuchang@genomics.org.cn)
 *
 *     This program is free software; you can redistribute it and/or
 *     modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2
 *     of the License, or (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program; if not, write to the Free Software
 *     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 *
 * =====================================================================================
 */

#include "slide.h"

char buf[1024] = "";

void sliding_core(chr_cytosine_t *chr_cytosine, const slide_opt_t *o)
{

	fprintf(stderr, "[tDMR] sliding windows ... \n");
	int n = chr_cytosine->n;
	char *chr = chr_cytosine->chr;

	slide_opt_t oo;
	oo.win = o->win;
	oo.slide = o->slide;
	oo.len = o->len;
	oo.nearby = o->nearby;
	oo.depth = o->depth;
	oo.pvalue = o->pvalue;

	int aa, bb, cc, dd;

	int i, nearby;
	double p = 0;
	cytosine_t *beg = chr_cytosine->c_list;
	cytosine_t *cyt;
	for (i = 0; i <= n - oo.win; i += oo.slide, beg += oo.slide) {
		aa = bb = cc = dd = 0;
		cyt = beg;
		while(cyt != beg + oo.win) {
			aa += cyt->a;
			bb += cyt->b;
			cc += cyt->c;
			dd += cyt->d;
			/*
			if ((cyt+1)->pos - cyt->pos + 1 > oo.nearby) {
				nearby = 0;
				break;
			}
			*/
			cyt++;
		}
		p = fisher2(aa, bb, cc, dd);
		cyt--;

		if (p <= oo.pvalue) { 
			fprintf(stdout, "%s\t%u\t%u\t%.7E\t%.2f\t%.2f\t%d\t%d\n", chr, beg->pos, cyt->pos, p, (aa+bb)?((float)aa/(float)(aa+bb)):(float)0, (cc+dd)?((float)cc/(float)(cc+dd)):(float)0, aa+bb, cc+dd);
		}
	}
}

void slide_main(const char *fn1, const char *fn2, const slide_opt_t *o)
{
	fprintf(stderr, "Fucking the tDMR sliding windows\n");
	chr_cytosine_t *chr_cytosine = load_cytosines(fn1, fn2, o->patt, 0);
	sliding_core(chr_cytosine, o);
	free(chr_cytosine->c_list);
	free(chr_cytosine);
}

slide_opt_t *ini_slide_opt(void)
{
	slide_opt_t *opt = (slide_opt_t *)calloc(1, sizeof(slide_opt_t));
	strcpy(opt->patt, "CG");
	opt->win = 5;
	opt->slide = 1;
	opt->len = 1000;
	opt->nearby = 200;
	opt->depth = 20;
	opt->pthreads = 1;
	opt->pvalue = 0.05;
	return opt;
}

int tdmr_slide(int argc, char *argv[])
{
	slide_opt_t *opt = ini_slide_opt();

	int c;
	while((c = getopt(argc, argv, "c:w:s:l:p:n:")) != -1) {
		switch (c) {
			case 'c':
				snprintf(opt->patt, PATT_LEN, "%s", optarg);
				break;
			case 'w': 
				opt->win = atoi(optarg); 
				break;
			case 's': 
				opt->slide = atoi(optarg); 
				break;
			case 'l': 
				opt->len = atoi(optarg); 
				break;
			case 'p': 
				opt->pvalue = atof(optarg); 
				break;
			case 'n':
				opt->nearby = atoi(optarg);
				break;
			default: 
				fprintf(stderr, "Unrecognized option %c\n", c);
				fprintf(stderr, "\n");
				return 1;
		}
	}

	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:        tdmr slide [OPT] <Tissue1 Cytosine> <Tissue2 Cytosine>\n");
		fprintf(stderr, "Options:      -c  STR     cytosine pattern (CG/CHH/CHG) [%s]\n", opt->patt);
		fprintf(stderr, "              -w  INT     # of Cytosine(CG/CHH/CHG) contains in one window [%d]\n", opt->win);
		fprintf(stderr, "              -s  INT     step size (# of C) for sliding windows [%d]\n", opt->slide);
		fprintf(stderr, "              -p  FLOAT   maximum p-value for significant tDMR [%.1le]\n", opt->pvalue);
		fprintf(stderr, "              -l  INT     maximum tDMR bp length [%d]\n", opt->len);
		fprintf(stderr, "              -n  INT     maximum distance between nearby Cytosine conteined in one tDMR [%d]\n", opt->nearby);
#ifdef PTHREADS
		fprintf(stderr, "              -t  INT     # of multi-threads [%d]\n", opt->pthreads);
#endif
		fprintf(stderr, "\n");
		return 1;
	}

	slide_main(argv[optind], argv[optind+1], opt);
	free(opt);
	return 1;
}
