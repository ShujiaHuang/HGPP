/*
 *    Filename:   extend.c
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

/*
   18-May-2010, 0.0.12

   - HYPO defined as 0.1 not 0.2, for detecting the tDMR in promotor region 

   12-May-2010, 0.0.11

   - Calculate the mC density

   13-Apr-2010, 0.0.10

   - Calculate the skewness and kurtosis

   6-Apr-2010, 0.0.9

   - Add depth cutoff for initial window

   - Not allowed too many uncovered C added

   - Add bounder test last 5 C should be significant

   3-Apr-2010, 0.0.8

   - Add bootstrap test for tDMR (test 0.05)

   31-Mar-2010, 0.0.7
 
   - Support multiple threads when doing extension
 
   30-Mar-2010, 0.0.6
 
   - Fix a minor bug when calculating the variance of methylation rate
 
   29-Mar-2010, 0.0.5
 
   - Allowed add mC that would change the trends
 
   - Change the output format [chr, beg, end, p-value, methy_rate1, methy_rate2,
     methy_rate1 var, methy_rate2 var, cov1, cov2, cov1 var, cov2 var, # of C,
     # of uncovered C1, # of uncovered C2]
 
   24-Mar-2010, 0.0.4
 
   - Build the extend as a module for tDMR detection
 
   - Change the input file type
 
   - Redefined the structure for Cytosine
 
   21-Mar-2010, 0.0.3
 
   - Getopt from command line
 
   - Add cutoff for nearby CpG distance when extend it into tDMR

*/

#include "extend.h"

#define online(x, i) {								\
	delta[i] = x - mean[i];							\
	mean[i] = mean[i] + delta[i] / n;				\
	M2[i] = M2[i] + delta[i] * (x - mean[i]);		\
}

int stat_tdmr(cytosine_t *ll, cytosine_t *rr, float *var, int *uncov1, int
		*uncov2, float *mC_density1, float *mC_density2, float *skewness, float *kurtosis)
{
	cytosine_t *p = ll;
	int u1, u2;
	u1 = u2 = 0;
	int c1, c2 = 0;
	int mC1, mC2;
	mC1 = mC2 = 0;
	int n = 0;
	float delta[4] = {0,0,0,0};
	float M2[4]    = {0,0,0,0};
	float mean[4]  = {0,0,0,0};
	float avg, moment2, moment3, moment4;
	avg = moment2 = moment3 = moment4 = 0;

	while(p <= rr){
		n++;
		c1 = (p->a + p->b);
		c2 = (p->c + p->d);
		mC1 += p->mCa;
		mC2 += p->mCc;
		if (!c1) u1++;
		if (!c2) u2++;
		online(c1, 0);
		online(c2, 1);
		if (p->mrate1 != -1) online((p->mrate1), 2);
		if (p->mrate2 != -1) online((p->mrate2), 3);
		if (p->mrate1 != -1 && p->mrate2 != -1) avg += fabsf(p->mrate1 - p->mrate2);
		p++;
	}
	avg /= n;
	p = ll;
	float tmp ;
	while(p <= rr) {
		tmp = fabsf(p->mrate1 - p->mrate2) - avg;
		moment2 += tmp * tmp;
		moment3 += tmp * tmp * tmp;
		moment4 += tmp * tmp * tmp * tmp;
		p++;
	}
	moment2 /= n;
	moment3 /= n;
	moment4 /= n;
	tmp = moment2 * moment2;
	*skewness = moment3 / sqrt(tmp * moment2);
	*kurtosis = moment4 / tmp;
	var[0] = M2[0] / (n-1);
	var[1] = M2[1] / (n-1);
	var[2] = M2[2] / (n-1);
	var[3] = M2[3] / (n-1);
	*uncov1 = u1;
	*uncov2 = u2;
	*mC_density1 = mC1 / (float) (n * 2);
	*mC_density2 = mC2 / (float) (n * 2);
	return n;
}

inline float backtrace(cytosine_t *wh, cytosine_t *l, cytosine_t *r, const int n, const int win, const double p) 
{
	if(n == win) return 0;
	cytosine_t *cp = l;
	int i, j;
	int a, b, c, d;
	while (cp != l) {
		
	}
}


inline float bootstrap(const int run, cytosine_t *l, cytosine_t *r, const int n, const double p)
{
	if (!run) return 1;
	cytosine_t *cp;
	int i, j;
	int a, b, c, d;
	int fail = 0;
	for(i = 0; i < BOOTSTRAP; ++i) {
		a = b = c = d = 0;
		for(j = 0; j < n; ++j) {
			cp = l + (rand() % n);
			a += cp->a; b += cp->b;
			c += cp->c; d += cp->d;
		}
		if (fisher2(a, b, c, d) >= p) fail++;
	}
	return (1 - (float)fail/(float)(BOOTSTRAP+1));
}

inline int bounder_test(cytosine_t *l, cytosine_t *r, const double pvalue, const int mdev)
{
	cytosine_t *p = l;
	int a, b, c, d;
	a = b = c = d = 0;
	while (p <= r) {
		a += p->a; b += p->b;
		c += p->c; d += p->d;
		p++;
	}
	float temp1 = (float) a / (float) (a + b);
	float temp2 = (float) c / (float) (c + d);
	if (fisher2(a, b, c, d) >= pvalue || (temp1 - temp2) * mdev < 0)
		return 0;
	else 
		return 1;
}

void extend_detect(const int tid, const chr_cytosine_t *chr_cytosine, const chr_tdmr_t *chr_tdmr, const ext_opt_t *o)
{
	srand(time(NULL));
	fprintf(stderr, "[tDMR] Running extend ...\n");
	int cytosine_n, tdmr_n;
	unsigned int beg, end, pre_beg, pre_end;
	cytosine_t *cytosine, *r, *l, *last;

	ext_opt_t oo;
	oo.win = o->win; oo.len = o->len; oo.pvalue = o->pvalue;
	oo.diff = o->diff; oo.depth = o->depth;
	oo.bounder = o->bounder - 1;
	oo.bdtest = o->bdtest;
	oo.btsp = o->btsp;
	const int len = o->len;
	const int nearby = o->nearby;
	const double pvalue = o->pvalue;

	tdmr_t *tdmr = chr_tdmr->tdmr_list;
	cytosine_n = chr_cytosine->n; tdmr_n = chr_tdmr->n;
	pre_beg = pre_end = 0;
	cytosine = r = chr_cytosine->c_list;
	last = cytosine + cytosine_n;
	const char *chr = chr_cytosine->chr;

	double rp, lp = 0;
	unsigned a, b, c, d = 0;
	unsigned ra, rb, rc, rd = 0;
	unsigned la, lb, lc, ld = 0;
	int mdev = 0;
	int cov1, cov2 = 0;
	int uncov1, uncov2 = 0;
	float var[4];
	float skewness, kurtosis;
	int i;
	cytosine_t *cur = cytosine;
/*
	tdmr_t *tdmr_exted = (tdmr_t *)calloc(tdmr_n, sizeof(tdmr_t));
	tdmr_t *exted = tdmr_exted;
	//*/
	for (i = 0; i < tdmr_n; ++i, ++tdmr) {

#ifdef PTHREADS
		if(o->nthreads > 1) {
			pthread_mutex_lock(&lock);
			tdmr_t *p = chr_tdmr->tdmr_list + i;
			if (p->tid < 0) {
				int j;
				for (j = 0; j < tdmr_n - i&& j < N_PER_THREAD; j+=1){
					(p+j)->tid = tid;
				}
			} else if (p->tid != tid) {
				pthread_mutex_unlock(&lock);
				continue;
			}
			pthread_mutex_unlock(&lock);
		}
#endif

		beg = tdmr->beg; end = tdmr->end;
		double pre_p = tdmr->p;
		mdev = (tdmr->mrate1 > tdmr->mrate2) ? 1 : -1;
		while (cur != last && beg > cur->pos) cur++;
		if (cur == last || beg != cur->pos) continue;
		l = cur;
		a = l->a; b = l->b;
		c = l->c; d = l->d;
		r = l + 1;
		while (r != last && end >= r->pos && (r-1)->pos + nearby >= r->pos) {
			a += r->a; b += r->b;
			c += r->c; d += r->d;
			r++;
		}
		la = ra = a;
		lb = rb = b;
		lc = rc = c;
		ld = rd = d;
		if(end != (r-1)->pos || a + b + c + d < oo.depth) continue;
		else r--;
		lp = rp = pre_p;
		float temp1, temp2;
		float mC_density1, mC_density2;
		mC_density1 = mC_density2 = 0;

		while (1) {
			/* pre version not allowed add Cytosine change the tDMR differetial trends
			   mdev * (l-1)->flag < 0
			   mdev * (r+1)->flag < 0
			   */
			if ((l == cytosine  || ((l-1)->pos + nearby < l->pos) || ((l-1)->pos + len < r->pos))) goto REXT;
			if ((r == last || (r->pos + nearby < (r+1)->pos) || (l->pos + len < (r+1)->pos))) goto LEXT;

			la = a + (l-1)->a; lb = b + (l-1)->b; lc = c + (l-1)->c; ld = d + (l-1)->d;
			temp1 = (float) la / (float) (la + lb);
			temp2 = (float) lc / (float) (lc + ld);
			if ((oo.bdtest && ! bounder_test(l, l + oo.bounder, pvalue, mdev)) ||
					(temp1 < temp2 && (mdev > 0 || temp1 * DIFF > temp2)) ||
					(temp2 < temp1 && (mdev < 0  ||  temp2 * DIFF > temp1))) {
				ra = a; rb = b; rc = c; rd = d;
				goto REXT;
			}
			ra = a + (r+1)->a; rb = b + (r+1)->b; rc = c + (r+1)->c; rd = d + (r+1)->d;
			temp1 = (float) ra / (float) (ra + rb);
			temp2 = (float) rc / (float) (rc + rd);
			if ((oo.bdtest && ! bounder_test(r - oo.bounder, r, pvalue, mdev) )
					|| (temp1 < temp2 && (mdev > 0 || temp1 * DIFF >= temp2)) ||
					(temp2 < temp1 && (mdev < 0 ||  temp2 * DIFF >= temp1))) {
				la = a; lb = b; lc = c; ld = d;
				goto LEXT;
			} 
			lp = fisher2(la, lb, lc, ld);
			rp = fisher2(ra, rb, rc, rd);
			if (rp >= pvalue && lp >= pvalue) {
				break;
			} else if (rp >= pvalue) {
				la = a; lb = b; lc = c; ld = d;
				pre_p = lp;
LEXT:
				while(l != cytosine && ((l-1)->pos + nearby >= l->pos) && (l-1)->pos + len >= r->pos) {
					la = a + (l-1)->a; lb = b + (l-1)->b; lc = c + (l-1)->c; ld = d + (l-1)->d;
					temp1 = (float) la / (float) (la + lb);
					temp2 = (float) lc / (float) (lc + ld);
					if ((oo.bdtest && ! bounder_test(l, l + oo.bounder, pvalue, mdev))
							|| (temp1 < temp2 && (mdev > 0 || temp1 * DIFF > temp2)) 
							|| (temp2 < temp1 && (mdev < 0 || temp2 * DIFF > temp1)) 
							|| ((lp = fisher2(la, lb, lc, ld)) >= pvalue)) {
						break;
					} else {
						--l;
						a = la; b = lb; c = lc; d = ld;
						pre_p = lp;
					}
				}
				break;
			} else if (lp >= pvalue) {
				ra = a; rb = b; rc = c; rd = d;
				pre_p = rp;
REXT:
				
				while(r != last && (r->pos + nearby >= (r+1)->pos) && l->pos + len >= (r+1)->pos) {
					ra = a + (r+1)->a; rb = b + (r+1)->b; rc = c + (r+1)->c; rd = d + (r+1)->d;
					temp1 = (float) ra / (float) (ra + rb);
					temp2 = (float) rc / (float) (rc + rd);
					if ((oo.bdtest && ! bounder_test(r - oo.bounder, r, pvalue, mdev)) || 
							(temp1 < temp2 && (mdev > 0 || temp1 * DIFF > temp2)) || 
							(temp2 < temp1 && (mdev < 0 || temp2 * DIFF > temp1)) || 
							((rp = fisher2(ra, rb, rc, rd)) >= pvalue)) {
						break;
					} else {
						++r;
						a = ra; b = rb; c = rc; d = rd;
						pre_p = rp;
					}
				}
				break;
			} else {
				if (rp < lp) {
					++r;
					pre_p = rp;
					a = ra; b = rb; c = rc; d = rd;
				} else {
					--l;
					pre_p = lp;
					a = la; b = lb; c = lc; d = ld;
				}
			}
		}
		int n = stat_tdmr(l, r, var, &uncov1, &uncov2, &mC_density1, &mC_density2, &skewness, &kurtosis);
		cov1 = a + b;
		cov2 = c + d;
		float bt = bootstrap(oo.btsp, l, r, n, 0.05);
		temp1 = (float) a / (float) cov1;
		temp2 = (float) c / (float) cov2;
		if ((temp1 < temp2 && temp2 >= HYPO && temp1 * oo.diff <= temp2) ||
				(temp1 > temp2 && temp1 >= HYPO && temp2 * oo.diff <= temp1))
			fprintf(stdout,
					"%s\t%u\t%u\t%.7E\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\t%.2f\t%.2f\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
					chr, l->pos, r->pos, pre_p, temp1,
					temp2, var[2], var[3], cov1, cov2, var[0],
					var[1], n, uncov1, uncov2, mC_density1, mC_density2, bt, skewness, kurtosis);
	}
}

#ifdef PTHREADS
typedef struct {
	int tid;
	chr_cytosine_t *chr_cytosine;
	chr_tdmr_t *chr_tdmr;
	ext_opt_t *opt;
} thread_aux_t;

static void *workers(void *thread_aux)
{
	thread_aux_t *aux = (thread_aux_t *)thread_aux;
	extend_detect(aux->tid, aux->chr_cytosine, aux->chr_tdmr, aux->opt);
}
#endif

void extend_main(const char *cfn1, const char *cfn2, const char *tdmr_file, ext_opt_t * const o)
{

	fprintf(stderr, "Fucking the extention for tDMR\n");
	chr_cytosine_t *chr_cytosine = load_cytosines(cfn1, cfn2, o->patt, o->mC_depth);


	FILE *tdmr_fp = fopen(tdmr_file, "r");
	chr_tdmr_t *chr_tdmr = (chr_tdmr_t *) calloc(1, sizeof(chr_tdmr_t));
	memcpy(chr_tdmr, chr_cytosine->chr, 256);
    tdmr_t *tdmr = chr_tdmr->tdmr_list;
	char chr[CHRN_LEN];

	int n, allocated;
	n = allocated = 0;
	fprintf(stderr, "[tDMR] loading candidate region ...\n");
	while(!feof(tdmr_fp)){
		if (n >= allocated) {
			allocated += 0xFFFF;
			chr_tdmr->tdmr_list = (tdmr_t *)realloc(chr_tdmr->tdmr_list, sizeof(tdmr_t) * allocated);
			tdmr = chr_tdmr->tdmr_list + n;
		}
		fscanf(tdmr_fp, "%s\t%u\t%u\t%lf\t%f\t%f\t%hu\t%hu\n",
				chr, &(tdmr->beg), &(tdmr->end), &(tdmr->p), &(tdmr->mrate1), &(tdmr->mrate2), &(tdmr->cov1), &(tdmr->cov2));
		tdmr->tid = -1;
		if (tdmr->p >= PVALUE) continue;
		if (tdmr->mrate2 < HYPO && tdmr->mrate1 < HYPO) continue ;
		if ((tdmr->mrate1 < tdmr->mrate2 && tdmr->mrate1 * DIFF <= tdmr->mrate2) || (tdmr->mrate1 > tdmr->mrate2 &&  tdmr->mrate2 * DIFF <= tdmr->mrate1)) {
			tdmr++; n++;
		}
	}
	fclose(tdmr_fp);
	chr_tdmr->n = n;

#ifndef PTHREADS
	extend_detect(0, chr_cytosine, chr_tdmr, o);
#else 
	if (o->nthreads > 1) {
		pthread_t *tid;
		pthread_attr_t attr;
		thread_aux_t *thread_aux; 
		int j;
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
		thread_aux = (thread_aux_t *)calloc(o->nthreads, sizeof(thread_aux_t));
		tid = (pthread_t*)calloc(o->nthreads, sizeof(pthread_t));
		for (j = 0; j < o->nthreads; ++j) {
			thread_aux[j].tid = j;
			thread_aux[j].chr_cytosine = chr_cytosine;
			thread_aux[j].chr_tdmr = chr_tdmr;
			thread_aux[j].opt = o;
			pthread_create(&tid[j], &attr, workers, thread_aux + j);
		}
		pthread_attr_destroy(&attr);
		for (j = 0; j < o->nthreads; ++j) pthread_join(tid[j], 0);
		free(thread_aux); free(tid);
	} else {
		extend_detect(0, chr_cytosine, chr_tdmr, o);
	}
#endif
	free(chr_cytosine->c_list); free(chr_cytosine);
	free(chr_tdmr->tdmr_list); free(chr_tdmr);
}

ext_opt_t *ini_ext_opt(void)
{
	ext_opt_t *opt = (ext_opt_t *) calloc(1, sizeof(ext_opt_t));
	strcpy(opt->patt, "CG");
	opt->win = 5;
	opt->pvalue = 0.05;
	opt->len = 10000;
	opt->nearby = 200;
	opt->diff = 2;
	opt->depth = 10;
	opt->mC_depth = 4;
	opt->bounder = 5;
	opt->bdtest = 1;
	opt->btsp = 0;
	opt->nthreads = 1;
	return opt;
}

int tdmr_extend(int argc, char *argv[])
{
	ext_opt_t *opt = ini_ext_opt();

	int c;
	while((c = getopt(argc, argv, "c:w:l:n:f:p:d:m:t:be")) != -1) {
		switch(c) {
			case 'c':
				snprintf(opt->patt, PATT_LEN, "%s", optarg);
				break;
			case 'w':
				opt->win = atoi(optarg);
				break;
			case 'l':
				opt->len = atoi(optarg);
				break;
			case 'n':
				opt->nearby = atoi(optarg);
				break;
			case 'f':
				opt->diff = atof(optarg);
				break;
			case 'p':
				opt->pvalue = atof(optarg);
				break;
			case 't':
				opt->nthreads = atoi(optarg);
				break;
			case 'd':
				opt->depth = atoi(optarg);
				break;
			case 'm':
				opt->mC_depth = atoi(optarg);
				break;
			case 'b':
				opt->btsp = 1;
				break;
			case 'e':
				opt->bdtest = 0;
				break;
			default: 
				fprintf(stderr, "Unrecognized option %c\n", c);
				return 1;
		}
	}

	if (optind + 3 > argc) {
		fprintf(stderr, " \n");
		fprintf(stderr, "Usage:        tdmr_ext [OPT] <in_cytosine1> <in_cytosine2> <in_tdmr> \n\n");
		fprintf(stderr, "Options:      -c  STR    cytosine pattern (CG/CHH/CHG) [%s]\n", opt->patt);
		fprintf(stderr, "              -l  INT    maximum tDMR length [%d]\n", opt->len);
		fprintf(stderr, "              -n  INT    maximum distance between nearby CG in one tDMR [%d]\n", opt->nearby);
		fprintf(stderr, "              -p  FLOAT  maximum p-value for significant tDMR [%.1le]\n", opt->pvalue);
		fprintf(stderr, "              -f  FLOAT  minimum differential folder [%d]\n", opt->diff);
		fprintf(stderr, "              -d  INT    minimum depth [%d]\n", opt->depth);
		fprintf(stderr, "              -m  INT    minimum depth for mC [%d]\n", opt->mC_depth);
		fprintf(stderr, "              -b         bootstrap for each significant tDMR (burn-in %d times) [NULL]\n", BOOTSTRAP);
		fprintf(stderr, "              -e         bounder test for last %d cytosines before a cytosine is added [DO]\n", opt->bounder);
		fprintf(stderr, "              -t  INT    # of multiple threads [%d]\n", opt->nthreads);
		fprintf(stderr, " \n");
		return 1;
	}
	extend_main(argv[optind], argv[optind+1], argv[optind+2], opt);
	free(opt);
	return 1;
}
