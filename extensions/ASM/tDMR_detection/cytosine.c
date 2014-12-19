/*
 *    Filename:   cytosine.c
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
   12-May-2010, 0.0.4

   - Stat. mC from binormal distribution, and also need depth cutoff

    3-May-2010, 0.0.3

   - Must all C, G repeat number <= 1.5

   12-Apr-2010, 0.0.2

   - Change the repeat number <= 1.5

 */

#include "cytosine.h"

chr_cytosine_t *load_cytosines(const char *cfn1, const char *cfn2, const char *patt, const int mC_depth)
{

	fprintf(stderr, "[tDMR] Loading Cytosines [%s] ...\n", patt);
	char line[1024] = "";
	gzFile fp1 = gzopen(cfn1, "r");
	gzFile fp2 = gzopen(cfn2, "r");
	if (fp1 == NULL || fp2 == NULL) {
		fprintf(stderr, "[tDMR] Can't open Files ...\n");
		exit(EXIT_FAILURE);
	}

	chr_cytosine_t *chr_cytosine = (chr_cytosine_t *) calloc(1, sizeof(chr_cytosine_t));
	int n, allocated;
	n = allocated = 0;
	cytosine_t *p = chr_cytosine->c_list;

	char chr[256] = "";
	unsigned int pos;
	char strain;
	char qpatt[3], seq[3];
	float rep;
	int uniq;
	depth_t methy, unmethy, bino;
	int cpg = patt[1] == 'G' ? 1 : 0;
	const int UNIQ = cpg ? 4 : 2;

	while (1) {
		if (gzgets(fp1, line, 1024) == 0) break;
		if (n >= allocated) {
			allocated += PER_ALLOC;
			chr_cytosine->c_list = (cytosine_t *) realloc(chr_cytosine->c_list, sizeof(cytosine_t) * allocated);
			p = chr_cytosine->c_list + n;
		}
		sscanf(line, "%s\t%u\t%c\t%s\t%s\t%f\t%hu\t%hu\t%hu\n", chr, &pos, &strain, qpatt, seq, &rep, &methy, &unmethy, &bino);
		uniq = (rep > REPEAT) ?  0 : 1;
		p->pos = pos;
		if (methy >= bino && methy + unmethy >= mC_depth) p->mCa = 1;
		else p->mCa = 0;
		if (strcmp(qpatt, patt) == 0) {
			if (!cpg) {
				p->a = methy; p->b = unmethy;
				if (gzgets(fp2, line, 1024) == 0) break;
				sscanf(line, "%s\t%u\t%c\t%s\t%s\t%f\t%hu\t%hu\t%hu\n", chr, &pos, &strain, qpatt, seq, &rep, &methy, &unmethy, &bino);
				uniq = (rep > REPEAT) ? (uniq - 1) : (uniq + 1);
				p->c = methy; p->d = unmethy;
				if (methy >= bino && methy + unmethy >= mC_depth) p->mCc = 1;
				else p->mCc = 0;
			} else if (cpg) {
				p->a = methy; p->b = unmethy;
				if (gzgets(fp1, line, 1024) == 0) break;
				sscanf(line, "%s\t%u\t%c\t%s\t%s\t%f\t%hu\t%hu\t%hu\n", chr, &pos, &strain, qpatt, seq, &rep, &methy, &unmethy, &bino);
				uniq = (rep > REPEAT) ? (uniq - 1) : (uniq + 1);
				p->a += methy; p->b += unmethy;
				if (methy >= bino && methy + unmethy >= mC_depth) p->mCa++;
				if (gzgets(fp2, line, 1024) == 0) break;
				sscanf(line, "%s\t%u\t%c\t%s\t%s\t%f\t%hu\t%hu\t%hu\n", chr, &pos, &strain, qpatt, seq, &rep, &methy, &unmethy, &bino);
				uniq = (rep > REPEAT) ? (uniq - 1) : (uniq + 1);
				p->c = methy; p->d = unmethy;
				if (methy >= bino && methy + unmethy >= mC_depth) p->mCc = 1;
				else p->mCc = 0;
				if (gzgets(fp2, line, 1024) == 0) break;
				sscanf(line, "%s\t%u\t%c\t%s\t%s\t%f\t%hu\t%hu\t%hu\n", chr, &pos, &strain, qpatt, seq, &rep, &methy, &unmethy, &bino);
				p->c += methy; p->d += unmethy;
				uniq = (rep > REPEAT) ? (uniq - 1) : (uniq + 1);
				if (methy >= bino && methy + unmethy >= mC_depth) p->mCc++;
			}
			if (uniq == UNIQ) {
				p->mrate1 = (p->a + p->b) ? (float) (p->a) / (float) (p->a + p->b) : -1;
				p->mrate2 = (p->c + p->d) ? (float) (p->c) / (float) (p->c + p->d) : -1;
				p->flag = (p->mrate1 > p->mrate2) ? 1 : ((p->mrate1 < p->mrate2) ? -1 : 0);
				p++; n++;
			}
		} else {
			if (gzgets(fp2, line, 1024) == 0) break;
		}
	}
	memcpy(chr_cytosine->chr, chr, CHRN_LEN);
	chr_cytosine->n = n;
	gzclose(fp1); gzclose(fp2);
	return chr_cytosine;
}
