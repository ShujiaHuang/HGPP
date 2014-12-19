/*
 * =====================================================================================
 *
 *     Filename:  main.c 
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "slide.h"
#include "extend.h"

int main(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Fuck the tDMR detection\n\n");
		fprintf(stderr, "Usage:         tdmr <prog> [opt] <infile>\n");
		fprintf(stderr, "Prog:          slide    using sliding window to calculate the significant regions\n");
		fprintf(stderr, "               extend   extend the significant regions\n");
		/*
		fprintf(stderr, "               all      run sliding windows then extend significant regions\n");
		//*/
		fprintf(stderr, "\n");
		exit(0);
	}
	if (strcmp(argv[1], "slide") == 0) tdmr_slide(argc-1, argv+1);
	else if (strcmp(argv[1], "extend") == 0) tdmr_extend(argc-1, argv+1);
	/*
	else if (strcmp(argv[1], "all") == 0) tmdr_all(argc-1, argv+1);
	//*/
	else {
		fprintf(stderr, "Unrecognized program identifier %s\n\n", argv[1]);
	}
	return 1;
}
