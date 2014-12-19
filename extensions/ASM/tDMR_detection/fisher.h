/* fisher exact test */
#ifndef FISHER_H
#define FISHER_H

/*
 *        drwan  not drawn
 * white    a       b      | n1_
 * black    c       d      | n2_
 * ------------------------------
 * total   n_1     n_2     | n
 *
 * a is an observation of random variable X, which confrom to hypergeometric distribution
 * hypergeometric distribution (without replacement)
 * p{X=a} = choose(n1_, a) * choose(n - n1_, n_1 - a) / choose(n, n_1)
 *
 * 
 * fisher0 calculate left tail probability: SUM(p(x), min <= x <= a)
 * fisher1 calculate right tail probability: SUM(p(x), a <= x <= max)
 * fisher2 calculate two tail probability: SUM(p(x), min <= x <= max && p(x) <= p(a))
 * fisher  calculate all above
 *
 * For more information, refer to http://en.wikipedia.org/wiki/Fisher's_exact_test
*/

#include <stdlib.h>
#include <math.h>
#include <float.h>

/*
 @abstract         Compute fisher p_value
 @param  a         Count of white balls drawn
 @param  b         Count of white balls remain
 @param  c         Count of black balls drawn
 @param  d         Count of black balls remain
 @param  p_left    Pointer to left-tail p value
 @param  p_right   Pointer to right-tail p value
 @param  p_twotail Pointer to twotail p value
*/
void fisher(int a, int b, int c, int d, double *p_left, double *p_right, double *p_twotail);


/*
 @abstract         Compute left-tail fisher p_value
 @param  a         Count of white balls drawn
 @param  b         Count of white balls remain
 @param  c         Count of black balls drawn
 @param  d         Count of black balls remain
 @return           Left-tail p value
*/
double fisher0(int a, int b, int c, int d);


/*
 @abstract         Compute right-tail fisher p_value
 @param  a         Count of white balls drawn
 @param  b         Count of white balls remain
 @param  c         Count of black balls drawn
 @param  d         Count of black balls remain
 @return           Right-tail p value
*/
double fisher1(int a, int b, int c, int d);


/*
 @abstract         Compute twotail fisher p_value
 @param  a         Count of white balls drawn
 @param  b         Count of white balls remain
 @param  c         Count of black balls drawn
 @param  d         Count of black balls remain
 @return           Right-tail p value
*/
double fisher2(int a, int b, int c, int d);

#endif

