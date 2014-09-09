/*
 * nwdp_v2.h
 *
 *  Created on: Jul 03, 2014
 *      Author: Stefan Seemann, seemann@rth.dk
 */

#ifndef NWDP_H_
#define NWDP_H_

#include <iostream>
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>

#define DEBUG 0

using namespace std;

extern float pnullSgl[];
extern float psubSgl[];
extern float pnullDbl[];
extern float psubDbl[];

extern float maxpnullSgl;
extern float maxpnullDbl;

extern std::map<const char, const int> nucIdx;

/* arguments */
extern float kappa;
extern float tau;
extern float alpha;
extern float beta;

/* sequence alignment */
const int a = 1;   /* Match */
const int b = 0;   /* Mismatch */
const int submatrix[ 4 ][ 4 ] = { { a, b, b, b },    /* substitution matrix */
                                  { b, a, b, b },
                                  { b, b, a, b },
                                  { b, b, b, a } } ;

extern float INFINITE;

typedef struct tagLocalHit {
	float similarity;		/* Similarity */
	int lstart_1;			/* Start position of LS_a */
	int lend_1;				/* End position of LS_a */
	int lstart_2;			/* Start position of LS_b */
	int lend_2;				/* End position of LS_b */
} LocalHit;

extern int readinput( istream & is, string & name, string & seq, vector<float> & prob );
extern void getunpaired( vector<float> & prob, int len, vector<float> & probSgl );
extern void getlogoddsDbl( vector<float> & probDbl, string seq, int len, float pnull );
extern void getlogoddsSgl( vector<float> & probSgl, string seq, int len, float pnull );
extern void reducematrix(vector<float> & prob, int len, int prec );
vector<string> &split(const string &s, char delim, vector<string> &elems);
vector<string> split(const string &s, char delim);

extern float nwdp( string seq_1, vector<float> & probSgl_1, int * idx_1_aln, int len_1, string seq_2, vector<float> & probSgl_2, int * idx_2_aln, int len_2, int & len_aln, float ** F, float ** Q, float ** P, char ** trF, char ** trQ, char ** trP );
extern float simdp( vector<float> & probDbl_1, int len_1, vector<float> & probDbl_2, int len_2, int * idx_1_aln, int * idx_2_aln, int len_aln );
extern void prob_backtracking( int * idx_1_subaln, int len_1, int * idx_2_subaln, int len_2, int & len_subaln, float ** F );
extern float simbp( string seq_1, vector<float> & probSgl_1, int len_1, string seq_2, vector<float> & probSgl_2, int len_2, int * idx_1_aln, int * idx_2_aln, int len_aln );
extern float gappenalty( int * idx_aln, int len_aln, int len );

extern void printalign(string & seq_1, int * idx_1_aln, string & seq_2, int * idx_2_aln, int len_aln );
extern void freeMatrix(float ** matrix, int leny);
extern void usage_dotaligner(char * program);

extern void nwdp_initF( float ** F, int L1, int L2 );
extern void nwdp_initF_affinegaps( float ** F, int L1, int L2, bool local );
extern void nwdp_initGap( float ** Q, int L1, int L2 );
extern float max3( float f1, float f2, float f3, char* ptr );
extern float max( float f1, float f2 );
template <typename T>
extern void print_matrixdp( T ** F, float * prob_1, int len_1, float * prob_2, int len_2 );
extern void affinegapcosts( int * idx_1_aln, int * idx_2_aln, int & len_aln, int & open, int & extended );
extern void nwdp_initTB( char ** traceback, int L1, int L2 );
extern void reverse( int * list, int len );
extern int statistical_sampling( unsigned int *sam_prob, int sam_len );

#endif /* NWDP_H_ */