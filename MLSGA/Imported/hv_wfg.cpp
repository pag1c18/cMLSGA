/*
 Copyright (C) 2010 Lyndon While, Lucas Bradstreet, Wesley Cox. 

This program is free software (software libre). You can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software Foundation; 
either version 2 of the License, or (at your option) any later version. 

This program is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE. See the GNU General Public License for more details. 
*/
/*This code has been modified by Przemyslaw Grudniewski for the puposes of the MLSGA-framework 2019*/

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
//#include <resource.h>
#include "hv_wfg.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



namespace HV
{

#define BEATS(x,y)   (x > y) 
#define WORSE(x,y)   (BEATS(y,x) ? (x) : (y)) 

	typedef double OBJECTIVE;

	typedef struct
	{
		OBJECTIVE *objectives;
	} POINT;

	typedef struct
	{
		int nPoints;
		int n;
		POINT *points;
	} FRONT;




	int n;     // the number of objectives 
	POINT ref; // the reference point 

	FRONT *fs;    // memory management stuff 
	int fr = 0;   // current depth 
	int maxm = 0; //number of points in the front
	int maxn = 0;	//number of objectives
	int safe;     // the number of points that don't need sorting
	int** nextp;
	int** prevp;
	int* firstp;
	int* lastp;
	int* psize;

	double totaltime;

	double hv(FRONT);

	int greater(const void *v1, const void *v2);
	int greaterabbrev(const void *v1, const void *v2);
	int dominates2way(POINT p, POINT q, int k);
	bool dominates1way(POINT p, POINT q, int k);
	void makeDominatedBit(FRONT ps, int p);
	double hv2(FRONT ps, int k);
	double inclhv(POINT p);
	double inclhv2(POINT p, POINT q);
	double inclhv3(POINT p, POINT q, POINT r);
	double inclhv4(POINT p, POINT q, POINT r, POINT s);
	double exclhv(FRONT ps, int p);
	double hv(FRONT ps);
	void Refine_PF(std::vector<std::vector<double>> &PF_temp, std::vector<double> &Real_PF_max_f);
}

int HV::greater(const void *v1, const void *v2)
// this sorts points worsening in the last objective
{
	POINT p = *(POINT*)v1;
	POINT q = *(POINT*)v2;
	for (int i = n - 1; i >= 0; i--) {
		if BEATS(p.objectives[i], q.objectives[i]) {
			return -1;
		}
		else if BEATS(q.objectives[i], p.objectives[i]) {
			return  1;
		}
	}
	return 0;
}

int HV::greaterabbrev(const void *v1, const void *v2)
// this sorts points worsening in the penultimate objective
{
	POINT p = *(POINT*)v1;
	POINT q = *(POINT*)v2;
	for (int i = n - 2; i >= 0; i--) {
		if BEATS(p.objectives[i], q.objectives[i]) {
			return -1;
		}
		else if BEATS(q.objectives[i], p.objectives[i]) {
			return  1;
		}
	}
	return 0;
}

int HV::dominates2way(POINT p, POINT q, int k)
// returns -1 if p dominates q, 1 if q dominates p, 2 if p == q, 0 otherwise 
// k is the highest index inspected 
{
	for (int i = k; i >= 0; i--) {
		if BEATS(p.objectives[i], q.objectives[i]) {
			for (int j = i - 1; j >= 0; j--) {
				if BEATS(q.objectives[j], p.objectives[j]) {
					return 0;
				}
			}
			return -1;
		}
		else if BEATS(q.objectives[i], p.objectives[i]) {
			for (int j = i - 1; j >= 0; j--) {
				if BEATS(p.objectives[j], q.objectives[j]) {
					return 0;
				}
			}
			return  1;
		}
	}
	return 2;
}

bool HV::dominates1way(POINT p, POINT q, int k)
// returns true if p dominates q or p == q, false otherwise 
// the assumption is that q doesn't dominate p 
// k is the highest index inspected 
{
	for (int i = k; i >= 0; i--) {
		if BEATS(q.objectives[i], p.objectives[i]) {
			return false;
		}
	}
	return true;
}

void HV::makeDominatedBit(FRONT ps, int p)
// creates the front ps[0 .. p-1] in fs[fr], with each point bounded by ps[p] and dominated points removed 
{
	int l = 0;
	int u = psize[fr] - 1;
	int index = lastp[fr];
	int removed = 0;
	for (int i = psize[fr] - 1; i >= 0; i--) {
		if (BEATS(ps.points[p].objectives[n - 1], ps.points[index].objectives[n - 1])) {
			fs[fr].points[u].objectives[n - 1] = ps.points[index].objectives[n - 1];
			for (int j = 0; j < n - 1; j++) {
				fs[fr].points[u].objectives[j] = WORSE(ps.points[p].objectives[j], ps.points[index].objectives[j]);
			}
			if (dominates1way(ps.points[p], ps.points[index], n - 2)) {
				int left = index - prevp[fr][index];
				int right = index + nextp[fr][index];
				if (left < 0) {
					firstp[fr] = right;
				}
				else {
					nextp[fr][left] = right - left;
				}
				if (right >= maxm) {
					lastp[fr] = left;
				}
				else {
					prevp[fr][right] = right - left;
				}
				removed++;
			}
			u--;
		}
		else {
			fs[fr].points[l].objectives[n - 1] = ps.points[p].objectives[n - 1];
			for (int j = 0; j < n - 1; j++) {
				fs[fr].points[l].objectives[j] = WORSE(ps.points[p].objectives[j], ps.points[index].objectives[j]);
			}
			l++;
		}
		index -= prevp[fr][index];
	}
	POINT t;
	// points below l are all equal in the last objective; points above l are all worse 
	// points below l can dominate each other, and we don't need to compare the last objective 
	// points above l cannot dominate points that start below l, and we don't need to compare the last objective 
	fs[fr].nPoints = 1;
	for (int i = 1; i < l; i++) {
		int j = 0;
		while (j < fs[fr].nPoints) {
			int k;
			switch (dominates2way(fs[fr].points[i], fs[fr].points[j], n - 2)) {
			case  0:
				j++;
				break;
			case -1: // AT THIS POINT WE KNOW THAT i CANNOT BE DOMINATED BY ANY OTHER PROMOTED POINT j 
				// SWAP i INTO j, AND 1-WAY DOM FOR THE REST OF THE js 
				t = fs[fr].points[j];
				fs[fr].points[j] = fs[fr].points[i];
				fs[fr].points[i] = t;
				while (j < fs[fr].nPoints - 1 && dominates1way(fs[fr].points[j], fs[fr].points[fs[fr].nPoints - 1], n - 1)) {
					fs[fr].nPoints--;
				}
				k = j + 1;
				while (k < fs[fr].nPoints) {
					if (dominates1way(fs[fr].points[j], fs[fr].points[k], n - 2)) {
						t = fs[fr].points[k];
						fs[fr].nPoints--;
						fs[fr].points[k] = fs[fr].points[fs[fr].nPoints];
						fs[fr].points[fs[fr].nPoints] = t;
					}
					else {
						k++;
					}
				}
			default:
				j = fs[fr].nPoints + 1;
			}
		}
		if (j == fs[fr].nPoints) {
			t = fs[fr].points[fs[fr].nPoints];
			fs[fr].points[fs[fr].nPoints] = fs[fr].points[i];
			fs[fr].points[i] = t;
			fs[fr].nPoints++;
		}
	}
	safe = WORSE(l, fs[fr].nPoints);
	for (int i = l; i < psize[fr]; i++) {
		int j = 0;
		while (j < safe) {
			if (dominates1way(fs[fr].points[j], fs[fr].points[i], n - 2)) {
				j = fs[fr].nPoints + 1;
			}
			else {
				j++;
			}
		}
		while (j < fs[fr].nPoints) {
			int k;
			switch (dominates2way(fs[fr].points[i], fs[fr].points[j], n - 1)) {
			case  0:
				j++;
				break;
			case -1: // AT THIS POINT WE KNOW THAT i CANNOT BE DOMINATED BY ANY OTHER PROMOTED POINT j 
				// SWAP i INTO j, AND 1-WAY DOM FOR THE REST OF THE js 
				t = fs[fr].points[j];
				fs[fr].points[j] = fs[fr].points[i];
				fs[fr].points[i] = t;
				while (j < fs[fr].nPoints - 1 && dominates1way(fs[fr].points[j], fs[fr].points[fs[fr].nPoints - 1], n - 1)) {
					fs[fr].nPoints--;
				}
				k = j + 1;
				while (k < fs[fr].nPoints) {
					if (dominates1way(fs[fr].points[j], fs[fr].points[k], n - 1)) {
						t = fs[fr].points[k];
						fs[fr].nPoints--;
						fs[fr].points[k] = fs[fr].points[fs[fr].nPoints];
						fs[fr].points[fs[fr].nPoints] = t;
					}
					else {
						k++;
					}
				}
			default:
				j = fs[fr].nPoints + 1;
			}
		}
		if (j == fs[fr].nPoints) {
			t = fs[fr].points[fs[fr].nPoints];
			fs[fr].points[fs[fr].nPoints] = fs[fr].points[i];
			fs[fr].points[i] = t;
			fs[fr].nPoints++;
		}
	}
	int last = lastp[fr];
	nextp[fr][last] = p - last;
	prevp[fr][p] = p - last;
	nextp[fr][p] = maxm - p;
	lastp[fr] = p;
	psize[fr] += 1 - removed;
	fr++;
}

double HV::hv2(FRONT ps, int k)
// returns the hypervolume of ps[0 .. k-1] in 2D 
// assumes that ps is sorted improving
{
	double volume = ps.points[0].objectives[0] * ps.points[0].objectives[1];
	for (int i = 1; i < k; i++) {
		volume += ps.points[i].objectives[1] * (ps.points[i].objectives[0] - ps.points[i - 1].objectives[0]);
	}
	return volume;
}

double HV::inclhv(POINT p)
// returns the inclusive hypervolume of p
{
	double volume = 1;
	for (int i = 0; i < n; i++) {
		volume *= p.objectives[i];
	}
	return volume;
}

double HV::inclhv2(POINT p, POINT q)
// returns the hypervolume of {p, q}
{
	double vp = 1;
	double vq = 1;
	double vpq = 1;
	for (int i = 0; i < n; i++) {
		vp *= p.objectives[i];
		vq *= q.objectives[i];
		vpq *= WORSE(p.objectives[i], q.objectives[i]);
	}
	return vp + vq - vpq;
}

double HV::inclhv3(POINT p, POINT q, POINT r)
// returns the hypervolume of {p, q, r}
{
	double vp = 1;
	double vq = 1;
	double vr = 1;
	double vpq = 1;
	double vpr = 1;
	double vqr = 1;
	double vpqr = 1;
	for (int i = 0; i < n; i++) {
		vp *= p.objectives[i];
		vq *= q.objectives[i];
		vr *= r.objectives[i];
		if (BEATS(p.objectives[i], q.objectives[i])) {
			if (BEATS(q.objectives[i], r.objectives[i])) {
				vpq *= q.objectives[i];
				vpr *= r.objectives[i];
				vqr *= r.objectives[i];
				vpqr *= r.objectives[i];
			}
			else {
				vpq *= q.objectives[i];
				vpr *= WORSE(p.objectives[i], r.objectives[i]);
				vqr *= q.objectives[i];
				vpqr *= q.objectives[i];
			}
		}
		else if (BEATS(p.objectives[i], r.objectives[i])) {
			vpq *= p.objectives[i];
			vpr *= r.objectives[i];
			vqr *= r.objectives[i];
			vpqr *= r.objectives[i];
		}
		else {
			vpq *= p.objectives[i];
			vpr *= p.objectives[i];
			vqr *= WORSE(q.objectives[i], r.objectives[i]);
			vpqr *= p.objectives[i];
		}

	}
	return vp + vq + vr - vpq - vpr - vqr + vpqr;
}

double HV::inclhv4(POINT p, POINT q, POINT r, POINT s)
// returns the hypervolume of {p, q, r, s}
{
	double vp = 1;
	double vq = 1;
	double vr = 1;
	double vs = 1;
	double vpq = 1;
	double vpr = 1;
	double vps = 1;
	double vqr = 1;
	double vqs = 1;
	double vrs = 1;
	double vpqr = 1;
	double vpqs = 1;
	double vprs = 1;
	double vqrs = 1;
	double vpqrs = 1;
	for (int i = 0; i < n; i++) {
		vp *= p.objectives[i];
		vq *= q.objectives[i];
		vr *= r.objectives[i];
		vs *= s.objectives[i];
		if (BEATS(p.objectives[i], q.objectives[i])) {
			if (BEATS(q.objectives[i], r.objectives[i])) {
				if (BEATS(r.objectives[i], s.objectives[i])) {
					vpq *= q.objectives[i];
					vpr *= r.objectives[i];
					vps *= s.objectives[i];
					vqr *= r.objectives[i];
					vqs *= s.objectives[i];
					vrs *= s.objectives[i];
					vpqr *= r.objectives[i];
					vpqs *= s.objectives[i];
					vprs *= s.objectives[i];
					vqrs *= s.objectives[i];
					vpqrs *= s.objectives[i];
				}
				else {
					OBJECTIVE z1 = WORSE(q.objectives[i], s.objectives[i]);
					vpq *= q.objectives[i];
					vpr *= r.objectives[i];
					vps *= WORSE(p.objectives[i], s.objectives[i]);
					vqr *= r.objectives[i];
					vqs *= z1;
					vrs *= r.objectives[i];
					vpqr *= r.objectives[i];
					vpqs *= z1;
					vprs *= r.objectives[i];
					vqrs *= r.objectives[i];
					vpqrs *= r.objectives[i];
				}
			}
			else if (BEATS(q.objectives[i], s.objectives[i])) {
				vpq *= q.objectives[i];
				vpr *= WORSE(p.objectives[i], r.objectives[i]);
				vps *= s.objectives[i];
				vqr *= q.objectives[i];
				vqs *= s.objectives[i];
				vrs *= s.objectives[i];
				vpqr *= q.objectives[i];
				vpqs *= s.objectives[i];
				vprs *= s.objectives[i];
				vqrs *= s.objectives[i];
				vpqrs *= s.objectives[i];
			}
			else {
				OBJECTIVE z1 = WORSE(p.objectives[i], r.objectives[i]);
				vpq *= q.objectives[i];
				vpr *= z1;
				vps *= WORSE(p.objectives[i], s.objectives[i]);
				vqr *= q.objectives[i];
				vqs *= q.objectives[i];
				vrs *= WORSE(r.objectives[i], s.objectives[i]);
				vpqr *= q.objectives[i];
				vpqs *= q.objectives[i];
				vprs *= WORSE(z1, s.objectives[i]);
				vqrs *= q.objectives[i];
				vpqrs *= q.objectives[i];
			}
		}
		else if (BEATS(q.objectives[i], r.objectives[i])) {
			if (BEATS(p.objectives[i], s.objectives[i])) {
				OBJECTIVE z1 = WORSE(p.objectives[i], r.objectives[i]);
				OBJECTIVE z2 = WORSE(r.objectives[i], s.objectives[i]);
				vpq *= p.objectives[i];
				vpr *= z1;
				vps *= s.objectives[i];
				vqr *= r.objectives[i];
				vqs *= s.objectives[i];
				vrs *= z2;
				vpqr *= z1;
				vpqs *= s.objectives[i];
				vprs *= z2;
				vqrs *= z2;
				vpqrs *= z2;
			}
			else {
				OBJECTIVE z1 = WORSE(p.objectives[i], r.objectives[i]);
				OBJECTIVE z2 = WORSE(r.objectives[i], s.objectives[i]);
				vpq *= p.objectives[i];
				vpr *= z1;
				vps *= p.objectives[i];
				vqr *= r.objectives[i];
				vqs *= WORSE(q.objectives[i], s.objectives[i]);
				vrs *= z2;
				vpqr *= z1;
				vpqs *= p.objectives[i];
				vprs *= z1;
				vqrs *= z2;
				vpqrs *= z1;
			}
		}
		else if (BEATS(p.objectives[i], s.objectives[i])) {
			vpq *= p.objectives[i];
			vpr *= p.objectives[i];
			vps *= s.objectives[i];
			vqr *= q.objectives[i];
			vqs *= s.objectives[i];
			vrs *= s.objectives[i];
			vpqr *= p.objectives[i];
			vpqs *= s.objectives[i];
			vprs *= s.objectives[i];
			vqrs *= s.objectives[i];
			vpqrs *= s.objectives[i];
		}
		else {
			OBJECTIVE z1 = WORSE(q.objectives[i], s.objectives[i]);
			vpq *= p.objectives[i];
			vpr *= p.objectives[i];
			vps *= p.objectives[i];
			vqr *= q.objectives[i];
			vqs *= z1;
			vrs *= WORSE(r.objectives[i], s.objectives[i]);
			vpqr *= p.objectives[i];
			vpqs *= p.objectives[i];
			vprs *= p.objectives[i];
			vqrs *= z1;
			vpqrs *= p.objectives[i];
		}
	}
	return vp + vq + vr + vs - vpq - vpr - vps - vqr - vqs - vrs + vpqr + vpqs + vprs + vqrs - vpqrs;
}

double HV::exclhv(FRONT ps, int p)
// returns the exclusive hypervolume of ps[p] relative to ps[0 .. p-1] 
{
	makeDominatedBit(ps, p);
	double volume = inclhv(ps.points[p]) - hv(fs[fr - 1]);
	fr--;
	return volume;
}

double HV::hv(FRONT ps)
// returns the hypervolume of ps[0 ..] 
{
	// process small fronts with the IEA 
	switch (ps.nPoints) {
	case 1:
		return inclhv(ps.points[0]);
	case 2:
		return inclhv2(ps.points[0], ps.points[1]);
	case 3:
		return inclhv3(ps.points[0], ps.points[1], ps.points[2]);
	case 4:
		return inclhv4(ps.points[0], ps.points[1], ps.points[2], ps.points[3]);
	}

	// these points need sorting 
	qsort(&ps.points[safe], ps.nPoints - safe, sizeof(POINT), greater);
	// n = 2 implies that safe = 0 
	if (n == 2) {
		return hv2(ps, ps.nPoints);
	}
	// these points don't NEED sorting, but it helps 
	qsort(ps.points, safe, sizeof(POINT), greaterabbrev);

	if (n == 3 && safe > 0) {
		double volume = ps.points[0].objectives[2] * hv2(ps, safe);
		n--;
		for (int i = 0; i < safe; i++) {
			nextp[fr][i] = 1;
			prevp[fr][i] = 1;
		}
		psize[fr] = safe;
		firstp[fr] = 0;
		lastp[fr] = safe - 1;
		nextp[fr][safe - 1] = maxm - safe + 1;
		for (int i = safe; i < ps.nPoints; i++) {
			// we can ditch dominated points here, but they will be ditched anyway in makeDominatedBit 
			volume += ps.points[i].objectives[n] * exclhv(ps, i);
		}
		n++;
		return volume;
	}
	else {
		double volume = inclhv4(ps.points[0], ps.points[1], ps.points[2], ps.points[3]);
		n--;
		for (int i = 0; i < 4; i++) {
			nextp[fr][i] = 1;
			prevp[fr][i] = 1;
		}
		psize[fr] = 4;
		firstp[fr] = 0;
		lastp[fr] = 3;
		nextp[fr][3] = maxm - 3;
		for (int i = 4; i < ps.nPoints; i++) {
			// we can ditch dominated points here, but they will be ditched anyway in makeDominatedBit 
			volume += ps.points[i].objectives[n] * exclhv(ps, i);
		}
		n++;
		return volume;
	}
}

double HV::HV_calc(pareto_front & PF, std::vector<std::vector<double>> & Real_PF)
// processes each front from the file 
{
	double HV_val = 0;		//output
	maxm = PF.Size_Show();		
	maxn = Real_PF[0].size();		
	if (maxm == 0)
		return 0;
	if (maxn == 0)
		abort();
	//evaluate the reference point
	//find the max values of each objective for the Real_PF
	std::vector<double> Real_PF_max_f(maxn, -1e30);

	for (int i = 0; i < Real_PF.size(); i++)
	{
		for (int j = 0; j < maxn; j++)
		{
			if (Real_PF[i][j] > Real_PF_max_f[j])
			{
				Real_PF_max_f[j] = Real_PF[i][j];
			}
		}
	}
	//increase the value of the ref points and assign them
	// initialise the reference point
	ref.objectives = (OBJECTIVE*)malloc(sizeof(OBJECTIVE) * maxn);
	double ref_norm_mod = 1;
	for (int i = 0; i < maxn; i++)
	{
		Real_PF_max_f[i] += 1;
		ref.objectives[i] = Real_PF_max_f[i];
		ref_norm_mod *= Real_PF_max_f[i];
	}

	//copy and refine the front - set the REAL_PF_max as the maximum values of all objectives and eliminate the points dominated
	std::vector<std::vector<double>> PF_temp;
	for (int i = 0; i < maxm; i++)
		PF_temp.push_back(PF.Indiv_Show(i).Fitness_Show());
	
	Refine_PF(PF_temp, Real_PF_max_f);
	if (PF_temp.empty())
		return 0.;
	else
		maxm = PF_temp.size();


	//FILECONTENTS *f = readFile(argv[1]);
	FRONT *f = (FRONT*)malloc(sizeof(FRONT));
	f->points = (POINT*)malloc(sizeof(POINT) * maxm);
	for (int j = 0; j < maxm; j++) 
	{
		f->points[j].objectives = (OBJECTIVE*)malloc(sizeof(OBJECTIVE) * (maxn));
	}

	//add your front to the *f

	for (int j = 0; j < maxm; j++)
	{
		std::vector<double> temp_fit = PF_temp[j];
		for (int i = 0; i < maxn; i++)
		{
			f->points[j].objectives[i] = temp_fit[i];
		}
	}
	f->n = maxn;
	f->nPoints = maxm;

	
	


	// allocate memory
	int maxdepth = maxn - 2;
	fs = (FRONT*)malloc(sizeof(FRONT) * maxdepth);
	for (int i = 0; i < maxdepth; i++) {
		fs[i].points = (POINT*)malloc(sizeof(POINT) * maxm);
		for (int j = 0; j < maxm; j++) {
			fs[i].points[j].objectives = (OBJECTIVE*)malloc(sizeof(OBJECTIVE) * (maxn - i - 1));
		}
	}
	nextp = (int**)malloc(sizeof(int*) * maxdepth);
	for (int i = 0; i < maxdepth; i++) {
		nextp[i] = (int*)malloc(sizeof(int) * maxm);
	}
	prevp = (int**)malloc(sizeof(int*) * maxdepth);
	for (int i = 0; i < maxdepth; i++) {
		prevp[i] = (int*)malloc(sizeof(int) * maxm);
	}
	firstp = (int*)malloc(sizeof(int) * maxdepth);
	lastp = (int*)malloc(sizeof(int) * maxdepth);
	psize = (int*)malloc(sizeof(int) * maxdepth);

	
	// modify the objective values relative to the reference point 

	for (int j = 0; j < f->nPoints; j++) {
		for (int k = 0; k < f->n; k++) {
			f->points[j].objectives[k] = fabs(f->points[j].objectives[k] - ref.objectives[k]);
		}
	}
	

	totaltime = 0;
	//for (int i = 0; i < f->nFronts; i++) {      
		//struct timeval tv1, tv2;
		//struct rusage ru_before, ru_after;
		//getrusage (RUSAGE_SELF, &ru_before);

	n = f->n;
	//safe = 0;
	//printf("hv(%d) = %1.10f\n", i+1, hv(f->fronts[i])); 
	HV_val = hv(f[0]);
	//getrusage (RUSAGE_SELF, &ru_after);
	//tv1 = ru_before.ru_utime;
	//tv2 = ru_after.ru_utime;
	//printf("Time: %f (s)\n", tv2.tv_sec + tv2.tv_usec * 1e-6 - tv1.tv_sec - tv1.tv_usec * 1e-6);
	//totaltime += tv2.tv_sec + tv2.tv_usec * 1e-6 - tv1.tv_sec - tv1.tv_usec * 1e-6;
//}
//printf("Total time = %f (s)\n", totaltime);
	
	for (int i = 0; i < maxdepth; i++) 
	{
	
		for (int j = 0; j < maxm; j++) 
		{
			free(fs[i].points[j].objectives);
		}
		free(fs[i].points);
	}
	for (long j = 0; j < maxm; j++)
	{
		free(f->points[j].objectives);
	}
	free(f->points);

	/*for (int i = 0; i < maxdepth; i++) 
	{
		free(nextp[i]);
	}
	for (int i = 0; i < maxdepth; i++) 
	{
		free(prevp[i]);
	}*/


	free(fs);    // memory management stuff 
	fr = 0;   // current depth 
	maxm = 0; // identify the biggest fronts in the file 
	maxn = 0;
	safe = 0;     // the number of points that don't need sorting
	free(nextp);
	free(prevp);


	free(f);
	free(firstp);
	free(lastp);
	free(psize);
	free(ref.objectives);

	//normalise the HV_val
	return HV_val/ref_norm_mod;
}


void HV::Refine_PF(std::vector<std::vector<double>> &PF_temp, std::vector<double> &Real_PF_max_f)
{
	short n_obj = Real_PF_max_f.size();
	//remove dominated points
	for (int i = 0; i < PF_temp.size(); i++)
	{
		//copy the fitness of individual
		std::vector<double> PF_temp_fit = PF_temp[i];
		short flag1 = 0, flag2 = 0;
		for (int k = 0; k < n_obj; k++)
		{
			//Copy ith fitness of each individual
			double fit_ind1;
			double fit_ind2;
			
			fit_ind1 = PF_temp_fit[k];		//ith fitness of 1st individual
			fit_ind2 = Real_PF_max_f[k];		//ith fitness of 2nd individual
			
			//Check which fitness is greater
			if (fit_ind1 < fit_ind2)
				flag1 = 1;
			else if (fit_ind1 > fit_ind2)
			{
				flag2 = 1;
				PF_temp[i][k] = fit_ind2;
			}
		}
		//Check which individual dominate
		if (flag1 == 0 && flag2 == 1)
		{
			PF_temp.erase(PF_temp.begin() + i);
			i--;
		}
	}
}
