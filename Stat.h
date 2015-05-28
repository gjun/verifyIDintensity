/*
 *   Author: Goo Jun (goo.jun@uth.tmc.edu)
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef __STAT_H__
#define __STAT_H__

#include <stdio.h>
#include <boost/math/tools/minima.hpp>
#include <boost/bind.hpp>

#include "common.h"
#include "Intensity.h"

class cStat
{
	public:
		virtual void Initialize(unsigned int, unsigned int, const char*) = 0;
		virtual void Estimate()  = 0;
		virtual double GetLLK(double)  = 0;
		virtual void WriteStat(const char*) = 0;
		virtual bool ReadStat(const char*) = 0;
		virtual void CheckCovs(double) = 0;
		virtual ~cStat() {};
		double minimizeLLK(unsigned int, double&, double&) ;

		double GetCallRate(unsigned int i) {return CallRate[i];};
		double GetMinPr(unsigned int);
		bool IsZeroOut(unsigned int j) {return bZeroOut[j];};
		bool ReadAbf(const char*);

	protected:
		double* CallRate;
		bool* bZeroOut;
		double *Pr;
		unsigned int n_marker, n_sample;
		unsigned int curSample;
		cIntensity Intensity;
};

double cStat::GetMinPr(unsigned int j)
{
	if (Pr[j*4+0]<Pr[j*4+2])
	{
		return (sqrt(Pr[j*4+0]) + 0.5*Pr[j*4+1]);
	}
	else
	{
		return (sqrt(Pr[j*4+2]) + 0.5*Pr[j*4+1]);
	}
}

double cStat::minimizeLLK(unsigned int i, double& alpha, double& LLK0)
{
	double LLK;
	double minLLK, minalpha = 0.0;
	double left, right;

	curSample = i;

	LLK0 = GetLLK(0);
	minLLK = LLK0;


	for(alpha=0.02;alpha<=0.5;alpha+=0.02)
	{
		LLK = GetLLK(alpha);
		if (minLLK>LLK)
		{
			minLLK = LLK;
			minalpha = alpha;
		}
	}


	if (minalpha < 0.01)
	{
		left = 0.0;
		right = 0.02;
	}
	else if (minalpha >0.49)
	{
		left = 0.48;
		right = 0.5;
	}
	else
	{
		left = minalpha-0.02;
		right = minalpha+0.02;
	}

	std::pair<double,double> r = boost::math::tools::brent_find_minima(boost::bind(&cStat::GetLLK,this,_1), left, right, 20);

	alpha = r.first;
	return r.second;
}


bool cStat::ReadAbf(const char* FileName)
{
	if (FileName != NULL)
	{
		FILE *fp = fopen(FileName,"r");
		char line[1024];
		char name[1024],sAbf[1024];
		fgets(line,1024,fp);
		printf("%s\n",line);
		for (unsigned j=0;j<n_marker;j++)
		{
			fscanf(fp, "%s %s", name, sAbf);

			//			printf("%s\tAbf:%s\n",name,sAbf);

			if	(strcmp(sAbf, "NA"))
			{
				float Abf = atof(sAbf);
				if (Abf>=0 && Abf <=1)
				{
					Pr[j*4+0] = (1-Abf) * (1-Abf);
					Pr[j*4+1] = (1-Abf) * Abf;
					Pr[j*4+2] = Abf*Abf;
				}
				else
				{
					std::cerr << "Abf not in [0,1]! \n";
				}
			}
			else
			{
				bZeroOut[j] = true;
			}
		}
	}
	return true;
}

#endif
