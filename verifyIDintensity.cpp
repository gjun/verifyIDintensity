/*
 *  Copyright (C) 2010-2014  Regents of the University of Michigan
 *
 *	Author: Goo Jun
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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <limits>
#include <boost/math/tools/minima.hpp>
#include <boost/bind.hpp>
#include <float.h>

#include "tclap/CmdLine.h"
#include "tclap/Arg.h"

const double pi=3.14159265;
double min_af;

double det(double *M)
{
	return (M[0]*M[3]-M[1]*M[2]);
}

double normpdf2(float nA, float nB, double m0, double m1, double* cov)
{
	double val,  D, prc[4];
	double y[2];
	y[0] = (double)nA;
	y[1] = (double)nB;

	D = det(cov);
	prc[0] = cov[3]/D;
	prc[3] = cov[0]/D;
	prc[1] = prc[2] = -1.0*cov[1]/D;

	val= -0.5*((y[0]-m0)*((y[0]-m0)*prc[0] + (y[1]-m1)*prc[1]) + (y[1]-m1)*((y[0]-m0)*prc[2]+(y[1]-m1)*prc[3]));
	val= 0.5*exp(val)/(pi*sqrt(D));

	return val;
}

double lognormpdf2(float nA, float nB, double m0, double m1, double* cov)
{
	double val,  D, prc[4];
	double y[2];
	y[0] = (double)nA;
	y[1] = (double)nB;


	D = det(cov);
	prc[0] = cov[3]/D;
	prc[3] = cov[0]/D;
	prc[1] = prc[2] = -1.0*cov[1]/D;

	val = -0.5*((y[0]-m0)*((y[0]-m0)*prc[0] + (y[1]-m1)*prc[1]) + (y[1]-m1)*((y[0]-m0)*prc[2]+(y[1]-m1)*prc[3]));
	val = val-log(2*pi)-0.5*log(D);

	return val;
}

class cIntensity
{
	public:
		void Initialize(const char*, unsigned int, unsigned int);
		void ReadData(unsigned int, unsigned int, float&, float&, float&, unsigned short&);
		~cIntensity();
	private:
		bool bInitialized;
		std::ifstream inFile;
		unsigned int pos_sample;
		unsigned int n_sample, n_marker;
		float *nA, *nB, *GC;
		unsigned short *gType;
};


void cIntensity::Initialize(const char* inFileName, unsigned int n_s, unsigned int n_m)
{
	inFile.open(inFileName, std::ios::in | std::ios::binary);
	if (!inFile.is_open())    
	{
		std::cerr << "Error opening file " << inFileName << "\n";
		exit (1);
	}

	n_sample = n_s;
	n_marker = n_m;
#ifdef DEBUG
	std::cerr << "Opening file " << inFileName << "\n";
#endif

	nA = new float[n_marker];
	nB = new float[n_marker];
	GC = new float[n_marker];
	gType = new unsigned short[n_marker];

	pos_sample = n_s+1;
	bInitialized = true;
}

void cIntensity::ReadData(unsigned int arg_i, unsigned int arg_j, float& arg_nA, float& arg_nB, float& arg_GC, unsigned short& arg_gType)
{
	unsigned short A, B;

	// Read intensity data of j-th marker of i-th sample from adpc.bin file

	// Read all intensity data of i-th sample at a time, for minimized disk I/O

	if (arg_i != pos_sample)
	{
#ifdef DEBUG
		std::cerr << "Reading intensity data of "<< arg_i << "-th sample from input file\n";
#endif
		std::streampos seek_pos = 16+(unsigned long)arg_i*(unsigned long)n_marker*18;
		inFile.seekg(seek_pos, std::ios::beg);

		for (unsigned int j=0;j<n_marker;j++)
		{
			try
			{
				inFile.read(reinterpret_cast <char *> (&A), sizeof(unsigned short));
				inFile.read(reinterpret_cast <char *> (&B), sizeof(unsigned short));
				inFile.read(reinterpret_cast <char *> (&nA[j]), sizeof(float));
				inFile.read(reinterpret_cast <char *> (&nB[j]), sizeof(float));
				inFile.read(reinterpret_cast <char *> (&GC[j]), sizeof(float));
				inFile.read(reinterpret_cast <char *> (&gType[j]), sizeof(unsigned short));
			}
			catch(int e)
			{
				std::cerr << "Error reading intensity data sample=" << arg_i << ", marker=" << j << " exception no " << e << "\n";
				exit (1);
			}
		}
		pos_sample = arg_i;
	}

	arg_nA = nA[arg_j];
	arg_nB = nB[arg_j];
	arg_GC = GC[arg_j];
	arg_gType = gType[arg_j];
}

cIntensity::~cIntensity()
{
	if (bInitialized)
	{
		if (inFile.is_open()) inFile.close();
		if (nA) delete[] nA;
		if (nB) delete[] nB;
		if (GC) delete[] GC;
		if (gType) delete[] gType;
	}
	bInitialized=false;
}

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
	curSample = i;
	double LLK;
	double minLLK, minalpha = 0.0;
	double left, right;

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

class cMarkerStat : public cStat
{
	double *meanA, *meanB, *covs;
	unsigned int *nums;

	public:
	void Initialize(unsigned int, unsigned int, const char*);
	void CheckCovs(double);
	void Estimate();
	double GetLLK(double);
	double GetCallRate(unsigned int);
	bool IsZeroOut(unsigned int);
	void WriteStat(const char*);
	bool ReadStat(const char*);
	~cMarkerStat();
};

cMarkerStat::~cMarkerStat()
{
	delete[] nums;
	delete[] Pr;
	delete[] CallRate;
	delete[] bZeroOut;

	delete[] meanA;
	delete[] meanB;
	delete[] covs;
}

double cMarkerStat::GetCallRate(unsigned int i)
{
	return CallRate[i];
}

bool cMarkerStat::IsZeroOut(unsigned int j)
{
	return bZeroOut[j];
}

void cMarkerStat::Initialize(unsigned int nS, unsigned int nM, const char* sIntensity)
{
	n_sample = nS;
	n_marker = nM;

	Pr = new double [4*n_marker];
	CallRate = new double [n_sample];
	bZeroOut = new bool [n_marker];

	nums = new unsigned [4*n_marker];
	meanA = new double [4*n_marker];
	meanB = new double [4*n_marker];
	covs = new double [4*4*n_marker];

	for (unsigned j=0;j<n_marker;j++)
	{
		bZeroOut[j] = false;
	}

	Intensity.Initialize(sIntensity, nS, nM);
}

bool cMarkerStat::ReadStat(const char* FileName)
{
	std::ifstream  InStatsFile(FileName, std::ios::in | std::ios::binary);
	if (!InStatsFile.is_open())
	{
		std::cout<<"Statistics file does not exist. Stats will be estimated from data." << std::endl;
		return false;
	}
	// If pre-calculated statistics exists, read
	InStatsFile.seekg(0, std::ios::beg);
	InStatsFile.read(reinterpret_cast <char *> (meanA), sizeof(double)*4*n_marker);
	InStatsFile.read(reinterpret_cast <char *> (meanB), sizeof(double)*4*n_marker);
	InStatsFile.read(reinterpret_cast <char *> (covs), sizeof(double)*4*4*n_marker);
	InStatsFile.read(reinterpret_cast <char *> (Pr), sizeof(double)*4*n_marker);
	InStatsFile.read(reinterpret_cast <char *> (nums), sizeof(unsigned)*4*n_marker);
	InStatsFile.read(reinterpret_cast <char *> (bZeroOut), sizeof(bool)*n_marker);
	InStatsFile.read(reinterpret_cast <char *> (CallRate), sizeof(double)*n_sample);

	InStatsFile.close();
	return true;
}

void cMarkerStat::WriteStat(const char* FileName)
{
	std::cout<<"Creating a new statistics file." << std::endl;
	//	Logger::gLogger->writeLog("Creating statistics file %s",FileName);

	std::ofstream  OutStatsFile(FileName, std::ios::out | std::ios::binary);
	if (!OutStatsFile.is_open())
	{
		//		Logger::gLogger->error("Error opening %s to write",sStatFile.c_str());
	}

	OutStatsFile.write(reinterpret_cast <char *> (meanA), sizeof(double)*4*n_marker);
	OutStatsFile.write(reinterpret_cast <char *> (meanB), sizeof(double)*4*n_marker);
	OutStatsFile.write(reinterpret_cast <char *> (covs), sizeof(double)*4*n_marker*4);
	OutStatsFile.write(reinterpret_cast <char *> (Pr), sizeof(double)*4*n_marker);
	OutStatsFile.write(reinterpret_cast <char *> (nums), sizeof(unsigned)*4*n_marker);
	OutStatsFile.write(reinterpret_cast <char *> (bZeroOut), sizeof(bool)*n_marker);
	OutStatsFile.write(reinterpret_cast <char *> (CallRate), sizeof(double)*n_sample);

	OutStatsFile.close();
}

void cMarkerStat::CheckCovs(double alphastep)
{
	for (unsigned j=0;j<n_marker;j++)
	{
		if (!bZeroOut[j])
		{
			unsigned n_alpha = (0.5/alphastep)+1;
			for (unsigned k=0;k<n_alpha;k++)
			{
				double alpha = alphastep*k;
				for(unsigned k1=0; k1<3; k1++)
				{
					for(unsigned k2=0; k2<3; k2++)
					{
						double mix_cov[4];

						for(unsigned l=0; l<4; l++)
						{
							mix_cov[l] = alpha*alpha*covs[(j*4+k1)*4+l] + (1-alpha)*(1-alpha)*covs[(j*4+k2)*4+l];
						}

						if (det(mix_cov)<1e-13)
						{
							bZeroOut[j] = true;
						}
					}
				}
			}
		}
	}
}


void cMarkerStat::Estimate()
{
	float nA, nB, GC;
	unsigned short gType;

	// If no pre-calculated statistics is given,
#ifdef DEBUG
	std::cerr << "Estimating per-marker statistics\n";
#endif
	double* sumA = new double[4*n_marker]();
	double* sumB = new double[4*n_marker]();

	double* sumAA = new double[4*n_marker]();
	double* sumBB = new double[4*n_marker]();
	double* sumAB = new double[4*n_marker]();

	for (unsigned i=0; i<n_sample; i++)
	{
		for (unsigned j=0; j<n_marker; j++)
		{
			Intensity.ReadData(i, j, nA, nB, GC, gType);

			if ((GC==0.0 && gType==3))
			{
				bZeroOut[j] = true;
			}
			else
			{
				bZeroOut[j] = (bZeroOut[j] || false);
			}
		}
	}

	unsigned int n_zeroout = 0;
	for (unsigned j=0;j<n_marker;j++)
	{
		if (bZeroOut[j]) 
		{
			n_zeroout++;
		}
	}

#ifdef DEBUG
	std::cerr << "Number of zero-ed out SNPs from the input file : "<< n_zeroout << "\n";
#endif

	unsigned int n_nonzero = n_marker - n_zeroout;

	for (unsigned i=0; i<n_sample; i++)
	{
		unsigned int callNum = 0;

		for(unsigned j=0; j<n_marker; j++)
		{
			try
			{
				Intensity.ReadData(i, j, nA, nB, GC, gType);
			}
			catch (int e)
			{
				std::cerr << "cMarkerStat::Estimate()  Error reading intensity data\n";
			}

			if (!bZeroOut[j] && gType<3 && !isnan(nA) && !isnan(nB))
			{
				callNum++;
			}
			else if (isnan(nA) || isnan(nB))
			{
				bZeroOut[j] = true;
			}
		}
		CallRate[i] = ((double)callNum) / ((double)n_nonzero);

		std::cerr << "Sample " << i << " Call Rate " << CallRate[i] << "\n";

		if (CallRate[i]>0.99) // Good sample
		{
			// Get Statistics (mean, cov) for each allele type
			for(unsigned j=0; j<n_marker; j++)
			{
				Intensity.ReadData(i, j, nA, nB, GC, gType);

				if (!(bZeroOut[j]) && gType<3)
				{
					sumA[j*4 + gType] += (double)nA;
					sumB[j*4 + gType] += (double)nB;

					sumAA[j*4 + gType] += (double)nA*(double)nA;
					sumBB[j*4 + gType] += (double)nB*(double)nB;
					sumAB[j*4 + gType] += (double)nA*(double)nB;

					nums[j*4 + gType] = nums[j*4 + gType]+1;

					if (isnan(sumA[j*4+gType]) || isnan(sumB[j*4+gType]))
					{
						std::cerr << "Error: Sum is NAN though previously checked \n";
						exit (1);
					}
				}
			}
		}
	}


#ifdef DEBUG
	std::cerr << "Estimating Mean and Cov.... \n";
#endif
	for(unsigned j=0; j<n_marker; j++)
	{
		if (!bZeroOut[j])
		{
			for(unsigned i=0; i<3; i++)
			{
				if (nums[j*4+i]>2)
				{
					meanA[j*4+i] = sumA[j*4+i]/nums[j*4+i];
					meanB[j*4+i] = sumB[j*4+i]/nums[j*4+i];

					covs[(j*4+i)*4+0] = sumAA[j*4+i]/(nums[j*4+i]) - (meanA[j*4+i]*meanA[j*4+i]);
					covs[(j*4+i)*4+1] = covs[(j*4+i)*4+2] = sumAB[j*4+i]/(nums[j*4+i]) - (meanA[j*4+i]*meanB[j*4+i]);
					covs[(j*4+i)*4+3] = sumBB[j*4+i]/(nums[j*4+i]) - (meanB[j*4+i]*meanB[j*4+i]);

				}
				else
				{
					bZeroOut[j] = true;
				}
				if (isnan(meanA[j*4+i]) || isnan(meanB[j*4+i]))
				{
					std::cerr << "Error: Mean is NAN though previously checked \n";
					exit (1);
				}
			}

			unsigned tot_num = nums[j*4+0]+nums[j*4+1]+nums[j*4+2];
			double Pr_A= ((double)(nums[j*4+0]) + (double)(nums[j*4+1])*0.5) / (double)tot_num;
			double Pr_B = ((double)(nums[j*4+2]) + (double)(nums[j*4+1])*0.5) / (double)tot_num;
			Pr[j*4+0] = Pr_A*Pr_A;
			Pr[j*4+1] = Pr_A*Pr_B;
			Pr[j*4+2] = Pr_B*Pr_B;
		}
	}

	delete [] sumA;
	delete [] sumB;
	delete [] sumAA;
	delete [] sumAB;
	delete [] sumBB;
}


double cMarkerStat::GetLLK(double alpha)
{
	double llk = 0;
	double lks = 0;
	float nA, nB, GC;
	unsigned short gType;

	unsigned int i=curSample;

	for(unsigned int j=0;j<n_marker;++j)
	{
		lks=0;
		if (GetMinPr(j)>min_af)
		{
			Intensity.ReadData(i, j, nA, nB, GC, gType);

			if (!bZeroOut[j] && !isnan(nA) && !isnan(nB))
			{
				for(unsigned k1=0; k1<3; k1++)
				{
					for(unsigned k2=0; k2<3; k2++)
					{
						double mA, mB, mix_cov[4];
						mA = meanA[j*4+k1]*alpha + meanA[j*4+k2]*(1-alpha);
						mB = meanB[j*4+k1]*alpha + meanB[j*4+k2]*(1-alpha);

						for(unsigned l=0; l<4; l++)
						{
							mix_cov[l] = alpha*alpha*covs[(j*4+k1)*4+l] + (1-alpha)*(1-alpha)*covs[(j*4+k2)*4+l];
						}

						double tmpval=normpdf2(nA, nB, mA, mB, mix_cov) * Pr[j*4+k1] * Pr[j*4+k2];
						if (isnan(tmpval))
						{
							fprintf(stderr,"tmpval is nan %f %f, %f %f %f, %f %f %f\n", nA, nB, meanA[j*4+0], meanA[j*4+1], meanA[j*4+2], mA, mB, det(mix_cov) );
						}
						lks+=tmpval;
					}
				}
				if (lks<1e-10)
				{
					double llk[9]={0.0}, max_llk=0-DBL_MAX;
					lks = 0;

					for(unsigned k1=0; k1<3; k1++)
					{
						for(unsigned k2=0; k2<3; k2++)
						{
							double mA, mB, mix_cov[4];
							mA = meanA[j*4+k1]*alpha + meanA[j*4+k2]*(1-alpha);
							mB = meanB[j*4+k1]*alpha + meanB[j*4+k2]*(1-alpha);

							for(unsigned l=0; l<4; l++)
							{
								mix_cov[l] = alpha*alpha*covs[(j*4+k1)*4+l] + (1-alpha)*(1-alpha)*covs[(j*4+k2)*4+l];
							}

							llk[k1*3+k2]=lognormpdf2(nA, nB, mA, mB, mix_cov) + log(Pr[j*4+k1]) + log(Pr[j*4+k2]);

							if (max_llk < llk[k1*3+k2])
							{
								max_llk = llk[k1*3+k2];
							}
						}
					}

					for (unsigned k=0;k<9;k++)
					{
						llk[k]  -= max_llk;
						lks += exp(llk[k]);
					}

					lks = log(lks) + max_llk;

					if (isnan(lks))
					{
						fprintf(stderr, "LLK is nan - j:  %u, x: %f %f, Pr: %f %f %f lks: %f max_llk:%f\n", j, nA, nB, Pr[j*4+0], Pr[j*4+1], Pr[j*4+2], log(lks), max_llk);
					}
				}
				else
				{
					double prev_lks = lks;
					lks = (double)log(lks);
					if (isnan(lks))
					{
						fprintf(stderr, "LLK is nan - j:  %u, x: %f %f, Pr: %f %f %f lks: %f %f \n", j, nA, nB, Pr[j*4+0], Pr[j*4+1], Pr[j*4+2], prev_lks, lks);
					}

				}
				llk+=lks;
			} // if !bZeroOut && !isnan NA&NB
		} // if >min_af
	} //for j

	return (0.-llk);
}


class cSampleStat : public cStat
{
	double *meanA, *meanB, *covs;
	unsigned int *nums;

	public:
	void Initialize(unsigned int, unsigned int, const char*);
	void Estimate();
	double GetLLK(double);
	bool IsZeroOut(unsigned int);
	void WriteStat(const char*);
	bool ReadStat(const char*);
	~cSampleStat();
	double GetCallRate(unsigned int);
	void CheckCovs(double);
};

void cSampleStat::Initialize(unsigned int n_s, unsigned int n_m, const char* sIntensity)
{
	n_marker = n_m;
	n_sample = n_s;

	Pr = new double[4*n_marker];
	bZeroOut = new bool[n_marker];
	CallRate = new double[n_sample];

	nums = new unsigned [4*n_sample];
	meanA = new double [4*n_sample];
	meanB = new double [4*n_sample];
	covs = new double [4*4*n_sample];

	for (unsigned j=0;j<n_marker;j++)
	{
		bZeroOut[j] = false;
	}

	Intensity.Initialize(sIntensity, n_s, n_m);
}

cSampleStat::~cSampleStat()
{
	delete[] Pr;
	delete[] bZeroOut;
	delete[] CallRate;

	delete[] nums;
	delete[] meanA;
	delete[] meanB;
	delete[] covs;
}

double cSampleStat::GetCallRate(unsigned int i)
{
	return CallRate[i];
}

void cSampleStat::CheckCovs(double alphastep)
{
	return;
}

bool cSampleStat::IsZeroOut(unsigned int j)
{
	return bZeroOut[j];
}

void cSampleStat::Estimate()
{
	float nA, nB, GC;
	unsigned short gType;

	double sumA[4] = {0.0}, sumB[4]= {0.0}, sumAA[4]={0.0}, sumAB[4]={0.0}, sumBB[4]={0.0};
	unsigned int* num_markers = new unsigned [4*n_marker];

	for (unsigned i=0; i<n_sample; i++)
	{
		for (unsigned j=0; j<n_marker; j++)
		{
			Intensity.ReadData(i, j, nA, nB, GC, gType);

			if ((GC==0.0 && gType==3))
			{
				bZeroOut[j] = true;
			}
			else
			{
				bZeroOut[j] = (bZeroOut[j] || false);
			}
		}
//		std::cerr << "reading sample " << i << std::endl;
	}

	unsigned int n_zeroout = 0;
	for (unsigned j=0;j<n_marker;j++)
	{
		if (bZeroOut[j]) 
		{
			n_zeroout++;
		}
	}
//	std::cerr << "n_zeroout: " << n_zeroout << std::endl;

	unsigned int n_nonzero = n_marker - n_zeroout;

	for (unsigned i=0; i<n_sample; i++)
	{
		unsigned int callNum = 0;

		for(unsigned j=0; j<n_marker; j++)
		{
			Intensity.ReadData(i, j, nA, nB, GC, gType);

			if (!bZeroOut[j] && gType<3 && !isnan(nA) && !isnan(nB))
			{
				callNum++;
			}
			else if (isnan(nA) || isnan(nB))
			{
				bZeroOut[j] = true;
			}
		}

		CallRate[i] = ((double)callNum) / ((double)n_nonzero);
		std::cerr << "Sample " << i << " Call Rate " << CallRate[i] << "\n";

		for(unsigned j=0; j<n_marker; j++)
		{
			Intensity.ReadData(i, j, nA, nB, GC, gType);

			if (!(bZeroOut[j]) && gType<3)
			{
				sumA[gType] += (double)nA;
				sumB[gType] += (double)nB;

				sumAA[gType] += (double)nA*(double)nA;
				sumBB[gType] += (double)nB*(double)nB;
				sumAB[gType] += (double)nA*(double)nB;

				nums[gType] = nums[gType]+1;
				num_markers[j*4+gType] = num_markers[j*4+gType]+1;
			}
			if (isnan(sumA[gType]) || isnan(sumB[gType]))
			{
				std::cerr << "Error: Sum is NAN though previously checked \n";
				exit (1);
			}
		}
		for(unsigned k=0; k<3; k++)
		{
			meanA[i*4+k] = sumA[k]/nums[k];
			meanB[i*4+k] = sumB[k]/nums[k];

			covs[i*4*4+k*4+0] = sumAA[k]/nums[k] - (meanA[i*4+k]*meanA[i*4+k]);
			covs[i*4*4+k*4+1] = covs[i*4*4+k*4+2] = sumAB[k]/(nums[k]) - (meanA[i*4+k]*meanB[i*4+k]);
			covs[i*4*4+k*4+3] = sumBB[k]/nums[k] - (meanB[i*4+k]*meanB[i*4+k]);

			if (isnan(meanA[i*4+k]) || isnan(meanB[i*4+k]))
			{
				std::cerr << "Error: Mean is NAN though previously checked \n";
				exit (1);
			}
			sumA[k]=0;
			sumB[k]=0;
			nums[k]=0;
			sumAA[k]=0;
			sumAB[k]=0;
			sumBB[k]=0;
		}
//		std::cerr << "Sample " << i  << " stat is estimated" << std::endl;
	} // for i < n_sample

	for(unsigned int j=0;j<n_marker;j++)
	{
		if (!bZeroOut[j])
		{
			unsigned tot_num = num_markers[j*4+0]+num_markers[j*4+1]+num_markers[j*4+2];
			double Pr_A= ((double)(num_markers[j*4+0]) + (double)(num_markers[j*4+1])*0.5) / (double)tot_num;
			double Pr_B = ((double)(num_markers[j*4+2]) + (double)(num_markers[j*4+1])*0.5) / (double)tot_num;
			Pr[j*4+0] = Pr_A*Pr_A;
			Pr[j*4+1] = Pr_A*Pr_B;
			Pr[j*4+2] = Pr_B*Pr_B;
		}
	}
	delete [] num_markers;
}

double cSampleStat::GetLLK(double alpha)
{
	double llk = 0;
	double lks = 0;
	float nA, nB, GC;
	unsigned short gType;
	unsigned int i = curSample;

	for(unsigned int j=0;j<n_marker;++j)
	{
		lks=0;
		if (GetMinPr(j) > min_af)
		{
			Intensity.ReadData(i, j, nA, nB, GC, gType);

			if (!bZeroOut[j] && !isnan(nA) && !isnan(nB))
			{
				for(unsigned k1=0; k1<3; k1++)
				{
					for(unsigned k2=0; k2<3; k2++)
					{
						double mA, mB, mix_cov[4];
						mA = meanA[i*4+k1]*alpha + meanA[i*4+k2]*(1-alpha);
						mB = meanB[i*4+k1]*alpha + meanB[i*4+k2]*(1-alpha);

						for(unsigned l=0; l<4; l++)
						{
							mix_cov[l] = alpha*alpha*covs[i*4*4+k1*4+l] + (1-alpha)*(1-alpha)*covs[i*4*4+k2*4+l];
						}

						double tmpval=normpdf2(nA, nB, mA, mB, mix_cov) * Pr[j*4+k1] * Pr[j*4+k2];
						lks+=tmpval;
					}
				}


				if (lks<1e-15)
				{
					double llk[9]={0.0}, max_llk=0-DBL_MAX;
					lks=0;

					for(unsigned k1=0; k1<3; k1++)
					{
						for(unsigned k2=0; k2<3; k2++)
						{
							double mA, mB, mix_cov[4];
							mA = meanA[i*4+k1]*alpha + meanA[i*4+k2]*(1-alpha);
							mB = meanB[i*4+k1]*alpha + meanB[i*4+k2]*(1-alpha);

							for(unsigned l=0; l<4; l++)
							{
								mix_cov[l] = alpha*alpha*covs[i*4*4+k1*4+l] + (1-alpha)*(1-alpha)*covs[i*4*4+k2*4+l];
							}
							llk[k1*3+k2]=lognormpdf2(nA, nB, mA, mB, mix_cov) + log(Pr[j*4+k1]) + log(Pr[j*4+k2]);

							if (max_llk < llk[k1*3+k2])
							{
								max_llk = llk[k1*3+k2];
							}
						}
					}
					for (unsigned k=0;k<9;k++)
					{
						llk[k]  -= max_llk;
						lks += exp(llk[k]);
					}

					lks = log(lks) + max_llk;

					if (isnan(lks))
					{
						fprintf(stderr, "LLK is nan - j:  %u, x: %f %f, Pr: %f %f %f lks: %f max_llk:%f\n", j, nA, nB, Pr[j*4+0], Pr[j*4+1], Pr[j*4+2], log(lks), max_llk);
					}
				}
				else
				{
					lks = (double)log(lks);
					if (isnan(lks))
					{
						fprintf(stderr, "LLK is nan - j:  %u, x: %f %f, Pr: %f %f %f lks: %f \n", j, nA, nB, Pr[j*4+0], Pr[j*4+1], Pr[j*4+2], log(lks));
					}

				}
				llk+=lks;
			} // if bZeroOut
		} // if (> min_af)
	} //for j
	return (0.-llk);
}

bool cSampleStat::ReadStat(const char *FileName)
{
	std::ifstream  InStatsFile(FileName, std::ios::in | std::ios::binary);
	if (!InStatsFile.is_open())
	{
		std::cout<<"Statistics file does not exist. Stats will be estimated from data." << std::endl;
		return false;
	}
	// If pre-calculated statistics exists, read
	InStatsFile.seekg(0, std::ios::beg);
	InStatsFile.read(reinterpret_cast <char *> (meanA), sizeof(double)*4*n_sample);
	InStatsFile.read(reinterpret_cast <char *> (meanB), sizeof(double)*4*n_sample);
	InStatsFile.read(reinterpret_cast <char *> (covs), sizeof(double)*4*4*n_sample);
	InStatsFile.read(reinterpret_cast <char *> (Pr), sizeof(double)*4*n_marker);
	InStatsFile.read(reinterpret_cast <char *> (nums), sizeof(unsigned)*4*n_sample);
	InStatsFile.read(reinterpret_cast <char *> (bZeroOut), sizeof(bool)*n_marker);
	InStatsFile.read(reinterpret_cast <char *> (CallRate), sizeof(double)*n_sample);

	InStatsFile.close();
	return true;
}

void cSampleStat::WriteStat(const char *FileName)
{
	std::cout<<"Creating a new statistics file " << FileName << std::endl;

	std::ofstream  OutStatsFile(FileName, std::ios::out | std::ios::binary);
	if (!OutStatsFile.is_open())
	{
		std::cerr << "Error opening " << FileName << " to write.\n";
		exit (1);
	}

	OutStatsFile.write(reinterpret_cast <char *> (meanA), sizeof(double)*n_sample*4);
	OutStatsFile.write(reinterpret_cast <char *> (meanB), sizeof(double)*n_sample*4);
	OutStatsFile.write(reinterpret_cast <char *> (covs), sizeof(double)*n_sample*4*4);
	OutStatsFile.write(reinterpret_cast <char *> (Pr), sizeof(double)*4*n_marker);
	OutStatsFile.write(reinterpret_cast <char *> (nums), sizeof(unsigned)*4*n_sample);
	OutStatsFile.write(reinterpret_cast <char *> (bZeroOut), sizeof(bool)*n_marker);
	OutStatsFile.write(reinterpret_cast <char *> (CallRate), sizeof(double)*n_sample);

	OutStatsFile.close();

}

int main(int argc, char** argv) 
{
	unsigned int n_marker, n_sample; 
	std::string sInFile, sStatFile, sRefFile, sAbfFile, sLogFile;

	bool bPerSample,bVerbose;
	cStat* pStat;

	try {
		TCLAP::CmdLine cmd("Command description message", ' ', "0.1");
		TCLAP::ValueArg<std::string> argIn("i","in","Input intensity (.adpc.bin) file",true,"","string");
		TCLAP::ValueArg<std::string> argStat("s","stat","Statistics file (created if not exist)",false,"","string");
		TCLAP::ValueArg<std::string> argAbf("b","abf","Allele frequency file (ABF), which is a plain text file with SNP_ID and Allele_B frequency",false,"","string");
		TCLAP::ValueArg<int> argNsample("n","number","Number of samples",true,0,"int");
		TCLAP::ValueArg<int> argNmarker("m","marker","Number of markers",true,196725,"int");
		TCLAP::ValueArg<float> argThreshold("t","threshold","Minimum allele frequency for likelihood estimation, default is 0.01", false,0.01,"float");
		TCLAP::SwitchArg switchPerSample("p","persample", "Do per-sample analysis, default is per-marker analysis",cmd,false);
		TCLAP::SwitchArg switchVerbose("v","verbose","Turn on verbose mode",cmd,false);

		cmd.add(argIn);
		cmd.add(argStat);
		cmd.add(argAbf);
		cmd.add(argNsample);
		cmd.add(argNmarker);
		//cmd.add(argAlphaStep);
		cmd.add(argThreshold);

		cmd.parse(argc, argv);

		sInFile = argIn.getValue();
		sStatFile = argStat.getValue();
		sAbfFile = argAbf.getValue();
		n_sample = argNsample.getValue();
		n_marker = argNmarker.getValue();
		//		alphastep = (double)argAlphaStep.getValue();
		min_af = (double)argThreshold.getValue();

		bPerSample = switchPerSample.getValue();
		bVerbose = switchVerbose.getValue();

	}
	catch (TCLAP::ArgException &e) {
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		abort();
	}

	if (bPerSample)
	{
		pStat = new cSampleStat;
	}
	else
	{
		pStat = new cMarkerStat;
	}
	pStat->Initialize(n_sample, n_marker, sInFile.c_str());

	// Try to read pre-stored stats from a file
	if (sStatFile.empty() || !pStat->ReadStat(sStatFile.c_str()))
	{
		try
		{
			pStat->Estimate();
		}
		catch (int e)
		{
			std::cerr << "Error while estimating statistics... \n";
			exit (1);
		}

		if (!sStatFile.empty())
		{
			pStat->WriteStat(sStatFile.c_str());
		}
	}
	if (!sAbfFile.empty())
	{
		pStat->ReadAbf(sAbfFile.c_str());
	}

	pStat->CheckCovs(0.01);

	std::cout << "ID\t%Mix\t\tLLK\t\tLLK0\t\n";
	std::cout << "-----------------------------------------------------------------\n";

	for(unsigned int i=0;i<n_sample;++i) {
		double alpha = 0;
		double minLLK = 1e15;
		double LLK0;
		minLLK = pStat->minimizeLLK(i, alpha, LLK0);

		std::cout << i << "\t" << alpha << "\t" << minLLK << "\t" << LLK0  << std::endl;
	}
	
	delete pStat;
	return 0;
}
