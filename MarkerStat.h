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


#ifndef __MARKERSTAT_H__
#define __MARKERSTAT_H__

#include "Stat.h"
#include "Intensity.h"

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


#endif
