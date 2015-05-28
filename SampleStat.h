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

#ifndef __SAMPLESTAT_H__
#define __SAMPLESTAT_H__

#include "Stat.h"
#include "Intensity.h"

class cSampleStat : public cStat
{
	double *meanA, *meanB, *covs;

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
	double nums[4] = {0.0};
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
	}

	unsigned int n_zeroout = 0;
	for (unsigned j=0;j<n_marker;j++)
	{
		if (bZeroOut[j]) 
		{
			n_zeroout++;
		}
	}

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
	OutStatsFile.write(reinterpret_cast <char *> (bZeroOut), sizeof(bool)*n_marker);
	OutStatsFile.write(reinterpret_cast <char *> (CallRate), sizeof(double)*n_sample);

	OutStatsFile.close();

}

#endif
