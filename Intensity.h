/*
 *  Author: Goo Jun
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

#ifndef __INTENSITY_H__
#define __INTENSITY_H__

#include <iostream>
#include <fstream>
#include <stdio.h>

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

//	std::cerr << "arg_i : " << arg_i << ", pos_sample : " << pos_sample << std::endl;

	if (arg_i != pos_sample)
	{
#ifdef DEBUG
		std::cerr << "Reading intensity data of "<< arg_i << "-th sample from input file\n";
#endif
		std::streampos seek_pos = 16+(unsigned long)arg_i*(unsigned long)n_marker*18;
		inFile.clear();
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

#endif
