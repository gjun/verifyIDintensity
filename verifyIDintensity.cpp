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

#include <string.h>
#include <iostream>

#include "tclap/CmdLine.h"
#include "tclap/Arg.h"

#include "common.h"
#include "MarkerStat.h"
#include "SampleStat.h"

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
