/*-----------------------------------------------------------------------*\
|																		  |
|	_   _           _            _____ __  __  ____  _  ________  		  |
|	| \ | |         | |          / ____|  \/  |/ __ \| |/ /  ____| 		  |
|	|  \| | ___   __| |_   _ ___| (___ | \  / | |  | | ' /| |__    		  |
|	| . ` |/ _ \ / _` | | | / __|\___ \| |\/| | |  | |  < |  __|   		  |
|	| |\  | (_) | (_| | |_| \__ \____) | |  | | |__| | . \| |____  		  |
|	|_| \_|\___/ \__,_|\__,_|___/_____/|_|  |_|\____/|_|\_\______|		  |                                                              |
|                                                                         |
|   Author: Matteo Mensi <matteo.mensi@mail.polimi.it>                    |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef NODUSSMOKE_RNM_UNITS_SPLITTER_H
#define NODUSSMOKE_RNM_UNITS_SPLITTER_H

// A class to evaluate splitters as units inside an RNM

#include "Units.h"

namespace NodusSMOKE
{
	class Units_Splitter : public Units
	{
	public:

		// CONSTRUCTOR (REQUIRES INLET STREAM, SPLIT RATIO, )
		Units_Splitter(UnitInfo unit_data);

		// SPLITTING OPERATION
		int Solve(std::vector<NodusSMOKE::StreamInfo> &streams_data_structure);

		/* Residual obtainer - RETURNS ERROR IF YOU CALL RESIDUALS FOR TRIVIAL UNITS */
		int GetResiduals(OpenSMOKE::OpenSMOKEVectorDouble &residuals, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure);

		/* Status printer - RETURNS ERROR */
		void PrintStatus(boost::filesystem::path output_folder, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure);
		
	private:
		unsigned int n_outlets_;
		std::vector<double> split_factor_;
	};

}

#include "splitters/Units_Splitter.hpp"

#endif /* NODUSSMOKE_RNM_UNITS_SPLITTER_H*/
