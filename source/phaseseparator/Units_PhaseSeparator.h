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

#ifndef NODUSSMOKE_UNITS_PHASESEPARATOR_H
#define NODUSSMOKE_UNITS_PHASESEPARATOR_H

// A class to evaluate phase splitters as units inside an RNM

#include "Units.h"

namespace NodusSMOKE
{
	class Units_PhaseSeparator : public Units
	{
		public:

			// CONSTRUCTOR (NOTICE THAT A THERMODYNAMICS FOR SOLIDS HAVE TO BE PROVIDED)
			Units_PhaseSeparator(OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN* thermo, NodusSMOKE::UnitInfo unit_data);
	
			// SPLITTING OPERATION
			int Solve(std::vector<NodusSMOKE::StreamInfo> &streams_data_structure);

			/* Residual obtainer - RETURNS ERROR IF YOU CALL RESIDUALS FOR TRIVIAL UNITS */
			int GetResiduals(OpenSMOKE::OpenSMOKEVectorDouble &residuals, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure);

			/* Printer - RETURNS ERROR */
			void PrintStatus(boost::filesystem::path output_folder, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure);

		private:
			int n_outlets_;
			std::vector<std::string> splitted_phase_;

			// Provided informations
			OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN& thermodynamicsMapSOLID;

			// local copy of inlet stream
			NodusSMOKE::StreamInfo StreamIn;
			// local working storage for outlet streams
			std::vector<NodusSMOKE::StreamInfo> StreamOut;
		
	};
}
#include "Units_PhaseSeparator.hpp"

#endif /* NODUSSMOKE_UNITS_PHASESEPARATOR_H*/
