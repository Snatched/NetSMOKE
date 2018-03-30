/*-----------------------------------------------------------------------*\
|																		  |
|			 _   _      _    _____ __  __  ____  _  ________         	  |
|			| \ | |    | |  / ____|  \/  |/ __ \| |/ /  ____|        	  |
|			|  \| | ___| |_| (___ | \  / | |  | | ' /| |__   			  |
|			| . ` |/ _ \ __|\___ \| |\/| | |  | |  < |  __|  		  	  |
|			| |\  |  __/ |_ ____) | |  | | |__| | . \| |____ 		 	  |
|			|_| \_|\___|\__|_____/|_|  |_|\____/|_|\_\______|		 	  |
|                                                                         |
|   Author: Matteo Mensi <matteo.mensi@mail.polimi.it>                    |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef NETSMOKE_UNITS_PHASESEPARATOR_H
#define NETSMOKE_UNITS_PHASESEPARATOR_H

// A class to evaluate phase splitters as units inside an RNM

#include "Units.h"

namespace NetSMOKE
{
	class Units_PhaseSeparator : public Units
	{
		public:

			// CONSTRUCTOR (NOTICE THAT A THERMODYNAMICS FOR SOLIDS HAVE TO BE PROVIDED)
			Units_PhaseSeparator(OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN* thermo, NetSMOKE::UnitInfo unit_data);
	
			// SPLITTING OPERATION
			int Solve(std::vector<NetSMOKE::StreamInfo> &streams_data_structure);

			/* Residual obtainer - RETURNS ERROR IF YOU CALL RESIDUALS FOR TRIVIAL UNITS */
			int GetResiduals(OpenSMOKE::OpenSMOKEVectorDouble &residuals, std::vector<NetSMOKE::StreamInfo> &streams_data_structure);

			/* Printer - RETURNS ERROR */
			void PrintStatus(boost::filesystem::path output_folder, std::vector<NetSMOKE::StreamInfo> &streams_data_structure);

		private:
			int n_outlets_;
			std::vector<std::string> splitted_phase_;

			// Provided informations
			OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN& thermodynamicsMapSOLID;

			// local copy of inlet stream
			NetSMOKE::StreamInfo StreamIn;
			// local working storage for outlet streams
			std::vector<NetSMOKE::StreamInfo> StreamOut;
		
	};
}
#include "Units_PhaseSeparator.hpp"

#endif /* NETSMOKE_UNITS_PHASESEPARATOR_H*/
