/*-----------------------------------------------------------------------*\
|																		  |
|	 _   _           _            _____ __  __  ____  _  ________  		  |
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

#ifndef NODUSSMOKE_UNITS_REACTORS_PFR_NONISOTHERMAL_H
#define NODUSSMOKE_UNITS_REACTORS_PFR_NONISOTHERMAL_H

// Parent class
#include "Units_Reactors_PFR.h"

namespace NodusSMOKE
{

	class Units_Reactors_PFR_NonIsothermal : public Units_Reactors_PFR {

	public:
		/* Default constructor */
		Units_Reactors_PFR_NonIsothermal(OpenSMOKE::ThermodynamicsMap_CHEMKIN *thermodynamicsMapXML, OpenSMOKE::KineticsMap_CHEMKIN *kineticsMapXML, NodusSMOKE::UnitInfo unit_data);

		/* Default deconstructor */
		~Units_Reactors_PFR_NonIsothermal();
		
		/* Solve */
		int Solve(std::vector<NodusSMOKE::StreamInfo> &streams_data_structure);

		/* Get Residuals */
		int GetResiduals(OpenSMOKE::OpenSMOKEVectorDouble &residuals, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure);

		/* Print reactor status */
		void PrintStatus(boost::filesystem::path output_folder, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure);
		
		/* Solve in the sequential method */
		int NonIterativeSolve(std::vector<NodusSMOKE::StreamInfo> &streams_data_structure);

		/* RTD aware solver */
		int RTD(OpenSMOKE::OpenSMOKEVectorDouble &residuals, const double t, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure);
	
	protected:
		// Maps
		OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap;
		OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap;
		
		// Heat exchange
		double UA_;

		//stuff for RTD
		NodusSMOKE::StreamInfo Previous_StreamIn;
		bool inert_in;
		bool inert_out;
		double t_start;
		int n_calls;
		
	};

}

#include "Units_Reactors_PFR_NonIsothermal.hpp"

#endif /* NODUSSMOKE_UNITS_REACTORS_PFR_NONISOTHERMAL_H */
