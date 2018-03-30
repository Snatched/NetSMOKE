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

#ifndef NETSMOKE_UNITS_MIXER_ADIABATIC_H
#define NETSMOKE_UNITS_MIXER_ADIABATIC_H

// A class to evaluate adiabatic mixers in a RNM

#include "Units.h"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

namespace NetSMOKE
{

	class Units_Mixer_Adiabatic : public Units
	{
	public: // FUNCTIONS

		// CONSTRUCTOR
		Units_Mixer_Adiabatic(OpenSMOKE::ThermodynamicsMap_CHEMKIN *thermodynamicsMapXML, NetSMOKE::UnitInfo unit_data);

		// MIXING OPERATION
		int Solve(std::vector<NetSMOKE::StreamInfo> &streams_data_structure);

		int GetResiduals(OpenSMOKE::OpenSMOKEVectorDouble &residuals, std::vector<NetSMOKE::StreamInfo> &streams_data_structure);

		void PrintStatus(boost::filesystem::path output_folder, std::vector<NetSMOKE::StreamInfo> &streams_data_structure);

	private: // INTERNAL VARIABLES

		// For species balance
		OpenSMOKE::OpenSMOKEVectorDouble InletComponentMassFlow;
		OpenSMOKE::OpenSMOKEVectorDouble omega_final;
		OpenSMOKE::OpenSMOKEVectorDouble mole_fraction;
		double mix_molecular_weight;
		// Energy balances
		double T_guess;
		double T_final;
		double H;
		double Cpstream;
		OpenSMOKE::OpenSMOKEVectorDouble HInlet;
		OpenSMOKE::OpenSMOKEVectorDouble xInlet;
		OpenSMOKE::OpenSMOKEVectorDouble TFirstGuess_Num;
		OpenSMOKE::OpenSMOKEVectorDouble TFirstGuess_Denum;

		// Indexes
		int n_species;
		int n_inlets;

		// Inlet structures
		OpenSMOKE::ThermodynamicsMap_CHEMKIN &thermodynamicsMap;
		NetSMOKE::StreamInfo StreamOut;

	};

}

#include "Units_Mixer_Adiabatic.hpp"

#endif /* NETSMOKE_UNITS_MIXER_ADIABATIC_H */
