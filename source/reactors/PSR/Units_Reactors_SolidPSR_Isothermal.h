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

#ifndef NETSMOKE_UNITS_REACTORS_SOLIDPSR_ISOTHERMAL_H
#define NETSMOKE_UNITS_REACTORS_SOLIDPSR_ISOTHERMAL_H

// Parent class
#include "Units_Reactors_PSR.h"
// Solid psr
#include "solidreactors/solidpsr/SolidPerfectlyStirredReactor_Isothermal_ConstantPressure.h"

namespace NetSMOKE
{

	class Units_Reactors_SolidPSR_Isothermal : public Units_Reactors_PSR {

	public:
		/* Default constructor */
		Units_Reactors_SolidPSR_Isothermal(OpenSMOKE::ThermodynamicsMap_CHEMKIN *thermodynamicsMapGasXML, OpenSMOKE::KineticsMap_CHEMKIN *kineticsMapGasXML, OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN *thermodynamicsMapSolidXML, OpenSMOKE::KineticsMap_Solid_CHEMKIN *kineticsMapSolidXML, double rho_solid, NetSMOKE::UnitInfo unit_data);

		/* Default destructor */
		~Units_Reactors_SolidPSR_Isothermal();

		/* Solve */
		int Solve(std::vector<NetSMOKE::StreamInfo> &streams_data_structure);

		/* Solve in the sequential solver */
		int NonIterativeSolve(std::vector<NetSMOKE::StreamInfo> &streams_data_structure);

		/* Get Residuals */
		int GetResiduals(OpenSMOKE::OpenSMOKEVectorDouble &residuals, std::vector<NetSMOKE::StreamInfo> &streams_data_structure);

		/* Print reactor status */
		void PrintStatus(boost::filesystem::path output_folder, std::vector<NetSMOKE::StreamInfo> &streams_data_structure);
		
		/* RTD - huge fuckery involved */
		int RTD(OpenSMOKE::OpenSMOKEVectorDouble &residuals, const double t, std::vector<NetSMOKE::StreamInfo> &streams_data_structure);
		
	protected:

		OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMapGas;
		OpenSMOKE::KineticsMap_CHEMKIN& kineticsMapGas;
		OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN& thermodynamicsMapSolid;
		OpenSMOKE::KineticsMap_Solid_CHEMKIN& kineticsMapSolid;

		// Static info for solver
		OpenSMOKE::SolidPerfectlyStirredReactor_Options* spsr_options;
		OpenSMOKE::ODE_Parameters* ode_parameters;

		// Solid density
		double rho_solid_;

		// masses and volumes
		double mass_gas;
		double mass_solid;
		double Vgas;
		double Vsolid;

		// mass flows
		double GasFlowIn;
		double SolidFlowIn;
		double GasFlowOut;
		double SolidFlowOut;

		// phase fractions
		double Fraction_Gas_In;
		double Fraction_Solid_In;
		double Fraction_Gas_Out;
		double Fraction_Solid_Out;

		// Production rate storage
		OpenSMOKE::OpenSMOKEVectorDouble Rgas_from_solid;
		OpenSMOKE::OpenSMOKEVectorDouble Rsolid;

		// Auxiliary vectors
		OpenSMOKE::OpenSMOKEVectorDouble Omega_previous_solid;
		OpenSMOKE::OpenSMOKEVectorDouble temp_omega_solid_out;


	};

}

#include "Units_Reactors_SolidPSR_Isothermal.hpp"

#endif /* NETSMOKE_UNITS_REACTORS_SOLIDPSR_ISOTHERMAL_H */
