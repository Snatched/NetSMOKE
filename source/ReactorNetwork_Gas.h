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

#ifndef NODUSSMOKE_REACTORNETWORK_GAS_H
#define NODUSSMOKE_REACTORNETWORK_GAS_H

// Include the Reactor Network Main Class
#include "ReactorNetwork.h"
#include "NodusSMOKE_Units"
#include "NodusSMOKE_Utilities"


	//!  A class for simulating gas phase reactor networks
	/*!
	The purpose of this class is to simulate reactor networks which have only gas phase streams.
	*/

namespace NodusSMOKE
{

	class ReactorNetwork_Gas : public ReactorNetwork 
	{


		public:

			/**
			*@brief Default constructor
			*@param thermodynamicsMap Thermodynamic interpreted XML file
			*@param kineticsMap Kinetics interpreted XML file
			*@param Stream structure containing the streams informations stored in a vector
			*@param FirstGuessY vector containing the first guess solution
			*@param InputFlowMatrix Matrix stream/units
			*@InputTripletsVector eigen sparse matrix defined as triplets for the unit/unit matrix construction
			*@input_flag_sequential_solver Boolean to skip some unnecessary step when the sequential serial solver is used
			*/
			ReactorNetwork_Gas(
				OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMap,
				OpenSMOKE::KineticsMap_CHEMKIN* kineticsMap,
				std::vector<NodusSMOKE::Units*> *DeviceMap,
				std::vector<NodusSMOKE::StreamInfo>* StreamMap,
				OpenSMOKE::OpenSMOKEVectorDouble FirstGuessY,
				std::vector<Eigen::Triplet<int>>* InputTripletsVector,
				std::vector<Eigen::Triplet<double>>* InputTripletsVector_UnitMatrix,
				bool input_flag_sequential_solver );


			// Distribute Mass Flow - Set mass flows in the network
			/**
			*@brief  Set mass flows in the network
			*/
			int DistributeMassFlows();

			//  Set values - used to change the values of the reactor structure by the solver
			/**
			*@brief  Assign the current solution vector to the network
			*@param itervalues current solution vector provided by the solver
			*/
			void SetValues(OpenSMOKE::OpenSMOKEVectorDouble itervalues);

			// Network equations in NLS - equations for all the species in all the Unit composing the NLS
			/**
			*@brief  Equation accessible to an NLS solver
			*@param y current solution vector
			*@param dy current residuals vector
			*/
			int NLSEquations(OpenSMOKE::OpenSMOKEVectorDouble &y, OpenSMOKE::OpenSMOKEVectorDouble &dy);

			// Network equations in ODE- equations for all the species in all the Unit composing the ODE system
			/**
			*@brief  Equations accessible to an ODE solver
			*@param t current time
			*@param y current solution vector
			*@param dy current residuals vector
			*/
			int ODEEquations(const double t, OpenSMOKE::OpenSMOKEVectorDouble &y, OpenSMOKE::OpenSMOKEVectorDouble &dy);

			// Network Solving NLS - Run the solver of the NLS
			/**
			*@brief  Constructs and runs the NLS solver
			*/
			void SolveNLS();

			// Network Solving ODE - Run the solver of the ODE system
			/**
			*@brief  Construct and runs the ODE solver
			*@param tf time of integration
			*/
			void SolveODE(const double tf);

			void HybridSolver(const double tf);

			// Set MWs and rho - assign to each stream a proper molecular weight, density of mix and mole fraction
			/**
			*@brief Assign to each stream a proper molecular weight, density of mix and mole fraction
			*/
			void SetMWsAndMoleFractions();

			// Print output file - generates a summary output containing all the streams and their composition
			/**
			*@brief Calls the fileprinting of the stream status
			*/
			void PrepareStreamSummary(std::ofstream& fOutput, const boost::filesystem::path output_file_ascii);
			void PrintStreamSummary(std::ostream& fOutput);

			// Test stuff rtd
			void SetTrack(std::string track_this);
			void PrintTrack(double t);
			void PrepareTrackHistory(boost::filesystem::path output_folder_track);

		private:

		    // needed to perform the massflow distribution
			std::vector<Eigen::Triplet<double>>& TripletsVector_UnitMatrix;

			unsigned int calls;

	};

}


#include "ReactorNetwork_Gas.hpp"

#endif /* NODUSSMOKE_REACTORNETWORK_GAS_H */
