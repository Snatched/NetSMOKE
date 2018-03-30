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

#ifndef NETSMOKE_ReactorNetwork_H
#define NETSMOKE_ReactorNetwork_H

// Include the Reactor Network Framework
#include "OpenSMOKEpp"  // OpenSMOKE definitions
#include "vector"  // std vector
#include "iostream" // input output stream
#include "NetSMOKE_Utilities"
#include "NetSMOKE_Units"



	//!  A class for simulating reactor networks
	/*!
	The purpose of this class is to simulate reactor networks. It provides a common interface to different
	reactor network models.
	*/

namespace NetSMOKE
{

	class ReactorNetwork

	{


	public:

		/**
		*@brief Default constructor
		*@param thermodynamicsMap Thermodynamic interpreted XML file
		*@param kineticsMap Kinetics interpreted XML file
		*@param DeviceMap vector of objects for each unit in the network
		*@param Stream structure containing the streams informations stored in a vector
		*@param InputFlowMatrix Matrix stream/units
		*@InputTripletsVector eigen sparse matrix defined as triplets for the unit/unit matrix construction
		*/
		ReactorNetwork(OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMap,
			OpenSMOKE::KineticsMap_CHEMKIN* kineticsMap,
			std::vector<NetSMOKE::Units*>* DeviceMap,
			std::vector<NetSMOKE::StreamInfo>* Stream,
			std::vector<Eigen::Triplet<int>>* InputTripletsVector);


		// COMMON FUNCTIONS //

			// Deconstructor - deletes the unit vector in the heap
			/**
			*@brief Clear the heap memory occupied by the device vector of objects
			*/
			~ReactorNetwork();
			

			// Flow Devices Processing - Process auxiliary units like mixers and splitters
			/**
			*@brief  Process auxiliary units like mixers and splitters
			*/
			int FlowDevicesProcessing();

			int FlowDevicesOrderedProcessing();
		
			//  Optimal Device Order - used to find the best order to perform auxiliary devices solving. Unnecessary but speed-up conversion.
			/**
			*@brief  Find optimal order of solving auxiliary devices
			*/
			void FindOptimalDeviceOrder();

			// Reactor processing - Used to evaluate the output of all reactors given the input
			/**
			*@brief  Used to evaluate the output of a reactor given the input
			*/
			int ReactorProcessing();
			int ReactorSequentialProcessing(); // For use in direct substitution solvers


			// Sparsity pattern - saved on the dedicated vectors rows_sparsity_pattern and cols_sparsity_pattern
			/**
			*@brief Used to analyze and get the sparsity pattern of the system
			*/
			void SparsityPattern();


			// Return Final Status - Stores the final, processed info of the Unit
			/**
			*@brief  Allows to extrapolate the results in terms of stream characteristics and solution vector to another scope, only at the final state
			*/
			void GetResults(OpenSMOKE::OpenSMOKEVectorDouble &solution, std::vector<NetSMOKE::StreamInfo> &StreamFin);

			// Print Stream On Video - Prints stream i characteristics on video using std iostream
			/**
			*@brief  Prints stream i characteristics on video using std iostream, at the state the solving is currently at
			*@param j index of the desired stream in the stream vector
			*/
			void PrintStreamOnVideo(int j);

			// Serial solver - Solves a network with no recycles sequentially
			/**
			*@brief  Sequential resoultion of units in series
			*/
			int SolveAsSeries();

			// Direct substitution method to make better 1st guess
			/**
			*@brief  Sequential resoultion of units in series
			*@param	 iter number of iterations to perform
			*/
			int SequentialSubstitution ( int iter );


			// Print Outlet Streams - Calls the console printing specifically for streams leaving the network. Uses PrintStreamOnVideo.
			/**
			*@brief  Calls the console printing specifically for streams leaving the network, at the state the solving is currently at. Usually called after solving.
			*/
			void PrintOutletStreams();

			// Print Final Reactors Status 
			/**
			*@brief  Perform OpenSMOKE++ print on file for each reactor
			*/
			void PrintFinalReactorsStatus(boost::filesystem::path output_folder);

			// Print iteration for ODE interface
			/**
			*@brief  Print function to display steps performed by the solver
			*@param y current residuals vector
			*/
			int PrintODE(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& dy);
			
			// Print iteration for NLS interface
			/**
			*@brief  Print function to display steps performed by the solver
			*@param dy current residuals vector
			*/
			int PrintNLS(const OpenSMOKE::OpenSMOKEVectorDouble& dy);

			/**
			*@brief WriteBackup saves current solution vector to a binary file
			*/
			void WriteBackup();

			/**
			*@brief InitializeFromBackUp gets first guess from binary file
			*/
			void InitializeFromBackUp(boost::filesystem::path backup_file);

			void CloseTrackFile();

			
			///////////////////////////////////////////////////////////////////////////////////////////////////

			//  Set values - used to change the values of the reactor structure by the solver
			/**
			*@brief  Assign the current solution vector to the network
			*@param itervalues current solution vector provided by the solver
			*/
			virtual void SetValues(OpenSMOKE::OpenSMOKEVectorDouble itervalues){};

			// Network equations in NLS - equations for all the species in all the Unit composing the NLS
			/**
			*@brief  Equation accessible to an NLS solver
			*@param y current solution vector
			*@param dy current residuals vector
			*/
			virtual int NLSEquations(OpenSMOKE::OpenSMOKEVectorDouble &y, OpenSMOKE::OpenSMOKEVectorDouble &dy) { return 0; };

			// Network equations in ODE- equations for all the species in all the Unit composing the ODE system
			/**
			*@brief  Equations accessible to an ODE solver
			*@param t current time
			*@param y current solution vector
			*@param dy current residuals vector
			*/
			virtual int ODEEquations(const double t, OpenSMOKE::OpenSMOKEVectorDouble &y, OpenSMOKE::OpenSMOKEVectorDouble &dy) { return 0; };

			// Network Solving NLS - Run the solver of the NLS
			/**
			*@brief  Constructs and runs the NLS solver
			*/
			virtual void SolveNLS(){};

			// Network Solving ODE - Run the solver of the ODE system
			/**
			*@brief  Construct and runs the ODE solver
			*@param tf time of integration
			*/
			virtual void SolveODE(const double tf){};

			// Print output file - generates a summary output containing all the streams and their composition
			/**
			*@brief Calls the fileprinting of the stream status
			*/
			virtual void PrintStreamSummary(std::ofstream& fOutput, const boost::filesystem::path output_file_ascii) {};

			// Test this rtd stuff omg omg
			virtual void PrintTrack(double t) {};
			virtual void SetMWsAndMoleFractions() {};
			

		protected:
			// Case switcher
			bool solid_case;
			bool use_sparse_solvers;
			
			// Indexes for searches inside relational matrixes
			unsigned int NS;
			unsigned int NSgas;
			unsigned int NSsolid;
			unsigned int Nunits;
			unsigned int Nstreams;
			unsigned int Nreactors;
			unsigned int Nreactors_HE;
			unsigned int Nreactors_Solid;
			unsigned int neq;

			OpenSMOKE::OpenSMOKEVectorDouble y0; // Initial condition/first guess
			OpenSMOKE::OpenSMOKEVectorDouble yf; // Final results 
			OpenSMOKE::OpenSMOKEVectorDouble residuals;
			
			// Kinetics and thermodynamics schemes for gas
			OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMapXMLinternal;
			OpenSMOKE::KineticsMap_CHEMKIN& kineticsMapXMLinternal;

			// Device objects
			// Device class objetcs are stored in a vector. Built reading UnitInfo data strcutures.
			std::vector<NetSMOKE::Units*> &Device;
			
			// Streams data strucutre
			std::vector<NetSMOKE::StreamInfo> &Stream;

			// Relational Matrixes
			Eigen::SparseMatrix<int> FlowMatrix; // stream/unit, compressed.
			std::vector<Eigen::Triplet<int>> &TripletsVector; //stream/unit, compressed triplets form.

			// Flow Devices optimal order
			bool ordered_devices;
			std::vector<int> FlowDevicesOptimalOrder;

			// Sequential solving vector
			std::vector<int> SolvingOrder;
			bool sequential_solver;

			// Sparsity pattern vectors
			std::vector<unsigned int> cols_sparsity_pattern;
			std::vector<unsigned int> rows_sparsity_pattern;

			// For printing
			std::vector<unsigned int> widths_of_output_species_;
			unsigned int NLSiter;
			unsigned int ODEiter;
			unsigned int n_step_video_NLS_;
			unsigned int counter_file_video_NLS_;
			unsigned int n_step_video_ODE_;
			unsigned int counter_file_video_ODE_;

			// Species tracking
			std::string tracked_species;
			bool tracking;
			std::vector<double> outlet_streams;
			std::ofstream fOutTrack;

			// Useful paths
			boost::filesystem::path working_folder;
			boost::filesystem::path backup_folder;

	};

}


#include "ReactorNetwork.hpp"

#endif /* NETSMOKE_ReactorNetwork_H */
