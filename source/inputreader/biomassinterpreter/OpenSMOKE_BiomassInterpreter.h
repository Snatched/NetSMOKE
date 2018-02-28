/*----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef OPENSMOKE_BIOMASSINTERPRETER_H
#define OPENSMOKE_BIOMASSINTERPRETER_H

// BIOMASS INTERPRETER

/* A class to evaluate CHEMKIN species composition of a biomass given its CHO elemental analysis */
#include "OpenSMOKEpp"

namespace OpenSMOKE
{
	class BiomassInterpreter {

		public:
			/* default constructor */
			BiomassInterpreter(){};

			/* Functions */
			void SetValues(double C_fraction, double H_fraction, double O_fraction, double others_fraction);
			void CheckValues();
			double LighMass() { return ((H / (1 - others_comp) - 0.054755)*(0.11028068) - (0.00211286)*(C / (1 - others_comp) - 0.5673579)) / ((0.11028068)*(0.01206341) - (0.03666363)*(0.00211286)); }
			double LigcMass() { return ((0.03666363)*(H / (1 - others_comp) - 0.054755) - (C / (1 - others_comp) - 0.5673579)*(0.01206341)) / ((0.00211286)*(0.03666363) - (0.01206341)*(0.11028068)); } // inserire formule
			double LigoMass() { return 1 - ((H / (1 - others_comp) - 0.054755)*(0.11028068) - (0.00211286)*(C / (1 - others_comp) - 0.5673579)) / ((0.11028068)*(0.01206341) - (0.03666363)*(0.00211286)) - ((0.03666363)*(H / (1 - others_comp) - 0.054755) - (C / (1 - others_comp) - 0.5673579)*(0.01206341)) / ((0.00211286)*(0.03666363) - (0.01206341)*(0.11028068)); }
			double LighMole() { return (((H / (1 - others_comp) - 0.054755)*(0.11028068) - (0.00211286)*(C / (1 - others_comp) - 0.5673579)) / ((0.11028068)*(0.01206341) - (0.03666363)*(0.00211286))) / (0.0025849*437.46026); }
			double LigcMole() { return (((0.03666363)*(H / (1 - others_comp) - 0.054755) - (C / (1 - others_comp) - 0.5673579)*(0.01206341)) / ((0.00211286)*(0.03666363) - (0.01206341)*(0.11028068))) / (0.0025849*301.31388); } // inserire formule
			double LigoMole() { return (1 - ((H / (1 - others_comp) - 0.054755)*(0.11028068) - (0.00211286)*(C / (1 - others_comp) - 0.5673579)) / ((0.11028068)*(0.01206341) - (0.03666363)*(0.00211286)) - ((0.03666363)*(H / (1 - others_comp) - 0.054755) - (C / (1 - others_comp) - 0.5673579)*(0.01206341)) / ((0.00211286)*(0.03666363) - (0.01206341)*(0.11028068))) / (0.0025849*423.39062); }

		private: 
			double C, H, O, others_comp;

	};

} // END OF NAMESPACE DECLARATION

#include "OpenSMOKE_BiomassInterpreter.hpp"

#endif /* OPENSMOKE_BIOMASSINTERPRETER_H */
