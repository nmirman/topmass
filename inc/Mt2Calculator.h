/*************************************************************************************************
MT2/MCT General Purpose Calculator

Usage:
	1)	#include "Mt2Calculator.h"

	2)	Declare instance of calculator with Mt2Calculator::Calculator my_calc_name;

	3)	If calculating MT2, you can set test masses for neutrino (210 subsystems) and W
		(220 and 221 subsystems) using my_calc_name.SetNeutrinoTestMass(mass) and
		my_calc_name.SetWTestMass(mass). These only need to be done once, not for each
		event.

	4)	For each event, you need to set the particles using
		my_calc_name.SetParticles(bjet1, bjet2, lepton1, lepton2, met);. Inputs are
		TLorentzVectors.

	5)	Get the desired MT2/MCT values using (p and c below refer to Konstantin's p and
		c):
		my_calc_name.GetMt2(p,c);
		my_calc_name.GetMt2Perp(p,c);
		my_calc_name.GetMct(p,c);
		my_calc_name.GetMctPerp(p,c);
		my_calc_name.GetBlInvariantMasses();

		The above all return doubles except GetBlInvariantMasses(), which returns a vector
		of doubles. Note that for the 220 subsytem and invariant mass, the values returned
		are chosen by combinatorics. If you want all the values to be returned, use the
		following:
		my_calc_name.GetAllMt2_220();
		my_calc_name.GetAllMt2Perp_220();
		my_calc_name.GetAllMct_220();
		my_calc_name.GetAllMctPerp_220();
		my_calc_name.GetAllBlInvariantMasses();

		These functions all return size two vectors of doubles except
		.GetAllBlInvariantMasses(), which returns two pairs of masses (2x2 vector).

	6)	Use my_calc_name.GetMt2Perp(p, c, phi) to calculate Mt2 Perp with phi set as the
		upstream direction
*************************************************************************************************/
#ifndef MT2_CALCULATOR_RES_H
#define MT2_CALCULATOR_RES_H
#include <vector>

#include "TLorentzVector.h"
#include "TMatrix.h"

typedef TMatrixT<double> TMatrixD;

namespace Mt2Calculator{
	class Calculator{
		public:
			Calculator();
			void SetNeutrinoTestMass(const double);
			void SetWTestMass(const double);
			void SetParticles(const TLorentzVector&, const TLorentzVector&,
					const TLorentzVector&, const TLorentzVector&,
					const TLorentzVector&);
			void SetUncorrectedJets(const TLorentzVector&, const TLorentzVector&);
			void SetResolutions(const double[], const double[], const double[],
					const double[], const double[], const double[],
					const TMatrixD&);//b1, b1u, b2, b2u, l1, l2, met

			double GetMt2(const int, const int);
			double GetMt2Perp(const int, const int);
			double GetMct(const int, const int);
			double GetMctPerp(const int, const int);

                        double GetMt2Grid(const int, const int, const int);
                        double GetMt2GridReturn(const int, const int, const int,
                                        double&, double&, double&, double& );
                        std::vector<double> ReturnNeutrinoPt();
                        double CalculateDelta(const TLorentzVector, const double, const double);
                        double Calculateb(const int, const double, const double, const double, const double);
                        double CalcTransverseMass(const TLorentzVector, const TLorentzVector, const long double);
                
                        std::vector<double> MaosReturn220(TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector&, double&, double&, double&, double&, double&, double&, double&, double&);
                        double MaosReturn210( TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector&, double&, double&, double&, double&, double&, double&, double&, double&);
                        
			std::vector<double> GetAllMt2_220();
			std::vector<double> GetAllMt2Perp_220();
			std::vector<double> GetAllMct_220();
			std::vector<double> GetAllMctPerp_220();

			std::vector<double> GetBlInvariantMasses();
			std::vector<double> GetBMasses();
			std::vector< std::vector<double> > GetAllBlInvariantMasses();

			double GetUpstream210PhiResolution();
			double GetUpstreamPhiResolution();
			TMatrixD GetObjCovarianceMatrix(TLorentzVector&, double, double);
			double GetMt2Perp221Resolution();
			double GetMt2Perp220Resolution();
			double GetMctPerp221Resolution();
			std::vector<double> GetBlInvariantMassResolutions();

			double GetMt2Perp(const int, const int, const double);
			double GetMctPerp(const int, const int, const double);

			TLorentzVector jet1, jet2, lep1, lep2, met;
			TLorentzVector jet1Uncor, jet2Uncor;
			long double neutrino_test_mass, W_test_mass;
			long double jet1Res[3], jet2Res[3], jet1UncorRes[3], jet2UncorRes[3];
			long double lep1Res[3], lep2Res[3];
			long double metRes[2][2];
			TMatrixD metResM;

		private:
			TLorentzVector CalcPerp(const TLorentzVector&, const TLorentzVector&);
			TLorentzVector CalcPerp(const TLorentzVector&, const long double);
			void MatrixMultiply(const long double[][2], const long double[][2],
					long double[][2]);

			double CalcMt2Perp(const TLorentzVector&, const TLorentzVector&,
					const TLorentzVector&, const long double);
			double CalcMt2Perp(const TLorentzVector&, const TLorentzVector&,
					const TLorentzVector&, const long double,
					const long double);
			double CalcMt2Numerically(const TLorentzVector&, const TLorentzVector&,
					const TLorentzVector&, const long double);
			double CalcMt2Analytically(const TLorentzVector&, const TLorentzVector&,
					const TLorentzVector&, const long double);

			double CalcMctPerp(const TLorentzVector&, const TLorentzVector&,
					const TLorentzVector&);
			double CalcMctPerp(const TLorentzVector&, const TLorentzVector&,
					const long double);
			double CalcMct(const TLorentzVector&, const TLorentzVector&);

                        double CalcMt2Grid(const TLorentzVector, const TLorentzVector, const TLorentzVector, const long double);       
                        double CalcMt2GridZooming(const TLorentzVector, const TLorentzVector, const TLorentzVector, const long double);
                        
                        double child1_Px_min, child1_Py_min, child2_Px_min, child2_Py_min;
                        
			TLorentzVector W1a();
			TLorentzVector W2a();
			TLorentzVector W1b();
			TLorentzVector W2b();

			double ResIM(const TLorentzVector&, const TLorentzVector&,
					const long double[]);

			long double PosSqrt(const long double);
			double PosSqrt(const double);
	};
}
#endif
