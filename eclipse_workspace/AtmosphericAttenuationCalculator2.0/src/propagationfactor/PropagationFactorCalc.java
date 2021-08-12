package propagationfactor;

import org.apache.commons.math3.complex.Complex;
public class PropagationFactorCalc {
	/** Calculate the Diffraction Factor and Reflection Factor. (dB)
	 * 
	 * @param EPSILON_R : Relative Dielectric Constant of a Surface
	 * @param SIGMA_E : Conductivity of Surface
	 * @param B: Polarization Parameter 0 = Horizontal Polarization, 1 = Vertical Polarization
	 * @param SIGMA_H : rms surface height deviation
	 * @param K_A : Effective Earth Radius
	 * @param LAMBDA : Wavelength of Transmission
	 * @param H_ONE : Height of Terminal One
	 * @param H_TWO : Height of Terminal Two
	 */


	 	
		//If any any point things become unclear please reference Barton's Book "Modern Radar System Analysis" TI 59 calculator 
		//These are the inputs to the calculator
		private double EPSILON_R,SIGMA_E,SIGMA_H,K_A,LAMBDA,H_ONE, H_TWO,R;
		private int B;
		//R_SUB_H exists to compare for line of sight
		private double R_H;
		//To calculate D_ONE, P and THETA must be calculated as well.
		private double D_ONE, P, THETA;
		//To calculate SIGMA_NAUGHT we require H_PRIME_SUB_ONE, H_PRIME_SUB_TWO, and PSI_NAUGHT
		private double DELTA_NAUGHT, PSI_NAUGHT, H_PRIME_ONE, H_PRIME_TWO;
		//To calculate DIFFRACTION FACTOR we need sub variables TOP, A_SUB_ONE, P_PRIME and BOTTOM to make it more legible
		private double P_PRIME, TOP, BOTTOM, A_ONE, F_S, S, G_ONE, G_TWO, H_C, H_MIN, H_L_ONE, H_L_TWO, F_D;
		//To calculate Reflection Factor
		private Complex RHO_NAUGHT, EPSILON_PRIME;
		private double RHO_S, RHO, F_R, ALPHA;
		//In Decibels
		private double reflectionFactorInDecibels, diffractionFactorInDecibels, convertedVariable;
		
		public PropagationFactorCalc(double EPSILON_R, double SIGMA_E, int B, double SIGMA_H, double K_A, double LAMBDA, double H_ONE, double H_TWO, double R){
			this.EPSILON_R = EPSILON_R;
			this.SIGMA_E = SIGMA_E;
			this.B = B;
			this.SIGMA_H = SIGMA_H;
			this.K_A = K_A;
			this.LAMBDA = LAMBDA;
			this.H_ONE = H_ONE;
			this.H_TWO = H_TWO;
			this.R = R;
			//These lines are used for testing will be removed when calculator is complete.
			}
			
		
		public void computePropagation() {
			calculateD_ONE();
			calculateDELTA_NAUGHT();
			calculateDiffractionFactor();
			if(F_D != 0) {
				diffractionFactorInDecibels = convertToDecibels(F_D);
				
			}
			calculateReflectionFactor();
			
			if(lineOfSightExists()) {
				reflectionFactorInDecibels = convertToDecibels(F_R);
			}else {
				this.reflectionFactorInDecibels = -200;
			}
		}
		public boolean lineOfSightExists() {
			R_H = Math.sqrt(2 * this.K_A) * (Math.sqrt(this.H_ONE) + Math.sqrt(this.H_TWO));
			
			if(R < R_H) {
				//line of sight exists
				return true;
			}else {
				return false;
			}
		
		}
		
		private double calculateD_ONE() {
			//calculate P
			P = Math.sqrt((((4 * K_A) * (H_ONE + H_TWO)) + Math.pow(R, 2)) / 3);
			//calculate THETA
			THETA = Math.acos(((2 * K_A) * (H_TWO - H_ONE)) * (R / Math.pow(P,3)));
			//with THETA and P; calculate D_ONE
			D_ONE = (R / 2) - (P * Math.cos((THETA + Math.PI)/3)); 
			return D_ONE;
		}
		
		
		private double calculateDELTA_NAUGHT() {
			H_PRIME_ONE = H_ONE - (Math.pow(D_ONE, 2)/(2*K_A));
			
			H_PRIME_TWO = (R - D_ONE) * (H_PRIME_ONE / D_ONE);
			PSI_NAUGHT = Math.atan((H_PRIME_ONE / D_ONE));
			//Something about no negative angles
			if(PSI_NAUGHT < 0) {
				PSI_NAUGHT = 0;
				
			}
			if(PSI_NAUGHT== 0) {
				DELTA_NAUGHT = 0;
			}else {
				DELTA_NAUGHT = (2 * H_PRIME_ONE * H_PRIME_TWO) / R;
				
			}
			return DELTA_NAUGHT;
		}
		private double calculateDiffractionFactor() {
			if(DELTA_NAUGHT >= LAMBDA/4) {
				F_D = 0;
				return F_D;
			}
			else {
				TOP = Math.sqrt(Math.pow((EPSILON_R -1) ,2) + Math.pow((60 * SIGMA_E * LAMBDA), 2));
				BOTTOM = Math.pow((Math.pow(EPSILON_R, 2) + Math.pow(60 * SIGMA_E * LAMBDA, 2)),B);
				
				P_PRIME = ((2 * Math.PI)/LAMBDA) * ((TOP)/(BOTTOM));
				
				A_ONE = 1 / (P_PRIME * R);
				
				//(A.12)
				S = (4.43 * Math.pow(10, -5)) * Math.pow(LAMBDA, -.3333);
				
				F_S = 2.507 * Math.pow((S * R), 1.5) * Math.exp((-1.607 * S * R));
				H_C = 30 * Math.pow(LAMBDA, .666);
				
				double SUB_CALCULATION_FOR_G_ONE = H_ONE/(2 * H_C);
				
				double SUB_TWO_CALCULATION_FOR_G_ONE =  0.948 * Math.sqrt(SUB_CALCULATION_FOR_G_ONE);
				
				double SUB_CALCULATION_FOR_G_TWO = H_TWO / (2 * H_C);
				
				double SUB_TWO_CALCULATION_FOR_G_TWO =  0.948 * Math.sqrt(SUB_CALCULATION_FOR_G_TWO);
				
				G_ONE = 0.1356 * Math.pow((SUB_CALCULATION_FOR_G_ONE),-0.904) * Math.pow(10,SUB_TWO_CALCULATION_FOR_G_ONE);
				G_TWO = 0.1356 * Math.pow((SUB_CALCULATION_FOR_G_TWO),-0.904) * Math.pow(10, SUB_TWO_CALCULATION_FOR_G_TWO);
				
				H_MIN = Math.sqrt((LAMBDA / (2 * Math.PI * P_PRIME)));
				
				H_L_ONE = (H_ONE * G_ONE) / H_MIN;
				
				H_L_TWO = (H_TWO * G_TWO) / H_MIN;
				
				F_D = 2 * A_ONE * F_S * H_L_ONE * H_L_TWO;
				
			}
			return F_D;
			
		}
		/*
		 * 
		 */
		private double calculateReflectionFactor() {
			//Complex Dielectric Constant EPSILON PRIME
			EPSILON_PRIME = new Complex(EPSILON_R, -60 * SIGMA_E * LAMBDA);
			//Temporary complex conversion of PSI_NAUGHT
			Complex PSI_NAUGHT_C = new Complex(PSI_NAUGHT);
			//vertical polarization
			if(B == 1) {
				Complex RHO_NAUGHT_T = (PSI_NAUGHT_C.sin().multiply(EPSILON_PRIME)).subtract((EPSILON_PRIME.subtract((PSI_NAUGHT_C.divide(EPSILON_PRIME.pow(B))).cos().pow(2))).sqrt());
				Complex RHO_NAUGHT_B = (PSI_NAUGHT_C.sin().multiply(EPSILON_PRIME)).add((EPSILON_PRIME.multiply(PSI_NAUGHT_C.divide(EPSILON_PRIME.pow(B))).sqrt()));
				RHO_NAUGHT = RHO_NAUGHT_T.divide(RHO_NAUGHT_B);
			}
			//horizontal polarization
			else {
				Complex RHO_NAUGHT_T = (PSI_NAUGHT_C.sin()).subtract((EPSILON_PRIME.subtract((PSI_NAUGHT_C.divide(EPSILON_PRIME.pow(B))).cos().pow(2))).sqrt());
				Complex RHO_NAUGHT_B = (PSI_NAUGHT_C.sin()).add((EPSILON_PRIME.subtract(PSI_NAUGHT_C.divide(EPSILON_PRIME.pow(B)).cos().pow(2)).sqrt()));
				RHO_NAUGHT = RHO_NAUGHT_T.divide(RHO_NAUGHT_B);
			}
			RHO_S = Math.sqrt(Math.exp(-(Math.pow((((4 * Math.PI * SIGMA_H) / LAMBDA)),2) * Math.pow(Math.sin(PSI_NAUGHT),2))));
			
			//RHO = Math.abs(RHO_NAUGHT) * RHO_S;
			RHO = RHO_NAUGHT.abs() * RHO_S;
			//if ANGLE < NOT LOW ANGLE;; then ALPHA = PI
			ALPHA = Math.PI;
			F_R = Math.sqrt(1 + Math.pow(RHO, 2) + 2 * RHO * Math.cos((ALPHA + 2 * Math.PI * (DELTA_NAUGHT/LAMBDA))));
			return F_R;
		}
		private double convertToDecibels(double input) {
			convertedVariable = 20 * Math.log10(input);
			return convertedVariable;
		}
		public double getDiffractionFactorInDecibels() {
			return this.diffractionFactorInDecibels;
		}
		public double getReflectionFactorInDecibels() {
			return this.reflectionFactorInDecibels;
		}
		//double EPSILON_R, double SIGMA_E, int B, double SIGMA_H, double K_A, double LAMBDA, double H_ONE, double H_TWO, double R
		public double getSigma_E() {
			return this.SIGMA_E;
		}
		public void setSigma_E(double Sigma_E) {
			this.SIGMA_E = Sigma_E;
		}
		public int getB() {
			return this.B;
		}
		public void setB(int B) {
			this.B = B;
		}
		
		public double getSigma_H() {
			return this.SIGMA_H;
		}
		public void setSigma_H(double Sigma_H) {
			this.SIGMA_H = Sigma_H;
		}
		public double getK_A() {
			return this.K_A;
		}
		public void setK_A(double K_A) {
			this.K_A = K_A;
		}
		public double getLambda() {
			return this.LAMBDA;
		}
		public void setLambda(double LAMBDA) {
			this.LAMBDA = LAMBDA;
		}
		public double getH_One() {
			return this.H_ONE;
		}
		public void setH_ONE(double H_One) {
			this.H_ONE = H_One;
		}
		public double getH_Two() {
			return this.H_TWO;
		}
		public void setH_Two(double H_Two) {
			this.H_TWO = H_Two;
		}
		public double getR() {
			return this.R;
		}
		public void setR(double R) {
			this.R = R;
		}
			


}
