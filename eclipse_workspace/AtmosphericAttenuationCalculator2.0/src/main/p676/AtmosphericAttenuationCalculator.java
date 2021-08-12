package main.p676;

import org.apache.commons.math3.util.FastMath;

public class AtmosphericAttenuationCalculator {

		private static final double tempF = 60;
		private static final double t = (5.0 / 9.0) * (tempF - 32);
		private static final double w = 7.5;
		private static final double p = 1013.25 + w * (t + 273.15) / 216.7;
		private static final double r_p = p / 1013.25;
		private static final double r_t = 288 / (273 + t);
		private static final double z_1 = getZ1();
		private static final double z_2 = getZ2();
		private static final double z_3 = getZ3();
		private static final double n_1 = getN1();
		private static final double n_2 = getN2();
		
		
		public static void main(String[] args) {
			double f = 0.5;
			double angle = Math.toRadians(10);
			double range_nmi = 100.0 / 1.852;
			double altitude = Math.sqrt(6370.0*6370.0 + range_nmi*1.852*range_nmi*1.852- 2*6370.0*range_nmi*1.852*FastMath.cos(angle+FastMath.PI/2.0)) - 6370.0;;
			
			long sumTime = 0;
			int i=0;
			
			for(i=0; i<250000; ++i) {
				long beginTime = System.nanoTime();
				double loss = calcTotalAttenuation(f, angle, range_nmi, altitude);
				sumTime += System.nanoTime() - beginTime;
			}
			
			double loss = calcTotalAttenuation(f, angle, range_nmi, altitude);
			System.out.println(loss);
			System.out.println(((double)sumTime / (double)i)/1.0e6);
		}
		
		public static double calcTotalAttenuation(double f, double angle, double range_nmi, double altitude) {
			double a_o = calcOxygenAttenuation(f);
			double h_o = calcOxygenH(f);
			double a_w = calcWaterVaporAttenuation(f);
			double h_w = calcWaterH(f);
			
			altitude = altitude * 1.852;
			
			if(angle > Math.toRadians(5)) {
				return 2 * (a_o * calcForAltGTE5(h_o, altitude) + a_w*calcForAltGTE5(h_w, altitude)) / FastMath.sin(angle);
			}
			else if(angle > 0) {
				double lt5 = a_o * calcForAltLT5(h_o, angle, altitude) + a_w * calcForAltLT5(h_w, angle, altitude);
				double atten_at_sea_level = (a_o+a_w) * range_nmi * 1.852;
				return 2 * FastMath.min(lt5,  atten_at_sea_level);
			}
			else {
				return 2 * (a_o + a_w) * range_nmi * 1.852;
			}
		}
		
		private static double calcForAltLT5(double h_ref, double phi_1, double h_2) {
			double h_1 = 0;
			
			// effective earth radius in 4/3 model
			double earthRadiusKm = 8500.0;
			
			double phi_2 = FastMath.acos((earthRadiusKm + h_1) / (earthRadiusKm + h_2) * FastMath.cos(phi_1));
			double x_1 = FastMath.tan(phi_1) * FastMath.sqrt((earthRadiusKm + h_1) / h_ref);
			double x_2 = FastMath.tan(phi_2) * FastMath.sqrt((earthRadiusKm + h_2) / h_ref);
			
			double h = FastMath.sqrt(h_ref) * (FastMath.sqrt(earthRadiusKm + h_1) * F(x_1) * FastMath.exp(-1*h_1/h_ref) / FastMath.cos(phi_1) - 
					FastMath.sqrt(earthRadiusKm + h_2) * F(x_2) * FastMath.exp(-1*h_2/h_ref) / FastMath.cos(phi_2));
			return h;
		}
		private static double F(double x) {
			return 1.0 / (0.661 * x + 0.339 * FastMath.sqrt(x*x+5.51));
		}
		
		
		private static double calcForAltGTE5(double h_ref, double height) {
			return h_ref * (1 - FastMath.exp(-1*height/h_ref));
		}
		private static double calcOxygenAttenuation(double f) {
			double t1 = 7.2 * FastMath.pow(r_t,  2.8) / (f*f + 0.34 * r_p*r_p * FastMath.pow(r_t, 1.6));
			double t2 = 0.62 * z_3 / (FastMath.pow((54-f), 1.6*z_1) + 0.83 * z_2);
			return (t1+t2) * f* f * (r_p*r_p) * 1E-3;
		}
		private static double calcWaterVaporAttenuation(double f) {
			double t1 = m(f, 3.98, 2.23, 22.235, 9.42) * g(f, 22);
			double t2 = m(f, 11.96, 0.7, 183.31, 11.14);
			double t3 = m(f, 0.081, 6.44, 321.226, 6.29);
			double t4 = m(f, 3.66, 1.6, 325.153, 9.22);
			double t5 = 25.37 * n_1 * FastMath.exp(1.09 * (1-r_t)) / FastMath.pow(f-380, 2);
			double t6 = 17.4 * n_1 * FastMath.exp(1.46 * (1-r_t)) / FastMath.pow(f-448, 2);
			double t7 = 844.6 * n_1 * FastMath.exp(0.17 * (1-r_t)) / FastMath.pow(f-557, 2) * g(f, 557);
			double t8 = 290 * n_1 * FastMath.exp(0.41 * (1-r_t)) / FastMath.pow(f-752, 2) * g(f, 752);
			double t9 = 8.33284E4 * n_2 * FastMath.exp(0.99 * (1-r_t)) / FastMath.pow(f - 1780, 2) * g(f, 1780);
			
			
			return (t1 +t2 +t3 + t4 + t5 + t6 + t7 + t8 + t9) * f * f * FastMath.pow(r_t, 2.5) * w * 1E-4;
		}
		
		private static double calcOxygenH(double f) {
			double t1 = 4.64 / (1+0.066*FastMath.pow(r_p, -2.3)) * FastMath.exp(-1*FastMath.pow((f-59.7)/(2.87+12.4 * FastMath.exp(-7.9*r_p)), 2));
			double t2 = 0.14 * FastMath.exp(2.12 * r_p) / (FastMath.pow(f - 118.75,  2) + 0.031 * FastMath.exp(2.2*r_p));
			double t3 = 0.0114 / (1+0.14*FastMath.pow(r_p, -2.6)) * f * (-0.0247 + 0.0001*f + 1.61E-6*f*f) / (1 - 0.0169*f + 4.1E-5*f*f + 3.2E-7*f*f*f);
			return FastMath.min((6.1/(1 +0.17*FastMath.pow(r_p, -1.1))) * (1+t1+t2+t3), 10.7*FastMath.pow(r_p, 0.3));
		}
		private static double calcWaterH(double f) {
			double sw = 1.013 / (1 + FastMath.exp(-8.6 * (r_p-0.57)));
			double t1 = (1.39 * sw) / (FastMath.pow(f - 22.235, 2) + 2.56 * sw);
			double t2 = (3.37 * sw) / (FastMath.pow(f - 183.31, 2) + 4.69 * sw);
			double t3 = (1.58 * sw) / (FastMath.pow(f - 325.1, 2) + 2.89 * sw);
			return 1.66 * (1 + t1 + t2 + t3);
		}
		private static double m(double f, double a, double b, double c, double d) {
			return a * n_1 * FastMath.exp(b*(1-r_t)) / (FastMath.pow(f - c, 2) + d*n_1*n_1);
		}
		private static double getN1() {
			return 0.955 * r_p * FastMath.pow(r_t, 0.68) + 0.006 * w;
		}
		private static double getN2() {
			return 0.735 * r_p * FastMath.pow(r_t, 0.5) + 0.0353 * FastMath.pow(r_t, 4) * w;
		}
		private static double g(double f, double f_i) {
			return 1 + FastMath.pow((f - f_i)/(f+f_i), 2);
		}
		private static double getZ1() {
			return phi(r_p, r_t, 0.0717, -1.8132, 0.0156, -1.6515);
		}
		private static double getZ2() {
			return phi(r_p, r_t, 0.5146, -4.63868, -0.1921, -5.7416);
		}
		private static double getZ3() {
			return phi(r_p, r_t, 0.3414, -6.5851, 0.2130, -8.5854);
		}
		private static double phi(double r_p, double r_t, double a, double b, double c, double d) {
			return FastMath.pow(r_p, a) * FastMath.pow(r_t, b) * FastMath.exp(c * (1-r_p) + d*(1-r_t));
		}
}