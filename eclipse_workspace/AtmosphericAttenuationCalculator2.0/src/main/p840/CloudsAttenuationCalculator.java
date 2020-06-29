package main.p840;

/** Calculate loss in db/km due to clouds and fog
 * 
 * @author Collin Thornton
 * @note Based on Rec. ITU-R P.840-8, sections 1-2
 */
public class CloudsAttenuationCalculator {
	
	/** Calculate loss in db/km due to clouds and fog
	 * 
	 * @param f Frequency. (GHz) 0 <= f <= 200
	 * @param T Temperature of cloud (F)
	 * @param M Density of water vapor in cloud (g/m^3)
	 * @return db/km due to clouds. Doubled for 2-way path loss
	 * @note This function is based on Rec. ITU-R P.840-8, sections 1-2
	 */
	public static double calculate(double f, double T, double M) {
		if(f < 0) throw new IllegalArgumentException("Cloud atten. freq must be between 0 and 200 GHz");
		if(f > 200) throw new IllegalArgumentException("Cloud atten. freq must be between 0 and 200 GHz");
	
		T = (5/9)*(T-32)+273.15; // Convert to Kelvin
		
		double theta = 300/T;
		double e_2 = 3.52;
		double e_0 = 77.66 + 103.3*(theta-1);
		double e_1 = 0.0671*e_0;
		
		double f_p = 20.20 - 146*(theta-1) + 316*Math.pow(theta-1, 2);
		double f_s = 39.8*f_p;
		
		double e_prime = (e_0-e_1)/(1+Math.pow(f/f_p, 2)) + (e_1-e_2)/(1+Math.pow(f/f_s, 2)) + e_2;
		double e_pprime = (f*(e_0-e_1))/(f_p*(1+Math.pow(f/f_p, 2))) + (f*(e_1-e_2))/(f_s*(1+Math.pow(f/f_s, 2)));
		double nu = (2+e_prime)/e_pprime;
		
		double K_i = (0.819*f)/(e_pprime*(1+Math.pow(nu, 2)));
		return K_i*M;		
	}
}
