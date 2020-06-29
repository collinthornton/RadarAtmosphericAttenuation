package main.p528;


/** Rec. ITU-R P.528-4 Annex II Section XII
 * 
 * @author Collin Thornton
 *
 */
public class EffectiveRayLength {
	public static double compute(double a, double T_e, RayOptics.RayOpticsData optics) {
		
		// <<Step 1>>
		double alpha = (Math.PI/2) + optics.theta_h1;
		double z_t = a + T_e;
		
		// <<Step 2>>
		if(optics.z_2 <= z_t) return optics.r_0;
		
		// <<Step 3>>
		if(z_t <= optics.z_1) {
			double z_c = optics.z_1 * Math.sin(alpha);
			if(z_t <= z_c) return 0;
			return 2*z_t*Math.sin(Math.acos(z_c/z_t));
		}
		
		// <<Step 4>>
		double A_q = Math.asin(optics.z_1 * Math.sin(alpha)/z_t);
		double A_e = Math.PI - alpha - A_q;
		
		if(A_e == 0) return z_t - optics.z_1;
		return (optics.z_1*Math.sin(A_e))/(Math.sin(A_q));
	}
}
