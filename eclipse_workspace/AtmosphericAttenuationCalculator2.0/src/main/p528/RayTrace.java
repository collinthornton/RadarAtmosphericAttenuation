package main.p528;

/** Rec. ITU-R P.528-4 Annex II Section V
 * 
 * @author Collin Thornton
 *
 */
public class RayTrace {
	public static class RayTraceData {
		public final double d_r;			// arc length to Earth 								(km)
		public final double theta_r; 		// incident angle of grazing ray at terminal 		(rad)
		
		RayTraceData(double d_r, double theta_r) {
			this.d_r = d_r;
			this.theta_r = theta_r;
		}
	}
	
	// Height of atmospheric layers (km)
	private static final double[] atmos_layers = { 
			0.000, 0.010, 0.020, 0.050, 
			0.100, 0.200 ,0.305, 0.500, 0.700, 
			1.000, 1.524, 2.000, 3.048, 5.000, 7.000, 
			10.00, 20.00, 30.48, 50.00, 70.00, 90.00, 
			110.0, 225.0, 350.0, 475.0 };
	
	
	/** Compute the arc length to smooth Earth horizon and incident angle of the grazing ray at terminal
	 * 
	 * @param h_r Real terminal height above ground (km)
	 */
	public static RayTraceData compute(double h_r) {
		
		if(h_r == 0.0) return new RayTraceData(0.0, 0.0);

		// <<Step 1>> Compute the scale factor
		double deltaN = -7.32 * Math.exp(0.005577*Constants.N_S);
		
		
		// <<Step 2>>
		double C_e = Math.log(Constants.N_S/(Constants.N_S+deltaN));
		
		
		// <<Step 3>> Ray trace through the atmosphere
		int i = 0;
		double theta_low = 0.0;
		double theta_high = theta_low;
		
		double r_high, N_high, n_high;
		double r_low = Constants.A_0 + atmos_layers[0];
		double N_low = Constants.N_S;// * Math.exp(-C_e*atmos_layers[i]);
		double n_low = 1.0 + (N_low*1.0e-6);
		double tau = 0.0;	// atmospheric bending per layer	(rad)

		while(atmos_layers[i] < h_r && i < atmos_layers.length-1) {
			if(atmos_layers[i+1] > h_r) {
				r_high = Constants.A_0 + h_r;
				N_high = Constants.N_S * Math.exp(-C_e*h_r);				
			} else {
				r_high = Constants.A_0 + atmos_layers[i+1];
				N_high = Constants.N_S * Math.exp(-C_e*atmos_layers[i+1]);
			}
			n_high = 1.0 + (N_high*1.0e-6);	

			theta_high = Math.acos((r_low/r_high)*(n_low/n_high)*Math.cos(theta_low));
			
			double A = (Math.log(n_high/n_low)/Math.log(r_high/r_low));
			tau += (theta_high - theta_low)*-A/(A+1.0);
			
			
			theta_low = theta_high;
			N_low = N_high;
			n_low = n_high;

			++i;
		}
		
		
		// <<Step 4>> Handle ray leaving atmosphere if applicable
		if(h_r > 475) {
			// Account for final incident angle
			theta_high = Math.acos(((Constants.A_0+475.0)*n_low)/(Constants.A_0+h_r)*Math.cos(theta_low));
			
			// Account for final bending angle
			double A = Math.log(1/n_low)/Math.log(Constants.A_0+h_r/Constants.A_0+475.0);
			tau += (theta_high - theta_low) * -A/(A+1.0);
		}
		
		
		// <<Step 5>> Compute incident angle
		double theta_r = theta_high;
		
		
		// <<Step 6>> Compute total bending angle (tau)
		
		
		
		// <<Step 7>> Compute arc distance
		double d_r = (theta_r + tau)*Constants.A_0;
		RayTraceData output = new RayTraceData(d_r, theta_r);
		
		return output;
	}
	
	
	/** Compute the arc length to smooth Earth horizon and incident angle of the grazing ray at terminal
	 * 
	 * @param h_r Real terminal height above ground (km)
	 * @param N_s Surface refractivity of Earth		(N-Units)
	 */
	public static RayTraceData compute(double h_r, double N_s) {
		double deltaN = -7.32 * Math.exp(0.005577*N_s);
		double C_e = Math.log(N_s/(N_s+deltaN));
		
		int i = 0;
		double theta_low = 0.0;
		double theta_high = theta_low;
		
		double r_low=0.0, N_low=0.0, n_low=0.0, r_high=0.0, N_high=0.0, n_high=0.0;
		double[] tau = new double[24];	// atmospheric bending per layer	(rad)
		
		while(atmos_layers[i] < h_r && i < atmos_layers.length-1) {
			r_low = Constants.A_0 + atmos_layers[i];
			N_low = N_s * Math.exp(-1*C_e*atmos_layers[i]);
			n_low = 1.0 + (N_low*Math.pow(10.0, -6));
			
			if(h_r > atmos_layers[i+1]) {
				r_high = Constants.A_0 + atmos_layers[i+1];
				N_high = N_s * Math.exp(-C_e*atmos_layers[i+1]);
			} else {
				r_high = Constants.A_0 + h_r;
				N_high = N_s * Math.exp(-C_e*h_r);
			}
			n_high = 1.0 + (N_high*Math.pow(10.0, -6));	
			
			theta_high = Math.acos((r_low*n_low)/(r_high*n_high)*Math.cos(theta_low));
			
			double A = (Math.log(n_high/n_low)/Math.log(r_high/r_low));
			tau[i] = (theta_high - theta_low)*(-1.0*A/(A+1.0));
			theta_low = theta_high;

			++i;
		}
		
		if(h_r > atmos_layers[atmos_layers.length-1]) {
			theta_high = Math.acos(((Constants.A_0+475.0)*n_low)/(Constants.A_0+h_r)*Math.cos(theta_low));
		}
		
		double tau_sum = 0.0;
		for(int j=0; j<i; ++j) tau_sum += tau[i];
		
		double theta_r = theta_high;
		double d_r = (theta_r + tau_sum)*Constants.A_0;
		RayTraceData output = new RayTraceData(d_r, theta_r);
		
		return output;
	}	
}
