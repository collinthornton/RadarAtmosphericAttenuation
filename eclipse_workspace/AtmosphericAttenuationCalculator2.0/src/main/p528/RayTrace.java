package main.p528;

/** Rec. ITU-R P.528-4 Annex II Section V
 * 
 * @author Collin Thornon
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
	private static final double[] atmos_layers = { 0, 0.01, 0.02, 0.05 ,0.10, 0.20 ,0.305, 0.50, 
			0.70, 1.00, 1.524, 2.00, 3.048, 5.00, 7.00, 10.00, 20.00, 30.48, 50.00, 70.00, 90.00, 110.00, 225.00, 350.00, 475.00 };
	
	private static double[] tau = new double[24];	// atmospheric bending per layer	(rad)
	
	/** Compute the arc length to smooth Earth horizon and incident angle of the grazing ray at terminal
	 * 
	 * @param h_r Real terminal height above ground (km)
	 */
	public static RayTraceData compute(double h_r) {
		
		// <<Step 1>> Compute the scale factor
		double deltaN = -7.32 * Math.exp(0.005577*Constants.N_S);
		
		
		// <<Step 2>>
		double C_e = Math.log10(Constants.N_S/(Constants.N_S+deltaN));
		
		
		// <<Step 3>> Ray trace through the atmosphere
		int i = 0;
		double theta_low = 0;
		double theta_high = theta_low;
		
		double r_high, N_high, n_high;
		double r_low = Constants.A_0 + atmos_layers[i];
		double N_low = Constants.N_S * Math.exp(-C_e*atmos_layers[i]);
		double n_low = 1 + (N_low*Math.pow(10, -6));
		
		while(atmos_layers[i] < h_r && i < atmos_layers.length-1) {
			if(atmos_layers[i+1] > h_r) {
				r_high = Constants.A_0 + h_r;
				N_high = Constants.N_S * Math.exp(-C_e*h_r);				
			} else {
				r_high = Constants.A_0 + atmos_layers[i+1];
				N_high = Constants.N_S * Math.exp(-C_e*atmos_layers[i+1]);
			}
			n_high = 1 + (N_high*Math.pow(10, -6));	
			
			theta_high = Math.acos((r_low*n_low)/(r_high*n_high)*Math.cos(theta_low));
			
			double A = (Math.log10(n_high/n_low)/Math.log10(r_high/r_low));
			tau[i] = (theta_high - theta_low)*(-A/(A+1));
			
			
			theta_low = theta_high;
			N_low = N_high;
			n_low = n_high;

			++i;
		}
		
		
		// <<Step 4>> Handle ray leaving atmosphere if applicable
		if(h_r > atmos_layers[atmos_layers.length-1]) {
			theta_high = Math.acos(((Constants.A_0+475)*n_low)/(Constants.A_0+h_r)*Math.cos(theta_low));
		}
		
		
		// <<Step 5>> Compute incident angle
		double theta_r = theta_high;
		
		
		// <<Step 6>> Compute total bending angle
		double tau_sum = 0;
		for(int j=0; j<i; ++j) tau_sum += tau[i];
		
		
		// <<Step 7>> Compute arc distance
		double d_r = (theta_r + tau_sum)*Constants.A_0;
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
		double C_e = Math.log10(N_s/(N_s+deltaN));
		
		int i = 0;
		double theta_low = 0;
		double theta_high = theta_low;
		
		double r_low=0, N_low=0, n_low=0, r_high=0, N_high=0, n_high=0;
		
		while(atmos_layers[i] < h_r && i < atmos_layers.length-1) {
			r_low = Constants.A_0 + atmos_layers[i];
			N_low = N_s * Math.exp(-1*C_e*atmos_layers[i]);
			n_low = 1 + (N_low*Math.pow(10, -6));
			
			if(h_r > atmos_layers[i+1]) {
				r_high = Constants.A_0 + atmos_layers[i+1];
				N_high = N_s * Math.exp(-1*C_e*atmos_layers[i+1]);
			} else {
				r_high = Constants.A_0 + h_r;
				N_high = N_s * Math.exp(-1*C_e*h_r);
			}
			n_high = 1 + (N_high*Math.pow(10, -6));	
			
			theta_high = Math.acos((r_low*n_low)/(r_high*n_high)*Math.cos(theta_low));
			
			double A = (Math.log10(n_high/n_low)/Math.log10(r_high/r_low));
			tau[i] = (theta_high - theta_low)*(-1*A/(A+1));
			theta_low = theta_high;

			++i;
		}
		
		if(h_r > atmos_layers[atmos_layers.length-1]) {
			theta_high = Math.acos(((Constants.A_0+475)*n_low)/(Constants.A_0+h_r)*Math.cos(theta_low));
		}
		
		double tau_sum = 0;
		for(int j=0; j<i; ++j) tau_sum += tau[i];
		
		double theta_r = theta_high;
		double d_r = (theta_r + tau_sum)*Constants.A_0;
		RayTraceData output = new RayTraceData(d_r, theta_r);
		
		return output;
	}	
}