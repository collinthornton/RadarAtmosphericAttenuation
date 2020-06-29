package main.p528;

/** Class to calculate the atmospheric attenuation of radar transmission. Rec. ITU-R P.528-4, Annex II
 * 
 * @author Collin Thornton
 * @note Assumes horizontal polarization
 * @note Based on Rec. ITU-R P.528-4, Annex II
 */
public class AtmosphericAttenuationCalculator {
	/** Calculate atmospheric attenuation in dB/km
	 * 
	 * @param f		Frequency 				(GHz)
	 * @param h_r1	Height of low terminal 	(km)
	 * @param h_r2	Height of high terminal (km)
	 * @param q		Time percentage			(0.00-1.00)
	 * @param d		Distance				(km)
	 * @return	dB/km
	 * @note Follows the step-by-step method of Annex II, Section III
	 */
	public static double compute(double f, double h_r1, double h_r2, double q, double d) {
		f = f*1000; // Convert to MHz to comply with Rec. ITU-R P.528-4
		
		Path path = new Path(h_r1, h_r2, f, q, d);
		
		
		// <<Step 1>> Compute the geometric properties of each terminal
		TerminalGeometry low_terminal 	= new TerminalGeometry(h_r1);
		TerminalGeometry high_terminal 	= new TerminalGeometry(h_r2);

		
		// <<Step 2>> Compute the maximum line-of-sight distance between the terminals
		path.d_ML = low_terminal.geom.d + high_terminal.geom.d;
		
				
		// <<Step 3>> Compute the smooth Earth diffraction loss
		double d_3 = path.d_ML + 0.5*Math.pow(Constants.A_E*Constants.A_E/path.atten.A_fs, 1/3);
		double d_4 = path.d_ML + 1.5*Math.pow(Constants.A_E*Constants.A_E/path.atten.A_fs, 1/3);
		
		double A_d_3 = SmoothEarthDiffraction.compute(low_terminal.geom, high_terminal.geom, path, d_3);
		double A_d_4 = SmoothEarthDiffraction.compute(low_terminal.geom, high_terminal.geom, path, d_4);
		
		double M_d = (A_d_4 - A_d_3) / (d_4 - d_3);
		path.atten.A_d0 = A_d_4 - M_d*d_4;
		
		path.atten.A_dML = M_d*path.d_ML + path.atten.A_d0;
		path.d_d = -(path.atten.A_d0/M_d);
		
		
		// <<Step 4>> Determine if path is in line-of-sight region or transhorizon
		if(path.input.d < path.d_ML) {
				LOSAttenuationCalculator los = new LOSAttenuationCalculator();
				double total_attenuation = los.compute(low_terminal.geom, high_terminal.geom, path);
		
				return total_attenuation;	
		}
		
		
		// HANDLE OTHER MODES HERE
		return -999999;
	}
	
	
	public static void main(String[] args) {
		double f = 5;		// Frequency 				(GHz)
		
		double h1 = 0.36;	// Height of lower terminal (km)
		double h2 = 0.50;	// Height of upper terminal (km)
		double q = 0.99;	// Time percentile			(0.00-1.00)
		double d = 60.0;		// Path distance			(km)
		
		
		System.out.println(compute(f, h1, h2, q, d));
	}
}