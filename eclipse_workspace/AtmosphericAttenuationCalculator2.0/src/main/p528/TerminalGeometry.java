package main.p528;

/** Rec. ITU-R P.528-4 Annex II Section IV
 * 
 * @author Collin Thornon
 *
 */
public class TerminalGeometry {
	public class Geom {
		public final double d;			// arc length to Earth 								(km)
		public final double theta;		// incident angle from terminal to Earth horizon 	(rad)
		public final double h;			// adjusted height of terminal above MSL			(km)
		public final double delta_h;	// terminal height correction						(km)
		
		Geom(double d, double theta, double h, double delta_h) {
			this.d = d;
			this.theta = theta;
			this.h = h;
			this.delta_h = delta_h;
		}
	}
	
	public Geom geom;
	
	TerminalGeometry(double h_r) { compute(h_r); }
	
	
	/** Compute arc length to Earth, incident angle from terminal to Earth, and adjust height of terminal above MSL
	 * 
	 * @param h_r Real terminal height above ground (km)
	 */
	public void compute(double h_r) {
		double d, theta, h, delta_h;
		
		
		// <<Step 1>> Utilize ray tracing to compute arc distance and incident angle
		RayTrace.RayTraceData trace = RayTrace.compute(h_r);
		
		double d_r = trace.d_r;
		double theta_r = trace.theta_r;
		
		
		// <<Step 2>> Compute effective height
		double phi = d_r/Constants.A_E;
		double h_e;
		
		if(phi <= 0.1) h_e = d_r*d_r/(2.0*Constants.A_E);
		else h_e = Constants.A_E/Math.cos(phi) - Constants.A_E;
		
		
		// <<Step 3>> Account for overestimation
		if(h_e <= h_r) {
			h = h_e;
			d = d_r;
		}
		else {
			h = h_r;
			d = Math.sqrt(2.0*Constants.A_E*h_r);
		}
		theta = theta_r;
		
		
		// <<Step 4>> Compute terminal height correction
		delta_h = h_r - h;
		
		
		// <<Step 5>> Perform corrections if necessary
		if(delta_h <= 0.0) {
			theta = Math.sqrt(2.0*h_r/Constants.A_E);
			d = Math.sqrt(2.0*h_r*Constants.A_E);
		}
		
		geom = new Geom(d, theta, h, delta_h);
	}
	
	
	public static void main(String[] args) {
		double h = 600.0;
		
		TerminalGeometry terminal = new TerminalGeometry(h);
		System.out.format("Arc length:\t\t%7.3f%n", 		terminal.geom.d);
		System.out.format("Incident angle:\t\t%7.3f%n", 	terminal.geom.theta);
		System.out.format("Adjusted height:\t%7.3f%n", 		terminal.geom.h);
		System.out.format("Height correction:\t%7.3f%n", 	terminal.geom.delta_h);
		return;
	}
}
