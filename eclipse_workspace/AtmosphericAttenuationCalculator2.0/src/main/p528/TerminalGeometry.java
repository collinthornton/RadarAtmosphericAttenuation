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
		
		RayTrace.RayTraceData trace = RayTrace.compute(h_r);
		
		double d_r = trace.d_r;
		double theta_r = trace.theta_r;
		
		double phi = d_r/Constants.A_E;
		double h_e;
		
		if(phi <= 0.1) h_e = Math.pow(d_r, 2)/(2*Constants.A_E);
		else h_e = Constants.A_E/Math.cos(phi) - Constants.A_E;
		
		if(h_e <= h_r) {
			h = h_e;
			d = d_r;
		}
		else {
			h = h_r;
			d = Math.sqrt(2*Constants.A_E*h_r);
		}
		theta = theta_r;
		
		delta_h = h_r - h;
		
		if(delta_h <= 0) {
			theta = Math.sqrt(2*h_r/Constants.A_E);
			d = Math.sqrt(2*h_r*Constants.A_E);
		}
		
		geom = new Geom(d, theta, h, delta_h);
	}
}
