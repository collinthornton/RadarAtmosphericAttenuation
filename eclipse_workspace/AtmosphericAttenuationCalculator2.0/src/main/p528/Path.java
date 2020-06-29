package main.p528;

public class Path {
	public class UserInput {
		final double h_r1, h_r2;		// terminal heights (1 is lower)
		final double f;				// frequency in MHz
		final double q;				// Time percentage
		final double d;				// distance (km)
		
		UserInput(double h1, double h2, double fr, double qi, double dd) {
			h_r1 = h1;
			h_r2 = h2;
			f = fr;
			q = qi;
			d = dd;
		}
	}


	public class Attenuation {
		/** Intercept smooth Earth diffraction loss */
		double A_d0;
		
		/** Smooth Earth diffraction loss at maximum line-of-sight */
		double A_dML; 
		
		/** Line-of-sight loss */
		double A_LOS;
		
		/** Loss due to atmospheric absorption (o2 and h2o vapor) */
		double A_a;
		
		/** Free space loss */
		double A_fs;
		
		/** Loss due to variability */
		double A_Y;
		
		/** Total loss */
		double A;
	}
	
	
	// These values must be initialized after construction
	
	/** Distance at which diffraction begins to affect ray */
	double d_0;
	
	/** Maximum line-of-sight distance */
	double d_ML; 
	
	/** Distance predicted to have 0 diffraction loss */
	double d_d;	
	
	/** Effective reflection coefficient */
	double R_Tg;
	
	
	UserInput input;
	Attenuation atten;
	
	/** Class to hold path parameters
	 *  
	 * @param h1	Actual height of low terminal	(km)
	 * @param h2	Actual height of high terminal	(km)
	 * @param f		Frequency	(MHz)
	 * @param q		Time percentage	(0.00-1.00)
	 * @param d		Path distance	(km)
	 */
	Path(double h1, double h2, double f, double q, double d) {
		input = new UserInput(h1, h2, f, q, d);
		atten = new Attenuation();
	}
}