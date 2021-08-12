package main.p528;


/** Rec. ITU-R P.528-4 Annex II Section X
 * 
 * @author Collin Thornton
 *
 */
public class SmoothEarthDiffraction {

	public static double compute(TerminalGeometry.Geom lt, TerminalGeometry.Geom ht, Path path, double d_0) {
		
		// <<Step 1>> Compute the normalized distances
		double x0 = 1.607*Math.pow(path.input.f, 1.0/3.0)*d_0;
		double x1 = 1.607*Math.pow(path.input.f, 1.0/3.0)*lt.d;
		double x2 = 1.607*Math.pow(path.input.f, 1.0/3.0)*ht.d;
		
		
		// <<Step 2>> Compute distance dependent term
		double G_x0 = 0.05751*x0 - 10.0*Math.log10(x0);
		double G_x1 = 0.05751*x1 - 10.0*Math.log10(x1);
		double G_x2 = 0.05751*x2 - 10.0*Math.log10(x2);
		
		
		// <<Step 3>>
		double y1 = 40.0*Math.log10(x1) - 117.0;
		double y2 = 40.0*Math.log10(x2) - 117.0;
		
		
		// <<Step 4>> Compute the height functions
		double F_x1, F_x2;
		
		if(x1 >= 2000) F_x1 = G_x1;
		else if(x1 > 200) {
			double W = 0.0134*x1*Math.exp(-0.005*x1);
			F_x1 = W*y1 + (1-W)*G_x1;
		}
		else F_x1 = y1;
		
		if(x2 >= 2000) F_x2 = G_x2;
		else if(x2 > 200) {
			double W = 0.0134*x2*Math.exp(-0.005*x2);
			F_x2 = W*y2 + (1-W)*G_x2;
		}
		else F_x2 = y2;
		
		
		
		return G_x0 - F_x1 - F_x2 - 20.0;
	}
}
