package main.p528;


/** Rec. ITU-R P.528-4 Annex II Section X
 * 
 * @author Collin Thornton
 *
 */
public class SmoothEarthDiffraction {

	public static double compute(TerminalGeometry.Geom lt, TerminalGeometry.Geom ht, Path path, double d_0) {
		
		// <<Step 1> Compute the normalized distances
		double[] x = { 1.607*Math.pow(path.input.f, 1/3)*d_0, 1.607*Math.pow(path.input.f, 1/3)*lt.d, 1.607*Math.pow(path.input.f, 1/3)*ht.d };
		
		
		// <<Step 2>> Compute distance dependent term
		double[] G_x = new double[3];
		for(int i=0; i<3; ++i) G_x[i] = 0.05751*x[i] - 10*Math.log10(x[i]);
		
		
		// <<Step 3>>
		double[] y = new double[2];
		for(int i=0; i<2; ++i) y[i] = 40*Math.log10(x[i+1]) - 117;
		
		
		// <<Step 4>> Compute the height functions
		double[] F_x = new double[2];
		for(int i=1; i<3; ++i) {
			if(x[i] >= 2000) F_x[i-1] = G_x[i];
			else if(x[i] >= 200) {
				double W = 0.0134*x[i]*Math.exp(-0.005*x[i]);
				F_x[i-1] = W*y[i-1] + (1-W)*G_x[i];
			}
			else F_x[i-1] = y[i-1];
		}
		
		return G_x[0] - F_x[0] - F_x[1] - 20;
	}
}