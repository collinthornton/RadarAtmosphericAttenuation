package main.p528;


/** Rec. ITU-R P.528-4 Annex II Section IX
 * 
 * @author Collin Thornton
 *
 */
public class GroundReflection {
	public static class GroundReflectionData {
		double R_g;			// Real part of reflection coefficient
		double phi_g;		// Imaginary part of reflection coefficient
	}
	public static GroundReflectionData compute(double psi, Path path) {
		
		// <<Step 1>>
		double X = 18000.0*Constants.SIGMA/path.input.f;
		double Y = Constants.E_R - Math.cos(psi)*Math.cos(psi);
		double T = Math.sqrt(Y*Y + X*X) + Y;
		
		double P = Math.sqrt(0.5*T);
		double Q = X/(2.0*P);
		double B = 1.0/(P*P + Q*Q);
		double A = (2.0*P)/(P*P + Q*Q);
		
		
		// <<Step 2>>
		GroundReflectionData data = new GroundReflectionData();
		data.R_g = Math.sqrt( (1.0 + B*Math.sin(psi)*Math.sin(psi) - A*Math.sin(psi)) / (1.0 + B*Math.sin(psi)*Math.sin(psi) + A*Math.sin(psi)) );
		
		data.phi_g = Math.atan2(-Q, Math.sin(psi)-P) - Math.atan2(Q, Math.sin(psi)+P);
		return data;
	}
}
