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
		double X = 18000*Constants.SIGMA/path.input.f;
		double Y = Constants.E_R - Math.pow(Math.cos(psi), 2);
		double T = Math.sqrt(Math.pow(Y, 2) + Math.pow(X, 2));
		
		double P = Math.sqrt(0.5*T);
		double Q = X/(2*P);
		double B = 1/(Math.pow(P, 2) + Math.pow(Q, 2));
		double A = (2*P)/(Math.pow(P, 2) + Math.pow(Q, 2));
		
		
		// <<Step 2>>
		GroundReflectionData data = new GroundReflectionData();
		data.R_g = Math.sqrt(((1+B*Math.pow(Math.sin(psi), 2))-A*Math.sin(psi)) / ((1+B*Math.pow(Math.sin(psi), 2))+A*Math.sin(psi)));
		data.phi_g = Math.atan2(-Q, Math.sin(psi)-P) - Math.atan2(Q, Math.sin(psi)+P);
		return data;
	}
}
