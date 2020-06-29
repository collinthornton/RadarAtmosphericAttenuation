package main.p528;

/** Rec. ITU-R P.528-4 Annex II Section VII
 * 
 * @author Collin Thornton
 *
 */
public class RayOptics {
	
	public static class RayOpticsData {
		public double deltar, d;
		public double theta_h1, theta_h2;		
		public double a_a, D1, D2, r_0, r_12;
		public double z_1, z_2;
	}
	
	public static RayOpticsData compute(double psi, Path path, TerminalGeometry.Geom lt, TerminalGeometry.Geom ht) {
		RayOpticsData optics = new RayOpticsData();
				
		double z = (Constants.A_0/Constants.A_E) - 1;		// (Eq. 62)
		double k_a = 1/(1+z*Math.cos(psi));
		optics.a_a = Constants.A_0*k_a;
		
		double deltah_a1 = lt.delta_h*(optics.a_a - Constants.A_0)/(Constants.A_E-Constants.A_0);
		double deltah_a2 = ht.delta_h*(optics.a_a - Constants.A_0)/(Constants.A_E-Constants.A_0);
		
		double H1 = path.input.h_r1 - deltah_a1, H2 = path.input.h_r2 - deltah_a2;
		
		optics.z_1 = optics.a_a + H1;
		optics.z_2 = optics.a_a + H2;
		double theta_1 = Math.acos(optics.a_a*Math.cos(psi)/optics.z_1) - psi;
		double theta_2 = Math.acos(optics.a_a*Math.cos(psi)/optics.z_2) - psi;
		optics.D1 = optics.z_1*Math.sin(theta_1);
		optics.D2 = optics.z_2*Math.sin(theta_2);
		double Hprime_1, Hprime_2;
		if(psi > 1.56) {
			Hprime_1 = H1;
			Hprime_2 = H2;
		} else {
			Hprime_1 = optics.D1*Math.tan(psi);
			Hprime_2 = optics.D2*Math.tan(psi);
		}
		
		// double deltaz = Math.abs(z_1-z_2); 				// (Eq. 71)
		
		optics.d = Math.max(optics.a_a*(theta_1+theta_2), 0);
		
		double alpha = Math.atan((Hprime_2-Hprime_1)/(optics.D1+optics.D2));
		optics.r_0 = (optics.D1+optics.D2)/Math.cos(alpha);
		optics.r_12 = (optics.D1+optics.D2)/Math.cos(psi);
		
		optics.deltar = 4*Hprime_1*Hprime_2/(optics.r_0+optics.r_12);
		optics.theta_h1 = alpha - theta_1;
		optics.theta_h2 = -1*(alpha + theta_2);
		return optics;
	}
}
