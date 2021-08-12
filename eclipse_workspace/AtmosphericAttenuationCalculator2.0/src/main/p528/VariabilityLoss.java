package main.p528;

public class VariabilityLoss {
	private static final double[][] table7 = {
			{ 0.1440, -40.0 },
			{ 0.8421, -25.0 },
			{ 1.5544, -20.0 },
			{ 2.0014, -18.0 },
			{ 2.5931, -16.0 },
			{ 3.3872, -14.0 },
			{ 4.4715, -12.0 },
			{ 5.9833, -10.0 },
			{ 8.1418, -8.0  },
			{ 11.0972, -6.0 },
			{ 14.2546, -4.0 },
			{ 16.4258, -2.0 },
			{ 18.0527,  0.0 },
			{ 18.0527,  2.0 },
			{ 18.3573,  4.0 },
			{ 18.3361,  6.0 },
			{ 18.3864, 20.0 }
	};
	
	public static double compute(double A_T, TerminalGeometry.Geom lt, TerminalGeometry.Geom ht, Path path, RayOptics.RayOpticsData optics, double r_ew) {
		
		// <<Step 1>> Compute input value from ray optics
		double f_theta_h;
		if(optics.theta_h1 <= 0.0) 		f_theta_h = 1.0;
		else if(optics.theta_h1 >= 1.0) f_theta_h = 0.0;
		else 	f_theta_h = Math.max(0.5 - 0.3183*Math.atan(20.0*Math.log10(32.0*optics.theta_h1)), 0.0);
		
		
		// <<Step 2>> Compute contribution of long term variability for time percentage Y_pi_99
		double Y_e_q = LongTermVariability.compute(lt, ht, path, path.input.q, f_theta_h, A_T);
		
		
		// <<Step 3>> Compute long term variability for time percentage Y_pi_99 = 0.5
		double Y_e_05 = LongTermVariability.compute(lt, ht, path, 0.5, f_theta_h, A_T);
		
		
		// <<Step 4>> Compute parameter of effects of tropospheric multipath
		double F_A_Y, F_delta_r;
		double lambda = 0.2997925/path.input.f;
		if(path.atten.A_Y <= 0.0) F_A_Y = 1.0;
		else if(path.atten.A_Y >= 9.0) F_A_Y = 0.1;
		else F_A_Y = (1.1 + 0.9*Math.cos(path.atten.A_Y*Math.PI/9.0))/2.0;
		
		if(optics.deltar >= lambda/2.0) F_delta_r = 1.0;
		else if(optics.deltar >= lambda/6.0) F_delta_r = 0.1;
		else F_delta_r = 0.5*(1.1 - 0.9*Math.cos((3*Math.PI/lambda)*(optics.deltar - lambda/6.0)));
		
		double R_s = path.R_Tg*F_delta_r*F_A_Y;
		
		double W_a;
		if(r_ew <= 0.0) W_a = 0.0001;
		else {
			double Y_pi_99 = 10.0 * Math.log10(path.input.f*r_ew*r_ew*r_ew) - 84.26;
			
			double K = table7[16][1];
			if(Y_pi_99 < table7[0][0]) K = table7[0][1];
			else {
				int i=1;
				while(i < 16 && Y_pi_99 > table7[i][0]) ++i;
				if(Y_pi_99 == table7[i][0]) 		K = table7[i][1];
				else if(Y_pi_99 < table7[i][0]) 	K = ((table7[i][1]-table7[i-1][1])*(Y_pi_99-table7[i-1][0]))/(table7[i][0]-table7[i-1][0]) + table7[i-1][1];
			}
			
			W_a = Math.pow(10.0, 0.1*K);
		}
		
		double W_R = R_s*R_s + 0.0001;
		double W = W_R + W_a;
		double K_LOS;
		
		if(W <= 0.0) 	K_LOS = -40.0;
		else		K_LOS = 10.0*Math.log10(W);
		K_LOS = Math.max(K_LOS, -40.0);
		
		
		// <<Step 5>> Compute contribution of tropospheric multipath
		double Y_pi = NakagamiRice.compute(K_LOS, path.input.q);
		
		
		// <<Step 6>> Sum the effects to get total variability contribution
		double Y_total_05 = Y_e_05 + 0.0;
		double Y = Math.sqrt((Y_e_q-Y_e_05)*(Y_e_q-Y_e_05) + Y_pi*Y_pi);
		double Y_total;
		
		if(path.input.q < 0.50) 	Y_total = Y_total_05 + Y;
		else				Y_total = Y_total_05 - Y;
		return Y_total;
	}
	
	
	public static void main(String[] args) {
		double h_r1 = 2.0, h_r2 = 2.2, d= 50.0;	// km
		double f = 5.0*1000.0;					// MHz
		double q = 0.95;					// %
				
		Path path = new Path(h_r1, h_r2, f, q, d);
		
		TerminalGeometry low_terminal 	= new TerminalGeometry(h_r1);
		TerminalGeometry high_terminal 	= new TerminalGeometry(h_r2);
		TerminalGeometry.Geom lt = low_terminal.geom;
		TerminalGeometry.Geom ht = high_terminal.geom;
		
		RayOptics.RayOpticsData optics = RayOptics.compute(0.8, path, lt, ht);
		
		
		double r_ew = 5.752349152657722;
		double A_LOS = 3.0;
		
		double Y_q = VariabilityLoss.compute(A_LOS, lt, ht, path, optics, r_ew);
		System.out.println(Y_q);
	}
}
