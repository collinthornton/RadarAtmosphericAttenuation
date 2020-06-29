package main.p528;

public class LongTermVariability {
	private static final double[][] table4 = {
			{ 0.01, 1.9507 },
			{ 0.02, 1.7166 },
			{ 0.05, 1.3265 },
			{ 0.10, 1.0000 }
	};
	private static final double[][] table5 = {
			{ 0.01, -5.00 },
			{ 0.02, -4.50 },
			{ 0.05, -3.70 },
			{ 0.10,  0.00 }
	};
	
	public static double compute(TerminalGeometry.Geom lt, TerminalGeometry.Geom ht, Path path, double q, double f_theta_h, double A_T) {
		
		// <<Step 1>> Compute smooth Earth horizon distance via ray tracing
		RayTrace.RayTraceData trace = RayTrace.compute(lt.h, 329);
		double d_Lq_1 = trace.d_r;
		trace = RayTrace.compute(ht.h, 329);
		double d_Lq_2 = trace.d_r;
		
		
		// <<Step 2>> Compute effective distance between terminals
		double d_qs = 60*Math.pow(100/path.input.f, 1/3);
		double d_Lq = d_Lq_1 + d_Lq_2;
		double d_q = d_Lq + d_qs;
		double d_e;
		
		if(path.input.d <= d_q) d_e = (130*path.input.d)/d_q;
		else d_e = 130 + path.input.d - d_q;
		
		
		// <<Step 3>>
		double g_01, g_09;
		if(path.input.f <= 1600) {
			g_01 = 0.21*Math.sin(5.22*Math.log10(path.input.f/200)) + 1.28;
			g_09 = 0.18*Math.sin(5.22*Math.log10(path.input.f/200)) + 1.23;
		}
		else {
			g_01 = 1.05;
			g_09 = 1.05;
		}
		
		
		// <<Step 4>>
		double f_2v  = 0.0 + (3.9-0.0)*Math.exp(-1.56e-11*Math.pow(d_e, 4.08));
		double f_2y1 = 5.4 + (10.0-5.4)*Math.exp(-1.57e-6*Math.pow(d_e, 2.31));
		double f_2y9 = 3.2 + (8.2 -3.2)*Math.exp(-3.75e-8*Math.pow(d_e, 2.88));
				
		double V = (1.59e-5*Math.pow(d_e, 2.32)-f_2v) * Math.exp(-2.77e-8*Math.pow(d_e, 3.25)) + f_2v;
		double Y_1 = (5.24e-4*Math.pow(d_e, 1.97) - f_2y1) * Math.exp(-4.70e-7*Math.pow(d_e, 2.90)) + f_2y1;
		double Y_9 = (2.93e-4*Math.pow(d_e, 2.00)-f_2y9) * Math.exp(-1.02e-7*Math.pow(d_e, 3.15)) + f_2y9;
		
		
		// <<Step 5>> Compute variability associated w/ long-term power fading
		double Y_q, Y, z_q, c_q;
		if(q == 0.50) Y_q = V;
		else if(q > 0.50) {
			double z_09 = inverseComplementaryCumulativeNormal(0.9);
			z_q = inverseComplementaryCumulativeNormal(q);
			c_q = z_q/z_09;
			Y = c_q*(-Y_9*g_09);
			Y_q = Y + V;
		} else if(q >= 0.10) {
			double z_01 = inverseComplementaryCumulativeNormal(0.1);
			z_q = inverseComplementaryCumulativeNormal(q);
			c_q = z_q/z_01;
			Y = c_q*(Y_1*g_01);
			Y_q = Y + V;
		} else {
			// Linearly interpolate from Table 4
			c_q = table4[3][1];
			if(q < table4[0][0]) c_q = table4[0][0];
			else {
				int i=1;
				while(q > table4[i][0] && i < 4) ++i;
				if(q == table4[i][0]) 		c_q = table4[i][1];
				else if(q < table4[i][0]) 	c_q = ((table4[i][1]-table4[i-1][1])*(q-table4[i-1][0]))/(table4[i][0]-table4[i-1][0]) + table4[i-1][1];
			}
			Y = c_q*(Y_1*g_01);
			Y_q = Y + V;
		}
			
		
		// <<Step 6>> Determine long term power fading for q = 0.10
		double Y_01 = (Y_1*g_01) + V;
		
		
		// <<Step 7>>
		double Y_el_q 	= f_theta_h * Y_q;
		double Y_el_01 	= f_theta_h * Y_01;
		
		
		// <<Step 8>>
		double A_YI = Y_el_01 - A_T - 3;
		path.atten.A_Y = Math.max(A_YI, 0);
		
		
		// <<Step 9>> Compute total variability loss if q >- 0.10
		if(q >= 0.10) return Y_el_q - path.atten.A_Y;
		
		
		// <<Step 10>> Apply corrections for q < 0.10
		double Y_temp = Y_el_q - path.atten.A_Y - A_T;
		
		
		// <<Step 11>> Linearly interpolate c_yq from q using table 5
		c_q = table5[3][1];
		if(q < table5[0][0]) c_q = table5[0][0];
		else {
			int i=1;
			while(q > table5[i][0] && i < 4) ++i;
			if(q == table5[i][0]) 		c_q = table5[i][1];
			else if(q < table5[i][0]) 	c_q = ((table5[i][1]-table5[i-1][1])*(q-table5[i-1][0]))/(table5[i][0]-table5[i-1][0]) + table5[i-1][1];
		}
		
		if(Y_temp > -c_q) 	return -c_q + A_T;
		else 				return Y_temp + A_T;
	}
	
	private static double inverseComplementaryCumulativeNormal(double q) {
	    final double C_0 = 2.515516;
	    final double C_1 = 0.802853;
	    final double C_2 = 0.010328;
	    final double D_1 = 1.432788;
	    final double D_2 = 0.189269;
	    final double D_3 = 0.001308;

	    double x = q;

	    if (q > 0.5)  x = 1.0 - x;

	    double T_x = Math.sqrt(-2.0 * Math.log(x));
	    double zeta_x = ((C_2 * T_x + C_1) * T_x + C_0) / (((D_3 * T_x + D_2) * T_x + D_1) * T_x + 1.0);
	    double Q_q = T_x - zeta_x;

	    if (q > 0.5)  Q_q = -Q_q;

	    return Q_q;
	}
}