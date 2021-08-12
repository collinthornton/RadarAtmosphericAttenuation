package main.p528;


/** Rec. ITU-R P.528-4 Annex II Section VIII
 * 
 * @author Collin Thornton
 *
 */
public class LOSPathLoss {
	public static double compute(double psi, double psi_limit, RayOptics.RayOpticsData optics, Path path) {
		double A_LOS;
		
		// <<Step 1>>
		if(path.input.d > path.d_0) return ((path.input.d-path.d_0)*(path.atten.A_dML-path.atten.A_d0)/(path.d_ML-path.d_0)) + path.atten.A_d0;
		
		
		// <<Step 2>>
		if(psi < psi_limit) return 0;
		
		
		// <<Step 3>> Determine wavelength 
		double lambda = 0.2997925 / path.input.f;
		
		// Compute reflection coefficient <<Step 4>>
		GroundReflection.GroundReflectionData gr = GroundReflection.compute(psi, path);
		
		
		// <<Step 5>> Compute divergence factor  
		// --> r_1 and r_2 and R_r were determined from line 39  of https://github.com/NTIA/p528/blob/master/src/GetPathLoss.cpp
		double r_1 = optics.D1 / Math.cos(psi);
		double r_2 = optics.D2 / Math.cos(psi);
		double R_r = (r_1 * r_2) / optics.r_12;
		
		double D_v = 1.0 / Math.sqrt(  1.0 + ((2.0*R_r*(1.0+(Math.sin(psi)*Math.sin(psi))))/(optics.a_a*Math.sin(psi))) + ((2.0*R_r)/optics.a_a)*((2.0*R_r)/optics.a_a)  );
		
		// <<Step 6>> Compute the ray length factor 
		double F_r = Math.min(optics.r_0/optics.r_12, 1.0);
		
		
		// <<Step 7>> Compute effective reflection coefficients 
		path.R_Tg = gr.R_g*D_v*F_r;
		double phi_Tg = ((2.0*Math.PI*optics.deltar/lambda) + gr.phi_g);
		
		
		// <<Step 8>> Compute the loss 
		double C = path.R_Tg*Math.cos(phi_Tg);
		double S = -path.R_Tg*Math.sin(phi_Tg);
		
		double R = C + S;
		double W_RL = Math.min(Math.abs(1.0+R), 1.0);
		double W_R0 = W_RL*W_RL + 0.0001;
		
		A_LOS = 10.0 * Math.log10(W_R0);
		
		//System.out.format("%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f%n", psi, psi_limit, optics.D1, optics.D2, optics.r_0, optics.r_12, optics.deltar, gr.phi_g, gr.R_g, A_LOS, path.R_Tg, C, S, R, W_RL, W_R0, phi_Tg, F_r, D_v);
		
		return A_LOS;
	}
	
	public static double compute(double psi, double psi_limit, RayOptics.RayOpticsData optics, Path path, double A_d0) {
		double A_LOS;
		
		// <<Step 1>>
		if(path.input.d > path.d_0) return ((path.input.d-path.d_0)*(path.atten.A_dML-A_d0)/(path.d_ML-path.d_0)) + A_d0;
		
		
		// <<Step 2>>
		if(psi < psi_limit) return 0;
		
		
		// <<Step 3>> Determine wavelength 
		double lambda = 0.2997925 / path.input.f;
		
		// Compute reflection coefficient <<Step 4>>
		GroundReflection.GroundReflectionData gr = GroundReflection.compute(psi, path);
		
		
		// <<Step 5>> Compute divergence factor  
		// --> r_1 and r_2 and R_r were determined from line 39  of https://github.com/NTIA/p528/blob/master/src/GetPathLoss.cpp
		double r_1 = optics.D1 / Math.cos(psi);
		double r_2 = optics.D2 / Math.cos(psi);
		double R_r = (r_1 * r_2) / optics.r_12;
		
		double D_v = 1.0 / Math.sqrt(  1.0 + ((2*R_r*(1.0+Math.sin(psi)*Math.sin(psi)))/(optics.a_a*Math.sin(psi))) + ((2.0*R_r)/optics.a_a)*((2.0*R_r)/optics.a_a)  );
		
		// <<Step 6>> Compute the ray length factor 
		double F_r = Math.min(optics.r_0/optics.r_12, 1.0);
		
		
		// <<Step 7> Compute effective reflection coefficients 
		path.R_Tg = gr.R_g*D_v*F_r;
		double phi_Tg = (2.0*Math.PI*optics.deltar/lambda) + gr.phi_g;
		
		
		// <<Step 8>> Compute the loss 
		double R = path.R_Tg*Math.cos(phi_Tg) - path.R_Tg*Math.sin(phi_Tg);
		double W_RL = Math.min(Math.abs(1.0+R), 1.0);
		double W_R0 = W_RL*W_RL + 0.0001;
		
		A_LOS = 10.0 * Math.log10(W_R0);
		
		//System.out.format("%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f%n", psi, psi_limit, optics.D1, optics.D2, optics.r_12, gr.phi_g, gr.R_g, A_LOS, R, W_RL, W_R0, phi_Tg, F_r, D_v);
		
		return A_LOS;
	}
}