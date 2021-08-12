package main.p528;

import java.io.File;
import java.io.PrintWriter;
import java.io.IOException;

/** Class to calculate the atmospheric attenuation of radar transmission. Rec. ITU-R P.528-4, Annex II
 * 
 * @author Collin Thornton
 * @note Assumes horizontal polarization
 * @note Based on Rec. ITU-R P.528-4, Annex II
 * @note References https://github.com/NTIA/p528
 */
public class AtmosphericAttenuationCalculator {
	/** Calculate atmospheric attenuation in dB/km
	 * 
	 * @param f		Frequency 				(GHz)
	 * @param h_r1	Height of low terminal 	(km)
	 * @param h_r2	Height of high terminal (km)
	 * @param q		Time percentage			(0.01-0.99)
	 * @param d		Distance				(km)
	 * @return	dB
	 * @note Follows the step-by-step method of Annex II, Section III
	 */
	public static Path compute(double f, double h_r1, double h_r2, double q, double d) {
		f = f*1000.0; // Convert to MHz to comply with Rec. ITU-R P.528-4
		
		Path path = new Path(h_r1, h_r2, f, q, d);
		
		
		// <<Step 1>> Compute the geometric properties of each terminal
		TerminalGeometry low_terminal 	= new TerminalGeometry(h_r1);
		TerminalGeometry high_terminal 	= new TerminalGeometry(h_r2);

		
		// <<Step 2>> Compute the maximum line-of-sight distance between the terminals
		path.d_ML = low_terminal.geom.d + high_terminal.geom.d;
		
				
		// <<Step 3>> Compute the smooth Earth diffraction loss
		double term = Math.pow(Constants.A_E*Constants.A_E/path.input.f, 1.0/3.0);
		double d_3 = path.d_ML + 0.5*term;
		double d_4 = path.d_ML + 1.5*term;
		
		double A_d_3 = SmoothEarthDiffraction.compute(low_terminal.geom, high_terminal.geom, path, d_3);
		double A_d_4 = SmoothEarthDiffraction.compute(low_terminal.geom, high_terminal.geom, path, d_4);
		
		double M_d = (A_d_4 - A_d_3) / (d_4 - d_3);
		path.atten.A_d0 = A_d_4 - M_d*d_4;
		
		path.atten.A_dML = M_d*path.d_ML + path.atten.A_d0;
		path.d_d = -(path.atten.A_d0/M_d);
		
		
		// <<Step 4>> Determine if path is in line-of-sight region or transhorizon
		if(path.input.d < path.d_ML) {
				LOSAttenuationCalculator los = new LOSAttenuationCalculator();
				double total_attenuation = los.compute(low_terminal.geom, high_terminal.geom, path);
				path.atten.A = total_attenuation;
				return path;	
		}
		
		path.atten.A = -999;
		// HANDLE OTHER MODES HERE
		return path;
	}
	
	
	public static void main(String[] args) {
		double f = 10;		// Frequency 				(GHz)
		
		double h1 = 0.300;	// Height of lower terminal (km)
		double h2 = 0.500;	// Height of upper terminal (km)
		double q = 0.99;	// Time percentile			(0.00-1.00)
		double d = 1.0;	// Path distance			(km)
		
		final double MAX_DISTANCE = 150;
		final double MAX_FREQ = 20;
		
		long sum_time = 0;
		
		try {
			File file = new File("output.txt");
			file.createNewFile();
		} catch (IOException e) {
			System.out.println("Error creating file");
		}
		
		
		try {
			PrintWriter file = new PrintWriter("output.txt");
			file.println("f\tA\tA_LOS\tA_fs\tA_Y\tA_a");
			
			System.out.println("f\tA\tA_LOS\tA_fs\tA_Y\tA_a");
			//System.out.println("psi\tpsi_limit\toptics.D1\toptics.D2\toptics.r_0\toptics.r_12\toptics.deltar\tgr.phi_g\tgr.R_g\tA_LOS\tpath.R_Tg\tC\tS\tR\tW_RL\tW_R0\tphi_Tg\tF_r\tD_v");
		
			for(double i=0.1; i<MAX_FREQ; i+=0.02) {
				long begin_time = System.currentTimeMillis();
				Path path = compute(i, h1, h2, q, d);
				System.out.format("%4.3f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%12.9f%n", i, path.atten.A, path.atten.A_LOS, path.atten.A_fs, path.atten.A_Y, path.atten.A_a);
				//file.format("%4.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f%n", i, path.atten.A, path.atten.A_LOS, path.atten.A_fs, path.atten.A_Y, path.atten.A_a);
				sum_time += System.currentTimeMillis() - begin_time;
			}
			
			double avg_time = (double)sum_time / (double)150*100-1;
			System.out.println("Average time over " + MAX_DISTANCE + " iterations: " + avg_time + " ms.");
			
			file.close();
		} catch (IOException e) {
			System.out.println("Error writing to file:");
			e.printStackTrace();
		}
	

		
	}
}
