package main;

import main.p528.AtmosphericAttenuationCalculator;
import main.p838.RainAttenuationCalculator;
import main.p840.CloudsAttenuationCalculator;

/** Class to calculate the total atmospheric attenuation of radar transmission
 * 
 * @author Collin Thornton
 * @note Assumes horizontal polarization
 */
public class TransmissionLossCalculator {
	/** Cacluate the total attenuation. (dB)
	 * 
	 * @param f Frequency (GHz)
	 * @param d One-way distance (km)
	 * @param gtx Gain of transmitter (dB)
	 * @param grx Gain of receiver (dB)
	 * @return total attenuation adjusted for 2 way loss 
	 */
	public static double calcTotalAttenuation(double f, double h_r1, double h_r2, double q, double d, double T, double cd,
			double rr, double theta, double gtx, double grx) {
		
		double total_attenuation = calcFreeSpaceLoss(f, d);
		total_attenuation += CloudsAttenuationCalculator.calculate(f, T, cd)*d;
		total_attenuation += RainAttenuationCalculator.calculate(f, rr, theta, 0)*d;
		total_attenuation += AtmosphericAttenuationCalculator.compute(f, h_r1, h_r2, q, cd);
		
		total_attenuation -= gtx;
		total_attenuation -= grx;
		total_attenuation += 3;	
		
		return total_attenuation;
	}
	
	/** Calculate the attenuation due to propagation through free space
	 * 
	 * @param f Frequency 				(GHz)
	 * @param d Distance 				(km)
	 * @return Attenuation 				(dB)
	 * @note Based on Rec. ITU-R P.528-4 Annex III, Step 8
	 */
	public static double calcFreeSpaceLoss(double f, double d) {
		return 32.45 + 20*Math.log10(f*1000*d);
	}
	
	
	
	
	
	
	
	
	public static void main(String[] args) {
		double f = 9.5;
		double d = 10;
		
		double T = 293.15;
		double M = 0.5;
		
		double rain_rate = 16, theta = 0, tau = 0;
				
		double fsl = calcFreeSpaceLoss(f, d);
		double ca = CloudsAttenuationCalculator.calculate(f, T, M);
		double ra = RainAttenuationCalculator.calculate(f, rain_rate, theta, tau);
		
		//double total_attenuation = calcTotalAttenuation(f, d, T, M, rain_rate, theta, gain_tx, gain_rx);
		
		System.out.println("Free-space attenuation (dB): " + fsl);
		System.out.println("Cloud attenuation (dB/km): " + ca);
		System.out.println("Rain attenuation (dB/km): " + ra);
		//System.out.println("Total attenuation (dB): " + total_attenuation);

	}	
}
