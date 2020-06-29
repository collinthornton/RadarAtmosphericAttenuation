package main.p528;


/** Rec. ITU-R P.528-4 Annex II Section XIV
 * 
 * @author Collin Thornton
 *
 */
public class AtmosphericAbsorptionRates {
	public static class AbsorptionData {
		public double gamma_oo;	// Oxygen absorption rate 		(dB/km)	
		public double gamma_ow;	// Water vapor absorption rate 	(dB/km)
	}
	
	private static double[][] table = {
			{ 100, 0.00019, 0 },
			{ 150, 0.00042, 0 },
			{ 205, 0.00070, 0 },
			{ 300, 0.00096, 0 },
			{ 325, 0.00130, 0 },
			{ 350, 0.00150, 0 },
			{ 400, 0.00180, 0 },
			{ 550, 0.00250, 0 },
			{ 700, 0.00300, 0 },
			{ 1000, 0.0042, 0 },
			{ 1520, 0.0050, 0 },
			{ 2000, 0.0070, 0 },
			{ 3000, 0.0088, 0 },
			{ 3400, 0.0092, 0.0001  },
			{ 4000, 0.0100, 0.00017 },
			{ 4900, 0.0110, 0.00340 },
			{ 8300, 0.0140, 0.00210 },
			{ 10200, 0.015, 0.00900 },
			{ 15000, 0.017, 0.02500 },
			{ 17000, 0.018, 0.04500 }
	};
	
	
	public static AbsorptionData compute(Path path) {
		AbsorptionData output = new AbsorptionData();
		// <<Step 1>> Find the correct rows of in the table
		int i = 1;
		while(i < 18 && path.input.f >= table[i][0]) { ++i; }
		double f_prime 	= table[i-1][0];
		double f_pprime = table[i][0];
		double gamma_oo_prime = table[i-1][1];
		double gamma_oo_pprime = table[i][1];
		double gamma_ow_prime = table[i-1][2];
		double gamma_ow_pprime = table[i][2];
		
		
		// <<Step 2>> Compute the interpolation scale factor
		double R = (Math.log10(path.input.f)-Math.log10(f_prime)) / (Math.log10(f_pprime)-Math.log10(f_prime));
		
		// <<Step 3>> Interpolate y_oo
		double X = R*(Math.log10(gamma_oo_pprime) - Math.log10(gamma_oo_prime)) + Math.log10(gamma_oo_prime);
		output.gamma_oo = Math.pow(10, X);
		
		// <<Step 4>> INterpolate y_ow
		if(path.input.f >= 3400) {
			double Y = R*(Math.log10(gamma_ow_pprime) - Math.log10(gamma_ow_prime)) + Math.log10(gamma_ow_prime);
			output.gamma_ow = Math.pow(10, Y);
			return output;
		}
		output.gamma_ow = 0;
		
		return output;
	}
}
