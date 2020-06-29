package main.p528;

import java.util.Arrays;
import java.util.Comparator;

/** Rec. ITU-R P.528-4 Annex II Section VI
 * 
 * @author Collin Thornton
 *
 */
public class LOSAttenuationCalculator {
	private class Tuple {
		double psi, deltar, d;
		
		Tuple(double psi, double deltar, double d) {
			this.psi = psi;
			this.deltar = deltar;
			this.d = d;
		}
	}
	private class SortByDeltaR implements Comparator<Tuple> {
		@Override
		public int compare(Tuple o1, Tuple o2) {
			if(o1.deltar > o2.deltar) 		return 1;
			else if(o1.deltar == o2.deltar) return 0;
			else 							return -1;
		}
		
	}
	
	//private double A;	// Basic transmission loss (dB)
	//private double K;	// variability coefficient 

	private static double[] rtab = { 0.06, 0.1, 1/9, 1/8, 1/7, 1/6, 1/5, 1/4, 1/3, 1/2 };
	private static double[] psitab = { 0.2, 0.5, 0.7, 1.0, 1.2, 1.5, 1.7, 2.0, 2.5, 3.0, 3.5, 
			4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 20.0, 45.30, 70.0, 80.0, 85.0, 88.0, 89.0 };
	
	public Tuple[] table = new Tuple[46];
	
	public double compute(TerminalGeometry.Geom lt, TerminalGeometry.Geom ht, Path path) {
		double lambda = 0.2997925/path.input.f;
		buildTable(lambda, path, lt, ht);
		
		double d_halflambda		= extrapolateDeltaRtoD(lambda/2);
		double psi_limit 		= extrapolateDtoPsi(d_halflambda);
		double d_sixthlambda	= extrapolateDeltaRtoD(lambda/6);
		
		
		// <<Step 6>> Determine path.d_0 
		if(path.input.d >= path.d_d || path.d_d >= path.d_ML) {
			if(path.input.d > d_sixthlambda || d_sixthlambda > path.d_ML) path.d_0 = lt.d;
			else path.d_0 = d_sixthlambda;
		}
		else if(path.d_d < d_sixthlambda && d_sixthlambda < path.d_ML) path.d_0 = d_sixthlambda;
		else path.d_0 = path.d_d;
		
		// <<Step 7>> Tune path.d_0 
		double temp_d = path.d_0 - 0.001;
		RayOptics.RayOpticsData optics;
		double psi;
		do {
			temp_d += 0.001;
			psi = extrapolateDtoPsi(temp_d);
			optics = RayOptics.compute(psi, path, lt, ht);
		} while(optics.d < path.d_0 && optics.d <= path.d_ML);
		path.d_0 = optics.d;
			
		
		// <<Step 8>> Compute line-of-sight loss at path.d_0 
		double psi_d0 = extrapolateDtoPsi(path.d_0);
		optics = RayOptics.compute(psi_d0, path, lt, ht);
		path.atten.A_d0 = LOSPathLoss.compute(psi_d0, psi_limit, optics, path);
		
		
		// <<Step 9>> Tune psi
		psi = extrapolateDtoPsi(path.input.d);
		optics = RayOptics.compute(psi, path, lt, ht);
		
		double delta = 0.01, error = optics.d - path.input.d;
		int los_iterations = 0;
		while(Math.abs(error) > 0.001 && los_iterations < Constants.LOS_ITERATIONS) {
			if(error > 0) {
				psi += delta;
				delta /= 2;
				psi -= delta;
			}
			else {
				psi -= delta;
			}
			optics = RayOptics.compute(psi, path, lt, ht);
			error = optics.d - path.input.d;
			++los_iterations;
		}
		//double d_r0 = optics.d;
		
		
		// <<Step 10> Compute line-of-sight loss
		path.atten.A_LOS = LOSPathLoss.compute(psi, psi_limit, optics, path);
		
		
		// <<Step 11> Compute atmospheric attenuation
		double r_eo = EffectiveRayLength.compute(Constants.A_E, Constants.T_EO, optics);
		double r_ew = EffectiveRayLength.compute(Constants.A_E, Constants.T_OW, optics);
		
		AtmosphericAbsorptionRates.AbsorptionData abs = AtmosphericAbsorptionRates.compute(path);
		path.atten.A_a = -abs.gamma_oo*r_eo - abs.gamma_ow*r_ew;
		
		
		// <<Step 12>> Compute free space loss
		path.atten.A_fs = computeFPL(optics, path, lt, ht);
		
		
		// <<Step 13>> Compute contribution of variability to loss
		double Y_q = VariabilityLoss.compute(path.atten.A_LOS, lt, ht, path, optics, r_ew);
		
		
		// <<Step 14>> Sum the components
		path.atten.A = path.atten.A_fs + path.atten.A_a + path.atten.A_LOS + Y_q;	
		
		return path.atten.A;
	}
	
	
	private double computeFPL(RayOptics.RayOpticsData optics, Path path, TerminalGeometry.Geom lt, TerminalGeometry.Geom ht) {
		
		//! SOMETHING IS INCORRECT WITH THETA_FS... EITHER optics.a_a (which would imply an error with psi) OR terminal.theta IS BROKEN
		//double theta_fs = (optics.a_a*(lt.theta+ht.theta))/Constants.A_0;
		double z_1 = Constants.A_0 + path.input.h_r1;
		double z_2 = Constants.A_0 + path.input.h_r1;
		
		//! This accounts for curvature of earth and height differences. Is more accurate
		//double r_fs = Math.max(Math.sqrt(Math.pow(z_2-z_1, 2) + 4*z_1*z_2*Math.pow(Math.sin(0.5*theta_fs),2)), z_2-z_1);
		
		
		//double A_fs = -32.45 - 20*Math.log10(path.input.f*r_fs);
		double A_fs = -32.45 - 20*Math.log10(path.input.f*path.input.d);
		
		// double theta_should_be = 2*Math.asin(Math.sqrt((36-Math.pow(z_2-z_1, 2))/(4*z_1*z_2)));
		
		return A_fs;
	}
	private void buildTable(double lambda, Path path, TerminalGeometry.Geom lt, TerminalGeometry.Geom ht) {
		table[0] = new Tuple(0.00, 0.00, lt.d + ht.d);

		for(int i=1; i<45; ++i) {
			double psi;
			if(i<11) 		psi = Math.asin((lambda*rtab[i-1])/(2*lt.h));
			else if(i < 21) psi = Math.asin((lambda*rtab[i-11])/(2*lt.d));
			else 			psi = psitab[i-21]*Math.PI/180;
			
			RayOptics.RayOpticsData optics = RayOptics.compute(psi, path, lt, ht);

			table[i] = new Tuple(psi, optics.deltar, optics.d);
		}
		table[45] = new Tuple(Math.PI/2, 2*lt.h, 0.00);
		
		Arrays.sort(table, new SortByDeltaR());		
	}
	
	private double extrapolateDeltaRtoD(double deltaR) {
		if(deltaR <= table[0].deltar) return table[0].d;
	
		int i=1;
		while(table[i].deltar < deltaR && i < 44) ++i;
		
		if(deltaR < table[i].deltar) return ((table[i].d-table[i-1].d)*(deltaR-table[i-1].deltar))/(table[i].deltar-table[i-1].deltar) + table[i-1].d;
		else if (deltaR == table[i].deltar) return table[i].d; 
		return table[45].d;
	}
	private double extrapolateDtoPsi(double D) {
		if(D >= table[0].d) return table[0].psi;
		
		int i=1;
		while(table[i].d > D && i < 44) ++i;
		
		if(D > table[i].d) return ((table[i].psi-table[i-1].psi)*(D-table[i-1].d))/(table[i].d-table[i-1].d) + table[i-1].psi;
		else if (D == table[i].d) return table[i].psi; 
		return table[45].psi;
	}	
	
	
	
	
	public static void main(String[] args) {	
		double f = 5000;  // MHz
		double lambda = 0.2997925/f;
		double q = 0.99;
		double d = 6;
		
		double h_r1 = 0.36;
		double h_r2 = 0.39;
				
		Path path = new Path(h_r1, h_r2, f, q, d);

		TerminalGeometry low_term = new TerminalGeometry(h_r1);
		TerminalGeometry high_term = new TerminalGeometry(h_r2);
		TerminalGeometry.Geom lt = low_term.geom;
		TerminalGeometry.Geom ht = high_term.geom;
		
		LOSAttenuationCalculator los = new LOSAttenuationCalculator();
		los.buildTable(lambda, path, lt, ht);	

		path.d_ML = lt.h + ht.h;
							
		los.compute(lt, ht, path);
		
		for(int i=0; i<46; ++i) {
			System.out.format("psi: %4.3e, dr: %4.3e, d: %4.3e%n", los.table[i].psi, los.table[i].deltar, los.table[i].d);
		}
		System.out.format("%4.3e%n", lambda);
		System.out.println(" ");
		
		double _d = los.extrapolateDeltaRtoD(lambda/2);
		double psi = los.extrapolateDtoPsi(_d);
		System.out.format("dr:\t%4.3e%n",	lambda/2);
		System.out.format("d:\t%4.3e%n",	_d);
		System.out.format("psi:\t%4.3e%n",  psi);
	}
}