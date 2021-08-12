package main.p528;

public class Constants {
	private Constants() { }
	
	/** surface refractivity 	(301 N-Units)
	 */
	public static double N_S = 301.0;
	
	/** actual Earth radius (6370 km)
	 */
	public static double A_0 = 6370.0;
	
	/** effective Earth radius (8493 km)
	 */
	public static double A_E = 8493.0;	
	
	/** effective thickness of oxygen absorbing layer (3.25 km)
	 */
	public static double T_EO = 3.25;
	
	/** effective thickness of water vapor absorbing layer (1.36 km)
	 */
	public static double T_OW = 1.36;
	
	/** relative dielectric constant (15)
	 */
	public static double E_R = 15.0;		
	
	/** conductivity (0.005 S/m)
	 */
	public static double SIGMA = 0.005;	
	
	public static int LOS_ITERATIONS = 25;
}
