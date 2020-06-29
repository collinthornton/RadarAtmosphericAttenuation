package main.p838;

/** Calculate loss in dB/km due to rain
 *
 * @author Collin Thornton
 * @note This class is based on Rec. ITU-R P.838-3, full paper
 */
public class RainAttenuationCalculator {

	private static final double[] kh_aj = { -5.33980, -0.35351, -0.23789, -0.94158 };
	private static final double[] kh_bj = { -0.10008, 1.26970, 0.86036, 0.64552 };
	private static final double[] kh_cj = { 1.13098, 0.45400, 0.15354, 0.16817 };
	private static final double kh_mk = -0.18961;
	private static final double kh_ck = 0.71147;
		
	private static final double[] ah_aj = { -0.14318, 0.29591, 0.32177, -5.37610, 16.1721 };
	private static final double[] ah_bj = { 1.82442, 0.77564, 0.63773, -0.96230, -3.29980 };
	private static final double[] ah_cj = { -0.55187, 0.19822, 0.13164, 1.47828, 3.43990 };
	private static final double ah_ma = 0.67849;
	private static final double ah_ca = -1.95537;
		
	private static final double[] kv_aj = { -3.80595, -3.44965, -0.39902, 0.50167 };
	private static final double[] kv_bj = { 0.56934, -0.22911, 0.73042, 1.07319 };
	private static final double[] kv_cj = { 0.81061, 0.51059, 0.11899, 0.27195 };
	private static final double kv_mk = -0.16398;
	private static final double kv_ck = 0.63297;
	
	private static final double[] av_aj = { -0.07771, 0.56727, -0.20238, -48.2991, 48.5883 };
	private static final double[] av_bj = { 2.33840, 0.95545, 1.14520, 0.791669, 0.791459 };
	private static final double[] av_cj = { -0.76284, 0.54039, 0.26809, 0.116226, 0.116479 };
	private static final double av_ma = -0.053739;
	private static final double av_ca = 0.83433;
		
	/** Calculate loss in dB/km due to rain
	 * 
	 * @param f Frequency (GHz) 1 <= f <= 1000
	 * @param rain_rate Rate of rain (mm/hr)
	 * @param theta Path elevation angle (deg)
	 * @param tau Polarization tilt angle (deg). 0 = horizontal. 90 = vertical. 45 = circular
	 * @return dB/km due to rain. Doubled for 2-way loss
	 * @note This function is based on Rec. ITU-R P.838-3, full paper
	 */	
	public static double calculate(double f, double rain_rate, double theta, double tau) {
		if(f < 1) throw new IllegalArgumentException("Rain atten. freq must be between 1 and 1000 GHz");
		if(f > 1000) throw new IllegalArgumentException("Rain atten. freq must be between 1 and 1000 GHz");
		
		double log_kh = 0, log_kv = 0;
		double ah = 0, av = 0;
		
		for(int i=0; i<5; ++i) {
			if(i<4) {
				log_kh += kh_aj[i]*Math.exp(-1*Math.pow((Math.log10(f)-kh_bj[i])/kh_cj[i], 2));
				log_kv += kv_aj[i]*Math.exp(-1*Math.pow((Math.log10(f)-kv_bj[i])/kv_cj[i], 2));
			}
			ah += ah_aj[i]*Math.exp(-1*Math.pow((Math.log10(f)-ah_bj[i])/ah_cj[i], 2));
			av += av_aj[i]*Math.exp(-1*Math.pow((Math.log10(f)-av_bj[i])/av_cj[i], 2));
		}
		log_kh += kh_mk*Math.log10(f)+kh_ck;
		log_kv += kv_mk*Math.log10(f)+kv_ck;
		ah += ah_ma*Math.log10(f)+ah_ca;
		av += av_ma*Math.log10(f)+av_ca;
		
		double kh = Math.pow(10, log_kh);
		double kv = Math.pow(10, log_kv);
		
		theta = Math.toRadians(theta);
		tau = Math.toRadians(tau);
		
		double k = (kh + kv + (kh-kv)*Math.pow(Math.cos(theta), 2) * Math.cos(2*tau)) / 2;		
		double a = (kh*ah + kv*av + (kh*ah - kv*av)*Math.pow(Math.cos(theta), 2)*Math.cos(2*tau)) / (2*k);
		
		return k*Math.pow(rain_rate, a);			
	}
}