package correlescence;

public class Particle1D {
	/**particle coordinate**/
	public double x;
	/** particle state 
	 * 0=pause, 
	 * 1=moving with speed**/
	public int nState;
	/**particle speed**/
	public double v;
	/** particle intensity **/
	public double dInt;
	/**particle state timer**/
	public double dStateTimer;
	/**particle lifetime timer**/
	public double dLTTimer;
	/** whether particle is rendered **/
	public boolean bVisible;
	

}
