package tests;

import correlescence.Correlescence2D;
import ij.IJ;
import ij.ImageJ;

public class Corr2DTest {
	public static void main( final String[] args )
	{
		new ImageJ();
		IJ.open("/home/eugene/Desktop/cc_test.tif");
		Correlescence2D cc = new Correlescence2D();
		cc.run(null);
	}
}
