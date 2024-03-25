package tests;

import correlescence.Correlescence2D;
import correlescence.Temporal_ICS;
import ij.IJ;
import ij.ImageJ;

public class TICSTest {
	public static void main( final String[] args )
	{
		new ImageJ();
		IJ.open("/home/eugene/Desktop/cc_test.tif");
		Temporal_ICS tics = new Temporal_ICS();
		tics.run(null);
	}
}
