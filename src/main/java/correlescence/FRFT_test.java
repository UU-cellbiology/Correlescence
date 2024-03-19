package correlescence;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

public class FRFT_test implements PlugIn 
{
	ImagePlus imp;

	@Override
	public void run(String arg) {
		// TODO Auto-generated method stub
		ImageProcessor ip;
		long startTime=0;
		long processingTime=0;
		// TODO Auto-generated method stub
		imp = IJ.getImage();		
		if (imp==null)
		{
		    IJ.noImage();
		    return;
		}		
		else if (imp.getType() != ImagePlus.GRAY8 && imp.getType() != ImagePlus.GRAY16 && imp.getType() != ImagePlus.GRAY32) 
		{
		    IJ.error("8, 16 or 32 bit greyscale image required");
		    return;
		}


		
		String imTitle = imp.getTitle();
		ip = imp.getProcessor();
		//float [][] frft;
		int originalWidth = ip.getWidth();
		startTime = System.nanoTime();
		//float [] realin=FFTtools.getImageRow(ip, 0, originalWidth);
		//frft=GeneralFFT.transformBluesteinFracCenterReal(realin, 0.33333333333333333333f);
		ImageStack ppdft = PolarDFT.PPDFT(ip,1,1);
		new ImagePlus("power_phase_",GeneralFFT.fft2Dpowerphase(ppdft)).show();
		new ImagePlus("complex_",ppdft).show();
		processingTime = System.nanoTime() - startTime;
		IJ.log("Processing time: " + String.format("%.2f",((double)Math.abs(processingTime))*0.000000001) + " s");
//		originalWidth = ip.getWidth();
	}

}
