package Correlescence;



import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.FHT;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;



public class FFTtoolsUI implements PlugIn 
{
	ImagePlus imp;
	boolean bPSInverseFFT;
	int nChoiceFFT2;
	
	@Override
	public void run(String arg) {
		
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
		
		if(arg.equals("periodicsmooth"))
		{
		
			if(!xDialogPS())
				return;
			startTime = System.nanoTime();
			ImageProcessor dupip = (FloatProcessor) ip.duplicate().convertToFloat();
			
			ImageStack fin = FFTtools.PeriodicSmooth(dupip, bPSInverseFFT);
			if(bPSInverseFFT)
				{new ImagePlus("ComplexFFT_smooth_"+imTitle,fin).show();}
			else
				{new ImagePlus("PeriodicNSmooth_"+imTitle,fin).show();}
		}
		if(arg.equals("2DFFT"))
		{
			if(!xDialogFFT2D())
				return;
			startTime = System.nanoTime();
			if(nChoiceFFT2==0)
			{ 
				new ImagePlus("Complex_FFT_of_"+imTitle,GeneralFFT.fft2Dtransform(ip, false)).show();
				/*FHT fht;
				ImageProcessor paddedip1=FFTtools.padzeros(ip);
				fht = new FHT(paddedip1);
				fht.transform();
				fht.setShowProgress(false);
				new ImagePlus("Complex_FFT_of_"+imTitle,fht.getComplexTransform()).show();
				IJ.log("FHT transform");*/
			}
			if(nChoiceFFT2==1)
			{
				if(imp.getStackSize()!=2)
				{
					IJ.log("Stack with two frames(slices) is required as an input.\nFirst image is real and second image is imaginary part of FFT.");
			    	return;
					
				}
				else
				{
					
					ImageStack stack_in = new ImageStack(ip.getWidth(),ip.getHeight());
					FloatProcessor real = (FloatProcessor)imp.getStack().getProcessor(1).duplicate().convertToFloat();
					GeneralFFT.GFFTswapQuadrantsBack(real);
					stack_in.addSlice(real);
					FloatProcessor imag = (FloatProcessor)imp.getStack().getProcessor(2).duplicate().convertToFloat();
					GeneralFFT.GFFTswapQuadrantsBack(imag);
					stack_in.addSlice(imag);
					new ImagePlus("InverseFFT_of_"+imTitle,GeneralFFT.fft2Dinverse(stack_in)).show();
				}
			}
			if(nChoiceFFT2==2)
			{
				startTime = System.nanoTime();
				new ImagePlus("PowerSpectrumNPhase_of_"+imTitle,GeneralFFT.fft2Dpowerphase(GeneralFFT.fft2Dtransform(ip, true))).show();
			}
		}
		processingTime = System.nanoTime() - startTime;
		IJ.log("Processing time: " + String.format("%.2f",((double)Math.abs(processingTime))*0.000000001) + " s");
	}
	/** 
	 * Dialog displaying options of periodic/smooth decomposition
	 * **/
	public boolean xDialogPS() 
	{
		int nChoice;
		GenericDialog xDialPS = new GenericDialog("Periodic-plus-smooth Decomposition");

		String [] itemsPS = new String [] {
				"Periodic-and-smooth decomposition","Complex FFT of smooth part"};

		xDialPS.addChoice("Return:", itemsPS, Prefs.get("Correlescence.PSdecompose", "Periodic-and-smooth decomposition"));

		xDialPS.setResizable(false);
		xDialPS.showDialog();
		if (xDialPS.wasCanceled())
	        return false;

		nChoice = xDialPS.getNextChoiceIndex();
		Prefs.set("Correlescence.PSdecompose", itemsPS[nChoice]);
		if (nChoice==0)
			{bPSInverseFFT=false;}
		else
			{bPSInverseFFT=true;}


		return true;
	}

	public boolean xDialogFFT2D() 
	{
		
		GenericDialog xDialFFT2 = new GenericDialog("2D FFT");

		String [] itemsFFT2 = new String [] {
				"Forward FFT","Inverse FFT", "Power spectrum and phase"};
		xDialFFT2.addMessage("2D FFT of image with non-square/non power of 2 dimensions");
		xDialFFT2.addChoice("Perform:", itemsFFT2, Prefs.get("Correlescence.FFT2", "Forward FFT"));

		xDialFFT2.setResizable(false);
		xDialFFT2.showDialog();
		if (xDialFFT2.wasCanceled())
	        return false;

		nChoiceFFT2 = xDialFFT2.getNextChoiceIndex();
		Prefs.set("Correlescence.FFT2", itemsFFT2[nChoiceFFT2]);



		return true;
	}

}
