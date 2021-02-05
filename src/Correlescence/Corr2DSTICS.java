package Correlescence;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.measure.Measurements;
import ij.plugin.PlugIn;
import ij.plugin.ZProjector;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;

public class Corr2DSTICS implements PlugIn {

	/** original stack **/
	ImagePlus imp;
	/** correlation calculation method **/
	int nCalcMethod;
	/** object calculating correlation **/
	ImCrossCorrelation x2D;
	/** object performing averaging**/
	ZProjector Zproj = new ZProjector();

	/** maximum delay in frames **/
	int nMaxUserDelay;
	/** image width in px**/
	int origW;
	/** image height in px**/
	int origH;
@Override
public void run(String arg) {
	

	
	int nStackSize;
	int i, nMaxDelay;
	ImageStack stSTICS;
	String sTitle;
	int nHalfStack;
	ImagePlus finalImp;

	
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
	
	nStackSize = imp.getStackSize();
	if(nStackSize ==1)
	{
	    IJ.log("Stack of at least 2 images is required, aborted.");
	    	return;
	}
	
	if(!xDialog2STICS())
		return;
	
	
	origW = imp.getWidth();
	origH =imp.getHeight();

	
    x2D = new ImCrossCorrelation();
    nHalfStack=(int)Math.floor((double)nStackSize*0.5);
	if(nMaxUserDelay==0)
		nMaxDelay = nHalfStack;
	else
	{
		if(nMaxUserDelay> nHalfStack)
		{
			
			nMaxDelay=nHalfStack;
			IJ.log("Maximum delay value is too high, reducing it to "+String.format(" %d frames", nHalfStack));
		}
		else
			nMaxDelay=nMaxUserDelay;
	}
	//stSTICS = new ImageStack(nCorrW,nCorrW);
	stSTICS = new ImageStack(origW,origH);
	
	for (i=0;i<=nMaxDelay;i++)
	{
		IJ.showStatus(String.format("Calculating average correlation for delay of %d frames..", i));
		stSTICS.addSlice(String.format("Delay %d frames", i), average2Dcorr(i,nStackSize));
		
	}
	sTitle=imp.getTitle();
	IJ.showStatus("2D STICS calculation done.");
	finalImp=new ImagePlus(sTitle+"_2D_STICS", stSTICS);
	//calibration
	Calibration cal = finalImp.getCalibration();
    cal.xOrigin=Math.round(origW*0.5);
    cal.yOrigin=Math.round(origH*0.5);
    
    finalImp.setCalibration(cal);
    finalImp.show();
	
}//end of run()

/**
 *  Wrapper calculating average 2D spatial correlation with specified delay in frames 
 */
ImageProcessor average2Dcorr(int nDelay, int nStackSize)
{
	ImageProcessor finalip, ip1,ip2;
	ImageStack crosscorrstack;
	int i,j;
	
	i=1; j=i+nDelay;
	crosscorrstack = new ImageStack(origW,origH);
	
	while(j<=nStackSize)
    {
		imp.setSliceWithoutUpdate(i);
    	ip1 = getFloatversion(imp);
    	imp.setSliceWithoutUpdate(j);
    	ip2 = getFloatversion(imp);

		//subtract average intensity value
    	//in general, something else can be here, like background value
		ip1.subtract(ImageStatistics.getStatistics(ip1,Measurements.MEAN,null).mean);
		ip2.subtract(ImageStatistics.getStatistics(ip2,Measurements.MEAN,null).mean);
				 
			
		if(nCalcMethod==0)
			{finalip = x2D.calcDirectCorrelationImage(ip1,ip2);}
		else
			{finalip = x2D.calcFFTCorrelationImage(ip1, ip2);}
				
		crosscorrstack.addSlice(null, finalip);
		i++; j++;
    }
	//do average projection
	Zproj.setImage(new ImagePlus("", crosscorrstack));
	Zproj.setMethod(ZProjector.AVG_METHOD);
	Zproj.doProjection();
	finalip=Zproj.getProjection().getProcessor();
	
	return finalip;
	
}

/** 
 * Dialog displaying options of 2D correlation
 * **/
public boolean xDialog2STICS() 
{
	GenericDialog x2DDial = new GenericDialog("Spatio-temporal image correlation 2D");

	String [] itemsC = new String [] {
			"Direct cross-correlation (slow)","FFT cross-correlation (fast)"};

	x2DDial.addNumericField("Max delay (zero=half of stack)", Prefs.get("Correlescence.2DSTICSMaxDelay", 0), 0, 4, " frames");
	x2DDial.addChoice("Calculation method:", itemsC, Prefs.get("Correlescence.2Dcorrmethod", "FFT cross-correlation (fast)"));
	
	x2DDial.setResizable(false);
	x2DDial.showDialog();
	if (x2DDial.wasCanceled())
        return false;

	nMaxUserDelay = (int)x2DDial.getNextNumber();
	Prefs.set("Correlescence.2DSTICSMaxDelay", nMaxUserDelay);
	nCalcMethod = x2DDial.getNextChoiceIndex();
	Prefs.set("Correlescence.2Dcorrmethod", itemsC[nCalcMethod]);

	return true;
}

/** Gets float copy of current imageprocessor */	
ImageProcessor getFloatversion(ImagePlus impx)

{
	if(imp.getType()!=ImagePlus.GRAY32)
	{
		return impx.getProcessor().convertToFloat();
	}else
	{
		return  impx.getProcessor().duplicate();
	}
}


}
