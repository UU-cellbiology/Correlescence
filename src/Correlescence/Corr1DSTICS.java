package Correlescence;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.plugin.PlugIn;

public class Corr1DSTICS implements PlugIn {

	/** original stack **/
	ImagePlus imp;
	/** correlation calculation method 
	 * 0 = direct correlation
	 * 1 = correlation using FFT **/
	int nCalcMethod;
	/** correlation normalization method 
	 * 0 = overlap area
	 * 1 = full area **/
	int nNormMethod;
	/** object calculating correlation **/
	imCC1D x1D = new imCC1D();
	/** maximum delay in frames 1D**/
	int nMaxUserDelay1D;
	
	
	@Override
	public void run(String arg) 
	{
		int nStackSize;
		int imHeight;
		String sTitle;
		ImagePlus finalImp;
	
		imp = IJ.getImage();		
		imHeight= imp.getHeight();
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
		if(nStackSize >1)
		{
		    IJ.log("Just one image required, aborted.");
		    	return;
		}
		
		if(!xDialog1STICS())
			return;
		
		sTitle=imp.getTitle();
		
		if(nMaxUserDelay1D==0)
			nMaxUserDelay1D = (int)(imHeight*0.5);
		else
		{
			if(nMaxUserDelay1D> imHeight)
			{			
				nMaxUserDelay1D=imHeight;
				IJ.log("Maximum delay value is too high, reducing it to "+String.format(" %d frames", nMaxUserDelay1D));
			}
		}
		//compose title
		if(nCalcMethod==0)
			sTitle=sTitle+"_1Dxt_direct";
		else
			sTitle=sTitle+"_1Dxt_fft";
		
		finalImp = new ImagePlus(sTitle,x1D.xCorrSpaceTime(imp.getProcessor(), nMaxUserDelay1D,nCalcMethod,nNormMethod));
		Calibration cal = finalImp.getCalibration();
	    cal.xOrigin=Math.round(imp.getWidth()*0.5)-1;
	    cal.yOrigin=0;
		finalImp.show();
			//x1D.xCorrSpaceTimeFFT();
	}
	
	/** 
	 * Dialog displaying options of 1D correlation
	 * **/
	public boolean xDialog1STICS() 
	{
		GenericDialog x1DDial = new GenericDialog("Spatio-temporal image correlation 1D");

		String [] itemsC = new String [] {
				"Direct cross-correlation (slow)","FFT cross-correlation (fast)"};
		String [] itemsN = new String [] {
				"Normalize by overlap area","Normalize by full area"};

		x1DDial.addNumericField("Max delay (zero=half of height)", Prefs.get("Correlescence.nMaxUserDelay1D", 0), 0, 4, " pixels");
		x1DDial.addChoice("Calculation method:", itemsC, Prefs.get("Correlescence.1Dcorrmethod", "FFT cross-correlation (fast)"));
		x1DDial.addChoice("Normalization:", itemsN, Prefs.get("Correlescence.1Dnormmethod", "Normalize by overlap area"));
		x1DDial.setResizable(false);
		x1DDial.showDialog();
		if (x1DDial.wasCanceled())
	        return false;

		nMaxUserDelay1D = (int)x1DDial.getNextNumber();
		Prefs.set("Correlescence.nMaxUserDelay1D", nMaxUserDelay1D);
		nCalcMethod = x1DDial.getNextChoiceIndex();
		Prefs.set("Correlescence.1Dcorrmethod", itemsC[nCalcMethod]);
		nNormMethod= x1DDial.getNextChoiceIndex();
		Prefs.set("Correlescence.1Dnormmethod", itemsN[nNormMethod]);

		return true;
	}
}
