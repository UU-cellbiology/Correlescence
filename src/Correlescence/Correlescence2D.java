package Correlescence;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;

public class Correlescence2D implements PlugIn {

	ImagePlus imp;
	
	/** delay between frames **/
	int nImNumber;
	/** way to calc correlation: consecutive or one fixed**/
	String sChoice;
	/** whether correct for the drift or not **/
	boolean bDrift;
	/** limit displacement for x and y**/
	boolean bDriftlimit;
	/** max displacement for x**/
	int nDriftlimitX;
	/** max displacement for y**/
	int nDriftlimitY;
	/** method to calculate xcorr: slow or fast **/
	int nCalcMethod;
	
	
	ResultsTable ptable;

@Override
public void run(String arg) {
	
	int nStackSize;
	ImageProcessor ip, ip1,ip2;
	int originalWidth;
    int originalHeight;

    int maxN;
	int nCorrW;
    int i,j;
    int [] xymax;
    String sTitle;
    String sImTitle;
    String sInTitle;
    ImageStack crosscorrstack;
    int [][] xydrifttable;
    boolean bAutoCorr =false;
    int nCurrentSlice;
	ImagePlus impCC;
	
    
	ImCrossCorrelation x2D = new ImCrossCorrelation();
	

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
	nCurrentSlice = imp.getCurrentSlice();
	
	if(nStackSize ==1)
	{
	    IJ.log("Only one image in stack, limiting choice to autocorrelation");
	    	bAutoCorr=true;
	}
	sInTitle = imp.getTitle();
	
	if(!x2Dialog(bAutoCorr))
		return;
	
	originalWidth = imp.getWidth();
    originalHeight = imp.getHeight();
	//halfW=(int)Math.round(originalWidth*0.5);
	//halfH=(int)Math.round(originalHeight*0.5);
	maxN = Math.max(originalWidth, originalHeight);
	nCorrW = 2;
    while(nCorrW<maxN) nCorrW *= 2;
    
    //crosscorrstack= new ImageStack(nCorrW,nCorrW);
    crosscorrstack= new ImageStack(originalWidth,originalHeight);
    ptable = ResultsTable.getResultsTable();
    
    xydrifttable = new int [nStackSize][2];
    sImTitle = "xcorr_";
    
    IJ.showStatus("Calculating correlation...");
    if(sChoice.equals("consecutive images"))
    {
    	sImTitle =sImTitle +sChoice;
    	sImTitle = sImTitle+"_delay_"+Integer.toString(nImNumber)+"_"+sInTitle;
    	if (bDrift)
    	{
    		xydrifttable[0][0]=0;
    		xydrifttable[0][1]=0;
    	}
    	ptable.reset();
	    i=1; j=i+nImNumber;
	    if(j>nStackSize)
	    {
	    	IJ.error("Frame delay larger than image stack size");
	    	return;
	    }
	    
	    
	  
	    while(j<=nStackSize)
	    {
	    	imp.setSliceWithoutUpdate(i);

	    	ip1 = getFloatversion(imp);
	    	imp.setSliceWithoutUpdate(j);
	    	ip2 = getFloatversion(imp);

			//subtract average intensity value
			ip1.subtract(ImageStatistics.getStatistics(ip1,Measurements.MEAN,null).mean);
			ip2.subtract(ImageStatistics.getStatistics(ip2,Measurements.MEAN,null).mean);
			 
		
			if(nCalcMethod==0)
				{ip = x2D.calcDirectCorrelationImage(ip1,ip2);}
			else
				{ip = x2D.calcFFTCorrelationImage(ip1, ip2);}
			
			sTitle = String.format("corr_%d_x_%d", i,j);
			//xymax = getmaxpositions(ip);
			xymax = getmaxpositionscenterlimit(ip);
			ptable.incrementCounter();
			ptable.addLabel(sTitle);
			ptable.addValue("Xmax_(px)",xymax[0]);	
			ptable.addValue("Ymax_(px)",xymax[1]);
			if (bDrift)
			{
				xydrifttable[i][0]=(int) (xymax[0] + xydrifttable[i-1][0]);
				xydrifttable[i][1]=(int) (xymax[1] + xydrifttable[i-1][1]);
				//reserve
				ptable.addValue("Xdrift_(px)",0);	
				ptable.addValue("Ydrift_(px)",0);
			}
			
			crosscorrstack.addSlice(null, ip);
			IJ.showProgress(j, nStackSize);
			i++; j++;
		   
	    }
	    IJ.showProgress(j, nStackSize);
    }
    //choice:
    //current image in stack and all others
    if(sChoice.equals("current image in stack and all others"))
    {
    	
    	i=imp.getCurrentSlice();
    	sImTitle =sImTitle +"vs_frame"+Integer.toString(i)+"_"+sInTitle;
    	ip1=getFloatversion(imp);
    	//subtract average intensity value
		ip1.subtract(ImageStatistics.getStatistics(ip1,Measurements.MEAN,null).mean);
    	j=1;
    	while(j<=nStackSize)
    	{
			imp.setSliceWithoutUpdate(j);
			ip2 = getFloatversion(imp);
			ip2.subtract(ImageStatistics.getStatistics(ip2,Measurements.MEAN,null).mean);
			if(nCalcMethod==0)
				{ip = x2D.calcDirectCorrelationImage(ip1,ip2);}
			else
				{ip = x2D.calcFFTCorrelationImage(ip1, ip2);}
		
			sTitle = String.format("corr_%d_x_%d", i,j);
			xymax = getmaxpositionscenterlimit(ip);
			//xymax = getmaxpositions(ip);
			ptable.incrementCounter();
			ptable.addLabel(sTitle);
			ptable.addValue("Xmax_(px)",xymax[0]);	
			ptable.addValue("Ymax_(px)",xymax[1]);
			if (bDrift)
			{
				xydrifttable[j-1][0]=(int) (xymax[0]);
				xydrifttable[j-1][1]=(int) (xymax[1]);
				//reserve
				ptable.addValue("Xdrift_(px)",0);	
				ptable.addValue("Ydrift_(px)",0);
			}

			
			crosscorrstack.addSlice(null, ip);
			IJ.showProgress(j, nStackSize);
			j++;
    		
    	}
    	IJ.showProgress(j, nStackSize);
    }
    //autocorrelation
    if(sChoice.equals("autocorrelation"))
    {	
    	sImTitle =sImTitle +"autocorrelation_"+sInTitle;
    
		//i=imp.getCurrentSlice();
		//ip1=getFloatversion(imp);
		//subtract average intensity value
		//ip1.subtract(ImageStatistics.getStatistics(ip1,Measurements.MEAN,null).mean);
		j=1;
		while(j<=nStackSize)
		{
			imp.setSliceWithoutUpdate(j);
			ip2 = getFloatversion(imp);
			ip2.subtract(ImageStatistics.getStatistics(ip2,Measurements.MEAN,null).mean);
			if(nCalcMethod==0)
				{ip = x2D.calcDirectCorrelationImage(ip2,ip2);}
			else
				{ip = x2D.calcFFTCorrelationImage(ip2, ip2);}
		
			sTitle = String.format("corr_%d_x_%d", j,j);
			xymax = getmaxpositionscenterlimit(ip);
			//xymax = getmaxpositions(ip);

			ptable.incrementCounter();
			ptable.addLabel(sTitle);
			ptable.addValue("Xmax_(px)",xymax[0]);	
			ptable.addValue("Ymax_(px)",xymax[1]);
			
			crosscorrstack.addSlice(null, ip);
			IJ.showProgress(j, nStackSize);
    		j++;
	    }
		IJ.showProgress(j, nStackSize);
    }

    IJ.showStatus("Calculating correlation...done.");
    impCC= new ImagePlus(sImTitle, crosscorrstack);
    Calibration cal = impCC.getCalibration();
    cal.xOrigin=Math.round(originalWidth*0.5);
    cal.yOrigin=Math.round(originalHeight*0.5);
    impCC.setCalibration(cal);
    impCC.show();
    IJ.resetMinAndMax();
    //IJ.run("Enhance Contrast", "saturated=0.35");

    
    
    if (bDrift)
	{
    	 if(sChoice.equals("consecutive images"))
    	 {
    		 //correct drift table with respect to the current slide
    		 double xbase,ybase;
    		 xbase = xydrifttable[nCurrentSlice-1][0];
    		 ybase = xydrifttable[nCurrentSlice-1][1];
    		 for (i=0;i<nStackSize;i++)
    		 {
    			 xydrifttable[i][0] -= xbase;
    			 xydrifttable[i][1] -= ybase;
    		 }
    	 }
    	//translating images
	    for (i=0;i<nStackSize;i++)
	    {
	    	imp.setSliceWithoutUpdate(i+1);
			imp.getProcessor().translate(xydrifttable[i][0], xydrifttable[i][1]);
			ptable.setValue("Xdrift_(px)", i, xydrifttable[i][0]);
			ptable.setValue("Ydrift_(px)", i, xydrifttable[i][1]);
			ptable.setValue("Frame", i,i+1);
			
	    }
	    //imp.show();
	}
    imp.setSlice(nCurrentSlice);
    //show results
    ptable.show("Results");

}//end of run()


/** 
 * Dialog displaying options of 2D correlation
 * **/
public boolean x2Dialog(boolean bAutoOn) 
{
	GenericDialog x2DDial = new GenericDialog("2D cross-correlation options");
	String [] items = new String [] {
			"consecutive images","current image in stack and all others","autocorrelation"};
	String [] itemsAuto = new String [] {"autocorrelation"};
	String [] itemsC = new String [] {
			"Direct cross-correlation (slow)","FFT cross-correlation (fast)"};

	if(bAutoOn)
	{	x2DDial.addRadioButtonGroup("Calculate 2D cross correlation between:", itemsAuto, 1, 1, "autocorrelation");}
	else
	{	x2DDial.addRadioButtonGroup("Calculate 2D cross correlation between:", items, 3, 1, Prefs.get("Correlescence.2Dcorr", "consecutive images"));}
	
	x2DDial.addNumericField("for consecutive, interval between images", Prefs.get("Correlescence.2Ddist", 1), 0, 4, "frames (or slices) ");
	x2DDial.addChoice("Calculation method:", itemsC, Prefs.get("Correlescence.2Dcorrmethod", "FFT cross-correlation (fast)"));

	x2DDial.addCheckbox("Correct drift (max of corr)?", Prefs.get("Correlescence.2Ddrift", false));
	x2DDial.addCheckbox("Limit max displacement?", Prefs.get("Correlescence.2Ddriftlimit", false));
	x2DDial.addNumericField("Max displacement X", Prefs.get("Correlescence.2DdriftXmax", 0), 0, 4, "px");
	x2DDial.addNumericField("Max displacement Y", Prefs.get("Correlescence.2DdriftYmax", 0), 0, 4, "px");
	
	x2DDial.setResizable(false);
	x2DDial.showDialog();
	if (x2DDial.wasCanceled())
        return false;

	sChoice=x2DDial.getNextRadioButton();
	if(!bAutoOn)
	{
		Prefs.set("Correlescence.2Dcorr", sChoice);
	}
	nImNumber = (int)x2DDial.getNextNumber();
	Prefs.set("Correlescence.2Ddist", nImNumber);
	nCalcMethod = x2DDial.getNextChoiceIndex();
	Prefs.set("Correlescence.2Dcorrmethod", itemsC[nCalcMethod]);
	bDrift = x2DDial.getNextBoolean();
	Prefs.set("Correlescence.2Ddrift", bDrift);
	bDriftlimit = x2DDial.getNextBoolean();
	Prefs.set("Correlescence.2Ddriftlimit", bDriftlimit);	
	
	nDriftlimitX = (int)x2DDial.getNextNumber();
	Prefs.set("Correlescence.2DdriftXmax", nDriftlimitX);
	nDriftlimitY = (int)x2DDial.getNextNumber();
	Prefs.set("Correlescence.2DdriftYmax", nDriftlimitY);
	
	//drift correction can be used only at frame step=1
	if(bDrift)
	{
		if(nImNumber!=1 && sChoice.equals("consecutive images"))
		{
			IJ.log("Drift correction is possible only with difference of 1 frame");
			IJ.log("Drift correction cancelled");
			bDrift=false;
		}
	}
	
	return true;
}

/** Gets x and y coordinates of maximum intensity pixel on image. */	
int [] getmaxpositions(ImageProcessor ipp)
{
	int [] results = new int [2];
	results[0]=0;
	results[1]=0;

	double s;
	double smax=Double.MIN_VALUE;
	
	for (int i=0;i<ipp.getWidth();i++)
	{
		for (int j=0;j<ipp.getHeight();j++)
		{
			s=ipp.get(i, j);	
			if (s>smax)
			{
				smax=s;
				results[0]=i;
				results[1]=j;
			}
		}
	}
	
	return results;		
}

/** Gets x and y coordinates of maximum intensity pixel on image (with respect to the center
 * and return them, taking into account the min and max  */	
int [] getmaxpositionscenterlimit(ImageProcessor ipp_in)
{
	int [] results = new int [2];
	results[0]=0;
	results[1]=0;
	int nWidth = ipp_in.getWidth();
	int nHeight = ipp_in.getHeight();
	int halfW=(int)Math.round(nWidth*0.5);
	int halfH=(int)Math.round(nHeight*0.5);
	
	ImageProcessor ipp;
	
	double s;
	double smax=Double.MIN_VALUE;
	
	if(bDriftlimit)
	{
		ipp=ipp_in.duplicate();
		ipp.setRoi(halfW-nDriftlimitX, halfH-nDriftlimitY, 2*nDriftlimitX, 2*nDriftlimitY);
		ipp=ipp.crop();
		halfW=nDriftlimitX;
		halfH=nDriftlimitY;
	}
	else
	{
		ipp=ipp_in;
	}
	
	for (int i=0;i<ipp.getWidth();i++)
	{
		for (int j=0;j<ipp.getHeight();j++)
		{
			s=ipp.get(i, j);	
			if (s>smax)
			{
				smax=s;
				results[0]=i;
				results[1]=j;
			}
		}
	}
	results[0]=results[0]-halfW;
	results[1]=results[1]-halfH;
	return results;		
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