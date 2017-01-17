package Correlescence;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
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
	/** method to calculate xcorr: slow or fast **/
	int nCalcMethod;
	ResultsTable ptable;

@Override
public void run(String arg) {
	
	int nStackSize;
	ImageProcessor ip, ip1,ip2;
	int originalWidth ;
    int originalHeight ;
	int maxN;
	int nCorrW;
    int i,j;
    int [] xymax;
    String sTitle;
    String sImTitle;
    ImageStack crosscorrstack;
    int [][] xydrifttable;
    boolean bAutoCorr =false;
	
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
	
	if(nStackSize ==1)
	{
	    IJ.log("Only one image in stack, limiting choice to autocorrelation");
	    	bAutoCorr=true;
	}
	
	if(!x2Dialog(bAutoCorr))
		return;
	
	originalWidth = imp.getWidth();
    originalHeight = imp.getHeight();
	maxN = Math.max(originalWidth, originalHeight);
	nCorrW = 2;
    while(nCorrW<maxN) nCorrW *= 2;
    
    crosscorrstack= new ImageStack(nCorrW,nCorrW);
    ptable = ResultsTable.getResultsTable();
    
    xydrifttable = new int [nStackSize][2];
    sImTitle = "xcorr_";
    
    if(sChoice.equals("consecutive images"))
    {
    	sImTitle =sImTitle +sChoice;
    	sImTitle = sImTitle+"_delay_"+Integer.toString(nImNumber);
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
	    	
	    //for (i=0;i<nStackSize-1;i++)
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
			xymax = getmaxpositions(ip);
			if (bDrift)
			{
				xydrifttable[i][0]=(int) (xymax[0]-nCorrW*0.5 + xydrifttable[i-1][0]);
				xydrifttable[i][1]=(int) (xymax[1]-nCorrW*0.5 + xydrifttable[i-1][1]);
			}
			ptable.incrementCounter();
			ptable.addLabel(sTitle);
			ptable.addValue("X",xymax[0]-nCorrW*0.5);	
			ptable.addValue("Y",xymax[1]-nCorrW*0.5);
			
			crosscorrstack.addSlice(null, ip);
		    i++; j++;
		    
	
	    }
    }
    //choice:
    //current image in stack and all others
    if(sChoice.equals("current image in stack and all others"))
    {
    	sImTitle =sImTitle +"vs_current_image";
    	i=imp.getCurrentSlice();
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
			xymax = getmaxpositions(ip);
			if (bDrift)
			{
				xydrifttable[j-1][0]=(int) (xymax[0]-nCorrW*0.5);
				xydrifttable[j-1][1]=(int) (xymax[1]-nCorrW*0.5);
			}
			ptable.incrementCounter();
			ptable.addLabel(sTitle);
			ptable.addValue("X",xymax[0]-nCorrW*0.5);	
			ptable.addValue("Y",xymax[1]-nCorrW*0.5);
			
			crosscorrstack.addSlice(null, ip);
    		j++;
    	}
    
    }
    //autocorrelation
    if(sChoice.equals("autocorrelation"))
    {	
    	sImTitle =sImTitle +"autocorrelation_";
    
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
			xymax = getmaxpositions(ip);
			if (bDrift)
			{
				xydrifttable[j-1][0]=(int) (xymax[0]-nCorrW*0.5);
				xydrifttable[j-1][1]=(int) (xymax[1]-nCorrW*0.5);
			}
			ptable.incrementCounter();
			ptable.addLabel(sTitle);
			ptable.addValue("X",xymax[0]-nCorrW*0.5);	
			ptable.addValue("Y",xymax[1]-nCorrW*0.5);
			
			crosscorrstack.addSlice(null, ip);
    		j++;
	    }
    }
    
    
    
    
    new ImagePlus(sImTitle, crosscorrstack).show();
    IJ.resetMinAndMax();
    //IJ.run("Enhance Contrast", "saturated=0.35");
    ptable.show("Results");
    
    
    if (bDrift)
	{
	    for (i=0;i<nStackSize;i++)
	    {
	    	imp.setSliceWithoutUpdate(i+1);
			imp.getProcessor().translate(xydrifttable[i][0], xydrifttable[i][1]);
	    }
	    imp.show();
	}
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
	
	x2DDial.addNumericField("for consecutive, distance between images", Prefs.get("Correlescence.2Ddist", 1), 0, 4, " ");
	x2DDial.addCheckbox("Correct drift (max of corr)?", Prefs.get("Correlescence.2Ddrift", false));
	x2DDial.addChoice("Calculation method:", itemsC, Prefs.get("Correlescence.2Dcorrmethod", "FFT cross-correlation (fast)"));
	
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
	bDrift = x2DDial.getNextBoolean();
	Prefs.set("Correlescence.2Ddrift", bDrift);
	nCalcMethod = x2DDial.getNextChoiceIndex();
	Prefs.set("Correlescence.2Dcorrmethod", itemsC[nCalcMethod]);
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