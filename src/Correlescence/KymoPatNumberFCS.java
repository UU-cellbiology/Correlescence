package Correlescence;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.plugin.filter.Analyzer;
import ij.process.ImageProcessor;

public class KymoPatNumberFCS implements PlugIn  
{

	ImagePlus imp;
	double dPatNBG;
	public ResultsTable ptable = Analyzer.getResultsTable(); 
	@Override
	public void run(String arg) {
		
		int nStackSize;
		int imWidth, imHeight;
		int i,j;
		double dTemp;
		double [] dIntensityT;
		double dIntensityM;
		double dIntensityVAR;
		ImageProcessor ip;
		

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
				//nCurrentSlice = imp.getCurrentSlice();
				
		if(nStackSize >1)
		{
		    IJ.log("Error! Only one frame is required");
		    	return;
		}		
		if(!patNDialog())
			return;
		
		imWidth=imp.getWidth();
		imHeight=imp.getHeight();
		ip=imp.getProcessor();
		
		dIntensityT = new double[imHeight];
	
		
		//calculating average intensity		
		for(j=0;j<imHeight;j++)
		{
			dTemp = 0;
			for(i=0;i<imWidth;i++)
			{
				dTemp+=ip.getf(i, j)-dPatNBG;
				
			}
			dIntensityT[j]=dTemp;
		}
		
		//calculating Mean
		dIntensityM=0;
		for(j=0;j<imHeight;j++)
			dIntensityM+=dIntensityT[j];
		dIntensityM=dIntensityM/(double)imHeight;
		//calculating variance 
		dIntensityVAR=0;
		for(j=0;j<imHeight;j++)
			dIntensityVAR+=Math.pow(dIntensityT[j]-dIntensityM,2);
		dIntensityVAR=dIntensityVAR/(double)imHeight;

		ptable.reset();
		ptable.incrementCounter();
		ptable.addValue("Background", dPatNBG);
		ptable.addValue("Average_Intensity", dIntensityM);
		ptable.addValue("Intensity_Variance", dIntensityVAR);
		ptable.addValue("Particle_Number", dIntensityM*dIntensityM/dIntensityVAR);
		ptable.show("Results");
		
	}
	/** 
	 * Dialog displaying options for background
	 * **/
	public boolean patNDialog() 
	{
		GenericDialog patNDialog = new GenericDialog("Background value");
	
		patNDialog.addNumericField("Average background intensity:", Prefs.get("Correlescence.dPatNBG", 0), 1, 6, " a.u.");
		
		patNDialog.setResizable(false);
		patNDialog.showDialog();
		if (patNDialog.wasCanceled())
	        return false;
		dPatNBG = patNDialog.getNextNumber();
		Prefs.set("Correlescence.dPatNBG", dPatNBG);
		return true;
		
	}

}
