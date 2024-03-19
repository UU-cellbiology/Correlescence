package correlescence;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

public class AngleletsFront implements PlugIn {

	ImagePlus imp;
	double startAngle;
	double endAngle;
	double stepAngle;
	double widthAngle;
	int dBGbox;
	Anglelets anglelet;
	boolean bShowSector;
	
	
	@Override
	public void run(String arg) {
		
		int nStackSize;
		String sTitle;
		double dStartA,dEndA;
		
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
		if(!angleletsDialog(arg))
			return;
		sTitle = imp.getTitle();
		anglelet=new Anglelets();
		
		if(arg.equals("sector"))
		{
			dStartA=checkPositive(360-endAngle);
			dEndA=checkPositive(360-startAngle);
			sTitle=sTitle+"sector_"+Double.toString(startAngle)+"_to_"+Double.toString(endAngle);
			new ImagePlus(sTitle,anglelet.decomposeSector(imp.getProcessor(), dStartA, dEndA, bShowSector)).show();
		}
		else
		{
			//dStartA=checkPositive(360-endAngle);
			//dEndA=checkPositive(360-startAngle);
			dStartA=startAngle;
			dEndA=endAngle;

			new ImagePlus(sTitle,anglelet.decomposeWedge(imp.getProcessor(), dStartA, dEndA,stepAngle, widthAngle, bShowSector, dBGbox)).show();
		}
		
		
		
	}
	/** 
	 * Dialog displaying options of AngleletsFront splitting
	 * **/
	public boolean angleletsDialog(String sOption) 
	{
		GenericDialog angleDialog = new GenericDialog("Split kymograph by angle");
		if(sOption.equals("sector"))
		{
			angleDialog.addNumericField("Start angle in FFT (0-360):", Prefs.get("Anglelets.startAngle", 0), 1, 5, " degrees");
			angleDialog.addNumericField("End angle in FFT (0-360):", Prefs.get("Anglelets.endAngle", 0), 1, 5, " degrees");
			angleDialog.addCheckbox("Show sector shape?", Prefs.get("Anglelets.bShowSector", false));
		}
		else
		{
			angleDialog.addNumericField("Start angle (0-180):", Prefs.get("Anglelets.startAngle", 0), 1, 5, " degrees");
			angleDialog.addNumericField("End angle: (0-180)", Prefs.get("Anglelets.endAngle", 0), 1, 5, " degrees");
			angleDialog.addNumericField("Angle step (0-180):", Prefs.get("Anglelets.stepAngle", 15), 1, 5, " degrees");
			angleDialog.addNumericField("Angle width (0-180):", Prefs.get("Anglelets.widthAngle", 15), 1, 5, " degrees");
			angleDialog.addCheckbox("Show FFT masks?", Prefs.get("Anglelets.bShowSector", false));
			angleDialog.addNumericField("Subtrackt background FFT (0=disable):", Prefs.get("Anglelets.dBGbox", 0), 0, 4, " pixels");
		}

		angleDialog.setResizable(false);
		angleDialog.showDialog();
		if (angleDialog.wasCanceled())
	        return false;
		
		//STUB, add range verification
		startAngle = angleDialog.getNextNumber();
		Prefs.set("Anglelets.startAngle", startAngle);
		endAngle = angleDialog.getNextNumber();
		Prefs.set("Anglelets.endAngle", endAngle);
		if(sOption.equals("sector"))
		{
			bShowSector = angleDialog.getNextBoolean();
			Prefs.set("Anglelets.bShowSector", bShowSector);
		}
		else
		{
			stepAngle = angleDialog.getNextNumber();
			Prefs.set("Anglelets.stepAngle", stepAngle);
			widthAngle = angleDialog.getNextNumber();
			Prefs.set("Anglelets.widthAngle", widthAngle);
			bShowSector = angleDialog.getNextBoolean();
			Prefs.set("Anglelets.bShowSector", bShowSector);
			dBGbox = (int)angleDialog.getNextNumber();
			Prefs.set("Anglelets.dBGbox", dBGbox);

		}
		return true;
	}
	
	/**
	 * function checks that the angle remain positive (from 0 to 360)
	 * **/
	double checkPositive(double angle)
	{
		if(angle<0)
	        return angle+360;
	    
	    if(angle>360)
	         return angle-360;
	    else
	    	return angle;
	}
}
