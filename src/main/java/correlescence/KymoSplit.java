package correlescence;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

public class KymoSplit implements PlugIn {

	ImagePlus imp;
	KymoSplitter kymspl= new KymoSplitter();
	double dMinKymoAngle;
	boolean bPaved;

	@Override
	public void run(String arg) {
		
		int nStackSize;
		String newTitle;
		
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
		
		nStackSize = imp.getStackSize();
		//nCurrentSlice = imp.getCurrentSlice();
		
		if(nStackSize >1)
		{
		    IJ.log("Error! Only one frame is required");
		    	return;
		}
		if(!kymoDialog())
			return;
		if (bPaved)
			{newTitle = imp.getTitle() + String.format("_splitted_minanlge_%d_wPaving", (int)dMinKymoAngle);}
		else
		{newTitle = imp.getTitle() + String.format("_splitted_minanlge_%d_noPaving", (int)dMinKymoAngle);}
		new ImagePlus(newTitle,kymspl.splitkymograph(imp.getProcessor(),dMinKymoAngle,bPaved)).show();
		//new ImagePlus("all", istack).show();
	}
	
	/** 
	 * Dialog displaying options of kymograph splitting
	 * **/
	public boolean kymoDialog() 
	{
		GenericDialog kymoDialog = new GenericDialog("Split kymograph by angle");
		double dSpeed;
		kymoDialog.addNumericField("Minimum speed:", Prefs.get("Correlescence.dSpeed", 0), 2, 4, " px/frame");
		kymoDialog.addCheckbox("Use paving? (longer, better quality)", Prefs.get("Correlescence.bPaved", true));
		kymoDialog.setResizable(false);
		kymoDialog.showDialog();
		if (kymoDialog.wasCanceled())
	        return false;
		
		//STUB, add verification
		dSpeed = kymoDialog.getNextNumber();
		if(dSpeed<0)
		{
			dSpeed=0;
			IJ.log("Minimum speed should be non-negative, resetting it to zero.");
		}
		Prefs.set("Correlescence.dSpeed", dSpeed);
		//calculate minAngle in degrees
		dMinKymoAngle = Math.atan(dSpeed)*180/Math.PI;
		//bPaving
		bPaved = kymoDialog.getNextBoolean();
		Prefs.set("Correlescence.bPaved", bPaved);
	
		return true;
	}

}
