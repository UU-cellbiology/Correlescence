package correlescence;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

public class MovieSplit  implements PlugIn {
	
	
	ImagePlus imp;
	KymoSplitter kymspl= new KymoSplitter();
	double dMinKymoAngle;
	boolean bPaved;
	
	@Override
	public void run(String arg) {

		
		int nStackSize;
		int imwidth;
		int imheight;
		int i;
		String newTitle;
		ImageStack isOriginalResliced;

		ImageStack isforward,isbackward,isstatic, issplit;
		
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
		
		if(!movieDialog())
			return;
		
		nStackSize = imp.getStackSize();
		imwidth  = imp.getWidth();
		imheight = imp.getHeight();
		isOriginalResliced = rectangularReslice(imp.getStack());

		
		isforward=new ImageStack(imwidth,nStackSize);
		isbackward=new ImageStack(imwidth,nStackSize);
		isstatic=new ImageStack(imwidth,nStackSize);
		IJ.showStatus("Decomposing...");
		for (i=1;i<=imheight;i++)
		{
			//imp.setSliceWithoutUpdate(i);
			//stub
			issplit = kymspl.splitkymograph(isOriginalResliced.getProcessor(i),dMinKymoAngle, bPaved);
			isforward.addSlice(issplit.getProcessor(1));
			isbackward.addSlice(issplit.getProcessor(2));
			isstatic.addSlice(issplit.getProcessor(3));
			IJ.showProgress(i, imheight);
		}
		IJ.showStatus("Decomposing...done.");
		if (bPaved)
		{newTitle = imp.getTitle() + String.format("_splitted_minanlge_%d_wPaving", (int)dMinKymoAngle);}
		else
		{newTitle = imp.getTitle() + String.format("_splitted_minanlge_%d_noPaving", (int)dMinKymoAngle);}
		
		new ImagePlus(newTitle + "_left_to_right",rectangularReslice(isforward)).show();
		IJ.run("Enhance Contrast", "saturated=0.35");
		new ImagePlus(newTitle +  "_right_to_left",rectangularReslice(isbackward)).show();
		IJ.run("Enhance Contrast", "saturated=0.35");
		new ImagePlus(newTitle +  "_static",rectangularReslice(isstatic)).show();
		IJ.run("Enhance Contrast", "saturated=0.35");
		
		//*/
	}
	
	/** 
	 * Dialog displaying options of kymograph splitting
	 * **/
	public boolean movieDialog() 
	{
		GenericDialog movieDialog = new GenericDialog("Split kymograph by angle");
		double dSpeed;
		movieDialog.addNumericField("Minimum speed:", Prefs.get("Correlescence.dSpeed", 0), 2, 4, " px/frame");
		movieDialog.addCheckbox("Use paving? (longer, better quality)", Prefs.get("Correlescence.bPaved", true));
		movieDialog.setResizable(false);
		movieDialog.showDialog();
		if (movieDialog.wasCanceled())
	        return false;
		
		//STUB, add verification
		dSpeed = movieDialog.getNextNumber();
		if(dSpeed<0)
		{
			dSpeed=0;
			IJ.log("Minimum speed should be non-negative, resetting it to zero.");
		}
		Prefs.set("Correlescence.dSpeed", dSpeed);
		//calculate minAngle in degrees
		dMinKymoAngle = Math.atan(dSpeed)*180/Math.PI;
		//bPaving
		bPaved = movieDialog.getNextBoolean();
		Prefs.set("Correlescence.bPaved", bPaved);
	
		return true;
	}
	
	public ImageStack rectangularReslice(ImageStack stacktoreslice)
	{
		ImageStack istack;
		int stackwidth;
		int stackheight;
		int stacksize;
		int dy;
		
		stackwidth = stacktoreslice.getWidth();
		stackheight = stacktoreslice.getHeight();
		stacksize=stacktoreslice.getSize();
		istack = new ImageStack(stackwidth, stacksize);
		for (dy=0; dy<stackheight; dy++)
		{
			istack.addSlice(getSlice(stacktoreslice,dy));
		}
		return istack;
		
	}
	//function returns horizontal reslice of stack at position y
	ImageProcessor getSlice(ImageStack ilstack, int y) 
	{
		//ImageStack stack = imp.getStack();
		ImageProcessor returnip;
		ImageProcessor ipx;
		int stackwidthslice;
		int stacksizeslice;	
		stackwidthslice=ilstack.getWidth();
		stacksizeslice=ilstack.getSize();	
		returnip = ilstack.getProcessor(1).createProcessor(stackwidthslice, stacksizeslice);
		//= new ImageProcessor(imwidth,nStackSize);

		int i,j;
		for(i=0;i<stacksizeslice;i++)
		{
			ipx = ilstack.getProcessor(i+1);
			for(j=0;j<stackwidthslice;j++)
				returnip.putPixelValue(j, i, ipx.getPixelValue(j, y));
		}
		
		return returnip;
		
	}

}
