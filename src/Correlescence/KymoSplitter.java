package Correlescence;

import java.awt.Color;
import java.awt.Rectangle;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.CanvasResizer;
import ij.plugin.ImageCalculator;
import ij.process.FHT;
import ij.process.ImageProcessor;

public class KymoSplitter {
	
	/** image height and width **/
	int imwidth, imheight;
	/** Fourier transform height and width/padded image size **/
	int nFTsize;
	int nPoints;
	
	public ImageStack splitkymograph(ImageProcessor ip, double angle, boolean bPaved)
	{
		ImageStack istack;
		ImageProcessor paddedip;
		ImageProcessor pavedip;
		ImageProcessor flips;
		//ImageProcessor transformip;
		//ImageProcessor forwardip;
		ImageProcessor tempft;
		ImageStack spectrum;
		ImageStack spectrumf;
		ImageStack invF;
		FHT fht;
		int bitDepth=ip.getBitDepth();
		int[][][] xyPoints=new int [1][1][1];//dummy
		//int[] yPoints;
		Roi roiFilt;
		int i,nFB;
		int nQuad;
		
		imwidth=ip.getWidth();
		imheight = ip.getHeight();
		istack = new ImageStack(imwidth, imheight);
		
		if (bPaved)
		{
			//create paved image
			CanvasResizer CR = new CanvasResizer();
			pavedip=CR.expandImage(ip, 3*imwidth, 3*imheight, imwidth, imheight);
			flips=ip.duplicate();
			flips.flipHorizontal();
			pavedip.insert(flips, 2*imwidth, imheight);
			pavedip.insert(flips, 0, imheight);
			flips.flipVertical();
			pavedip.insert(flips, 0, 0);
			pavedip.insert(flips, 2*imwidth, 0);
			pavedip.insert(flips, 0, 2*imheight);
			pavedip.insert(flips, 2*imwidth, 2*imheight);
			flips=ip.duplicate();
			flips.flipVertical();
			pavedip.insert(flips, imwidth, 0);
			pavedip.insert(flips, imwidth, 2*imheight);
			//new ImagePlus("zzz", pavedip).show();
			//pavedip.
			paddedip=FFTtools.padzeros(pavedip);
		}
		else
		{
			paddedip=FFTtools.padzeros(ip);
		}
		//paddedip=ImCrossCorrelation.padzeros(ip);
		nFTsize=paddedip.getWidth();
		//do fft
	    fht = new FHT(paddedip);
	    fht.setShowProgress(false);
	    fht.transform();
	    
	    
	    spectrum = fht.getComplexTransform();
	    if(angle<=45)
	    {nPoints=4;}
	    else
	    {nPoints=5;}
	    
	    //doing forward and backward transforms
	    
	    for (nFB=0;nFB<2;nFB++)
	    {
		  
		    spectrumf=spectrum.duplicate();
		    
		    if(nFB==0)
		    	//forward filtering
		    	{xyPoints=getRoixyPointsforward(angle);}
		    else
		    	//backward filtering
		    	{xyPoints=getRoixyPointsbackward(xyPoints);}
		    
		    for (nQuad=0;nQuad<2;nQuad++)
		    {
		    	roiFilt = new PolygonRoi(xyPoints[nQuad][0],xyPoints[nQuad][1],nPoints,Roi.POLYGON);
		    	roiFilt.setFillColor(Color.black);
		    	for(i=1;i<3;i++)
		    	{
			    	tempft=spectrumf.getProcessor(i);
			    	tempft.setColor(0);
			    	tempft.fill(roiFilt);
			    	spectrumf.setProcessor(tempft, i);
				}
		    }
		    //new ImagePlus("zzz", spectrumf.getProcessor(1)).show();
		    
		    invF = FFTtools.doComplexInverseTransform(spectrumf,0, 0);
		    
		    //convert back
		    tempft = invF.getProcessor(1);
		    //new ImagePlus("zzz", tempft).show();
		    tempft=setBitDepth(bitDepth, tempft);
		    if(bPaved)
		    {tempft.setRoi(new Rectangle(imwidth,imheight,imwidth,imheight));}
		    else
		    {tempft.setRoi(new Rectangle(0,0,imwidth,imheight));}
		    tempft=tempft.crop();
		   
		    if(nFB==0)
		    	{istack.addSlice("forward",tempft);}
		    else
		    	{istack.addSlice("backward",tempft);}
		    //new ImagePlus("forward", forwardip).show();		    
		   // ImageStack ct1 = fht.ge.getComplexTransform();
			
		 
	    }
	    
	    //get static image by subtracting "forward" and "backward" from original
	    ImageCalculator Icalc =  new ImageCalculator();
	    //Icalc.run('subtract', istack.getProcessor(1), istack.getProcessor(2))
	    ImagePlus tempimp = new ImagePlus("static",ip.duplicate());
	    Icalc.run("subtract", tempimp, new ImagePlus("minforward",istack.getProcessor(1)));
	    Icalc.run("subtract", tempimp, new ImagePlus("minbackward",istack.getProcessor(2)));
	    //tempimp.show();
	    //tempimp = Icalc.run("subtract", new ImagePlus("subtracted",tempimp.getProcessor()), new ImagePlus("beckward",istack.getProcessor(2)));
	    istack.addSlice("static",tempimp.getProcessor());
	    
	    //new ImagePlus("all", istack).show();
	    
	    return istack;		
	}
	
	public ImageProcessor setBitDepth(int bitDepth, ImageProcessor ipx)
	{
		 switch (bitDepth) 
		 {
         	case 8: ipx = ipx.convertToByte(false); break;
         	case 16: ipx = ipx.convertToShort(false); break;
         	case 32: break;
		 }
		 return ipx;
	}
	
	int[][][] getRoixyPointsforward(double angle)
	{
		int[][][] xyPoints;
		int i,j;
		int addedshift;
		/**
		 * first index corresponds to quadrant (left top or right bottom for forward, for example)
		 * second index corresponds to x or y coordinate of point 0=x and 1=y
		 * third index corresponds to the consecutive number of point in polygon
		 * **/
		xyPoints= new int[2][2][nPoints];
		
		
		//left top corner polygon
		//left top point
		xyPoints[0][0][0]=0;
		xyPoints[0][1][0]=0;
		
		//right top point
		//xyPoints[0][0][1]=nFTsize/2-1;
		xyPoints[0][0][1]=nFTsize/2+1;
		xyPoints[0][1][1]=0;

		
		//bottom right point
		xyPoints[0][0][2]=nFTsize/2;
		xyPoints[0][1][2]=nFTsize/2+1;
		//xyPoints[0][0][2]=nFTsize/2+1;
		//xyPoints[0][1][2]=nFTsize/2+4;
		
		
		//bottom left points (one or two, depends on angle)
		
		//angle less than 45 degrees
		if(nPoints==4)
		{
			xyPoints[0][0][3]=0;
			//now calculate y
			addedshift = (int)Math.round(0.5*(double)nFTsize*Math.tan(angle*Math.PI/180.0));
			if(addedshift>1)
				xyPoints[0][1][3]=nFTsize/2+1+addedshift;
			else
				//xyPoints[0][1][3]=nFTsize/2+4;
				xyPoints[0][1][3]=nFTsize/2+1;
		}
		//angle more than 45 degrees
		else
		{
				xyPoints[0][1][3]=nFTsize;
				//now calculate x
				addedshift = (int)Math.round(0.5*(double)nFTsize*Math.tan((90-angle)*Math.PI/180.0));
				xyPoints[0][0][3]=nFTsize/2-addedshift;
				
				//bottom left corner extra point
				xyPoints[0][0][4]=0;
				xyPoints[0][1][4]=nFTsize;
				
				
		}
		//xyPoints[0][1][3]=nFTsize/2+1;
		

		//right bottom corner, flip it	

		for (i=0;i<2;i++)
				for (j=0;j<nPoints;j++)
				{
					xyPoints[1][i][j]=nFTsize-xyPoints[0][i][j];
				}
		return xyPoints;
	}
	int[][][] getRoixyPointsbackward(int[][][] xyPointsforw)
	{
		int[][][] xyPoints;
		int i,j;
		
		
		
		xyPoints= new int[2][2][nPoints];
		
		//y is not changing		
		for (i=0;i<2;i++)
			for (j=0;j<nPoints;j++)
				xyPoints[i][1][j]=xyPointsforw[i][1][j];
		//x is flipped
		for (i=0;i<2;i++)
			for (j=0;j<nPoints;j++)
				xyPoints[i][0][j]=nFTsize-xyPointsforw[i][0][j];
		/*
		//left top corner
		xyPoints[0][0][0]=0;
		xyPoints[0][1][0]=0;
		xyPoints[0][0][1]=nFTsize/2;
		xyPoints[0][1][1]=0;
		xyPoints[0][0][2]=nFTsize/2;
		xyPoints[0][1][2]=nFTsize/2+1;
		xyPoints[0][0][3]=0;
		xyPoints[0][1][3]=nFTsize/2+1;

		//right bottom corner		
		
		xyPoints[1][0][0]=nFTsize/2;
		xyPoints[1][1][0]=nFTsize/2-1;
		xyPoints[1][0][1]=nFTsize;
		xyPoints[1][1][1]=nFTsize/2-1;
		xyPoints[1][0][2]=nFTsize;
		xyPoints[1][1][2]=nFTsize;
		xyPoints[1][0][3]=nFTsize/2;
		xyPoints[1][1][3]=nFTsize;
		*/
		return xyPoints;
	}
	

}
