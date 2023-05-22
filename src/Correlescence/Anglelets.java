package Correlescence;

import java.awt.Color;
import java.awt.Polygon;
import java.awt.Rectangle;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.process.Blitter;
import ij.process.ByteProcessor;
import ij.process.FHT;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class Anglelets {
	
	/** Fourier transform height and width/padded image size **/
	int nFTsize;
	FHT fht;
	/** angles of array corners **/
	double [] arraycorners; 
	/** angles of array corners **/
	int [][] cornerscoords; 
	ImageStack imMask;
	
	
	//empty constructor
	public Anglelets()
	{
		arraycorners = new double[4];
		arraycorners[0]=45;
		arraycorners[1]=45+90;
		arraycorners[2]=45+2*90;
		arraycorners[3]=45+3*90;
		
		
	}
	/** function calculating coordinates of corresponding corners
	 * @return **/
	public void calcCornerCoordinates()
	{
		cornerscoords = new int[4][2];
		//int i;
	
		/*
		
		for(i=0;i<4;i++)
		{
			cornerscoords[i][0]=(int)(nFTsize*0.5*Math.sqrt(2)*Math.cos(arraycorners[i]*Math.PI/180)+0.5*nFTsize);
			cornerscoords[i][1]=(int)(nFTsize*0.5*Math.sqrt(2)*Math.sin(arraycorners[i]*Math.PI/180)+0.5*nFTsize);
		}*/
		
		//manually more precise
		cornerscoords[0][0]=nFTsize;
		cornerscoords[0][1]=nFTsize;

		cornerscoords[1][0]=0;
		cornerscoords[1][1]=nFTsize;

		cornerscoords[2][0]=0;
		cornerscoords[2][1]=0;

		cornerscoords[3][0]=nFTsize;
		cornerscoords[3][1]=0;

		return;
		
	}
	
	public ImageStack decomposeSector(ImageProcessor ip, double startAngle, double endAngle, boolean bShow)
	{
		ImageStack spectrum;
		ImageStack finstack;
		ImageProcessor paddedip;
		int imWidth,imHeight;
		ImageProcessor tempip;
		
		imWidth=ip.getWidth();
		imHeight=ip.getHeight();
		
		//initialize final stack
		finstack = new ImageStack(imWidth,imHeight);
		
		paddedip=FFTtools.padzeros(ip);
		nFTsize=paddedip.getWidth();
		calcCornerCoordinates();
		//do fft
	    fht = new FHT(paddedip);
	    fht.setShowProgress(false);
	    fht.transform();
	    
	    
	    spectrum = fht.getComplexTransform();
	    //new ImagePlus("before",spectrum).show();
	    tempip=getSector(spectrum, startAngle, endAngle,bShow);
	    tempip.setRoi(new Rectangle(0,0,imWidth,imHeight));
	    tempip=tempip.crop();
	    finstack.addSlice(tempip);
		return finstack;
		
	}
	
	public ImageStack decomposeWedge(ImageProcessor ip, double startAngle, double endAngle, double stepAngle, double widthAngle, boolean bShow, int dBGBox)
	{
		ImageStack spectrum;
		ImageStack finstack;
		ImageProcessor paddedip;
		int imWidth,imHeight;
		double dAngle;
		double dSectorStart,dSectorEnd;
		ImageProcessor tempip;
		
		
		imWidth=ip.getWidth();
		imHeight=ip.getHeight();
		
		//initialize final stack
		finstack = new ImageStack(imWidth,imHeight);
		
		paddedip=FFTtools.padzeros(ip);
		nFTsize=paddedip.getWidth();
		calcCornerCoordinates();
		//do fft
	    fht = new FHT(paddedip);
	    fht.setShowProgress(false);
	    fht.transform();
	    
	    if(bShow)
	    {
	    	imMask=new ImageStack(nFTsize,nFTsize);
	    }
	    spectrum = fht.getComplexTransform();
	    getStraightWedgeFilter(startAngle,startAngle+stepAngle);
	    //new ImagePlus("before",spectrum).show();
	    for(dAngle=startAngle;dAngle<=endAngle;dAngle+=stepAngle)
	    {
	    	dSectorStart =checkPositive(dAngle - widthAngle*0.5);
	    	dSectorEnd =checkPositive(dAngle + widthAngle*0.5);
	    	tempip=getWedge(spectrum, dSectorStart, dSectorEnd,bShow,dBGBox);
	    	tempip.setRoi(new Rectangle(0,0,imWidth,imHeight));
	    	tempip=tempip.crop();
	    	finstack.addSlice(tempip);
	    	IJ.showProgress((int)dAngle, (int)(endAngle-startAngle));
	    }
	    IJ.showProgress(2, 2);
	    if(bShow)
	    {
	    	new ImagePlus("FFT Masks",imMask).show();
	    }
		return finstack;
		
	}
	
	/** Gets image spectrum and returns inverse FFT only from the wedge (symmetric with respect to the center of image)
	 * located at between angles startA and endA (taking care of rectangular FFT shape)
	 * **/
	public ImageProcessor getWedge(ImageStack ft, double startA, double endA, boolean bShow, int dBGbox)
	{
		ImageProcessor finalip=null;
		ImageProcessor tempip;
		ImageProcessor tempft;
		ImageStack invF=null;
		ImageStack filtSpectrum = new ImageStack(nFTsize,nFTsize);
		PolygonRoi roiFilt = null;
		PolygonRoi roiMirr = null;
		Roi BGroi =null;
		Rectangle r=null;
		int i;
	
	
		roiFilt=getSectorROI(startA, endA);
		//copy selected ROI only to new (filtered) spectrum stack
		r=roiFilt.getBounds();

		//add empty (zeros) spectrum stack
		filtSpectrum.addSlice(new FloatProcessor(nFTsize,nFTsize));
		filtSpectrum.addSlice(new FloatProcessor(nFTsize,nFTsize));
		//upper wedge part
		for (i=1;i<3;i++)
		{
			tempip=filtSpectrum.getProcessor(i);
			tempft=ft.getProcessor(i);	
			tempft.setRoi(roiFilt);
			
			tempip.snapshot();
			tempip.copyBits(tempft.crop(), r.x, r.y, Blitter.COPY);
			ImageProcessor mask = roiFilt.getMask();
			tempip.setMask(mask);
			tempip.setRoi(roiFilt.getBounds());
			tempip.reset(tempip.getMask());
		}
		//lower wedge part, mirrored with respect to center
		roiMirr=getSectorMirrored(roiFilt);
		r=roiMirr.getBounds();
		for (i=1;i<3;i++)
		{
			tempip=filtSpectrum.getProcessor(i);
			tempft=ft.getProcessor(i);	
			tempft.setRoi(roiMirr);
			
			tempip.snapshot();
			tempip.copyBits(tempft.crop(), r.x, r.y, Blitter.COPY);
			ImageProcessor mask = roiMirr.getMask();
			tempip.setMask(mask);
			tempip.setRoi(roiMirr.getBounds());
			tempip.reset(tempip.getMask());
		}	
		
		//BG subtraction
		//*
		//int nWidth=20;
		if (dBGbox>0)
		{
			BGroi=new Roi((nFTsize-dBGbox)*0.5,(nFTsize-dBGbox)*0.5,dBGbox,dBGbox);
			BGroi.setFillColor(Color.black);
			for(i=1;i<3;i++)
	    	{
				tempip=filtSpectrum.getProcessor(i);
				tempip.setColor(0);
				tempip.fill(BGroi);
				filtSpectrum.setProcessor(tempip, i);
			}
		}
		//*/
		//new ImagePlus("fdf",roiFilt.getMask()).show();
		//new ImagePlus("after",filtSpectrum).show();
		invF=  FFTtools.doComplexInverseTransform(filtSpectrum,0, 0);
		
		if(bShow)
		{
			tempip = new ByteProcessor(nFTsize,nFTsize);
			roiFilt.setFillColor(Color.white);
			tempip.setColor(255);
			tempip.fill(roiFilt);
			roiMirr.setFillColor(Color.white);
			tempip.setColor(255);
			tempip.fill(roiMirr);
			//new ImagePlus("SectorShape",tempip).show();
			if (dBGbox>0)
			{
				BGroi.setFillColor(Color.black);
				tempip.setColor(0);
				tempip.fill(BGroi);			
			}
			
			imMask.addSlice(tempip);
		}
	    
	    //convert back
		finalip = invF.getProcessor(1);
		return finalip;
	}	
	
	/** Returns FFT filter image (nFTsize nFTsize) in the shape of a wedge (symmetric with respect to the center of image)
	 * located at between angles startA and endA (taking care of rectangular FFT shape)
	 * **/
	public FloatProcessor getStraightWedgeFilter(double startA, double endA)//, boolean bShow, int dBGbox)
	{
		FloatProcessor ipFilter = new FloatProcessor(nFTsize, nFTsize);
		PolygonRoi roiFilt = null;
		PolygonRoi roiMirr = null;
		
		roiFilt=getSectorROI(startA, endA);
		//roiFilt.setFillColor(1);
		ipFilter.setColor(1);
		ipFilter.fill(roiFilt);	
		//lower wedge part, mirrored with respect to center
		roiMirr=getSectorMirrored(roiFilt);
		ipFilter.fill(roiMirr);
		ipFilter.smooth();
		
		new ImagePlus("filt",ipFilter).show();
		
		return ipFilter;
	}
	/** Gets image spectrum and returns inverse FFT only from the sector
	 * located at between angles startA and endA (taking care of rectangular FFT shape)
	 * **/
	public ImageProcessor getSector(ImageStack ft, double startA, double endA, boolean bShow)
	{
		ImageProcessor finalip=null;
		ImageProcessor tempip;
		ImageProcessor tempft;
		ImageStack invF=null;
		ImageStack filtSpectrum = new ImageStack(nFTsize,nFTsize);
		Roi roiFilt = null;
		Rectangle r=null;
		int i;
	
	
		roiFilt=getSectorROI(startA, endA);
		//copy selected ROI only to new (filtered) spectrum stack
		r=roiFilt.getBounds();

		//add empty (zeros) spectrum stack
		filtSpectrum.addSlice(new FloatProcessor(nFTsize,nFTsize));
		filtSpectrum.addSlice(new FloatProcessor(nFTsize,nFTsize));
		//upper wedge
		for (i=1;i<3;i++)
		{
			tempip=filtSpectrum.getProcessor(i);
			tempft=ft.getProcessor(i);	
			tempft.setRoi(roiFilt);
			
			tempip.snapshot();
			tempip.copyBits(tempft.crop(), r.x, r.y, Blitter.COPY);
			ImageProcessor mask = roiFilt.getMask();
			tempip.setMask(mask);
			tempip.setRoi(roiFilt.getBounds());
			tempip.reset(tempip.getMask());
		}
		//lower wedge, mirrored with respect to center
		//roiFilt=getWedgeBottom(roiFilt);
		
		
		
		//new ImagePlus("fdf",roiFilt.getMask()).show();
		//new ImagePlus("after",filtSpectrum).show();
		invF=  FFTtools.doComplexInverseTransform(filtSpectrum,0, 0);
	    
	    //convert back
		finalip = invF.getProcessor(1);
		//new ImagePlus("restored",finalip).show();
		
		//show chosen sector (for debug purposes)
		if(bShow)
		{
			tempip = new ByteProcessor(nFTsize,nFTsize);
			roiFilt.setFillColor(Color.white);
			tempip.setColor(255);
			tempip.fill(roiFilt);
			new ImagePlus("SectorShape",tempip).show();
		}
		return finalip;
		
		
	}
	
	/**returns polygon ROI mirrored with respect to image (FFT) center**/
	public PolygonRoi getSectorMirrored(PolygonRoi roiIn)
	{
		PolygonRoi fin=null;
		Polygon poly=null;
		int nPoints,i;

		//get polygon points
		poly=roiIn.getPolygon();
		nPoints = poly.npoints;
		for(i=0;i<nPoints;i++)
		{
			poly.xpoints[i]= nFTsize-poly.xpoints[i];
			poly.ypoints[i]= nFTsize-poly.ypoints[i];
			
		}
		fin = new PolygonRoi(poly,Roi.POLYGON);
		return fin;
		
	}
	
	/**returns polygon roi covering sector from startA till endA**/
	public PolygonRoi getSectorROI(double startA, double endA)
	{
		PolygonRoi roi = null;
		int[][] xyPoints;
		
		int [] nCornersN;
		int [] xyPoint;
		int i;
		int ind;
		
		nCornersN=getCornerPointsN(startA, endA);
		//number of points in polygon:
		//number of corners
		//+center + angles boundaries
		xyPoints= new int[2][nCornersN[0]+3];
		//center
		xyPoints[0][0]=(int) (nFTsize*0.5);
		xyPoints[1][0]=(int) (nFTsize*0.5);
		//beginning angle ray
		xyPoint=GetBorderCoordinates(startA);
		xyPoints[0][1]=xyPoint[0];
		xyPoints[1][1]=xyPoint[1];
		
		//corners
		for(i=0;i<nCornersN[0];i++)
		{
			ind = i+nCornersN[1];
			if (ind>3)
				ind=ind-4;
			xyPoints[0][i+2]=cornerscoords[ind][0];
			xyPoints[1][i+2]=cornerscoords[ind][1];
		}
		
		//last angle
		xyPoint=GetBorderCoordinates(endA);
		xyPoints[0][nCornersN[0]+2]=xyPoint[0];
		xyPoints[1][nCornersN[0]+2]=xyPoint[1];
		
		//make 
		roi=new PolygonRoi(xyPoints[0],xyPoints[1],nCornersN[0]+3,Roi.POLYGON);
		return roi;
	}
	
	/**returns x and y coordinates of intersection of perimeter 
	 * of the image (squared) 
	 * and ray from the middle of the image emanating at the dAngle
	 * (positive, in degrees, from 0 to 360)
	 * **/
	public int [] GetBorderCoordinates(double dAngle)
	{
		int [] xy = new int[2];
		double dHalf=nFTsize*0.5;
		
		//probably can be done more elegantly with rotation,
		//but well, it works :)
		
		//bottom
		if(dAngle>=arraycorners[0]&&dAngle<=arraycorners[1])
		{
			xy[0]=(int)Math.round(dHalf/Math.tan(dAngle*Math.PI/180)+dHalf);
			xy[1]=(int)nFTsize;
			return xy;
		}
		
		//left
		if(dAngle>=arraycorners[1]&&dAngle<=arraycorners[2])
		{
			xy[0]=0;
			xy[1]=(int)Math.round(dHalf/Math.tan((dAngle-90)*Math.PI/180)+dHalf);
			return xy;
		}
		//top
		if(dAngle>=arraycorners[2]&&dAngle<=arraycorners[3])
		{
			xy[0]=(int)Math.round(dHalf/Math.tan((180-dAngle)*Math.PI/180)+dHalf);
			xy[1]=0;
			return xy;
		}
		//right
		else
		{
			xy[0]=nFTsize;
			if(dAngle<=arraycorners[0])
				xy[1]=(int)Math.round(dHalf*Math.tan(dAngle*Math.PI/180)+dHalf);
			else
				xy[1]=(int)Math.round(dHalf*Math.tan((dAngle-360)*Math.PI/180)+dHalf);
			return xy;
		}
		

		
	}
	
	/**returns image corner points number (array val [0])  
	 * between startA angle and endA angle (in degrees, only positive values)
	 * and index of the most close to startA corner (array val [1])
	 * **/
	public int [] getCornerPointsN(double startA, double endA)
	{
		int [] returnarr=new int[2];
		double [] arraycornrot = new double [4];
		boolean [] included = new boolean [4];
		double threshA;
		boolean bInv=false;
		int nCornerNumber=0;
		int i;
		double dMax=400;
		int dMaxInd=0;
		
		double sA, eA;
		
		//check if 
		if(startA>endA)
		{
			sA=endA;
			eA=startA;
			bInv=true;
		}
		else
		{
			sA=startA;
			eA=endA;
		}
		
		threshA=eA-sA;
		//rotate corners by angle A
		for(i=0;i<4;i++)
		{
			arraycornrot[i]=checkPositive(arraycorners[i]-sA);
		}
		//check if angles are within sector defined by angles
		for(i=0;i<4;i++)
		{
			if(arraycornrot[i]<=threshA)
			{
				included[i]=true;
				if(bInv)
					included[i]=false;
			}
			else
			{
				included[i]=false;
				if(bInv)
					included[i]=true;
			}
		}
		//find closest corner among included
		//and count all included ones
		for(i=0;i<4;i++)
		{
			if(included[i])
			{
				nCornerNumber++;
				threshA=checkPositive(arraycorners[i]-startA);
				if(threshA<dMax)
				{
					dMax=threshA;
					dMaxInd=i;
				}
			}
		}
		
		returnarr[0]=nCornerNumber;
		returnarr[1]=dMaxInd;
		return returnarr;
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
