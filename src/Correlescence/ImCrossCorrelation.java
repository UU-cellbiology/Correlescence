package Correlescence;


import java.awt.Rectangle;

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FHT;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;


public class ImCrossCorrelation {

	public ImageProcessor calcFFTCorrelationImage(ImageProcessor ip1, ImageProcessor ip2)
	{
	
		
		int i,j;
		float dRefRe,dRefIm, dWarpRe, dWarpIm, dVal;
		ImageProcessor paddedip1, paddedip2, normIP; 
		//ip1 and ip2 assumed to be same size
		int origW =ip1.getWidth();
		int origH =ip1.getHeight();
		int nPaddedS;
		
		int imheightFFT, imwidthFFT;
	    FHT fht1,fht2;
	   
	    paddedip1=FFTtools.padzeros(ip1);
	    paddedip2=FFTtools.padzeros(ip2);
	    //since now width and height is the same, just one measurement
	    nPaddedS=paddedip1.getWidth();
	    
	    fht1 = new FHT(paddedip1);
	    fht2 = new FHT(paddedip2);
	    fht1.setShowProgress(false);
	    fht2.setShowProgress(false);
	    
	    fht1.transform();
	    fht2.transform();
	    ImageStack ct1 = fht1.getComplexTransform();
	    ImageStack ct2 = fht2.getComplexTransform();
	    ImageStack istemp;
		
	    FloatProcessor [] is_refFFT = new FloatProcessor[2];
		FloatProcessor [] is_warpFFT = new FloatProcessor[2];
		FloatProcessor [] is_FFTMult = new FloatProcessor[2];
		
		ImageProcessor crossCorrIm;
		
		imheightFFT=ct1.getHeight();
		imwidthFFT =ct1.getWidth();
		for(i=0;i<2;i++)
		{
			is_refFFT[i]=(FloatProcessor) ct1.getProcessor(i+1);
		}
		for(i=0;i<2;i++)
		{
			is_warpFFT[i]=(FloatProcessor) ct2.getProcessor(i+1);
			is_FFTMult[i] = new FloatProcessor(imwidthFFT,imheightFFT);
		}
			
		//calculate convolution between two FFT
		for(i=0;i<imwidthFFT;i++)
			for(j=0;j<imheightFFT;j++)
			{
				dRefRe=is_refFFT[0].getf(i,j);
				dRefIm=is_refFFT[1].getf(i,j);
				
				dWarpRe=is_warpFFT[0].getf(i,j);
				//conjugate
				dWarpIm=(float) ((-1.0)*is_warpFFT[1].getf(i,j));
		
				
				is_FFTMult[0].setf(i,j,dRefRe*dWarpRe-dWarpIm*dRefIm);
				is_FFTMult[1].setf(i,j,dRefRe*dWarpIm+dWarpRe*dRefIm);
				
				
			}

		istemp = new ImageStack(imwidthFFT, imheightFFT);
		
		istemp.addSlice("Real", is_FFTMult[0]);
		istemp.addSlice("Imaginary", is_FFTMult[1]);
		
		ct1 = FFTtools.doComplexInverseTransform(istemp,0, 0);
		
		crossCorrIm = ((FloatProcessor) ct1.getProcessor(1)).duplicate();

		fht1.swapQuadrants(crossCorrIm);
		
		
		// normalize
		
		normIP = calcNormCorrelationCoeff(paddedip1,paddedip2);
		//new ImagePlus("normale", normIP).show();
		for(i=0;i<imwidthFFT;i++)
			for(j=0;j<imheightFFT;j++)
			{
				dVal=crossCorrIm.getf(i,j)/normIP.getf(i,j);
				crossCorrIm.setf(i,j,dVal);
			}
			
		//new ImagePlus("xcorr", crossCorrIm.duplicate()).show();
		
		//crop to original size
		crossCorrIm.setRoi(new Rectangle((int)((nPaddedS-origW)*0.5),(int)((nPaddedS-origH)*0.5),origW,origH));
		
	    return crossCorrIm.crop();
		
	}

	/** Function calculates slow (direct) cross correlation between to images
	 * by shifting them and calculating result
	 * **/
	public ImageProcessor calcDirectCorrelationImage(ImageProcessor ipp1, ImageProcessor ipp2)
	{
		
		ImageProcessor ip1,ip2;

		//to have the same size as FFT cross correlation for comparison
		ip1=FFTtools.padzeros(ipp1);
		ip2=FFTtools.padzeros(ipp2);
		
		int origW = ipp1.getWidth();
	    int origH = ipp1.getHeight();
		int nCorrW = ip1.getWidth();
		
		
	    double dCC;//,dCC1,dCC2;
	    int i,j,m,n;
	    double val1,val2;
	    //double mean1,mean2;
	    
	    //mean1 = ImageStatistics.getStatistics(ip1,Measurements.MEAN,null).mean;
	    //mean2 = ImageStatistics.getStatistics(ip2,Measurements.MEAN,null).mean;
	    
	    FloatProcessor crosscorr= new FloatProcessor (nCorrW,nCorrW);
	    
	    int dx,dy;
	    int dMaxx=(int)Math.round(nCorrW*0.5-1);
	    int dMaxy=(int)Math.round(nCorrW*0.5-1);
	    
	    //correlation shifts
	    for (dx=-dMaxx;dx<=dMaxx;dx++)
	    	for (dy=-dMaxy;dy<=dMaxy;dy++)
	    	{
	    		//now calculate correlation value for these shifts
	    		dCC=0;
				//dCC1=0;
				//dCC2=0;
				for(i=0;i<nCorrW;i++)
					for(j=0;j<nCorrW;j++)
					{
						m=i+dx;
						n=j+dy;
						if(m>-1 &&m<nCorrW&& n>-1 &&n<nCorrW)
						{
							val1=ip1.getf(m,n);//-mean1;
							val2=ip2.getf(i,j);//-mean2;
							dCC+=val1*val2;
							//dCC1+=val1*val1;
							//dCC2+=val2*val2;
						}
					}
						
				//if(dCC1!=0.0 && dCC2!=0.0)
					//crosscorr.setf(dx+(int)(nCorrW*0.5),dy+(int)(nCorrW*0.5),(float)(dCC/(Math.sqrt(dCC1)*Math.sqrt(dCC2))));
					crosscorr.setf(dx+(int)(nCorrW*0.5),dy+(int)(nCorrW*0.5),(float)(dCC));
				//else
					//crosscorr.setf(dx+(int)(nCorrW*0.5),dy+(int)(nCorrW*0.5),0);
	    	}
	    
	    //normalization	    
	    ImageProcessor normIP;
	    float dVal;
	    normIP = calcNormCorrelationCoeff(ip1,ip2);
		for(i=0;i<nCorrW;i++)
			for(j=0;j<nCorrW;j++)
			{
				dVal=crosscorr.getf(i,j)/normIP.getf(i,j);
				crosscorr.setf(i,j,dVal);
			}
		
		//crop to original size
		crosscorr.setRoi(new Rectangle((int)((nCorrW-origW)*0.5),(int)((nCorrW-origH)*0.5),origW,origH));
			    
	    return crosscorr.crop();
	}
	
	
	/** Function calculates normalization coefficient (denominator) 
	 *  for correlation
	 * **/
	public ImageProcessor calcNormCorrelationCoeff(ImageProcessor ip1, ImageProcessor ip2)
	{
		FloatProcessor fpNorm, fpIntegr1, fpIntegr2;
		//FloatProcessor fpCheck1, fpCheck2;
		int imSize = ip1.getWidth(); //width is assumed to be the same as height
		fpNorm = new FloatProcessor(imSize,imSize);
		fpIntegr1 = new FloatProcessor(imSize+1,imSize+1);
		fpIntegr2 = new FloatProcessor(imSize+1,imSize+1);
		//fpCheck1 = new FloatProcessor(imSize,imSize);
		//fpCheck2 = new FloatProcessor(imSize,imSize);
		float dVal1=0, dVal2=0;
		int i,j,dx,dy;
		int imHalfSize;
		
		
		// let's calculate squared integral running sum images
		for (i=0;i<imSize;i++)
			for (j=0;j<imSize;j++)
			{
				dVal1 = ip1.getf(i, j);
				dVal1 = dVal1*dVal1 + fpIntegr1.getf(i,j+1)+ fpIntegr1.getf(i+1,j)-fpIntegr1.getf(i,j);
				fpIntegr1.setf(i+1, j+1, dVal1);
				
				dVal2 = ip2.getf(i, j);
				dVal2 = dVal2*dVal2 + fpIntegr2.getf(i,j+1)+ fpIntegr2.getf(i+1,j)-fpIntegr2.getf(i,j);
				fpIntegr2.setf(i+1, j+1, dVal2);				
			}
		//new ImagePlus("xcorrcheck1", fpIntegr1.duplicate()).show();
		//now let's calculate normalization values for each image
		imHalfSize = (int)(imSize*0.5);
		for(dx=(-imHalfSize);dx<imHalfSize;dx++)
			for(dy=(-imHalfSize);dy<imHalfSize;dy++)
			{
				//for both positive values (bottom right square)
				if(dx>=0 && dy>=0)
				{
					 dVal1=fpIntegr1.getf(imSize,imSize)-fpIntegr1.getf(imSize,dy)-fpIntegr1.getf(dx,imSize)+fpIntegr1.getf(dx,dy);
				     dVal2=fpIntegr2.getf(imSize-dx,imSize-dy);
				}
				//for both negative values (top left square)
				if(dx<=0 && dy<=0)
				{
					dVal1=fpIntegr1.getf(imSize+dx,imSize+dy); 
					dVal2=fpIntegr2.getf(imSize,imSize)-fpIntegr2.getf(imSize,-dy)-fpIntegr2.getf(-dx,imSize)+fpIntegr2.getf(-dx,-dy);
				}
				//for x positive and y negative (bottom left square)
				if(dx>0 && dy<0)
				{       
					dVal1=fpIntegr1.getf(imSize,imSize+dy)-fpIntegr1.getf(dx,imSize+dy);
				    dVal2=fpIntegr2.getf(imSize-dx,imSize)-fpIntegr2.getf(imSize-dx,-dy);
				}
				//for x negative and y positive (top right square)
				if(dx<0 && dy>0)
				{
					dVal1=fpIntegr1.getf(imSize+dx,imSize)-fpIntegr1.getf(imSize+dx,dy);
					dVal2=fpIntegr2.getf(imSize,imSize-dy)-fpIntegr2.getf(-dx,imSize-dy);
				} 
				
				fpNorm.setf(dx+imHalfSize,dy+imHalfSize, (float)((Math.sqrt(dVal1*dVal2))));
				//fpCheck1.setf(dx+imHalfSize,dy+imHalfSize, dVal1);
				//fpCheck2.setf(dx+imHalfSize,dy+imHalfSize, dVal2);
				
			}
		//new ImagePlus("xcorrcheck1", fpCheck1.duplicate()).show();
		//new ImagePlus("xcorrcheck2", fpCheck2.duplicate()).show();
		return fpNorm;
	}
	
	/** Function calculates shift in X and Y between images using
	 * cross correlation calculated through FFT transform
	 * */
	public int [] calcShiftFFTCorrelation(ImageProcessor ip1, ImageProcessor ip2)
	{
		int [] xyshift;
		
		ImageProcessor crossCorrIm = calcFFTCorrelationImage(ip1, ip2);
		xyshift = getmaxpositions(crossCorrIm);
		
		return xyshift;
		
	}

	/** Function calculates shift in X and Y between images using
	 * cross correlation calculated through FFT transform
	 * */
	public double [] calcShiftFFTCorrelationDouble(ImageProcessor ip1, ImageProcessor ip2)
	{
		int [] xyshift;
		int i;
		double [] xyshiftd;
		
		ImageProcessor crossCorrIm = calcFFTCorrelationImage(ip1, ip2);
		xyshift = getmaxpositions(crossCorrIm);
		xyshiftd= new double[xyshift.length];
		for(i=0;i<xyshift.length;i++)
		{
			xyshiftd[i]=(double)xyshift[i];
		}
		
		return xyshiftd;
		
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

}
