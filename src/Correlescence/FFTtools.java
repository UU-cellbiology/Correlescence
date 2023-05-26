package Correlescence;

import ij.ImagePlus;
import ij.ImageStack;
import ij.Undo;
import ij.process.FHT;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.StackProcessor;

public class FFTtools {

	/** 
	 * **. Function resizes input ImageProcessor to the closest 2^N size 
	 *  (both in width and height, original image remains at x=0,y=0) 
	 *  and pads this extension with zeros
	 * @param ip = input ImageProcessor
	 * @return padded ImageProcessor
	 */
	  
    public static ImageProcessor padzeros(ImageProcessor ip) 
    {
        int originalWidth = ip.getWidth();
        int originalHeight = ip.getHeight();
        int maxN = Math.max(originalWidth, originalHeight);
        int i = 2;
        while(i<maxN) i *= 2;
        if (i==maxN && originalWidth==originalHeight) 
        {
        
            return ip;
        }
        maxN = i;
   
        ImageProcessor ip2 = ip.createProcessor(maxN, maxN);
        ip2.setValue(0);
        ip2.fill();
        ip2.insert(ip, 0, 0);
  
        Undo.reset();
        //new ImagePlus("padded", ip2.duplicate()).show();
        return ip2;
    }
	/** 
	 * **. Function resizes input ImageProcessor to the closest 2^N size 
	 *  (both in width and height, original image remains at x=0,y=0) 
	 *  and pads this extension with zeros
	 * @param ip = input ImageProcessor
	 * @param value = value to pad
	 * @return padded ImageProcessor
	 */
    
    public static ImageProcessor padvalue(ImageProcessor ip, double value) 
    {
        int originalWidth = ip.getWidth();
        int originalHeight = ip.getHeight();
        int maxN = Math.max(originalWidth, originalHeight);
        int i = 2;
        while(i<maxN) i *= 2;
        if (i==maxN && originalWidth==originalHeight) 
        {
        
            return ip;
        }
        maxN = i;
   
        ImageProcessor ip2 = ip.createProcessor(maxN, maxN);
        ip2.setValue(value);
        ip2.fill();
        
        ip2.insert(ip, 0, 0);
        //ip2.insert(ip, (int)Math.floor((maxN-originalWidth)*0.5), (int)Math.floor((maxN-originalHeight)*0.5));
  
        Undo.reset();
        //new ImagePlus("padded", ip2.duplicate()).show();
        return ip2;
    }
    
	/** 
	 * **. Function resizes input ImageProcessor to the closest 2^N size 
	 *  (both in width and height, original image remains at x=0,y=0) 
	 *  and mirrors/flips this extension with bottom/right parts of the image
	 * @param ip = input ImageProcessor
	 * @param value = value to pad
	 * @return padded ImageProcessor
	 */  
    
    public static ImageProcessor padmirror(ImageProcessor ip) 
    {
        int originalWidth = ip.getWidth();
        int originalHeight = ip.getHeight();
        ImageProcessor t1;
        ImageProcessor ip_in;
        boolean bHmoreW = false;
        boolean bAddMore = false;
        int maxN = Math.max(originalWidth, originalHeight);
        int i = 2;
        while(i<maxN) i *= 2;
        if (i==maxN && originalWidth==originalHeight) 
        {    
            return ip;
        }
        maxN = i;
   
        ImageProcessor ip2 = ip.createProcessor(maxN, maxN);
        
        if(originalHeight<=originalWidth)
        {
        	ip_in=ip;
        }
        else
        {
        	ip_in=ip.rotateRight();
        	ip_in.flipHorizontal();
        	bHmoreW=true;
        	i=originalWidth;
        	originalWidth=originalHeight;
        	originalHeight=i;
        }
        ip2.insert(ip_in, 0, 0);
        //right part
        int dW=maxN-originalWidth;
        if (dW>0)
        {
        	ip.setRoi(originalWidth-dW-1,0,dW,originalHeight);
        	t1=ip_in.crop();
        	t1.flipHorizontal();
        	ip2.insert(t1, originalWidth, 0);
        }
        //bottom part
        int dH=maxN-originalHeight;
        if(dH>originalHeight)
        {
        	dH=originalHeight;
        	bAddMore=true;
        }
        if(dH>0)
        {
        	ip_in.setRoi(0,originalHeight-dH-1,originalWidth,dH);
        	t1=ip_in.crop();
        	t1.flipVertical();
        	ip2.insert(t1, 0, originalHeight);
        	
        }
        
        //bottom right part
        if (dW>0 && dH>0)
        {
        	ip_in.setRoi(originalWidth-dW-1,originalHeight-dH-1,dW,dH);
        	t1=ip_in.crop();
        	t1.flipVertical();
        	t1.flipHorizontal();
        	ip2.insert(t1, originalWidth, originalHeight);
        	
        }
        //flip some more to fill the image
        if(bAddMore)
        {
        	int nFlips = (int)Math.floor(maxN/originalHeight);
        	int nLeftover= maxN-nFlips*originalHeight;
        	if(nFlips>1)
        	{
    			//new image, full width
    			ip2.setRoi(0,originalHeight-1,maxN,originalHeight);
    			t1=ip2.crop();
        		for (i=2;i<(nFlips+1);i++)
        		{
        			t1.flipVertical();
        			ip2.insert(t1, 0, originalHeight*i-1);
        		
        		}
        	}
        	if (nLeftover>0)
        	{
        		ip2.setRoi(0,originalHeight*nFlips-1-nLeftover,maxN,nLeftover);
        		t1=ip2.crop();
        		t1.flipVertical();
    			ip2.insert(t1, 0, originalHeight*nFlips-1);
        	}
        }
        
        if(bHmoreW)
        {
        	ip2.flipHorizontal();
        	ip2=ip2.rotateLeft();
        }
        
        return ip2;
    }
    
    
    public static ImageStack unpad(ImageStack stack, int width, int height) 
    {
      
 
        if (width==0 || height==0 || (width==stack.getWidth()&&height==stack.getHeight()))
            return stack;
        StackProcessor sp = new StackProcessor(stack, null);
        ImageStack stack2 = sp.crop(0, 0, width, height);
        return stack2;
    }
    
    /**	Swap quadrants 1 and 3 and 2 and 4 of the specified ImageProcessor 
	so the power spectrum origin is at the center of the image.
	<pre>
	    2 1
	    3 4
	</pre>
*/
    public static void swapQuadrants(ImageStack stack) {
        FHT fht = new FHT(new FloatProcessor(1, 1));
        fht.setShowProgress(false);
        for (int i=1; i<=stack.getSize(); i++)
            fht.swapQuadrants(stack.getProcessor(i));
    }
    
    /** Complex to Complex Inverse Fourier Transform
     *   @author Joachim Wesner
     */
     public static void c2c2DFFT(float[] rein, float[] imin, int maxN, float[] reout, float[] imout) 
     {
             FHT fht = new FHT(new FloatProcessor(maxN,maxN));
             fht.setShowProgress(false);
             float[] fhtpixels = (float[])fht.getPixels();
             // Real part of inverse transform
             for (int iy = 0; iy < maxN; iy++)
                   cplxFHT(iy, maxN, rein, imin, false, fhtpixels);
             fht.inverseTransform();
             // Save intermediate result, so we can do a "in-place" transform
             float[] hlp = new float[maxN*maxN];
             System.arraycopy(fhtpixels, 0, hlp, 0, maxN*maxN);
             // Imaginary part of inverse transform
             for (int iy = 0; iy < maxN; iy++)
                   cplxFHT(iy, maxN, rein, imin, true, fhtpixels);
             fht.inverseTransform();
             System.arraycopy(hlp, 0, reout, 0, maxN*maxN);
             System.arraycopy(fhtpixels, 0, imout, 0, maxN*maxN);
       }
     
     /** Build FHT input for equivalent inverse FFT
     *   @author Joachim Wesner
     */
     public static void cplxFHT(int row, int maxN, float[] re, float[] im, boolean reim, float[] fht) 
     {
             int base = row*maxN;
             int offs = ((maxN-row)%maxN) * maxN;
             if (!reim) {
                   for (int c=0; c<maxN; c++) {
                         int l =  offs + (maxN-c)%maxN;
                         fht[base+c] = ((re[base+c]+re[l]) - (im[base+c]-im[l]))*0.5f;
                   }
             } else {
                   for (int c=0; c<maxN; c++) {
                         int l = offs + (maxN-c)%maxN;
                         fht[base+c] = ((im[base+c]+im[l]) + (re[base+c]-re[l]))*0.5f;
                   }
             }
       }
     
     public static ImageStack doComplexInverseTransform(ImageStack stack, int width, int height ) 
     {
    
         int maxN = stack.getWidth();
         FFTtools.swapQuadrants(stack);
         float[] rein = (float[])stack.getPixels(1);
         float[] imin = (float[])stack.getPixels(2);
         float[] reout= new float[maxN*maxN];
         float[] imout = new float[maxN*maxN];
         c2c2DFFT(rein, imin, maxN, reout, imout);
         ImageStack stack2 = new ImageStack(maxN, maxN);
         FFTtools.swapQuadrants(stack);
         stack2.addSlice("Real", reout);
         stack2.addSlice("Imaginary", imout);
         stack2 = unpad(stack2,width,height);
         /*String name = WindowManager.getUniqueName(imp.getTitle().substring(10));
         ImagePlus imp2 = new ImagePlus(name, stack2);
         imp2.getProcessor().resetMinAndMax();
         imp2.show();*/
         return stack2;
     }
     

 	/**
 	 *  This module implements the periodic-plus-smooth decomposition of an
 	 *  image introduced by L. Moisan [J. Math. Imag. Vision 39(2), 161-179,
 	 *  doi:10.1007/s10851-010-0227-1]. 
 	 * **/
 	public static ImageStack PeriodicSmooth(ImageProcessor ip, boolean bInverseFFT)
 	{
 		int imW,imH;
 		int i,j,k;

 		ImageProcessor smoothIP=null;
 		
 		imW=ip.getWidth();
 		imH=ip.getHeight();
 		
 		ImageProcessor vip = ip.createProcessor(imW,imH);
 		vip = vip.convertToFloat();
 		float diff;
 		for (i=0;i<imW;i++)
 		{
 			diff = ip.getf(i, 0) - ip.getf(i, imH-1);
 			vip.putPixelValue(i, 0, diff);
 			vip.putPixelValue(i, imH-1, (-1)*diff);
 		}
		for (j=0;j<imH;j++)
		{
			diff = ip.getf(0, j) - ip.getf(imW-1, j);
			
			vip.putPixelValue(0, j, vip.getf(0, j)+diff);
			vip.putPixelValue(imW-1, j, vip.getf(imW-1, j)-diff);
		}
 		
 		//new ImagePlus("boundaries",vip).show();
 
		ImageStack ct_dft = GeneralFFT.fft2Dtransform(vip, false);
 	
		//new ImagePlus("FFT_VIP",ct_dft).show();
 		ImageStack ct_inv = null;
 		
 		FloatProcessor k_dft =k_dft (imW,  imH);
 		//new ImagePlus("cosines",k_dft).show();
 		ImageProcessor complexFFT;
 		for (k=1;k<3;k++)
 		{
 			complexFFT=ct_dft.getProcessor(k);
 			for(i=0;i<imW; i++)
 				for (j=0;j<imH;j++)
 				{				
 					complexFFT.putPixelValue(i, j, complexFFT.getPixelValue(i, j)*k_dft.getf(i,j));
 				}
 			complexFFT.putPixelValue(imW/2, imW/2, 0.);
 		}
 		
 		
 		//return complex FFT of smooth part 
 		if(bInverseFFT)
 			{return ct_dft;}
 		//return original image minus smooth part and smooth part
 		else
 		{
 			smoothIP=GeneralFFT.fft2Dinverse(ct_dft);
 			ImageProcessor ip_subtract = smoothIP.createProcessor(imW, imH);
 			ct_inv = new ImageStack(imW,imH);
 			//original minus smooth
 			for(i=0;i<imW; i++)
 				for (j=0;j<imH;j++)
 				{ip_subtract.putPixelValue(i, j, ip.getf(i,j)-smoothIP.getf(i,j));}
 			
 			ct_inv.addSlice(ip_subtract);
 			//add smooth
 			ct_inv.addSlice(smoothIP);
 			return ct_inv;
 		}
 			
 		
 		
 	}
 	/**
 	 * function generates 2D image of cosines, part of PeriodicSmooth process
 	 * **/
 	
 	public static FloatProcessor  k_dft (int imW, int imH)
 	{
 		FloatProcessor result = new FloatProcessor(imW, imH);
 		
 		double [] m =cos_array(imW);
 		double [] n =cos_array(imH);

 		for (int i=0;i<imW;i++)
 			for (int j=0;j<imH;j++)
 			{
 				//result.setf(i, j, (float)(2.0*(m[i]+n[j]-2.0)));
 				result.setf(i, j, (float)(0.5/(2.0 - m[i]-n[j])));
 			}
 		result.setf(0, 0, 1.f);
 		//FFTswapQuadrants(result);
 		return result;		
 	}
 	/**
 	 * function generates frequency values for 1dFFT
 	 * and multiplies them by cosines, part of PeriodicSmooth process
 	 * @param n
 	 * @return
 	 */
 	
 	public static double [] cos_array(int n)
 	{
 		double [] cos_n = new double[n];
 		int i;
 		double val,V;
 		int N;
 		
 		val=2.*Math.PI/(double)n;
 		N=Math.floorDiv(n-1, 2)+1;
 		V=Math.floorDiv(n, 2);
 		for (i=0;i<N;i++)
 		{
 			cos_n[i]=(double)i;
 		}
 		for (i=N;i<n;i++)
 		{
 			cos_n[i]=V;
 			V=V-1;
 		}
 		
 		for(i=0;i<n;i++)
 			cos_n[i]=Math.cos(cos_n[i]*val);
 		
 		return cos_n;
 	}
 	

 	

 	
 	/**
 	 * function that swaps quadrants on one imageprocessor
 	 * assumes square image as an input
 	 * **/
 	public static void FFTswapQuadrants(ImageProcessor ip) {
 		//IJ.log("swap");
  		ImageProcessor t1, t2;
 		int size = ip.getWidth()/2;
 		ip.setRoi(size,0,size,size);
 		t1 = ip.crop();
   		ip.setRoi(0,size,size,size);
 		t2 = ip.crop();
 		ip.insert(t1,0,size);
 		ip.insert(t2,size,0);
 		ip.setRoi(0,0,size,size);
 		t1 = ip.crop();
   		ip.setRoi(size,size,size,size);
 		t2 = ip.crop();
 		ip.insert(t1,size,size);
 		ip.insert(t2,0,0);
 		ip.resetRoi();
 	}
 	
	/** function returns float array corresponding to the image 
	 * @param ip image processor
	 * @param nRow row number
	 * @param nWidth image width
	 * **/
	public static float [] getImageCol(ImageProcessor ip, int nCol, int nHeight)
	{
		float [] data = new float [nHeight];

		for(int j=0;j<nHeight;j++)
			data[j]= ip.getf(nCol,j);
		return data;
	}
	/** function returns float array corresponding to the image 
	 * @param ip image processor
	 * @param nRow row number
	 * @param nWidth image width
	 * **/
	public static float [] getImageRow(ImageProcessor ip, int nRow, int nWidth)
	{
		float [] data = new float [nWidth];

		for(int i=0;i<nWidth;i++)
			data[i]= ip.getf(i,nRow);
		return data;
	}
	/** function returns float array corresponding to the image 
	 * @param ip image processor
	 * @param nRow row number
	 * @param nWidth image width
	 * **/
	public static void getImageRowData(ImageProcessor ip, int nRow, int nWidth, float [] data)
	{
		//float [] data = new float [nWidth];

		for(int i=0;i<nWidth;i++)
			data[i]= ip.getf(i,nRow);
		//return data;
	}
	public static void putImageCol(ImageProcessor ip, int nCol, int nHeight, float[] data)
	{
	
		for(int j=0;j<nHeight;j++)
			ip.putPixelValue(nCol, j, data[j]);
		return ;
	}

	public static void putImageRow(ImageProcessor ip, int nRow, int nWidth, float[] data)
	{
	
		for(int i=0;i<nWidth;i++)
			ip.putPixelValue(i, nRow, data[i]);
		return ;
	}
	
}
