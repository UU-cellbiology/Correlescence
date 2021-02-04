package Correlescence;


import ij.IJ;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class imCC1D {
	
	/** Function calculates cross-correlation 
	 *  on kymograph image, assuming image y axis is time and x-axis is x (displacement).
	 *  
	 * **/
	public FloatProcessor xCorrSpaceTime(ImageProcessor ipp1, int nMaxDelay,int nCalcMethod, int nNormMethod)
	{
		
		int imWidth, imHeight;
		FloatProcessor resultip;
		float [] fCurrCC;
		float [] fAverCC;
		float [] row1;
		float [] row2;
		int  resHeight;
		int dt,t,i;
	
		

		imWidth=ipp1.getWidth();
		imHeight=ipp1.getHeight();
		
		resHeight=nMaxDelay;
	
		resultip = new FloatProcessor(imWidth, resHeight);
		
		for(dt=0;dt<resHeight;dt++)
		{
			IJ.showProgress(dt, resHeight-1);

			fAverCC=new float[imWidth];
			fCurrCC=new float[imWidth];
			for (t=0;t<(imHeight-dt);t++)
			{
				row1=getImageRow(ipp1, t, imWidth);
				row2=getImageRow(ipp1, t+dt, imWidth);
				subtractAverage(row1);
				subtractAverage(row2);
				//direct cross-correlation
				if(nCalcMethod==0)
				{
					fCurrCC = calcCCRowsDirect(row1,row2);
				}
				//FFT cross-correlation
				else
				{
					fCurrCC = calcCCRowsFFT(prepareForFFT(row1),prepareForFFT(row2));
					fCurrCC = trimFFTCC(fCurrCC, imWidth);
				}
				//normalize
				normalizeCC(row1,row2,fCurrCC,nNormMethod);
				
				for(i=0;i<imWidth;i++)
				{
					fAverCC[i]+=fCurrCC[i];
				}
			}
			//averaging
			for(i=0;i<imWidth;i++)
			{
				fAverCC[i]/=(float)(imHeight-dt);
				resultip.setf(i, dt,fAverCC[i]);
			}
			
			
		}
		
		return resultip;
		
	}
	/** Function returns directly calculated cross-correlation (non-normalized!) 
	 * of row1 and row2 float arrays of the same size**/
	public float[] calcCCRowsDirect(float [] row1, float [] row2)
	{
		float dCCtop;
		int nWidth = row1.length;
		int i,dx;
		int dxhalf=(int)Math.floor(0.5*(nWidth-1));
		float [] dCCfin = new float [nWidth];
		int imin,imax;

		int nCount;
		

		for(nCount=0;nCount<nWidth;nCount++)
		{
			//shift
			dx=nCount-dxhalf;
			 
			 if(dx>=0)
			 {
				 imin=0;
	             imax=(nWidth-dx);
			 }
	         else
	         {
	             imin=-dx;
	             imax=nWidth;
	         }
			

			 dCCtop=0;

			 for(i=imin;i<imax;i++)
			 {
				 dCCtop += row1[i]*row2[i+dx];
			 }		
			 if(Math.abs(dCCtop)>0.0)
			 {
				 dCCfin[nCount]=(float) (dCCtop);
			 }

				 
		}
		return dCCfin;
	}

	/** function calculates non-normalized cross-correlation between float arrays row1 and row2,
	 * requires their size is equal and length is power of 2 **/
	public float[] calcCCRowsFFT(float [] row1, float [] row2)
	{

		float [] fCCproper;
		float [] fCrossCorr;
		fft1d fftcalc = new fft1d();
		int i,n;
		

		fCrossCorr = fftcalc.correl(row1, row2);
		//reorder output, removing wraparound
		n=fCrossCorr.length;

		fCCproper=new float [n];
		int nhalf=n/2;

		for(i=0;i<(nhalf);i++)
		{
			fCCproper[nhalf-i-1]=fCrossCorr[i];
			fCCproper[nhalf+i]=fCrossCorr[n-i-1];
			
		}
		
		return fCCproper; 
		
	}
	/** function pads float array with zeros to the closest power of 2
	 * ? add windowing ?
	 * **/
	public float [] prepareForFFT(float [] datain)
	{
		
		int n = datain.length;
		int nShift;

		int size = 2;
		while (size<n) size *= 2; // find power of 2 where the data fit
		float[] dataout = new float[size]; // leave the original data untouched, work on a copy
		if (n==size)
			return datain;
		else
		{ 
			//position at center
			nShift=(int)Math.floor((size-n)*0.5);
			System.arraycopy(datain, 0, dataout, nShift, n); // pad to 2^n-size
		}
		return dataout;
	}

	/** function returns float array corresponding to the image 
	 * @param ip image processor
	 * @param nRow row number
	 * @param nWidth image width
	 * **/
	public float [] getImageRow(ImageProcessor ip, int nRow, int nWidth)
	{
		float [] data = new float [nWidth];

		for(int i=0;i<nWidth;i++)
			data[i]= ip.getf(i,nRow);
		return data;
	}
	/** function subtracts average value from provided array
	 * **/
	public void subtractAverage(float [] data)
	{
		int n=data.length;
		float fAver=0;
		int i;
		for(i=0;i<n;i++)
		{
			fAver+=data[i];
		}
		fAver=fAver/(float)n;
		for(i=0;i<n;i++)
		{
			data[i]-=fAver;
		}
		
		return;
	}
	/** function trims float array coming from FFT (power of two)
	 *  to the image width
	 * **/
	public float [] trimFFTCC(float [] xcorr, int imWidth)
	{
		float [] trimmed = new float [imWidth];
		int n=xcorr.length;
		

		System.arraycopy(xcorr, (int)((n-imWidth)*0.5), trimmed, 0, imWidth);
		return trimmed;
		
	}
	/** function that normalizes xcorr correlation array using values
	 * from row1 and row2 and
	 * nNormMethod=0 - overlap area
	 * nNormMethod=1 - full area
	 * **/
	public void normalizeCC(float [] row1, float [] row2, float [] xcorr, int nNormMethod)
	{
		int n=row1.length;
		int i,dx;
		float fnorm1=0;
		float fnorm2=0;
		float [] fInt1;
		float [] fInt2;
		int dxhalf=(int)Math.floor(0.5*(n-1));
		
		if(nNormMethod==1)
		{
			// FULL AREA
			for (i=0;i<n;i++)
			{
				fnorm1+=row1[i]*row1[i];
				fnorm2+=row2[i]*row2[i];
			}
			fnorm1=(float)(1.0/(Math.sqrt(fnorm1*fnorm2)));
			for (i=0;i<n;i++)
			{
				xcorr[i]*=fnorm1; 
			}
		}
		else
		{
			// OVERLAP AREA
			//making integral tables (well, rows)
			fInt1 = new float[n+1];			
			fInt2 = new float[n+1];
			for(i=0;i<n;i++)
			{
				fInt1[i+1] =fInt1[i]+row1[i]*row1[i];
				fInt2[i+1] =fInt2[i]+row2[i]*row2[i];				
			}
			for(i=0;i<n;i++)
			{
				//shift
				dx=i-dxhalf;
				/*if(i==255)
				{
					n++;
					n--;
				}*/
				if(dx>=0)
				{
					fnorm1=fInt1[n-dx];
					fnorm2=fInt2[n]-fInt2[dx];
				}
				else
				{					
					fnorm2=fInt2[n+dx];
					fnorm1=fInt1[n]-fInt1[-dx];										
				}
				fnorm1=(float)(1.0/(Math.sqrt(fnorm1*fnorm2)));
				xcorr[i]*=fnorm1; 
			}
			
		}
		return;
	}
}
