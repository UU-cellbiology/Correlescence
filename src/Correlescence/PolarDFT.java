package Correlescence;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/** Based on the paper "Fast and accurate Polar Fourier transform"
* by A.Averbuch, R.R.Coifman, D.L.Donoho, M.Elad, M.Israeli
* and adapted from the Matlab code available at
* https://elad.cs.technion.ac.il/wp-content/uploads/2018/02/PolarLab.zip
* see prof. Elad page https://elad.cs.technion.ac.il/
 * **/
public class PolarDFT {
	/* 
	 * Computes the Fractional discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. 
	 * y[n]=sum_{k=0}^{N-1} x(k)*exp(-i*2*pi*k*n*alpha)  n=0,1,2, ... ,N-1
	 * 
	 * So that for alpha=1/N we get the regular FFT, and for alpha=-1/N we get the regular IFFT. 
	 */
	public static void transformBluesteinFrac(float[] real, float[] imag, float alpha) {
		// Find a power-of-2 convolution length m such that m >= n * 2 + 1
		int n = real.length;
		if (n != imag.length)
			throw new IllegalArgumentException("Mismatched lengths");
		if (n >= 0x20000000)
			throw new IllegalArgumentException("Array too large");
		//int m = Integer.highestOneBit(n) * 4;
		double valJ;
		// Trigonometric tables
		double[] cosTable = new double[2*n];
		double[] sinTable = new double[2*n];
		for (int i = 0; i < n; i++) {
			//int j = i * i ; 
			valJ=i*i*alpha*Math.PI;
			//cosTable[i] = Math.cos(Math.PI * j * alpha / n);
			//sinTable[i] = Math.sin(Math.PI * j * alpha / n);
			cosTable[i] = Math.cos(valJ);
			sinTable[i] = Math.sin(valJ);
		}
		for (int i = -n; i < 0; i++) {
			valJ=i*i*alpha*Math.PI;
			//cosTable[i] = Math.cos(Math.PI * j * alpha / n);
			//sinTable[i] = Math.sin(Math.PI * j * alpha / n);
			cosTable[i+2*n] = Math.cos(valJ);
			sinTable[i+2*n] = Math.sin(valJ);
		}
		
		// Temporary vectors and preprocessing
		//float[] areal = new float[m];
		//float[] aimag = new float[m];
		float[] areal = new float[2*n];
		float[] aimag = new float[2*n];
		for (int i = 0; i < n; i++) {
			areal[i] = (float)(real[i] * cosTable[i] + imag[i] * sinTable[i]);
			aimag[i] = (float)(-real[i] * sinTable[i] + imag[i] * cosTable[i]);
		}
		//float[] breal = new float[m];
		//float[] bimag = new float[m];
		float[] breal = new float[2*n];
		float[] bimag = new float[2*n];
		breal[0] = (float)cosTable[0];
		bimag[0] = (float)sinTable[0];
		for (int i = 1; i < 2*n; i++) {
			//breal[i] = breal[m - i] = (float)cosTable[i];
			//bimag[i] = bimag[m - i] = (float)sinTable[i];
			breal[i] = breal[2*n - i] = (float)cosTable[i];
			bimag[i] = bimag[2*n - i] = (float)sinTable[i];

		}
		
		// Convolution
		//float[] creal = new float[m];
		//float[] cimag = new float[m];
		float[] creal = new float[2*n];
		float[] cimag = new float[2*n];

		GeneralFFT.convolve(areal, aimag, breal, bimag, creal, cimag);
		
		// Postprocessing
		for (int i = 0; i < n; i++) {
			real[i] = (float) (creal[i] * cosTable[i] + cimag[i] * sinTable[i]);
			imag[i] = (float)(-creal[i] * sinTable[i] + cimag[i] * cosTable[i]);
		}
	}
	
	/* 
	 * Computes the Fractional discrete Fourier transform (DFT) of the given complex vector, 
	 * returning result as float array [0] real and [1] imag.
	 * The vector can have any length. 
	 * y[n]=sum_{k=0}^{N-1} x(k)*exp(-i*2*pi*k*n*alpha)  n=0,1,2, ... ,N-1
	 * 
	 * So that for alpha=1/N we get the regular FFT, and for alpha=-1/N we get the regular IFFT. 
	 */
	public static float[][] transformBluesteinFracReal(float[] realin, float alpha) {
		int n = realin.length;
		float[] real = new float [n];
		float[] imag = new float [n];
		float[][] fin = new float [2][n];
		System.arraycopy(realin, 0, real, 0, n);
		transformBluesteinFrac(real, imag, alpha);	
		fin[0]=real;
		fin[1]=imag;
		return fin;
	}
	
	
	/* 
	 * Computes the Fractional discrete Fourier transform (DFT) of the given complex vector, 
	 * returning result as float array [0] real and [1] imag.
	 * The vector can have any length. 
	 *          y[n]=sum_{k=0}^{N-1} x(k)*exp(-i*2*pi*k*n*alpha)  n=-N/2, ... ,N/2-1
	 *       
	 */
	public static void transformBluesteinFracCenter(float[] real, float[] imag, float alpha) {
		// Find a power-of-2 convolution length m such that m >= n * 2 + 1
		double cosx, sinx;
		int n = real.length;
		if (n != imag.length)
			throw new IllegalArgumentException("Mismatched lengths");
		if (n >= 0x20000000)
			throw new IllegalArgumentException("Array too large");
		//int m = Integer.highestOneBit(n) * 4;
		
		//float[] real = new float [n];
		//float[] imag = new float [n];
		//float[][] fin = new float [2][n];
		//System.arraycopy(realin, 0, real, 0, n);
		
		//correction
		double valJ;
		for (int i = 0; i < n; i++) {
			valJ=i*n*alpha*Math.PI;
			//int j = (int)((long)i * i % (n * 2));  // This is more accurate than j = i * i
			//cosTable[i] = Math.cos(Math.PI * j * alpha / n);
			//sinTable[i] = Math.sin(Math.PI * j * alpha / n);
			cosx=Math.cos(valJ);
			sinx=Math.sin(valJ);
			valJ=real[i];
			real[i] = (float) (real[i]*cosx-imag[i]*sinx);
			imag[i] = (float) (valJ*sinx+imag[i]*cosx);
		}	
		
		//double valJ;
		// Trigonometric tables
		double[] cosTable = new double[2*n];
		double[] sinTable = new double[2*n];
		for (int i = 0; i < n; i++) {

			valJ=i*i*alpha*Math.PI;
			cosTable[i] = Math.cos(valJ);
			sinTable[i] = Math.sin(valJ);
		}
		for (int i = -n; i < 0; i++) {
			valJ=i*i*alpha*Math.PI;
			cosTable[i+2*n] = Math.cos(valJ);
			sinTable[i+2*n] = Math.sin(valJ);
		}
		
		// Temporary vectors and preprocessing

		float[] areal = new float[2*n];
		float[] aimag = new float[2*n];
		for (int i = 0; i < n; i++) {
			areal[i] = (float)(real[i] * cosTable[i] + imag[i] * sinTable[i]);
			aimag[i] = (float)(-real[i] * sinTable[i] + imag[i] * cosTable[i]);
		}

		float[] breal = new float[2*n];
		float[] bimag = new float[2*n];
		breal[0] = (float)cosTable[0];
		bimag[0] = (float)sinTable[0];
		for (int i = 1; i < 2*n; i++) 
		{
			breal[i] = breal[2*n - i] = (float)cosTable[i];
			bimag[i] = bimag[2*n - i] = (float)sinTable[i];

		}
		
		// Convolution

		float[] creal = new float[2*n];
		float[] cimag = new float[2*n];

		GeneralFFT.convolve(areal, aimag, breal, bimag, creal, cimag);
		
		// Postprocessing
		for (int i = 0; i < n; i++) {
			real[i] = (float) (creal[i] * cosTable[i] + cimag[i] * sinTable[i]);
			imag[i] = (float)(-creal[i] * sinTable[i] + cimag[i] * cosTable[i]);
		}
		
	}
	/* 
	 * Computes the Fractional discrete Fourier transform (DFT) of the given complex vector, 
	 * returning result as float array [0] real and [1] imag.
	 * The vector can have any length. 
	 *          y[n]=sum_{k=0}^{N-1} x(k)*exp(-i*2*pi*k*n*alpha)  n=-N/2, ... ,N/2-1
	 *       
	 */
	public static float[][] transformBluesteinFracCenterReal(float[] realin, float alpha) {
		// Find a power-of-2 convolution length m such that m >= n * 2 + 1
		int n = realin.length;
		if (n >= 0x20000000)
			throw new IllegalArgumentException("Array too large");
		//int m = Integer.highestOneBit(n) * 4;
		
		float[] real = new float [n];
		float[] imag = new float [n];
		float[][] fin = new float [2][n];
		//System.arraycopy(realin, 0, real, 0, n);
		
		//correction
		double valJ;
		for (int i = 0; i < n; i++) {
			valJ=i*n*alpha*Math.PI;
			//int j = (int)((long)i * i % (n * 2));  // This is more accurate than j = i * i
			//cosTable[i] = Math.cos(Math.PI * j * alpha / n);
			//sinTable[i] = Math.sin(Math.PI * j * alpha / n);
			real[i] = (float) (realin[i]*Math.cos(valJ));
			imag[i] = (float) (realin[i]*Math.sin(valJ));
		}	
		
		//double valJ;
		// Trigonometric tables
		double[] cosTable = new double[2*n];
		double[] sinTable = new double[2*n];
		for (int i = 0; i < n; i++) {
			//int j = i * i ; 
			valJ=i*i*alpha*Math.PI;
			//cosTable[i] = Math.cos(Math.PI * j * alpha / n);
			//sinTable[i] = Math.sin(Math.PI * j * alpha / n);
			cosTable[i] = Math.cos(valJ);
			sinTable[i] = Math.sin(valJ);
		}
		for (int i = -n; i < 0; i++) {
			valJ=i*i*alpha*Math.PI;
			//cosTable[i] = Math.cos(Math.PI * j * alpha / n);
			//sinTable[i] = Math.sin(Math.PI * j * alpha / n);
			cosTable[i+2*n] = Math.cos(valJ);
			sinTable[i+2*n] = Math.sin(valJ);
		}
		
		// Temporary vectors and preprocessing
		//float[] areal = new float[m];
		//float[] aimag = new float[m];
		float[] areal = new float[2*n];
		float[] aimag = new float[2*n];
		for (int i = 0; i < n; i++) {
			areal[i] = (float)(real[i] * cosTable[i] + imag[i] * sinTable[i]);
			aimag[i] = (float)(-real[i] * sinTable[i] + imag[i] * cosTable[i]);
		}
		//float[] breal = new float[m];
		//float[] bimag = new float[m];
		float[] breal = new float[2*n];
		float[] bimag = new float[2*n];
		breal[0] = (float)cosTable[0];
		bimag[0] = (float)sinTable[0];
		for (int i = 1; i < 2*n; i++) {
			//breal[i] = breal[m - i] = (float)cosTable[i];
			//bimag[i] = bimag[m - i] = (float)sinTable[i];
			breal[i] = breal[2*n - i] = (float)cosTable[i];
			bimag[i] = bimag[2*n - i] = (float)sinTable[i];

		}
		
		// Convolution
		//float[] creal = new float[m];
		//float[] cimag = new float[m];
		float[] creal = new float[2*n];
		float[] cimag = new float[2*n];

		GeneralFFT.convolve(areal, aimag, breal, bimag, creal, cimag);
		
		// Postprocessing
		for (int i = 0; i < n; i++) {
			real[i] = (float) (creal[i] * cosTable[i] + cimag[i] * sinTable[i]);
			imag[i] = (float)(-creal[i] * sinTable[i] + cimag[i] * cosTable[i]);
		}
		
		fin[0]=real;
		fin[1]=imag;
		return fin;
	}
	
	
	/** This function performs a Recto (Pseudo) Polar transform on a Imageprocessor ip
	 * If input is N*N, the output will have 2NS1 by 2NS2 output values. 
	 * If the input is not square, it is squared before the transform.
	 * Also, the size is increased so as to get even number of rows and columns.
	 * s1 - oversampling factor along the rays
	 * s2 - oversampling factor along the slopes (e.g. angles)
	 * Output - Imagestack corresponding to complex [0]real [1]imag
	 *  matrix (theta,r) with dimensions  (2*S1*N)*(2*S2*N)  **/
	public static ImageStack PPDFT(ImageProcessor ip, int s1, int s2)
	{
		int n1, n2;
		int N;
		int i,k;
		int ll;
		double fDiv,fDiv2;
		double alpha, cosx, sinx, tmp;
		ImageProcessor ip_pad;

		//ImageStack isFin;
		FloatProcessor [] ipArrFin = new FloatProcessor[2];
		FloatProcessor ipTemp;
		//ImageStack isTemp;
		FloatProcessor [] ipArrTemp = new FloatProcessor[2];
		FloatProcessor [] ipArrTemp2 = new FloatProcessor[2];
		//ImageStack isTemp2;
		float [][] lltemp;

		
		//Stage 1: Checking input size and defining the sizes of the input/output arrays
		n1=ip.getWidth();
		n2=ip.getHeight();
		N=(int) (Math.ceil(Math.max(n1,n2)*0.5)*2);
		ip_pad = new FloatProcessor(N,N);
		ip_pad.insert(ip,(int)(N/2-Math.floor(n1*0.5)), (int)(N/2-Math.floor(n2*0.5)));
		//new ImagePlus("padded", ip_pad).show();
		//isFin = new ImageStack(2*s1*N,2*s2*N); 
		for(k=0;k<2;k++)
		{
			ipArrFin[k]=new FloatProcessor(2*s1*N,2*s2*N);
			//isFin.addSlice(new FloatProcessor(2*s1*N,2*s2*N));
		}
		///////////////////////////////////////
		// Stage 2: Constructing quadrant 1 and 3
		////////////////////////////////////
		ipTemp=new FloatProcessor(N,(s1*2)*N);
		ipTemp.insert(ip_pad, 0, 0);
		//isTemp = GeneralFFT.fft2DtransformCols(ipTemp);
		ipArrTemp=GeneralFFT.fft2DtransformColsArr(ipTemp);
		for(k=0;k<2;k++)
			//GeneralFFT.swapHalfs(isTemp.getProcessor(k+1), 1);
			GeneralFFT.swapHalfs(ipArrTemp[k], 1);
		//isTemp2= new ImageStack(N*(s2),(s1*2)*N);
		for(k=0;k<2;k++)
		{
			ipTemp=new FloatProcessor(N*(s2),(s1*2)*N);
			ipTemp.insert(ipArrTemp[k], 0, 0);
			ipArrTemp2[k]=ipTemp;
			//isTemp2.addSlice(ipTemp);
		}
		lltemp = new float [2][];
		fDiv = 1./(double)(s1*s2*N*N);
		
		for(ll=0;ll<2*N*s1;ll++)
		{
			for(k=0;k<2;k++)
				{lltemp[k]=FFTtools.getImageRow(ipArrTemp2[k],ll,N*(s2));}
			alpha=(double)(ll-N*s1)*fDiv;
			transformBluesteinFracCenter(lltemp[0], lltemp[1], (float)alpha);
			for(k=0;k<2;k++)
			{
				for(i=0;i<N*(s2);i++)
				{
					ipArrFin[k].setf(i, ll, lltemp[k][N*(s2)-i-1]);
				}
				
			}
			IJ.showProgress(ll, 4*N*s1-1);
		}
		//////////////////////////////////////////////
		// Stage 3: Constructing quadrant 2 and 4
		///////////////////////////////////////////////
		ipTemp=new FloatProcessor((s1*2)*N,N);
		ipTemp.insert(ip_pad, 0, 0);
		//isTemp = GeneralFFT.fft2DtransformRows(ipTemp);
		ipArrTemp= GeneralFFT.fft2DtransformRowsArr(ipTemp);
		for(k=0;k<2;k++)
			//GeneralFFT.swapHalfs(isTemp.getProcessor(k+1), 2);
			GeneralFFT.swapHalfs(ipArrTemp[k], 2);
		//transpose without conjugation
		//isTemp2= new ImageStack(N,(s1*2)*N);
		for(k=0;k<2;k++)
		{
			//isTemp2.addSlice(transpose(isTemp.getProcessor(k+1)));
			ipArrTemp2[k]=transpose(ipArrTemp[k]);
			
		}
		//extend
		//isTemp = new ImageStack(N*s2,(s1*2)*N);
		for(k=0;k<2;k++)
		{
			ipTemp=new FloatProcessor(N*s2,(s1*2)*N);
			ipTemp.insert(ipArrTemp2[k], 0, 0);
			ipArrTemp[k]=ipTemp;
		}
		lltemp = new float [2][];
		fDiv = 1./(double)(s1*s2*N*N);
		fDiv2 = 2.0*Math.PI*(N*s2*0.5-1);
		for(ll=0;ll<2*N*s1;ll++)
		{
			for(k=0;k<2;k++)
				{lltemp[k]=FFTtools.getImageRow(ipArrTemp[k],ll,N*(s2));}
			alpha=(double)(ll-N*s1)*fDiv;
			for(i=0;i<s2*N;i++)
			{
				cosx=Math.cos(fDiv2*i*alpha);
				sinx=Math.sin(fDiv2*i*alpha);
				tmp=lltemp[0][i];
				lltemp[0][i]=(float) (lltemp[0][i]*cosx-lltemp[1][i]*sinx);
				lltemp[1][i]=(float) (tmp*sinx+lltemp[1][i]*cosx);				
			}
			transformBluesteinFrac(lltemp[0], lltemp[1], (float)alpha);
				
			for(k=0;k<2;k++)
			{
				for(i=0;i<N*s2;i++)
				{
					ipArrFin[k].setf(i+N*s2, ll, lltemp[k][i]);
				}
				
			}
			IJ.showProgress(ll+2*N*s1, 4*N*s1-1);
		}
		
		ImageStack isFin = new ImageStack(2*s1*N,2*s2*N); 
		for(k=0;k<2;k++)
		{

			isFin.addSlice(ipArrFin[k]);
		}
		return isFin;
		
	}
	

	
	public static FloatProcessor transpose(FloatProcessor ipin)
	{
		int imW, imH;
		int i,j;
		imW= ipin.getWidth();
 		imH= ipin.getHeight();
		FloatProcessor iFin= new FloatProcessor(imH,imW);
		
		for(i=0;i<imW;i++)
			for(j=0;j<imH;j++)
			{
				iFin.setf(j,i, ipin.getf(i,j));
			}
		
		return iFin;
	}
}
