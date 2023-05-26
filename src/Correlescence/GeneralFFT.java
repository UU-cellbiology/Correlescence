package Correlescence;
/* 
 * Free FFT and convolution (Java)
 * 
 * Copyright (c) 2020 Project Nayuki. (MIT License)
 * https://www.nayuki.io/page/free-small-fft-in-multiple-languages
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 * - The above copyright notice and this permission notice shall be included in
 *   all copies or substantial portions of the Software.
 * - The Software is provided "as is", without warranty of any kind, express or
 *   implied, including but not limited to the warranties of merchantability,
 *   fitness for a particular purpose and noninfringement. In no event shall the
 *   authors or copyright holders be liable for any claim, damages or other
 *   liability, whether in an action of contract, tort or otherwise, arising from,
 *   out of or in connection with the Software or the use or other dealings in the
 *   Software.
 */

import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/** this class/code (referenced above) was modified and extended here, thank you, Nayuki! 
 * 
 * mostly changed precision  to float (since input images almost never double)
 * **/

public final class GeneralFFT {
	
	// pre-computed tables of cos and sin 
	// for speed performance
	// Z-chirp cos/sin in Bluestein part
	double [] sinXB;
	double [] cosXB;
	//cos/sin of power of 2 for Radix2
	// this is removed now since using Radix2 from NR
	//double [] sinXR;
	//double [] cosXR;
	
	
	//helpful stuff
	float[] yreal;
	float[] yimag;
	

	/* 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This is a wrapper function.
	 */
	public static void transform(float[] real, float[] imag) {
		int n = real.length;
		if (n != imag.length)
			throw new IllegalArgumentException("Mismatched lengths");
		if (n == 0)
			return;
		else if ((n & (n - 1)) == 0)  // Is power of 2
			transformRadix2(real, imag);
		else  // More complicated algorithm for arbitrary sizes
			transformBluestein(real, imag);
	}
	
	/* 
	 * Computes the discrete Fourier transform (DFT) of the given real vector, returning complex array [0]=real and [1]=imag.
	 * The vector can have any length. This is a wrapper function.
	 */
	public static float[][] transformR(float[] realin) {
		int n = realin.length;
		
		float[][] fin = new float [2][n];
		System.arraycopy(realin, 0, fin[0], 0, n);
		transform(fin[0],fin[1]);
	
		return fin;
		
	}
	
	/* 
	 * Computes the discrete Fourier transform (DFT) of the given real vector, returning complex array [0]=real and [1]=imag.
	 * The vector can have any length. This is a wrapper function.
	 * uses pre-computed values of cos and sin
	 */
	public float[][] prC_transformR(float[] realin) {
		int n = realin.length;
		
		float[] real = new float [n];
		float[] imag = new float [n];
		float[][] fin = new float [2][n];
		System.arraycopy(realin, 0, real, 0, n);
		prC_transform(real,imag);
		fin[0]=real;
		fin[1]=imag;
		return fin;
		
	}
	/* 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This is a wrapper function.
	 */
	public void prC_transform(float[] real, float[] imag) {
		int n = real.length;
		if (n != imag.length)
			throw new IllegalArgumentException("Mismatched lengths");
		if (n == 0)
			return;
		else if ((n & (n - 1)) == 0)  // Is power of 2
			prC_transformRadix2(real, imag);
		else  // More complicated algorithm for arbitrary sizes
			prC_transformBluestein(real, imag);
	}
	
	/* 
	 * Computes the inverse discrete Fourier transform (IDFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This is a wrapper function. This transform does not perform scaling, so the inverse is not a true inverse.
	 */
	public static void inverseTransform(float[] real, float[] imag) {
		int n = real.length;
		for(int i=0;i<n;i++)
		{
			real[i]/=n;
			imag[i]/=n;
		}
		transform(imag, real);
	}
	
	/* 
	 * Computes the inverse discrete Fourier transform (IDFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This is a wrapper function. This transform does not perform scaling, so the inverse is not a true inverse.
	 */
	public void prC_inverseTransform(float[] real, float[] imag) 
	{
		int n = real.length;
		for(int i=0;i<n;i++)
		{
			real[i]/=n;
			imag[i]/=n;
		}
		
		prC_transform(imag, real);
	}

	
	/* 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector's length must be a power of 2. Uses the Cooley-Tukey decimation-in-time radix-2 algorithm.
	 */
	public static void transformRadix2(float[] real, float[] imag) {
		// Length variables
		int n = real.length;
		if (n != imag.length)
			throw new IllegalArgumentException("Mismatched lengths");
		int levels = 31 - Integer.numberOfLeadingZeros(n);  // Equal to floor(log2(n))
		if (1 << levels != n)
			throw new IllegalArgumentException("Length is not a power of 2");
		
		// Trigonometric tables
		double[] cosTable = new double[n / 2];
		double[] sinTable = new double[n / 2];
		double valI;
		for (int i = 0; i < n / 2; i++) {
			valI = 2.0 * i* Math.PI/ n;
			cosTable[i] = Math.cos(valI);
			sinTable[i] = Math.sin(valI);
		}
		
		// Bit-reversed addressing permutation
		for (int i = 0; i < n; i++) {
			int j = Integer.reverse(i) >>> (32 - levels);
			if (j > i) {
				float temp = real[i];
				real[i] = real[j];
				real[j] = temp;
				temp = imag[i];
				imag[i] = imag[j];
				imag[j] = temp;
			}
		}
		
		// Cooley-Tukey decimation-in-time radix-2 FFT
		for (int size = 2; size <= n; size *= 2) {
			int halfsize = size / 2;
			int tablestep = n / size;
			for (int i = 0; i < n; i += size) {
				for (int j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
					int l = j + halfsize;
					double tpre =  real[l] * cosTable[k] + imag[l] * sinTable[k];
					double tpim = -real[l] * sinTable[k] + imag[l] * cosTable[k];
					real[l] = real[j] - (float)tpre;
					imag[l] = imag[j] - (float)tpim;
					real[j] += tpre;
					imag[j] += tpim;
				}
			}
			if (size == n)  // Prevent overflow in 'size *= 2'
				break;
		}
	}
	
	
	public void NR_fourier(float[] real, float[] imag)
	{
		int n = real.length;
		int i;
		float [] datain = new float [2*n];
		for (i=0;i<n;i++)
		{
			datain[i*2]=real[i];
			datain[i*2+1]=imag[i];
		}
		fft1d.four1(datain, 1);
		for (i=0;i<n;i++)
		{
			real[i]=datain[i*2];
			imag[i]=datain[i*2+1];
		}
	}
	/* 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector's length must be a power of 2. Uses the Cooley-Tukey decimation-in-time radix-2 algorithm.
	 */
	public void prC_transformRadix2(float[] real, float[] imag) {
		
		NR_fourier(real, imag);
		/*
		// Length variables
		int n = real.length;
		if (n != imag.length)
			throw new IllegalArgumentException("Mismatched lengths");
		int levels = 31 - Integer.numberOfLeadingZeros(n);  // Equal to floor(log2(n))
		if (1 << levels != n)
			throw new IllegalArgumentException("Length is not a power of 2");
		

		// Bit-reversed addressing permutation
		for (int i = 0; i < n; i++) {
			int j = Integer.reverse(i) >>> (32 - levels);
			if (j > i) {
				float temp = real[i];
				real[i] = real[j];
				real[j] = temp;
				temp = imag[i];
				imag[i] = imag[j];
				imag[j] = temp;
			}
		}
		
		// Cooley-Tukey decimation-in-time radix-2 FFT
		for (int size = 2; size <= n; size *= 2) {
			int halfsize = size / 2;
			int tablestep = n / size;
			for (int i = 0; i < n; i += size) {
				for (int j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
					int l = j + halfsize;
					double tpre =  real[l] * cosXR[k] + imag[l] * sinXR[k];
					double tpim = -real[l] * sinXR[k] + imag[l] * cosXR[k];
					real[l] = real[j] - (float)tpre;
					imag[l] = imag[j] - (float)tpim;
					real[j] += tpre;
					imag[j] += tpim;
				}
			}
			if (size == n)  // Prevent overflow in 'size *= 2'
				break;
		}
		*/
	}
	
	/* 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This requires the convolution function, which in turn requires the radix-2 FFT function.
	 * Uses Bluestein's chirp z-transform algorithm.
	 */
	public static void transformBluestein(float[] real, float[] imag) {
		// Find a power-of-2 convolution length m such that m >= n * 2 + 1
		int n = real.length;
		if (n != imag.length)
			throw new IllegalArgumentException("Mismatched lengths");
		if (n >= 0x20000000)
			throw new IllegalArgumentException("Array too large");
		int m = Integer.highestOneBit(n) * 4;
		
		// Trigonometric tables
		double[] cosTable = new double[n];
		double[] sinTable = new double[n];
		double valI;
		for (int i = 0; i < n; i++) {
			//int j = (int)((long)i * i % (n * 2));  // This is more accurate than j = i * i
			valI= ((int)((long)i * i % (n * 2)))*Math.PI / n;
			cosTable[i] = Math.cos(valI);
			sinTable[i] = Math.sin(valI);
		}
		
		// Temporary vectors and preprocessing
		float[] areal = new float[m];
		float[] aimag = new float[m];
		for (int i = 0; i < n; i++) {
			areal[i] = (float)(real[i] * cosTable[i] + imag[i] * sinTable[i]);
			aimag[i] = (float)(-real[i] * sinTable[i] + imag[i] * cosTable[i]);
		}
		float[] breal = new float[m];
		float[] bimag = new float[m];
		breal[0] = (float)cosTable[0];
		bimag[0] = (float)sinTable[0];
		for (int i = 1; i < n; i++) {
			breal[i] = breal[m - i] = (float)cosTable[i];
			bimag[i] = bimag[m - i] = (float)sinTable[i];
		}
		
		// Convolution
		float[] creal = new float[m];
		float[] cimag = new float[m];
		convolve(areal, aimag, breal, bimag, creal, cimag);
		
		// Postprocessing
		for (int i = 0; i < n; i++) {
			real[i] = (float) (creal[i] * cosTable[i] + cimag[i] * sinTable[i]);
			imag[i] = (float)(-creal[i] * sinTable[i] + cimag[i] * cosTable[i]);
		}
	}
	
	/* 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This requires the convolution function, which in turn requires the radix-2 FFT function.
	 * Uses Bluestein's chirp z-transform algorithm.
	 */
	public void prC_transformBluestein(float[] real, float[] imag) {
		// Find a power-of-2 convolution length m such that m >= n * 2 + 1
		int n = real.length;
		if (n != imag.length)
			throw new IllegalArgumentException("Mismatched lengths");
		if (n >= 0x20000000)
			throw new IllegalArgumentException("Array too large");
		int m = Integer.highestOneBit(n) * 4;
		
		
		// Temporary vectors and preprocessing
		float[] areal = new float[m];
		float[] aimag = new float[m];
		for (int i = 0; i < n; i++) {
			areal[i] = (float)(real[i] * cosXB[i] + imag[i] * sinXB[i]);
			aimag[i] = (float)(-real[i] * sinXB[i] + imag[i] * cosXB[i]);
		}

		// Convolution
		float[] creal = new float[m];
		float[] cimag = new float[m];
		prC_convolve(areal, aimag, creal, cimag);
		
		// Postprocessing
		for (int i = 0; i < n; i++) {
			real[i] = (float) (creal[i] * cosXB[i] + cimag[i] * sinXB[i]);
			imag[i] = (float)(-creal[i] * sinXB[i] + cimag[i] * cosXB[i]);
		}
	}

	
	/* 
	 * Computes the circular convolution of the given real vectors. Each vector's length must be the same.
	 */
	public static void convolve(float[] xvec, float[] yvec, float[] outvec) {
		int n = xvec.length;
		if (n != yvec.length || n != outvec.length)
			throw new IllegalArgumentException("Mismatched lengths");
		convolve(xvec, new float[n], yvec, new float[n], outvec, new float[n]);
	}
	
	
	/* 
	 * Computes the circular convolution of the given complex vectors. Each vector's length must be the same.
	 */
	public static void convolve(float[] xreal, float[] ximag,
			float[] yreal, float[] yimag, float[] outreal, float[] outimag) {
		
		int n = xreal.length;
		if (n != ximag.length || n != yreal.length || n != yimag.length
				|| n != outreal.length || n != outimag.length)
			throw new IllegalArgumentException("Mismatched lengths");
		
		xreal = xreal.clone();
		ximag = ximag.clone();
		yreal = yreal.clone();
		yimag = yimag.clone();
		transform(xreal, ximag);
		transform(yreal, yimag);
		
		for (int i = 0; i < n; i++) {
			double temp = xreal[i] * yreal[i] - ximag[i] * yimag[i];
			ximag[i] = ximag[i] * yreal[i] + xreal[i] * yimag[i];
			xreal[i] = (float)temp;
		}
		//inverse non-scaled transform
		transform(ximag,xreal);
		
		for (int i = 0; i < n; i++) {  // Scaling (because this FFT implementation omits it)
			outreal[i] = xreal[i] / n;
			outimag[i] = ximag[i] / n;
		}
	}
	
	/* 
	 * Computes the circular convolution of the given complex vectors. Each vector's length must be the same.
	 */
	//public void prC_convolve(float[] xreal, float[] ximag,
	//		float[] yreal, float[] yimag, float[] outreal, float[] outimag) 
	public void prC_convolve(float[] xreal, float[] ximag, float[] outreal, float[] outimag) 
	{
		
		int n = xreal.length;
		if (n != ximag.length || n != yreal.length || n != yimag.length
				|| n != outreal.length || n != outimag.length)
			throw new IllegalArgumentException("Mismatched lengths");
		
		xreal = xreal.clone();
		ximag = ximag.clone();
		//yreal = yreal.clone();
		//yimag = yimag.clone();
		prC_transform(xreal, ximag);
		//prC_transform(yreal, yimag);
		
		for (int i = 0; i < n; i++) {
			double temp = xreal[i] * yreal[i] - ximag[i] * yimag[i];
			ximag[i] = ximag[i] * yreal[i] + xreal[i] * yimag[i];
			xreal[i] = (float)temp;
		}
		//inverse non scaled transform
		prC_transform(ximag, xreal);
		//prC_inverseTransform(xreal, ximag);
		
		for (int i = 0; i < n; i++) {  // Scaling (because this FFT implementation omits it)
			outreal[i] = xreal[i] / n;
			outimag[i] = ximag[i] / n;
		}
	}
	
	/** Function returns complex 2D Fourier transform of an arbitrary size image.
	 *  Returned ImageStack has real values at first image and imaginary at the second
	 *  if bSwapQuadrants is true, it will do that (for a visual representation)
	 * **/
	public static ImageStack fft2Dtransform(ImageProcessor ip, boolean bSwapQuadrants)
	{
		int imW=ip.getWidth();
		int imH=ip.getHeight();
		int i,j,k;
		GeneralFFT Gfft = new GeneralFFT();
		
		FloatProcessor dupip = (FloatProcessor) ip.duplicate().convertToFloat();
		FloatProcessor realip = new FloatProcessor(imW,imH);
		
		float [][][] fftc = new float[2][imW][imH];
		float [][] fft1d;
		//float [] datain=null;
		//FFT in rows
		Gfft.preComputeCosSin(imW);
		float [] datain = new float[imW];
		for(j=0;j<imH;j++)
		{
			//fft1d=GeneralFFT.transformR(FFTtools.getImageRow(dupip, j, imW));
			//fft1d=Gfft.prC_transformR(FFTtools.getImageRow(dupip, j, imW));
			FFTtools.getImageRowData(dupip, j, imW, datain);
			//fft1d = GeneralFFT.transformR(datain);
			fft1d = Gfft.prC_transformR(datain);
			for(i=0;i<imW;i++)
				for(k=0;k<2;k++)
				{
					fftc[k][i][j]=fft1d[k][i];
				}
		}
		//FFT in columns
		Gfft.preComputeCosSin(imH);
		for(i=0;i<imW;i++)
		{
			Gfft.prC_transform(fftc[0][i], fftc[1][i]);
			//GeneralFFT.transform(fftc[0][i], fftc[1][i]);
		}
	
		ImageStack fft2d = new ImageStack(imW,imH);
		realip = new FloatProcessor(fftc[0]);
		if(bSwapQuadrants)
			{GFFTswapQuadrants(realip);}
		fft2d.addSlice(realip);
		//dupip.setFloatArray(fftc[1]);
		dupip = new FloatProcessor(fftc[1]);
		if(bSwapQuadrants)
			{GFFTswapQuadrants(dupip);}
		fft2d.addSlice(dupip);
		
		return fft2d;
		
	}
	
	/** Function returns complex Fourier transform of an arbitrary size image per each row.
	 *  Returned ImageStack has real values at first image and imaginary at the second
	 * **/
	public static ImageStack fft2DtransformCols(ImageProcessor ip)
	{
		int imW=ip.getWidth();
		int imH=ip.getHeight();
		int i,j,k;
		
		GeneralFFT Gfft = new GeneralFFT();
		
		FloatProcessor dupip = (FloatProcessor) ip.duplicate().convertToFloat();
		FloatProcessor realip = new FloatProcessor(imW,imH);
		
		float [][][] fftc = new float[2][imW][imH];
		float [][] fft1d;
		Gfft.preComputeCosSin(imH);
		//FFT in cols
		for(i=0;i<imW;i++)
		{
			fft1d=Gfft.prC_transformR(FFTtools.getImageCol(dupip, i, imH));
			for(j=0;j<imH;j++)
				for(k=0;k<2;k++)
				{
					fftc[k][i][j]=fft1d[k][j];
				}
		}
		ImageStack fft2drows = new ImageStack(imW,imH);
		realip.setFloatArray(fftc[0]);
		fft2drows.addSlice(realip);
		dupip.setFloatArray(fftc[1]);
		fft2drows.addSlice(dupip);
		
		return fft2drows;
	}
	
	/** Function returns complex Fourier transform of an arbitrary size image per each row.
	 *  Returned ImageStack has real values at first image and imaginary at the second
	 * **/
	public static FloatProcessor [] fft2DtransformColsArr(ImageProcessor ip)
	{
		int imW=ip.getWidth();
		int imH=ip.getHeight();
		int i,j,k;
		
		GeneralFFT Gfft = new GeneralFFT();
		FloatProcessor [] twoip = new FloatProcessor [2];
		
		//FloatProcessor realip = new FloatProcessor(imW,imH);		
		//FloatProcessor dupip = (FloatProcessor) ip.duplicate().convertToFloat();

		twoip[0]=new FloatProcessor(imW,imH);
		twoip[1]=(FloatProcessor) ip.duplicate().convertToFloat();
				
		float [][][] fftc = new float[2][imW][imH];
		float [][] fft1d;
		Gfft.preComputeCosSin(imH);
		//FFT in cols
		for(i=0;i<imW;i++)
		{
			fft1d=Gfft.prC_transformR(FFTtools.getImageCol(twoip[1], i, imH));
			for(j=0;j<imH;j++)
				for(k=0;k<2;k++)
				{
					fftc[k][i][j]=fft1d[k][j];
				}
		}
		//ImageStack fft2drows = new ImageStack(imW,imH);
		twoip[0].setFloatArray(fftc[0]);
		//fft2drows.addSlice(realip);
		twoip[1].setFloatArray(fftc[1]);
		//fft2drows.addSlice(dupip);
		
		return twoip;
	}
	/** Function returns complex Fourier transform of an arbitrary size image per each column.
	 *  Returned ImageStack has real values at first image and imaginary at the second
	 * **/
	public static ImageStack fft2DtransformRows(ImageProcessor ip)
	{
		int imW=ip.getWidth();
		int imH=ip.getHeight();
		int i,j,k;
		
		GeneralFFT Gfft = new GeneralFFT();
		
		FloatProcessor dupip = (FloatProcessor) ip.duplicate().convertToFloat();
		FloatProcessor realip = new FloatProcessor(imW,imH);
		
		float [][][] fftc = new float[2][imW][imH];
		float [][] fft1d;
		
		
		Gfft.preComputeCosSin(imW);
		//FFT in rows
		for(j=0;j<imH;j++)
		{
			fft1d=Gfft.prC_transformR(FFTtools.getImageRow(dupip, j, imW));
			for(i=0;i<imW;i++)
				for(k=0;k<2;k++)
				{
					fftc[k][i][j]=fft1d[k][i];
				}
		}
		ImageStack fft2drows = new ImageStack(imW,imH);
		realip.setFloatArray(fftc[0]);
		fft2drows.addSlice(realip);
		dupip.setFloatArray(fftc[1]);
		fft2drows.addSlice(dupip);
		
		return fft2drows;
	}
	/** Function returns complex Fourier transform of an arbitrary size image per each column.
	 *  Returned ImageStack has real values at first image and imaginary at the second
	 * **/
	public static FloatProcessor [] fft2DtransformRowsArr(ImageProcessor ip)
	{
		int imW=ip.getWidth();
		int imH=ip.getHeight();
		int i,j,k;
		
		GeneralFFT Gfft = new GeneralFFT();
		FloatProcessor [] twoip = new FloatProcessor [2];
		
		twoip[0]=new FloatProcessor(imW,imH);
		twoip[1]=(FloatProcessor) ip.duplicate().convertToFloat();
		
		//FloatProcessor dupip = (FloatProcessor) ip.duplicate().convertToFloat();
		//FloatProcessor realip = new FloatProcessor(imW,imH);
		
		float [][][] fftc = new float[2][imW][imH];
		float [][] fft1d;
		
		
		Gfft.preComputeCosSin(imW);
		//FFT in rows
		for(j=0;j<imH;j++)
		{
			fft1d=Gfft.prC_transformR(FFTtools.getImageRow(twoip[1], j, imW));
			for(i=0;i<imW;i++)
				for(k=0;k<2;k++)
				{
					fftc[k][i][j]=fft1d[k][i];
				}
		}
		//ImageStack fft2drows = new ImageStack(imW,imH);
		twoip[0].setFloatArray(fftc[0]);
		//fft2drows.addSlice(realip);
		twoip[1].setFloatArray(fftc[1]);
		//fft2drows.addSlice(dupip);
		
		return twoip;
	}
	/**
	 * Function performs inverse fft2D transform given stack with
	 * image 1= real and image 2 complex 2D fft
	 * **/
	public static ImageProcessor fft2Dinverse(ImageStack istack)
	{
		int imW=istack.getWidth();
		int imH=istack.getHeight();
		int i,j,k;
		FloatProcessor realip = new FloatProcessor(imW, imH);
		float [][][] fftc = new float[2][imW][imH];
		float [][] fft1d;
		GeneralFFT Gfft = new GeneralFFT();
		ImageProcessor real = null;
		ImageProcessor imag = null;
		real=istack.getProcessor(1);
		imag=istack.getProcessor(2);
		Gfft.preComputeCosSin(imW);
		//FFT in rows
		for(j=0;j<imH;j++)
		{
			fft1d= new float [2][imW];
			real.getRow(0, j, fft1d[0], imW);
			imag.getRow(0, j, fft1d[1], imW);
			
			//GeneralFFT.inverseTransform(fft1d[0],fft1d[1]);
			Gfft.prC_inverseTransform(fft1d[0],fft1d[1]);
			for(i=0;i<imW;i++)
				for(k=0;k<2;k++)
				{
					fftc[k][i][j]=fft1d[k][i];
				}
		}
		//FFT in columns
		Gfft.preComputeCosSin(imH);
		for(i=0;i<imW;i++)
		{
			Gfft.prC_inverseTransform(fftc[0][i],fftc[1][i]);			
		}

		realip.setFloatArray(fftc[0]);
		return realip;		
		
	}
	
	/**
 	 * function that swaps quadrants on one ImageProcessor
 	 * it does not assume that it has square size
 	 * So it is not really symmetric and requires GFFTswapQuadrantsBack for inversion
 	 * **/
 	public static void GFFTswapQuadrants(ImageProcessor ip) {
 		//IJ.log("swap");
 		int imW, imH;
  		int imWR, imWL, imHT, imHB;
 		ImageProcessor q1, q2, q3, q4;

  		imW= ip.getWidth();
 		imH= ip.getHeight();
 		
 		imWR = Math.floorDiv(imW, 2);
 		imWL = imW-imWR;
 		imHB = Math.floorDiv(imH, 2);
 		imHT = imH-imHB;
 		
 		//left top quadrant (1)
 		ip.setRoi(0,0,imWL,imHT);
 		q1=ip.crop();
 		//bottom right quadrant (4)
 		ip.setRoi(imWL,imHT,imWR,imHB);
 		q4=ip.crop();
 		//top right quadrant (2)
 		ip.setRoi(imWL,0,imWR,imHT);
 		q2=ip.crop();
 		//bottom left quadrant (3)
 		ip.setRoi(0,imHT,imWL,imHB);
 		q3=ip.crop();
 		
 		//(1)<-->(4) and (2)<-->(3)
 		ip.insert(q4,0,0);
 		ip.insert(q3,imWR,0);
 		ip.insert(q2,0,imHB);
 		ip.insert(q1,imWR,imHB);
 		
 		ip.resetRoi();
 	}
 	
	/**
 	 * function that swaps quadrants back (after GFFTswapQuadrants) on one ImageProcessor
 	 * it does not assume that it has square size
 	 * **/
 	public static void GFFTswapQuadrantsBack(ImageProcessor ip) {
 		//IJ.log("swap");
 		int imW, imH;
  		int imWR, imWL, imHT, imHB;
 		ImageProcessor q1, q2, q3, q4;

  		imW= ip.getWidth();
 		imH= ip.getHeight();
 		
 		imWL = Math.floorDiv(imW, 2);
 		imWR = imW-imWL;
 		imHT = Math.floorDiv(imH, 2);
 		imHB = imH-imHT;
 		
 		//left top quadrant (1)
 		ip.setRoi(0,0,imWL,imHT);
 		q1=ip.crop();
 		//bottom right quadrant (4)
 		ip.setRoi(imWL,imHT,imWR,imHB);
 		q4=ip.crop();
 		//top right quadrant (2)
 		ip.setRoi(imWL,0,imWR,imHT);
 		q2=ip.crop();
 		//bottom left quadrant (3)
 		ip.setRoi(0,imHT,imWL,imHB);
 		q3=ip.crop();
 		
 		//(1)<-->(4) and (2)<-->(3)
 		ip.insert(q4,0,0);
 		ip.insert(q3,imWR,0);
 		ip.insert(q2,0,imHB);
 		ip.insert(q1,imWR,imHB);
 		
 		ip.resetRoi();
 	}
 	
 	/**
 	 * function that swaps top and bottom halves (dim=1)
 	 * or left and right halfs (dim=2) (similar to Matlabs fftshift) 
 	 * of ip ImageProcessor
 	 * it does not assume that its size is even
 	 * **/
 	public static void swapHalfs(ImageProcessor ip, int dim) 
 	{
 		int imW, imH;
 		int nHalf;
 		imW= ip.getWidth();
 		imH= ip.getHeight();
 		ImageProcessor h1,h2;
 		if(dim==1)
 		{
 			nHalf = Math.floorDiv(imH, 2);
 			ip.setRoi(0,0,imW,nHalf);
 			h1=ip.crop();
 			ip.setRoi(0,nHalf,imW,nHalf);
 			h2=ip.crop();
 			ip.insert(h2, 0, 0);
 			ip.insert(h1, 0, nHalf);
 		}
 		else
 		{
 			nHalf = Math.floorDiv(imW, 2);
 			ip.setRoi(0,0,nHalf,imH);
 			h1=ip.crop();
 			ip.setRoi(nHalf,0,nHalf,imH);
 			h2=ip.crop();
 			ip.insert(h2, 0, 0);
 			ip.insert(h1, nHalf, 0); 			
 		}
 	}
 	/**
 	 * Calculates power spectrum/phase stack from an input stack
 	 * containing complex FFT in 2D
 	 * **/
 	public static ImageStack fft2Dpowerphase(ImageStack istack)
	{
 		int imW=istack.getWidth();
		int imH=istack.getHeight();
		ImageStack fin = new ImageStack(imW,imH);
		int i,j;
		double freqpow;
		double phaseval;
		float realv,imv;
		
		ImageProcessor realip = istack.getProcessor(1);
		ImageProcessor imagip = istack.getProcessor(2);
		FloatProcessor powerip = new FloatProcessor(imW,imH);
		FloatProcessor phaseip = new FloatProcessor(imW,imH);
		
		for(i=0;i<imW;i++)
			for(j=0;j<imH;j++)
			{
				realv = realip.getf(i, j);
				imv = imagip.getf(i, j);
				freqpow = Math.pow(realv, 2)+Math.pow(imv, 2);
				freqpow = Math.sqrt(freqpow);
				freqpow = Math.log(freqpow);
				phaseval = Math.atan2(imv, realv);
				powerip.setf(i, j, (float)freqpow);
				phaseip.setf(i, j, (float)phaseval);
				
			}
		fin.addSlice(powerip);
		fin.addSlice(phaseip);
		return fin;

	}
 	/** Function pre-computes sin and cos values
 	 * if FFT will happen on the same length arrays,
 	 * to speed up calculations **/
 	public void preComputeCosSin(int N)
 	{
 		int i;
		// Trigonometric tables
		cosXB = new double[N];
		sinXB = new double[N];
		double valI;
		for (i = 0; i < N; i++) 
		{
			valI= ((int)((long)i * i % (N * 2)))*Math.PI / N;
			cosXB[i] = Math.cos(valI);
			sinXB[i] = Math.sin(valI);
		}
		int m = Integer.highestOneBit(N) * 4;
		// Trigonometric tables
	/*	cosXR = new double[m / 2];
		sinXR = new double[m / 2];
		for (i = 0; i < m / 2; i++) 
		{
			valI = 2.0 * i* Math.PI/ m;
			cosXR[i] = Math.cos(valI);
			sinXR[i] = Math.sin(valI);
		}*/
		
		// Temporary vectors and preprocessing

		yreal = new float[m];
		yimag = new float[m];
		yreal[0] = (float)cosXB[0];
		yimag[0] = (float)sinXB[0];
		for (i = 1; i < N; i++) {
			yreal[i] = yreal[m - i] = (float)cosXB[i];
			yimag[i] = yimag[m - i] = (float)sinXB[i];
		}
		transform(yreal,yimag);
 	}
 	
 	/*
	public static void main( String... args) throws Exception
	{
		
		
		float [] testArr = new float[]{20f,21f,22f,23f,24f,23f,22f,21f,45f};
		int size =testArr.length;
		float[][] FFTtr = GeneralFFT.transformR(testArr);
		
		
		//for(int i=0;i<FFTtr[0].length;i++)
		//{
		//	FFTtr[0][i]/=size;
		//	FFTtr[1][i]/=size;
		//}
		
		
		GeneralFFT.inverseTransform(FFTtr[0], FFTtr[1]);
		for(int i=0;i<FFTtr[0].length;i++)
		{
			System.out.println(FFTtr[0][i]);
		}
		System.out.print("done");
	}
*/
}