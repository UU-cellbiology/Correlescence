package Correlescence;

/** class performing 1D DFT (Discrete Fourier Transform) using "Numerical recipes" code
 * 	+ cross correlation addition
 * **/
public class fft1d {
	
	float [] fFreqPhase = null;

	
	//empty constructor
	public fft1d()
	{
		
	}
	
	/** Calculates DFT transform of complex input 
	 *  check "Numerical Recipes for explanation  
	 *  @param data supposed to be 0,2,4,... real part, 1,3,5... imaginary part
	 *  @param isign specifies whether it is forward or inverse transform
	 * **/
	public static void four1(float[] data, int isign)
	{
		int nn, mmax, m, j, istep, i;
		float wtemp, wr, wpr, wpi, wi, theta, tempr, tempi, tempf;
		
		nn = data.length;
		
		int n=nn>>1;
		j=1;
		for(i=1;i<nn;i+=2)
		{
			if(j>i)
			{
				tempf=data[j-1];
				data[j-1]=data[i-1];
				data[i-1]=tempf;
				tempf=data[j];
				data[j]=data[i];
				data[i]=tempf;
			}
			m=n;
			while (m>=2 && j>m)
			{
				j-=m;
				m>>=1;
			}
			j+=m;
		}
		//Danielson-Lanczos routine section
		mmax=2;
		while(nn>mmax)
		{
			istep=mmax<<1;
			theta = (float) (isign*(2.0*Math.PI/mmax));
			wtemp = (float) Math.sin(0.5*theta);
			wpr=-2.0f*wtemp*wtemp;
			wpi=(float)Math.sin(theta);
			wr=1.0f;
			wi=0.0f;
			for(m=1;m<mmax;m+=2)
			{
				for(i=m;i<=nn;i+=istep)
				{
					j=i+mmax;
					tempr=wr*data[j-1]-wi*data[j];
					tempi=wr*data[j]+wi*data[j-1];
					data[j-1]=data[i-1] - tempr;
					data[j]=data[i]-tempi;
					data[i-1]+=tempr;
					data[i]+=tempi;
				}
				wtemp=wr;
				wr=wtemp*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			mmax=istep;
			
		}

	}
	
	/** returns discrete Fourier transform of real float array datain
	 *  (taken from "Numerical recipes", check it for details)
	 *  assumes that input has length of power of 2!
	 * **/
	public float[] realft(float[] datain, int isign)
	{
		int i,i1,i2,i3,i4;
		float c1=0.5f,c2,h1r,h1i,h2r,h2i,wr,wi,wpr,wpi,wtemp;
		float theta;

		
		int n = datain.length;
		//int size = 2;
		//while (size<n) size *= 2; // find power of 2 where the data fit
		float[] data = new float[n]; // leave the original data untouched, work on a copy
		System.arraycopy(datain, 0, data, 0, n); // pad to 2^n-size with zeros
		//n=size;
		
		theta=(float) (Math.PI/(double)(n>>1));
		
		if (isign == 1) 
		{
			c2 = -0.5f;
			four1(data,1);
		} 
		else 
		{
			c2=0.5f;
			theta = -theta;
		}
		wtemp=(float) Math.sin(0.5*theta);
		wpr = -2.0f*wtemp*wtemp;
		wpi=(float) Math.sin(theta);
		wr=1.0f+wpr;
		wi=wpi;
		for (i=1;i<(n>>2);i++) 
		{
			i2=1+(i1=i+i);
			i4=1+(i3=n-i1);
			h1r=c1*(data[i1]+data[i3]);
			h1i=c1*(data[i2]-data[i4]);
			h2r= -c2*(data[i2]+data[i4]);
			h2i=c2*(data[i1]-data[i3]);
			data[i1]=h1r+wr*h2r-wi*h2i;
			data[i2]=h1i+wr*h2i+wi*h2r;
			data[i3]=h1r-wr*h2r+wi*h2i;
			data[i4]= -h1i+wr*h2i+wi*h2r;
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		if (isign == 1) 
		{
			data[0] = (h1r=data[0])+data[1];
			data[1] = h1r-data[1];
		} else 
		{
			data[0]=c1*((h1r=data[0])+data[1]);
			data[1]=c1*(h1r-data[1]);
			four1(data,-1);
		}
		return data;
	}
	
	/**
	 * function calculating cross-correlation between two signals
	 * results are returned in wraparound order, i.e. correlations
	 * at increasingly negative lags are in ans[n-1] on down to ans[n/2],
	 * while correlations at increasingly positive lags are 
	 * in ans[0] (zero lag) on up to ans[n/2-1].
	 * if data1 lags data2 (i.e. is shifter to the right).
	 * then ans will show peak at positive lags
	 * **/
	public float[] correl(float[] data1, float[] data2)
	{
		float [] ans;
		float [] temp;
		float tmp;
		int n=data1.length;
		int i, no2;
		
		ans=realft(data1,1);
		temp=realft(data2,1);
				
		no2=n>>1;
		
		for(i=2;i<n;i+=2)
		{
			tmp=ans[i];
			ans[i]=(ans[i]*temp[i]+ans[i+1]*temp[i+1])/no2;
			ans[i+1]=(ans[i+1]*temp[i]-tmp*temp[i+1])/no2;
		}
		ans[0]=ans[0]*temp[0]/no2;
		ans[1]=ans[1]*temp[1]/no2;
		
		ans=realft(ans,-1);
		
		return ans;
	
	}
	
	/**
	 * function calculating auto-correlation between two signals
	 * results are returned in wraparound order, i.e. correlations
	 * at increasingly negative lags are in ans[n-1] on down to ans[n/2],
	 * while correlations at increasingly positive lags are 
	 * in ans[0] (zero lag) on up to ans[n/2-1].
	 * if (nCalcFreqPhase == true),calculates dominant frequency (with max amplitude)
	 * and phase and stores it in the corresponding member variables
	 * **/
	public float[] autocorrel(float[] data1, boolean nCalcFreqPhase)
	{
		float [] ans;

		int n=data1.length;
		float [] invPrep = new float[n];
		int i, no2;
		float maxAmp =Float.MIN_VALUE;
		int indMax=0;
		
		ans=realft(data1,1);				
		no2=n>>1;
		
		for(i=2;i<n;i+=2)
		{
			invPrep[i]=(ans[i]*ans[i]+ans[i+1]*ans[i+1])/no2;
			if(maxAmp<invPrep[i]&& i<no2)
			{
				maxAmp=invPrep[i];
				indMax=i;
			}
			//ans[i]=tmp;
			//ans[i+1]=(ans[i+1]*tmp-tmp*ans[i+1])/no2;
			//ans[i+1]=0;
		}
		invPrep[0]=ans[0]*ans[0]/no2;
		invPrep[1]=ans[1]*ans[1]/no2;
		//if(invPrep[0]+invPrep[1]>maxAmp)
		//	indMax=0;
		
		invPrep=realft(invPrep,-1);
		
		if(nCalcFreqPhase)
		{
			fFreqPhase= new float [2];
			fFreqPhase[0]= ((float)(indMax)*0.5f/(float)n);
			fFreqPhase[1] = (float)Math.atan2(ans[indMax+1], ans[indMax]);
		}
		
		return invPrep;
	
	}

}
