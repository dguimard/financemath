package com.guimard.finance.math;



public class EmpiricalDistribution extends hep.aida.bin.DynamicBin1D{
    
    
    int nSamples;       // number of samples.
   


    public EmpiricalDistribution(RandomVariable X, int t, int N) 
    {
        super();
        for(int j=0;j<N;j++)add(X.getValue(t));
        nSamples=N;
        
    }


    public int get_nSamples(){ return nSamples; }
    

    public void set_nSamples(int N){ nSamples=N; }

}