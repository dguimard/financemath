package com.guimard.finance.product;


 


public abstract class Asset{
        
    /** time steps to horizon */
    int T;         
    
    /** size of time step */
    double dt;     // size of time step
    
    /** asset price at time t=0 */
    double S_0;    
    
    /** constant dividend yield */
    double q;      
    
    /** number of times the Z-increments are reused through sign changes.
     */
    int nSignChange;
    
    /** flag indicating wether the volatility is deterministic
     */
    boolean volatilityIsDeterministic;
    
    /** array of Z-increments driving the price path */
    double[] Z;
    
    /** price path of the riskfree bond */
    double[] B;
    
    /** the discounted asset price path S[t]=S^B(t*dt) */
    double[] S;
    


    public Asset(int T, double dt, double S_0, double q, int nSignChange)
    {
        this.T=T; this.dt=dt; this.S_0=S_0;
        this.q=q; this.nSignChange=nSignChange;
        volatilityIsDeterministic=false;              // default
        
        //allocate path arrays
        Z=new double[T];
        B=new double[T+1];
        S=new double[T+1];
        
        //initialize paths
        B[0]=1;
        S[0]=S_0;
        
    } // end constructor
          
    
    

    public int get_T(){ return T; }
    

    public double get_dt(){ return dt; }

    public double get_S_0(){ return S_0; }
    
    /** <p>Constant dividend yield.</p>
     */
    public double get_q(){ return q; }
    
    /** <p>Number of times the <a href=#Zinc>Z-increments<a> are reused through 
     * random sign changes before new increments are generated.</p>
     */
    public int get_nSignChange(){ return nSignChange; } 
    
    /** <p>True if the volatility of the asset is deterministic,
     *  false otherwise.</p>
     */
    public boolean get_volatilityIsDeterministic()
    { return volatilityIsDeterministic; }
    
    /** <p>Reference to the array Z[ ] of standard normal increments driving
     *  the current asset price path.</p>
     */
    public double[] get_Z(){ return Z; }
    
    /** <p>Reference to the array B[ ] containing the riskfree bond.</p>
     */
    public double[] get_B(){ return B; }
    
    /** <p>Reference to the array S[ ] containing the <i>discounted</i>
     *  asset price path.</p>
     */
    public double[] get_S(){ return S; } 
  
    
        
    /** <p>Reduces price by future dividends.</p>
     *
     * @param t Current time.
     */
    public double dividendReductionFactor(int t)
    {
        return Math.exp(-q*(T-t)*dt);  }
    
    
    /** <p>Forward price at horizon T.</p>
     *
     * @param t Current time.
     */
    public double forwardPrice(int t)
    {
        return S[t]*B[T]*dividendReductionFactor(t);  }
    
    
    /** <p><code>sigma(t)*sqrt(dt)</code>, where sigma(t) is the volatility of 
     * the asset.</p>
     *
     * <p>Polymorphic link to subclasses implementing this method.
     * Default implementation is an error message.</p>
     *
     * @param t Current time.
     */
    public double get_sigmaSqrtdt(int t)
    { 
        System.err.println("sigma*sqrt(dt) undefined in present context");
        System.exit(0);
        return 0;
     }
   
    
    /** <p>sqrt(\int_t^T sigma^2(u)du), where sigma(u) is the volatility 
     * of the asset. For a constant volatility asset this is simply
     * sigma*sqrt(tau), where tau is time to the horizon.</p>
     *
     * <p>Polymorphic link to subclasses implementing this method.
     * Default implementation is an error message.</p>
     *
     * @param t Current time.
     */
    public double Sigma(int t)          
    { 
        System.err.println("Sigma(t) undefined in present context");
        System.exit(0);
        return 0;
     }                                     
  
    /** <p>Sets up a path simulation (t=0) or a simulation of 
     * branches of an existing path (t>0, conditional expectations).
     * Default: nothing to do.</p>
     *
     * @param t Time of branching.
     */
    public void simulationInit(int t){ } 


    public abstract void timeStep(int whichProbability, int t);
    
    

    public void timeStep(int whichProbability, int t, int s)
    {
        for(int u=t;u<s;u++)timeStep(whichProbability,u);  } 



    public void newPathBranch(int whichProbability, int t)
    {
        for(int u=t;u<T;u++)timeStep(whichProbability,u);  }
    
    
    /** <p>New independent path of riskfree bond and discounted asset
     *  (driven by new independent <a href=#Zinc>Z-increments<a>).</p> 
     *
     *  <p>Not implemented as newPathBranch(whichProbability,0) since this 
     *  method may be overridden in subclasses so as to generate groups of 
     *  dependent paths.</p>
     *
     * @param whichProbability Probability for simulation (market/risk neutral).
     */
    public void newPath(int whichProbability)
    {
        for(int t=0;t<T;t++)timeStep(whichProbability,t);  }
      


   public int pathSegment(int whichProbability, int t, Trigger trg)
   {
      int s=t;
      do{ timeStep(whichProbability,s);
          s++; }
      while( (s<T)&&!trg.isTriggered(t,s) );
      
      return s;
   }
    
    
} // end asset

