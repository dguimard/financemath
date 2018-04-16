package com.guimard.finance.product;


public abstract class Asset {


    int T;

    /**
     * size of time step
     */
    double dt;     // size of time step

    /**
     * asset price at time t=0
     */
    double S_0;

    double q;


    int nSignChange;


    boolean volatilityIsDeterministic;


    double[] Z;

    double[] B;


    double[] S;


    public Asset(int T, double dt, double S_0, double q, int nSignChange) {
        this.T = T;
        this.dt = dt;
        this.S_0 = S_0;
        this.q = q;
        this.nSignChange = nSignChange;
        volatilityIsDeterministic = false;              // default

        //allocate path arrays
        Z = new double[T];
        B = new double[T + 1];
        S = new double[T + 1];

        //initialize paths
        B[0] = 1;
        S[0] = S_0;

    } // end constructor


    public int get_T() {
        return T;
    }


    public double get_dt() {
        return dt;
    }

    public double get_S_0() {
        return S_0;
    }

    public double get_q() {
        return q;
    }

    public int get_nSignChange() {
        return nSignChange;
    }

    public boolean get_volatilityIsDeterministic() {
        return volatilityIsDeterministic;
    }

    public double[] get_Z() {
        return Z;
    }

    public double[] get_B() {
        return B;
    }

    public double[] get_S() {
        return S;
    }

    public double dividendReductionFactor(int t) {
        return Math.exp(-q * (T - t) * dt);
    }


    public double forwardPrice(int t) {
        return S[t] * B[T] * dividendReductionFactor(t);
    }


    public double get_sigmaSqrtdt(int t) {
        System.err.println("sigma*sqrt(dt) undefined in present context");
        System.exit(0);
        return 0;
    }


    public double Sigma(int t) {
        System.err.println("Sigma(t) undefined in present context");
        System.exit(0);
        return 0;
    }


    public void simulationInit(int t) {
    }


    public abstract void timeStep(int whichProbability, int t);


    public void timeStep(int whichProbability, int t, int s) {
        for (int u = t; u < s; u++) timeStep(whichProbability, u);
    }


    public void newPathBranch(int whichProbability, int t) {
        for (int u = t; u < T; u++) timeStep(whichProbability, u);
    }


    public void newPath(int whichProbability) {
        for (int t = 0; t < T; t++) timeStep(whichProbability, t);
    }


    public int pathSegment(int whichProbability, int t, Trigger trg) {
        int s = t;
        do {
            timeStep(whichProbability, s);
            s++;
        }
        while ((s < T) && !trg.isTriggered(t, s));

        return s;
    }


}