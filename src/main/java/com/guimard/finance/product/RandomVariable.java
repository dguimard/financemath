package com.guimard.finance.product;

public abstract class RandomVariable {

    public abstract double getValue(int t);

    public double  conditionalExpectation(int t, int N)
    {
        double sum_X=0;
        for(int n=0;n<N;n++)sum_X+=getValue(t);

        return sum_X/N;
    }

    public double expectation(int N)
    {
        return conditionalExpectation(0,N);
    }

    public abstract double conditionalExpectation(int t, double precision, double confidence, int nSignChange);
}
