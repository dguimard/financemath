package com.guimard.finance.product;


import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class Call extends PathIndependentOption {
    private static final Logger logger = LoggerFactory.getLogger(Call.class);
    private double strikePrice;

    public Call(double K, Asset asset)
    {
        super(asset,"Call");
        this.strikePrice=K;

        hasAnalyticPrice=false;
        C[0]=controlledDiscountedMonteCarloPrice(40000);

    }

    public Call(double K, Asset asset, int forget)
    {
        super(asset,"Call");
        this.strikePrice=K;

    }


    public double get_K(){ return strikePrice; }

    public double currentDiscountedPayoff()
    {
        double x=S[T]-strikePrice/B[T];
        return (x>0)? x:0;
    }
}
