package com.guimard.finance.product;


public abstract class Trigger {

    int T;   

    public Trigger(int T)  //can only be called from concrete subclasses
    {
        this.T = T;
    }


    public abstract boolean isTriggered(int t, int s);


    public int nextTime(int t) {
        int s = t;
        while (!isTriggered(t, s)) s++;
        return s;
    }


}