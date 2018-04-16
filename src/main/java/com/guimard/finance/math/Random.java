package com.guimard.finance.math;

import cern.jet.random.Uniform;
import cern.jet.random.engine.MersenneTwister;

public class Random {
    public static final Uniform
            uniform_1 = new Uniform(new MersenneTwister(113));
    
    public static final Uniform
            uniform_2 = new Uniform(new MersenneTwister(2113));
    
    static double Y = 0;         //second normal deviate from Box-Muller
    static boolean New = true;   //flag indicating wether a new normal deviate

   
    public static double U1() {
        return uniform_1.nextDouble();
    }

    public static double U2() {
        return uniform_2.nextDouble();
    }

    
    public static int Sign() {
        double x = U1();
        return (x > 0.5) ? 1 : -1;
    }


    public static int Sign(double p) {
        double x = U1();
        return (x < p) ? 1 : -1;
    }

    public static double STN() {
        return FinMath.N_Inverse(uniform_1.nextDouble());

    } 
} 