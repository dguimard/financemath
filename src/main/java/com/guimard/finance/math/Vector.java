package com.guimard.finance.math;

public class Vector {


    public static void setZero(double[] X) {
        for (int j = 0; j < X.length; j++) X[j] = 0;
    }


    public static double[] add(double[] X, double[] Y) {
        for (int j = 0; j < X.length; j++) X[j] += Y[j];
        return X;
    }


    public static double[] subtract(double[] X, double[] Y) {
        for (int j = 0; j < X.length; j++) X[j] -= Y[j];
        return X;
    }

    public static double[] mult(double[] X, double[] Y) {
        for (int j = 0; j < X.length; j++) X[j] *= Y[j];
        return X;

    }


    public static double[] scalarMult(double f, double[] X) {
        for (int j = 0; j < X.length; j++) X[j] *= f;
        return X;

    }


    public static double[] sqrt(double[] X) {
        for (int j = 0; j < X.length; j++) X[j] = Math.sqrt(X[j]);
        return X;
    }

    public static void print(double[] X) {

        for (int j = 0; j < X.length; j++) System.out.println(X[j]);


    }


}
