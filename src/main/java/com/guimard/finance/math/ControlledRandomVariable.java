package com.guimard.finance.math;

public abstract class ControlledRandomVariable extends RandomVariable {


    public static final int nBeta = 2000;

    private double beta_at_time_zero;

    public ControlledRandomVariable() {
    }

    public abstract double[] getControlledValue(int t);

    public abstract double getControlVariateMean(int t);

    public double getValue(int t) {
        return getControlledValue(t)[0];
    }


    public void controlVariateMeanTest() {
        String str = "Control variate mean:\n";
        System.out.println(str);

        // analytic control variate mean
        double cv_mean = getControlVariateMean(0);
        cv_mean = 1.0 * Math.round(cv_mean * 10000) / 10000;

        str = "analytic: " + cv_mean + "\n";
        System.out.println(str);

        // the control variate as a random variable
        RandomVariable CVM = new RandomVariable() {

            public double getValue(int t) {
                return getControlledValue(t)[1];
            }

        }; // end CVM

        // Monte Carlo control variate mean
        for (int i = 0; i < 10; i++) {
            double cvmc_mean = CVM.expectation(50000);
            cvmc_mean = 1.0 * Math.round(cvmc_mean * 10000) / 10000;

            str = "Monte Carlo: " + cvmc_mean;
            System.out.println(str);
        } // end for i

    } // end controlVariateMeanTest


    public double betaCoefficient(int t, int N) {

        double sum_X = 0,
                sum_Y = 0,
                sum_XX = 0,
                sum_YY = 0,
                sum_XY = 0;

        for (int n = 0; n < N; n++) {
            double[] d = getControlledValue(t);
            double x = d[0],
                    y = d[1];

            sum_X += x;
            sum_Y += y;
            sum_XX += x * x;
            sum_XY += x * y;
        }

        return (N * sum_XY - sum_X * sum_Y) / (N * sum_XX - sum_X * sum_X);

    } //end betaCoefficient


    public double betaCoefficient(int N) {
        return betaCoefficient(0, N);
    }


    public double correlationWithControlVariate(int t, int N) {
        //allocate X (this) as a local 2 dimensional random vector
        RandomVector X = new RandomVector(2) {

            public double[] getValue(int t) {
                return getControlledValue(t);
            }
        }; //end definition of X

        return X.conditionalCorrelation(0, 1, t, N);

    } //end correlationWithControlVariate


    public double correlationWithControlVariate(int N) {
        return correlationWithControlVariate(0, N);
    }


    public RandomVariable controlled_X(final int t) {

        final double
                beta = betaCoefficient(t, nBeta),      // beta coefficient
                mean_y = getControlVariateMean(t);    // control variate mean E_t(Y)


        // random variable Xc=X-beta(t)(Y-E_t(Y))
        return new RandomVariable() {

            // random draw from controlled_X
            // avoid name collision with of time s with time t!!!
            public double getValue(int s) {
                double[] xc = getControlledValue(t);
                double x = xc[0],
                        y = xc[1];

                return x - beta * (y - mean_y);

            } //end getValue

        }; // end new

    } // end controlled_X


    public double conditionalExpectation(int t, int N) {
        return controlled_X(t).expectation(N);
    }


    public double expectation(int N) {
        return conditionalExpectation(0, N);
    }


    public double
    conditionalExpectation(int t, int N, int m) {
        return controlled_X(t).expectation(N, m);
    }


    public double
    expectation(int N, int m) {
        return conditionalExpectation(0, N, m);
    }


    public double
    conditionalExpectation(int t, double precision, double confidence) {
        return controlled_X(t).expectation(precision, confidence);
    }


    public double expectation(double precision, double confidence) {
        return conditionalExpectation(0, precision, confidence);
    }


    public double conditionalExpectation
            (int t, double precision, double confidence, int sampleGroupSize) {
        return controlled_X(t).
                expectation(precision, confidence, sampleGroupSize);
    }


    public double expectation
            (double precision, double confidence, int sampleGroupSize) {
        return conditionalExpectation(0, precision, confidence, sampleGroupSize);
    }


}