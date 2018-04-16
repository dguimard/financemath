package com.guimard.finance.math;


public abstract class RandomVariable {


    boolean hasAnalyticMean;


    boolean hasConditionalAnalyticMean;


    boolean hasAnalyticVariance;


    boolean hasConditionalAnalyticVariance;


    boolean hasAnalyticMoment;


    boolean hasConditionalAnalyticMoment;


    boolean hasAnalyticCentralMoment;


    boolean hasConditionalAnalyticCentralMoment;
    EmpiricalDistribution empiricalDist;
    boolean empiricalDistributionIsInitialized;

    public RandomVariable() {
        hasAnalyticMean = false;
        hasConditionalAnalyticMean = false;
        hasAnalyticVariance = false;
        hasConditionalAnalyticVariance = false;
        hasAnalyticMoment = false;
        hasConditionalAnalyticMoment = false;
        hasAnalyticCentralMoment = false;
        hasConditionalAnalyticCentralMoment = false;
        empiricalDistributionIsInitialized = false;
    }


    public RandomVariable
            (boolean hasAnalyticMean,
             boolean hasConditionalAnalyticMean,
             boolean hasAnalyticVariance,
             boolean hasConditionalAnalyticVariance,
             boolean hasAnalyticMoment,
             boolean hasConditionalAnalyticMoment,
             boolean hasAnalyticCentralMoment,
             boolean hasConditionalAnalyticCentralMoment) {
        this.hasAnalyticMean = hasAnalyticMean;
        this.hasConditionalAnalyticMean = hasConditionalAnalyticMean;
        this.hasAnalyticVariance = hasAnalyticVariance;
        this.hasConditionalAnalyticVariance = hasConditionalAnalyticVariance;
        this.hasAnalyticMoment = hasAnalyticMoment;
        this.hasConditionalAnalyticMoment = hasConditionalAnalyticMoment;
        this.hasAnalyticCentralMoment = hasAnalyticCentralMoment;
        this.hasConditionalAnalyticCentralMoment =
                hasConditionalAnalyticCentralMoment;
        empiricalDistributionIsInitialized = false;
    }

    public boolean get_hasAnalyticMean() {
        return hasAnalyticMean;
    }

    public boolean get_hasConditionalAnalyticMean() {
        return hasConditionalAnalyticMean;
    }

    public boolean get_hasAnalyticVariance() {
        return hasAnalyticVariance;
    }

    public boolean get_hasConditionalAnalyticVariance() {
        return hasConditionalAnalyticVariance;
    }

    public boolean get_hasAnalyticMoment() {
        return hasAnalyticMoment;
    }

    public boolean get_hasConditionalAnalyticMoment() {
        return hasConditionalAnalyticMoment;
    }

    public boolean get_hasAnalyticCentralMoment() {
        return hasAnalyticCentralMoment;
    }

    public boolean get_hasConditionalAnalyticCentralMoment() {
        return hasConditionalAnalyticCentralMoment;
    }

    public void setHasAnalyticMean(boolean value) {
        hasAnalyticMean = value;
    }

    public void setHasAnalyticVariance(boolean value) {
        hasAnalyticVariance = value;
    }

    public abstract double getValue(int t);

    public double analyticConditionalMean(int t) {
        String message =
                "RandomVariable.analyticConditionalMean()\n:" +
                        "No analytic formula for the conditional mean implemented,\n" +
                        "use Monte Carlo instead. Exiting...";

        //to log it properly
        return 0;  // for the compiler
    }

    public double analyticMean() {
        String message =
                "RandomVariable.analyticMean()\n:" +
                        "No analytic formula for the mean implemented,\n" +
                        "use Monte Carlo instead. Exiting...";
        //to log it properly

        return 0; // for the compiler
    }

    public double analyticConditionalVariance(int t) {
        String message =
                "RandomVariable.analyticConditionalVariance()\n:" +
                        "No analytic formula for the conditional variance implemented,\n" +
                        "use Monte Carlo instead. Exiting...";

        return 0; // for the compiler
    }

    public double analyticVariance() {
        String message =
                "RandomVariable.analyticVariance()\n:" +
                        "No analytic formula for the variance implemented,\n" +
                        "use Monte Carlo instead. Exiting...";

        //to log it properly
        return 0; // for the compiler
    }

    public double analyticConditionalMoment(int t, int n) {
        String message =
                "RandomVariable.analyticConditionalMoment()\n:" +
                        "No analytic formula for the conditional moments implemented,\n" +
                        "use Monte Carlo instead. Exiting...";
        System.out.println(message);
        System.exit(0);

        return 0; // for the compiler


    }

    public double analyticMoment(int n) {
        String message =
                "RandomVariable.analyticMean()\n:" +
                        "No analytic formula for the moments implemented,\n" +
                        "use Monte Carlo instead. Exiting...";
        System.out.println(message);
        System.exit(0);

        return 0; // for the compiler

    }

    public double analyticConditionalCentraMoment(int t, int n) {
        String message =
                "RandomVariable.analyticConditionalCentralMoment()\n:" +
                        "No analytic formula for the conditional central moments implemented,\n" +
                        "use Monte Carlo instead. Exiting...";
        System.out.println(message);
        System.exit(0);

        return 0; // for the compiler

    }

    public double analyticCentralMoment(int n) {
        String message =
                "RandomVariable.analyticCentralMoment()\n:" +
                        "No analytic formula for the central moments implemented,\n" +
                        "use Monte Carlo instead. Exiting...";
        System.out.println(message);
        System.exit(0);

        return 0; // for the compiler
    }

    public double
    conditionalExpectation(int t, int N) {
        double sum_X = 0;
        for (int n = 0; n < N; n++) sum_X += getValue(t);

        return sum_X / N;
    }


    public double expectation(int N) {
        return conditionalExpectation(0, N);
    }

    public double
    conditionalExpectation(int t, int N, boolean report) {
        double sum_X = 0;
        for (int n = 0; n < N; n++) {

            sum_X += getValue(t);
            if (n % (N / 30) == 0) System.out.print("*");
        }

        return sum_X / N;
    }


    public double expectation(int N, boolean report) {
        return conditionalExpectation(0, N, report);
    }

    public double[]
    conditionalMeanAndStandardDeviation(int t, int N) {
        double sum_X = 0,                 //X_1+X_2+...+X_n
                sum_X_square = 0,          //X_1^2+X_2^2+...+X_n^2
                mean_X,                  //mean
                mean_X_square,           //mean of the squares
                variance_X,              //variance
                sigma_X;                 //standard deviation


        for (int n = 0; n < N; n++) {
            double x = getValue(t);
            sum_X += x;
            sum_X_square += x * x;
        }

        mean_X = sum_X / N;
        mean_X_square = sum_X_square / N;
        variance_X = mean_X_square - mean_X * mean_X;
        sigma_X = Math.sqrt(variance_X);

        double[] results = {mean_X, sigma_X};
        return results;
    }

    public double[] meanAndStandardDeviation(int N) {
        return conditionalMeanAndStandardDeviation(0, N);
    }

    public double[]
    conditionalMeanAndStandardDeviation(int t, int N, int sampleGroupSize) {
        int nSampleGroups = N / sampleGroupSize;  //note integer division

        double group_sum,           //sum of sample values in a group
                group_mean,          //mean of the current sample group
                sum = 0,               //sum of group means
                sumSquares = 0,        //sum of squares of group means
                mean,                //mean
                meanSquares,         //mean of squares of group means
                variance,            //current variance of group means
                sigma;               //current standard deviation of group means


        for (int n = 0; n < nSampleGroups; n++) {
            //compute the mean over the next sample group
            group_sum = 0;
            for (int k = 0; k < sampleGroupSize; k++) group_sum += getValue(t);
            group_mean = group_sum / sampleGroupSize;

            //update sum of group means, group mean squares
            sum += group_mean;
            sumSquares += group_mean * group_mean;
        }

        mean = sum / nSampleGroups;
        meanSquares = sumSquares / nSampleGroups;
        variance = meanSquares - mean * mean;
        sigma = Math.sqrt(variance);

        double[] results = {mean, sigma};
        return results;

    } //end meanAndStandardDeviation

    public double[]
    meanAndStandardDeviation(int N, int sampleGroupSize) {
        return conditionalMeanAndStandardDeviation(0, N, sampleGroupSize);
    }

    public double
    conditionalExpectation(int t, int N, int m) {

        double sum = 0;
        for (int n = 0; n < N; n++) {
            //progress report every m samples

            sum += getValue(t);
        }

        return sum / N;
    }

    public double expectation(int N, int m) {
        return conditionalExpectation(0, N, m);
    }

    public double[] conditionalMeanAndStandardDeviation
            (int t, int N, int m, int sampleGroupSize) {
        int nSampleGroups = N / sampleGroupSize;  //note integer division

        double group_sum,           //sum of sample values in a group
                group_mean,          //mean of the current sample group
                sum = 0,               //sum of group means
                sumSquares = 0,        //sum of squares of group means
                mean,                //mean
                meanSquares,         //mean of squares of group means
                variance,            //current variance of group means
                sigma;               //current standard deviation of group means

        //initialize progress report

        long before = System.currentTimeMillis();

        for (int n = 0; n < nSampleGroups; n++) {

            group_sum = 0;
            for (int k = 0; k < sampleGroupSize; k++) group_sum += getValue(t);
            group_mean = group_sum / sampleGroupSize;

            //update sum of group means, group mean squares
            sum += group_mean;
            sumSquares += group_mean * group_mean;
        }

        mean = sum / nSampleGroups;
        meanSquares = sumSquares / nSampleGroups;
        variance = meanSquares - mean * mean;
        sigma = Math.sqrt(variance);

        double[] results = {mean, sigma};
        return results;

    } //end meanAndStandardDeviation

    public double[] meanAndStandardDeviation
            (int N, int m, int sampleGroupSize) {
        return conditionalMeanAndStandardDeviation
                (0, N, m, sampleGroupSize);
    }

    public double
    conditionalExpectation(int t, double precision, double confidence) {
        double sum = 0,               //X_1+X_2+...+X_n
                sumSquares = 0,        //X_1^2+X_2^2+...+X_n^2
                mean,                //current mean
                meanSquares,         //current mean of the squares
                variance,            //current sample variance
                sigma;               //current standard deviation
        //sigma(X,N)

        int n = 0;
        boolean done = false;


        while ((n < 100) || (!done)) {
            n++;
            double x = getValue(t);      //x=X_n
            sum += x;
            sumSquares += x * x;
            //check for termination every 100 samples
            if (n % 100 == 99) {
                mean = sum / n;
                meanSquares = sumSquares / n;
                variance = meanSquares - mean * mean;
                sigma = Math.sqrt(variance);

                double f = precision * Math.sqrt(n) / sigma;
                done = (2 * FinMath.N(f) > confidence);

                if (n == 1000000) done = true;
            } //end if

        } //end while

        return sum / n;

    }

    public double  expectation(double precision, double confidence) {
        return conditionalExpectation(0, precision, confidence);
    }

    public double conditionalExpectation (int t, double precision, double confidence, int sampleGroupSize) {
        double group_sum,           //sum of sample values in a group
                group_mean,          //mean of the current sample group
                sum = 0,               //sum of group means
                sumSquares = 0,        //sum of squares of group means
                mean,                //mean
                meanSquares,         //mean of squares of group means
                variance,            //current variance of group means
                sigma;               //current standard deviation of group means

        int n = 0;
        boolean done = false;


        while ((n < 100) || (!done)) {
            n++;
            //compute the average over the next sample group
            group_sum = 0;
            for (int k = 0; k < sampleGroupSize; k++) group_sum += getValue(t);
            group_mean = group_sum / sampleGroupSize;

            sum += group_mean;
            sumSquares += group_mean * group_mean;

            //check for termination every 100 sample groups
            if (n % 100 == 99) {
                mean = sum / n;
                meanSquares = sumSquares / n;
                variance = meanSquares - mean * mean;
                sigma = Math.sqrt(variance);

                double f = precision * Math.sqrt(n) / sigma;
                done = (2 * FinMath.N(f) > confidence);

                if (n == 1000000) done = true;
            } //end if

        } //end while

        return sum / n;

    } //end getExpectation

    public double expectation (double precision, double confidence, int sampleGroupSize) {
        return conditionalExpectation(0, precision, confidence, sampleGroupSize);
    }

    public double conditionalVariance(int t, int N) {
        double sum_X = 0,
                sum_X_square = 0,
                mean_X,
                mean_X_square;

        for (int n = 0; n < N; n++) {
            double x = getValue(t);

            sum_X += x;
            sum_X_square += x * x;
        }

        mean_X = sum_X / N;
        mean_X_square = sum_X_square / N;

        return mean_X_square - mean_X * mean_X;

    }

    public double variance(int N) {
        return conditionalVariance(0, N);
    }

    public double conditionalMoment(final int t, final int n, final int N) {
        /** allocate the random variable X^n
         *  conditioned on information available at time t.
         */
        RandVariable Xn = new RandVariable() {

            //the value is X^n
            public double nextValue(int t) {
                double x = getValue(t),
                        y = x;
                for (int j = 0; j < n - 1; j++) x *= y;   //will have to do for x^n
                return x;
            }

        }; //end Xn

        return Xn.conditionalMean(t, N);

    }

    double moment(final int n, final int N) {
        return conditionalMoment(0, n, N);
    }

    public RandVariable centered_X(final int t, final int n, final int N) {
        double m = 0;
        if (hasConditionalAnalyticMean) m = analyticConditionalMean(t);
        else m = conditionalExpectation(t, N);

        // make final for use in inner class
        final double mean_X = m;

        return new RandVariable() {

            public double nextValue(int s) {
                double x = getValue(t) - mean_X, y = x;
                for (int i = 0; i < n - 1; i++) x *= y;

                return x;
            } // end nextValue

        }; // end return new

    } // end centered_X

    double conditionalCentralMoment(final int t, final int n, final int N) {
        return centered_X(t, n, N).mean(N);

    }//end conditionalCentralMoment

    double centralMoment(final int n, final int N) {
        double mu = 0;  //the mean

        if (hasAnalyticMean) mu = analyticMean();
        else mu = expectation(N);

        RandVariable Xcn = new RandVariable(mu) {

            public double nextValue(int t) {
                double x = getValue(t) - mean,
                        y = x;
                for (int j = 0; j < n - 1; j++) x *= y;   //will have to do for (x-mu)^n
                return x;
            }
        };//end Xcn

        return Xcn.mean(N);

    }//end moment

    public boolean get_empiricalDistributionIsInitialized() {
        return empiricalDistributionIsInitialized;
    }


    public void initEmpiricalDistribution(int N) {
        if (!empiricalDistributionIsInitialized)
            empiricalDist = conditionalEmpiricalDistribution(0, N);

        empiricalDistributionIsInitialized = true;
    }


    public void fillSampleSet(int N) {
        if (!empiricalDistributionIsInitialized) initEmpiricalDistribution(N);
        else {
            int currentSize = empiricalDist.get_nSamples(),
                    moreSamples = N - currentSize;
            if (moreSamples > 0) {
                for (int i = 0; i < moreSamples; i++) empiricalDist.add(getValue(0));
                empiricalDist.set_nSamples(N);
            } // end if

        } // end else

    }

    public EmpiricalDistribution conditionalEmpiricalDistribution(int t, int N) {
        return new EmpiricalDistribution(this, t, N);
    }


    public EmpiricalDistribution empiricalDistribution(int N) {
        fillSampleSet(N);
        return empiricalDist;
    }


    public double quantile(double phi, int N) {
        fillSampleSet(N);
        return empiricalDist.quantile(phi);
    }


    public double cumulativeDistributionFunction(double x, int N) {
        fillSampleSet(N);
        return empiricalDist.quantileInverse(x);
    }


    public RandomVariable scale(final double lambda) {
        final RandomVariable X = this;
        return new RandomVariable() {

            public double getValue(int t) {

                return lambda * X.getValue(t);

            }
        }; // end return new

    } // end plus


    public RandomVariable plus(final RandomVariable Y) {
        final RandomVariable X = this;
        return new RandomVariable() {

            public double getValue(int t) {

                return X.getValue(t) + Y.getValue(t);

            }
        }; // end return new

    }

    public RandomVariable minus(final RandomVariable Y) {
        final RandomVariable X = this;
        return new RandomVariable() {

            public double getValue(int t) {

                return X.getValue(t) - Y.getValue(t);

            }
        }; // end return new

    }

    public RandomVariable mult(final RandomVariable Y) {
        final RandomVariable X = this;
        return new RandomVariable() {

            public double getValue(int t) {

                return X.getValue(t) * Y.getValue(t);

            }
        }; // end return new

    }

    public RandomVariable div(final RandomVariable Y) {
        final RandomVariable X = this;
        return new RandomVariable() {

            public double getValue(int t) {

                return X.getValue(t) / Y.getValue(t);

            }
        }; // end return new

    }
}