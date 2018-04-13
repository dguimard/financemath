package com.guimard.finance.product;


import com.guimard.finance.math.ControlledRandomVariable;

public abstract class PathIndependentOption extends Option {


    public PathIndependentOption(Asset asset, String name) {
        super(asset, name);
    }


    public RandomVariable discountedPayoff() {
        RandomVariable discounted_payoff = new RandomVariable() {

            public double getValue(int t) {
                underlying.timeStep(Flag.RISK_NEUTRAL_PROBABILITY, t, T);
                return currentDiscountedPayoff();
            }

            @Override
            public double conditionalExpectation(int t, double precision, double confidence, int nSignChange) {
                return 0;
            }
        }; // end discounted_payoff

        return discounted_payoff;

    } //end DiscountedPayoff



    public ControlledRandomVariable controlledDiscountedPayoff() {
        ControlledRandomVariable
                controlled_discounted_payoff = new ControlledRandomVariable() {

            /** We use the discounted underlying at expiration as the
             * control variate.
             */
            public double[] getControlledValue(int t) {
                underlying.timeStep(Flag.RISK_NEUTRAL_PROBABILITY, t, T);
                double x = currentDiscountedPayoff(),
                        cv = S[T];                 // control variate

                double[] value_control_variate_pair = {x, cv};
                return value_control_variate_pair;

            } //end getControlledValue


            /** This is simply S^B(t)exp(-q(T-t)).
             *  See book, chapter 3, section 1, equation 3.15.
             */
            public double getControlVariateMean(int t) {
                double[] S = underlying.get_S();
                double q = underlying.get_q();

                return S[t] * Math.exp(-q * (T - t) * dt);
            }

        }; // end controlled_discounted_payoff

        return controlled_discounted_payoff;

    } // end ControlledDiscountedPayoff


    public double minimumVarianceDelta
            (int whichProbability, int t, int nPath) {
        underlying.simulationInit(t);         //sets pathCounter to zero
        double numSum = 0, denomSum = 0;

        if (hasAnalyticPrice()) {
            // compute as E_t(DeltaS(t)DeltaC(t))/E_t(DeltaS(t)^2)

            C[t] = discountedAnalyticPrice(t);
            for (int n = 0; n < nPath; n++) {

                double deltaSt, deltaCt;
                underlying.timeStep(whichProbability, t);
                deltaSt = S[t + 1] - S[t];
                C[t + 1] = discountedAnalyticPrice(t + 1);
                deltaCt = C[t + 1] - C[t];

                numSum += deltaSt * deltaCt;
                denomSum += deltaSt * deltaSt;
            } // end for n

            return numSum / denomSum;
        } // end if

        else if (whichProbability == RISK_NEUTRAL_PROBABILITY) {
            // compute as E_t(DeltaS(t)h)/E_t(DeltaS(t)^2)

            for (int n = 0; n < nPath; n++) {

                double deltaSt, h;
                underlying.timeStep(whichProbability, t);
                deltaSt = S[t + 1] - S[t];
                underlying.timeStep(whichProbability, t + 1, T);
                h = currentDiscountedPayoff();

                numSum += deltaSt * h;
                denomSum += deltaSt * deltaSt;
            } // end for n

            return numSum / denomSum;
        } // end else if

        else { // whichProbability==MARKET_PROBABILITY
            // compute gain as E_t(DeltaS(t)DeltaC(t))/E_t(DeltaS(t)^2)
            // where C[t], C[t+1] are computed using path branching

            C[t] = discountedMonteCarloPrice(t, nPath);
            for (int n = 0; n < nPath; n++) {

                double deltaSt, deltaCt;
                underlying.timeStep(whichProbability, t);
                deltaSt = S[t + 1] - S[t];
                C[t + 1] = discountedMonteCarloPrice(t + 1, nPath);
                deltaCt = C[t + 1] - C[t];

                numSum += deltaSt * deltaCt;
                denomSum += deltaSt * deltaSt;
            } // end for n

            return numSum / denomSum;
        } // end else

    } // end minimumVarianceDelta


    public double minimumVarianceDelta
            (int whichProbability, int t, int nPath, Trigger rebalance) {
        underlying.simulationInit(t);         //sets pathCounter to zero

        if (hasAnalyticPrice()) {
            // compute as E_t(DeltaS(t)DeltaC(t))/E_t(DeltaS(t)^2)

            C[t] = discountedAnalyticPrice(t);
            double numSum = 0, denomSum = 0;
            for (int n = 0; n < nPath; n++) {

                double deltaSt, deltaCt;
                //path computation to time s of next hedge trade
                int s = underlying.pathSegment(whichProbability, t, rebalance);
                deltaSt = S[s] - S[t];
                C[s] = discountedAnalyticPrice(s);
                deltaCt = C[s] - C[t];

                numSum += deltaSt * deltaCt;
                denomSum += deltaSt * deltaSt;
            } // end for n

            return numSum / denomSum;
        } // end if

        else if (whichProbability == RISK_NEUTRAL_PROBABILITY) {
            // compute as E_t(DeltaS(t)h)/E_t(DeltaS(t)^2)

            double numSum = 0, denomSum = 0;
            for (int n = 0; n < nPath; n++) {

                double deltaSt, h;
                //path computation to time s of next hedge trade
                int s = underlying.pathSegment(whichProbability, t, rebalance);
                deltaSt = S[s] - S[t];
                underlying.timeStep(whichProbability, s, T);
                h = currentDiscountedPayoff();

                numSum += deltaSt * h;
                denomSum += deltaSt * deltaSt;
            } // end for n

            return numSum / denomSum;
        } // end else if

        else { // whichProbability==MARKET_PROBABILITY
            // compute gain as E_t(DeltaS(t)DeltaC(t))/E_t(DeltaS(t)^2)
            // where C[t], C[t+1] are computed using path branching

            C[t] = discountedMonteCarloPrice(t, nPath);
            double numSum = 0, denomSum = 0;
            for (int n = 0; n < nPath; n++) {

                double deltaSt, deltaCt;
                //path computation to time s of next hedge trade
                int s = underlying.pathSegment(whichProbability, t, rebalance);
                deltaSt = S[s] - S[t];
                C[s] = discountedMonteCarloPrice(s, nPath);
                deltaCt = C[s] - C[t];

                numSum += deltaSt * deltaCt;
                denomSum += deltaSt * deltaSt;
            } // end for n

            return numSum / denomSum;
        }
    }

    // to improve

}