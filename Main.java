//import org.knowm.xchart.SwingWrapper;
//import org.knowm.xchart.XYChart;
//import org.knowm.xchart.XYChartBuilder;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

class Harmonic {
  private int countOfHarmonics;
  private int limitFrequency;
  private int countOfDescreteCalls;
  private double[][] signalsForAllHarmonics;
  private double[] signalsForResultingHarmonic;

  public Harmonic(int countOfHarmonics, int limitFrequency, int countOfDescreteCalls) {
    this.countOfHarmonics = countOfHarmonics;
    this.limitFrequency = limitFrequency;
    this.countOfDescreteCalls = countOfDescreteCalls;
    this.signalsForAllHarmonics = new double[countOfHarmonics][countOfDescreteCalls];
    this.signalsForResultingHarmonic = new double[countOfDescreteCalls];
  }

  public double[] calculateSignalsForResultingHarmonic() {
    Random r = new Random();
    double A = r.nextDouble();
    double fi = r.nextDouble() * Math.PI;
    for (int i = 0; i < countOfHarmonics; i++) {
      for (int j = 0; j < countOfDescreteCalls; j++) {
        signalsForResultingHarmonic[j] += A * Math.sin(1. * limitFrequency * (i + 1) / countOfHarmonics * j + fi);
      }
    }
    return signalsForResultingHarmonic;
  }

  public double[] calculateSignalsForResultingHarmonic(int n) {
    Random r = new Random();
    double A = r.nextDouble();
    double fi = r.nextDouble() * Math.PI;
    double[] res;
    res = new double[n];
    for (int i = 0; i < countOfHarmonics; i++) {
      for (int j = 0; j < n; j++) {
        res[j] += A * Math.sin(1. * limitFrequency * (i + 1) / countOfHarmonics * j + fi);
      }
    }
    return res;
  }
  
  public double calculateMathExpectation(double[] tmp) {
    double sum = 0;
    for (double signal : tmp) {
      sum += signal;
    }
    return sum / tmp.length;
  }

  public double calculateMathExpectation() {
    double sum = 0;
    for (double signal : getSignalsForResultingHarmonic()) {
      sum += signal;
    }
    return sum / getSignalsForResultingHarmonic().length;
  }

  public double calculateDispersion() {
    double sum = 0;
    double mathExpectation = calculateMathExpectation();
    for (double signal : getSignalsForResultingHarmonic()) {
      sum += Math.pow(signal - mathExpectation, 2);
    }
    return sum / (getSignalsForResultingHarmonic().length - 1);
  }

  public double calculateDispersion(double[] tmp) {
    double sum = 0;
    double mathExpectation = calculateMathExpectation(tmp);
    for (double signal : tmp) {
      sum += Math.pow(signal - mathExpectation, 2);
    }
    return sum / (tmp.length - 1);
  }

  public double[] calculateCorrelation() {
    return calculateCorrelationWithOtherFunc(this);
  }

  public double[] calculateCorrelationWithOtherFunc(Harmonic otherHarmonic) {
    double[] correlation_arr = new double[getCountOfDescreteCalls() / 2];
    double mathExp = calculateMathExpectation();
    double mathExp2 = otherHarmonic.calculateMathExpectation();

    for (int tau = 0; tau < getCountOfDescreteCalls() / 2; tau++) {
      double correlation = 0;
      for (int t = 0; t < getCountOfDescreteCalls() / 2; t++) {
        correlation += (getSignalsForResultingHarmonic()[t] - mathExp)
            * (otherHarmonic.getSignalsForResultingHarmonic()[t + tau] - mathExp2);
      }
      correlation_arr[tau] = correlation / (getCountOfDescreteCalls() - 1);
    }
    return correlation_arr;
  }

  public double[] calculateDFT(double[] signalsOfResultingHarmonic){
        int N = signalsOfResultingHarmonic.length;

        Map<Integer, Double> coefficients = getCoefficientsMap(N);

        double[] dft_real = new double[N];
        double[] dft_imaginary = new double[N];
        double[] dft_final = new double[N];

        for (int p = 0; p < N; p++) {
            for (int k = 0; k < N; k++) {
                dft_real[p] += signalsOfResultingHarmonic[k] * Math.cos(coefficients.get((p*k) % N));
                dft_imaginary[p] += signalsOfResultingHarmonic[k] * Math.sin(coefficients.get((p*k) % N));
            }
            dft_final[p] = Math.sqrt(Math.pow(dft_real[p],2) + Math.pow(dft_imaginary[p],2));
        }
        return dft_final;
    }

  public Map<Integer,Double> getCoefficientsMap(int N){
        Map<Integer, Double> coefficients = new HashMap<>();
        for (int p = 0; p < N; p++) {
            for (int k = 0; k < N; k++) {
                coefficients.putIfAbsent((p*k) % N, 2 * Math.PI * p * k / N);
            }
        }
        return  coefficients;
  }

  public int getCountOfHarmonics() {
    return countOfHarmonics;
  }

  public void setCountOfHarmonics(int countOfHarmonics) {
    this.countOfHarmonics = countOfHarmonics;
  }

  public int getLimitFrequency() {
    return limitFrequency;
  }

  public void setLimitFrequency(int limitFrequency) {
    this.limitFrequency = limitFrequency;
  }

  public int getCountOfDescreteCalls() {
    return countOfDescreteCalls;
  }

  public void setCountOfDescreteCalls(int countOfDescreteCalls) {
    this.countOfDescreteCalls = countOfDescreteCalls;
  }

  public double[][] getSignalsForAllHarmonics() {
    return signalsForAllHarmonics;
  }

  public void setSignalsForAllHarmonics(double[][] signalsForAllHarmonics) {
    this.signalsForAllHarmonics = signalsForAllHarmonics;
  }

  public double[] getSignalsForResultingHarmonic() {
    return signalsForResultingHarmonic;
  }

  public void setSignalsForResultingHarmonic(double[] signalsForResultingHarmonic) {
    this.signalsForResultingHarmonic = signalsForResultingHarmonic;
  }
}

class FastFourier extends Harmonic {
    public FastFourier(int countOfHarmonics, int limitFrequency, int countOfDescreteCalls) {
        super(countOfHarmonics, limitFrequency, countOfDescreteCalls);
    }

    public double[] calculateFFT(double[] signals){
        int N = signals.length;
        double[] fft = new double[N];
        double fourierReal_1 = 0, fourierReal_2 = 0, fourierImaginary_1 = 0, fourierImaginary_2 = 0;

        for (int p = 0; p < N; p++) {
            for (int k = 0; k < N/2 - 1; k++) {
                fourierImaginary_1 += signals[2 * k] * Math.sin((4 * Math.PI * p * k) / N);
                fourierReal_1 += signals[2 * k] * Math.cos((4 * Math.PI * p * k) / N);
                fourierImaginary_2 += signals[2 * k + 1] * Math.sin((2 * Math.PI * p * (k * 2 + 1)) / N);
                fourierReal_2 += signals[2 * k + 1] * Math.cos((2 * Math.PI * p * (k * 2 + 1)) / N);
            }
            fft[p] = Math.sqrt(Math.pow(fourierReal_1 + fourierReal_2, 2) + 
                Math.pow(fourierImaginary_1 + fourierImaginary_2, 2));
        }
        return fft;
    }
}

public class Main {
    public static void main(String[] args) {
        FastFourier harmonic = new FastFourier(10, 1500, 256);

        double[] count = new double[harmonic.getCountOfDescreteCalls()];
        for (int k = 0; k < count.length; k++) {
            count[k] = k;
        }

        //XYChart chart = new XYChartBuilder().width(600).height(400).title("x(t)").xAxisTitle("t").yAxisTitle("x").build();
        double[] signals = harmonic.calculateSignalsForResultingHarmonic();
        //chart.addSeries("Fourier Function", count, harmonic.calculateFFT(signals));

        long[] DFTTime = new long[100];
        long[] FFTTime = new long[100];

        for (int i1 = 0; i1 < 100; i1++) {
            long time = System.currentTimeMillis();
            harmonic.calculateFFT(signals);
            FFTTime[i1] = System.currentTimeMillis() - time;

            time = System.currentTimeMillis();
            harmonic.calculateDFT(signals);
            DFTTime[i1] = System.currentTimeMillis() - time;
        }
        int tmp = 0;
        int tmp1 = 0;
        for (int j = 1; j < FFTTime.length; j++) {
          tmp += FFTTime[0];
          tmp1 += DFTTime[0];
        }

        tmp /= FFTTime.length;
        tmp1 /= DFTTime.length;

        System.out.println("Середній час виконання FFT: " + tmp);
        System.out.println("Середній час виконання DFT: " + tmp1);

        //new SwingWrapper<>(chart).displayChart();
    }
}