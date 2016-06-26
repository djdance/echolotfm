package ru.echolotfm.fingerprint;

import android.content.Context;
import android.database.Cursor;
import android.database.sqlite.SQLiteDatabase;
import android.database.sqlite.SQLiteException;
import android.graphics.*;
import android.os.Environment;
import android.widget.ImageView;
import com.musicg.dsp.FastFourierTransform;
import android.media.AudioFormat;
import android.media.AudioRecord;
import android.media.MediaRecorder.AudioSource;
import android.util.Log;

import com.musicg.dsp.WindowFunction;
import com.sun.media.sound.FFT;

import java.io.File;

public class ExtAudioRecorder {
    String TAG="djd";
    public ImageView pole;
    ImageView[] poleN;
    public Bitmap poleBM,poleBM1;
    Bitmap[] poleBMn;
    public Canvas poleC;
    Canvas[] poleCN;

    public Paint poleP,poleP1;
    public int recordMode=0;
    public int recordedFPs=0;
    byte[] Bark;
    int ww,hh;
    int numFrames=1;
    double averAmp, lastAverageAmp =0;
    double maxIntens, minIntens, averIntens;//,sumFullIntensAmp;
    SQLiteDatabase db;
    //private FingerprintProperties fingerprintProperties=FingerprintProperties.getInstance();
    int mycolor;
    double perceptualDiff;
    double minValidAmp = 0.00000000001;
    double[][][] barkArray;
    int lowFreq,highFreq;
    FastFourierTransform fft;
    float barkToFrequnits;
    long minA1plusB1 =99999;
    int currentFPprobability,oldFPprobability=0;

    float[] equal;
    float[] diff;
    float[] maxEqual;
    float[] totalMaximum;
    int[] slidingMaximumIndex;
    float[] slidingMaximum;
    int[] maxIndex;
    int recognized=-1;

    int[][][] originalFP=null;
    int[][] needleFP=null;
    int[] originalDuration;
    float[] poleK;
    int limitNeedle=30;
    Rect srcRect=null, destRect=null;
    double[][] signals;
    int needleIndex=0;

    public enum State {INITIALIZING, READY, RECORDING, ERROR, STOPPED};
    private AudioRecord     audioRecorder = null;
    private State          	state;

    private int                      bufferSize;
    byte[] bufferFinregprintB;
    long oldA1,oldB1,oldR1,maxA1,maxB1,maxR1,minA1,minB1,minR1;
    public static double period =  0.0625;//0.25 0.125;//0.0625;//0.03125;//0.0625;//0.125; // N раз в сек, СТЕПЕНЬ ДВОЙКИ
    private int fftSampleSize=4096*4;//1024;      // number of sample in fft, the value needed to be a number to power of 2
    private double[][] spectrogram; // relative spectrogram
    private double[][] absoluteSpectrogram; // absolute spectrogram
    private int fftLength;   // number of y-axis unit
    private double unitFrequency;   // frequency per y-axis unit

    private int                      framePeriod;
    private byte[]                   buffer;
    final static int sampleRate = 8000;//22050;//;

	public static ExtAudioRecorder getInstanse() {
		ExtAudioRecorder result = null;
        int i=0;
        result = new ExtAudioRecorder(	true,
                                        AudioSource.MIC,
                                        sampleRate,
                                        AudioFormat.CHANNEL_IN_MONO,
                                        //AudioFormat.CHANNEL_CONFIGURATION_MONO,
                                        AudioFormat.ENCODING_PCM_16BIT);
		return result;
	}

    public ExtAudioRecorder(boolean uncompressed, int audioSource, int sampleRate, int channelConfig, int audioFormat){
        try{
            buffer = new byte[(int) Math.round(4096*4*period)];//framePeriod*bSamples/8*nChannels];
            fftSampleSize=buffer.length/2;
            framePeriod =(int) Math.round(4096*2*period);//bufferSize / ( 2 * bSamples * nChannels / 8 );

            bufferSize = AudioRecord.getMinBufferSize(sampleRate, channelConfig, audioFormat);
            //Log.d(TAG,"MinBufferSize("+sampleRate+","+channelConfig+","+audioFormat+")="+bufferSize); //4096  LG  , 640 daxian
            if (bufferSize == AudioRecord.ERROR || bufferSize == AudioRecord.ERROR_BAD_VALUE || bufferSize == AudioRecord.ERROR_INVALID_OPERATION)
                throw new Exception("AudioRecord initialization failed! bad buffer size "+bufferSize);
            if (bufferSize<4096) bufferSize=4096;
            // Set frame period and timer interval accordingly
            //Log.d(TAG,"framePeriod="+framePeriod);
            audioRecorder = new AudioRecord(audioSource, sampleRate, channelConfig, audioFormat, bufferSize*10);
            if (audioRecorder.getState() != AudioRecord.STATE_INITIALIZED)
                throw new Exception("AudioRecord initialization failed! state="+audioRecorder.getState());
            audioRecorder.setRecordPositionUpdateListener(updateListener);
            audioRecorder.setPositionNotificationPeriod(framePeriod);
            state = State.INITIALIZING;
            state = State.READY;
        } catch (Exception e){
            if (e.getMessage() != null)
                Log.e(ExtAudioRecorder.class.getName(), e.getMessage());
            else
                Log.e(ExtAudioRecorder.class.getName(), "Unknown error occured while initializing recording");
            state = State.ERROR;
        }
    }

	public State getState() {
		return state;
	}

    private AudioRecord.OnRecordPositionUpdateListener updateListener = new AudioRecord.OnRecordPositionUpdateListener()
	{
        @Override
        public void onMarkerReached(AudioRecord audioRecord) {
            //To change body of implemented methods use File | Settings | File Templates.
        }

        public void onPeriodicNotification(AudioRecord recorder) {
            ///*
            poleBM1=poleBM.copy(Bitmap.Config.RGB_565, false);
            if (srcRect==null) {
                srcRect = new Rect(0, 0, poleBM.getWidth(), poleBM.getHeight());
                destRect = new Rect(srcRect);
                destRect.offset(numFrames, 0);
            }
            poleC.drawBitmap(poleBM1, srcRect, destRect, null);
            poleP.setColor(Color.BLACK);
            for (int i=0;i<numFrames;i++)
                poleC.drawLine(i, 0, i, hh, poleP);
                // */

			audioRecorder.read(buffer, 0, buffer.length); // Fill buffer

            // get spectrogram's data
            //short[] amplitudes=buffer;
            final int numSamples = buffer.length;
            numFrames=numSamples/(fftSampleSize*2);  // WTF? numFrames:=1!! todo!
            //Log.d(TAG, "numSamples="+numSamples+", fftSampleSize="+fftSampleSize+", numFrames= "+numFrames);       //2
            if (numFrames>0) {

                WindowFunction window = new WindowFunction();
                window.setWindowType("Hamming");//Rectangular");
                double[] win=window.generate(fftSampleSize);

                signals = new double[numFrames][];
                int x = 0, y = 0, x1 = 0, y1 = 0;
                poleP.setColor(Color.GREEN);
                for (int f = 0; f < numFrames; f++) {
                    averAmp = 0.0d;
                    signals[f] = new double[fftSampleSize];
                    int startSample = f * fftSampleSize*2;
                    //Log.d(TAG,"startSample="+startSample+". startSample+n*2+1(n=fftSampleSize-1)="+(startSample+(fftSampleSize-1)*2+1));
                    //Log.d(TAG,"buffer.length="+buffer.length+", startSample+n*2="+(startSample+fftSampleSize*2));
                    for (int n = 0; n < fftSampleSize; n++) {
                        signals[f][n] = getShort(buffer[startSample+n*2], buffer[startSample+n*2+1])*win[n];
                        //signals[f][n] = buffer[startSample+n*2+1];
                        //averAmp += Math.abs(Math.round(signals[f][n]));
                        averAmp += Math.abs(signals[f][n]);
                    }
                    //рисуем внизу амплитуду сигнала (скользящее)
                    averAmp = Math.min(100,0.1*averAmp/ fftSampleSize);
                    poleC.drawLine(numFrames - f - 1, (float) (hh-averAmp), numFrames - f , (float) (hh - lastAverageAmp), poleP);
                    lastAverageAmp= (lastAverageAmp+averAmp)/2;
                }
                //Log.d(TAG,"lastAverageAmp="+lastAverageAmp);

                absoluteSpectrogram = new double[numFrames][];
                for (int f = 0; f < numFrames; f++) {
                    FFT fft = new FFT(signals[f].length / 2, -1);
                    fft.transform(signals[f]);
                    absoluteSpectrogram[f]=new double[signals[f].length/2];
                    final int len=absoluteSpectrogram[f].length;
                    for(int i = 0; i < len; i++) {
                        absoluteSpectrogram[f][len-1-i] = Math.sqrt(signals[f][i*2] * signals[f][i*2] + signals[f][i*2 + 1] * signals[f][i*2 + 1]);
                    }

                }
                fftLength = absoluteSpectrogram[0].length; //256 для 8000
                barkToFrequnits=300.0f/(float) fftLength;
                unitFrequency = (double) sampleRate / 2 / fftLength; //15 для 8000
                lowFreq= (int) (250/(unitFrequency));//начнем анализ с 250Гц = 16 для 8000
                highFreq= (int) (3500/(unitFrequency));//закончим на 3500Гц =224 для 8000
                //Log.d(TAG,"fftLength="+fftLength+", unitFrequency="+unitFrequency+"MHz, lowFreq="+lowFreq+", highFreq="+highFreq+",signals.len="+signals[0].length);
                spectrogram = new double[numFrames][fftLength];
                barkArray = new double[spectrogram.length][17][2];
                for (int f = 0; f < numFrames; f++) {
                    maxIntens = Double.MIN_VALUE;
                    minIntens = Double.MAX_VALUE;
                    averIntens =0;
                    for (int j = lowFreq;j < highFreq; j++) {
                        averIntens+=absoluteSpectrogram[f][j];
                        if (absoluteSpectrogram[f][j] > maxIntens)
                            maxIntens = absoluteSpectrogram[f][j];
                        else if (absoluteSpectrogram[f][j] < minIntens)
                            minIntens = absoluteSpectrogram[f][j];
                    }
                    averIntens=averIntens/ fftLength;
                    // avoiding divided by zero
                    if (minIntens == 0)
                        minIntens = minValidAmp;
                    perceptualDiff = Math.log10(maxIntens / minIntens);  // perceptual difference

                    //количество ярких частот. Чем больше, тем шумнее.
                    //sumFullIntensAmp=1.0;
                    for (int j = lowFreq; j < highFreq; j++) {
                        if (absoluteSpectrogram[f][j] < minValidAmp || maxIntens<=0)
                            spectrogram[f][j] = 0;
                        else {
                            //фильтр 1 - максимайзер
                            spectrogram[f][j] = (Math.log10(absoluteSpectrogram[f][j] / minIntens)) / perceptualDiff;    // 0..1
                            //варианты максимайзера
                            //spectrogram[f][j] = (absoluteSpectrogram[f][j]-minIntens)/maxIntens;
                            //spectrogram[f][j] = absoluteSpectrogram[f][j]/maxIntens;

                            //фильтр 2 - гейт
                            if (spectrogram[f][j]<0.95)
                                spectrogram[f][j]=0;
                            //sumFullIntensAmp += spectrogram[f][j];//(spectrogram[f][j]>0.8?1:0);
                        }
                    }
                    //Log.d(TAG,"sumFullIntensAmp="+sumFullIntensAmp);
                    //sumFullIntensAmp=Math.pow(1.0/(sumFullIntensAmp-1),2);//1.....0, 1 - ok, 0 - шумно


                    //BARK
                    //for (int ii = 0; ii < 17; ii++)
                    //    for (int jj = 0; jj < 2; jj++)
                    //        barkArray[f][ii][jj]=0;
                    for (int j = lowFreq; j < highFreq; j++) {
                        if (Bark[(int) (j*barkToFrequnits)] >= 0 && Bark[(int) (j*barkToFrequnits)] <= 16) {
                            barkArray[f][Bark[(int) (j*barkToFrequnits)]][0] += spectrogram[f][j];
                            barkArray[f][Bark[(int) (j*barkToFrequnits)]][1]++;
                        }
                    }

                    ///*
                    for (int j = lowFreq; j < highFreq; j++) {
                        if (maxIntens>0) {
                            //рисуем наверху интенсивности сигнала
                            mycolor=(int) (Math.min(100.0, 0.2 * 100.0 * absoluteSpectrogram[f][j] / averIntens));
                            poleP.setColor(Color.argb(255, (int) (Math.min(255.0, 2.0 * 255.0 * absoluteSpectrogram[f][j] / maxIntens)), mycolor, mycolor));
                            poleC.drawPoint(numFrames - f - 1, j - lowFreq, poleP);
                            //рисуем наверху-2 нормализованные интенсивности сигнала
                            //mycolor=spectrogram[f][j]>0.8? (int) (255 * spectrogram[f][j]) :0;
                            mycolor=255;//(int) (255 * spectrogram[f][j]);
                            poleP.setColor(Color.argb(255, mycolor, mycolor, mycolor));
                            if (spectrogram[f][j]>0) {
                                poleC.drawPoint(numFrames - f - 1, j - lowFreq, poleP);
                            }
                            //рисуем внизу количество ярких частот (шум)
                            //poleP.setColor(Color.WHITE);
                            //poleC.drawLine(numFrames - f - 1, hh, numFrames - f - 1, hh - Math.round(sumFullIntensAmp * 30.0), poleP);
                        }
                    }// */
                    //Log.d(TAG, "fftLength="+fftLength+", minIntens= "+ minIntens +", maxIntens="+ maxIntens+", averIntens="+averIntens);       //7..10000
                    //Log.d(TAG, "y1="+y1+", minIntens= "+ minIntens +", maxIntens="+ maxIntens+", averIntens="+averIntens);       //7..10000
                    for (int j = 0; j < 17; j++) {
                        if (barkArray[f][j][1] > 0)
                            barkArray[f][j][0] = barkArray[f][j][0] / barkArray[f][j][1];
                        //рисуем наверху-3 барк-шкалу
                        //mycolor = (int) (Math.min(255,barkArray[f][j][0] * 255*2));
                        //poleP.setColor(Color.argb(255, mycolor, mycolor, (int) (mycolor*0.7f)));
                        //poleC.drawLine(numFrames - f - 1, highFreq-lowFreq+j*5, numFrames - f - 1, highFreq-lowFreq+(j+1)*5,poleP);
                    }


                    //Алгоритм №1
                    //LS method
                    //подобие среднего вектора ax+b, r - отклонение, разброс, добротность
                    double X_ = 0, Y_ = 0, YX_ = 0, X2_ = 0, Y2_ = 0, zz = 0, zz1 = 0, zr = 0, b, a, r;
                    int a1;
                    int b1;
                    int r1;
                    for (int j = 2; j <= 16; j++) {
                        zz = barkArray[f][j][0];
                        X_ = X_ + (j - 1);
                        Y_ = Y_ + zz;
                        YX_ = YX_ + ((j - 1) * zz);
                        X2_ = X2_ + ((j - 1) * (j - 1));
                        Y2_ = Y2_ + (zz * zz);
                    }
                    zz1 = (15 * X2_) - (X_ * X_);
                    zr = (15 * YX_) - (Y_ * X_);
                    if (zz1 != 0) {
                        b = ((Y_ * X2_) - (YX_ * X_)) / zz1;
                        a = zr / zz1;
                    } else {
                        b = 0;
                        a = 0;
                    }
                    zz1 = Math.sqrt(zz1 * ((15 * Y2_) - (Y_ * Y_)));
                    if (zz1 != 0) {
                        r = zr / zz1;
                    } else {
                        r = 0;
                    }
                    //legacy
                    a1 = (int) (220 + Math.round(a * 2000));
                    b1 = (int) (50 + Math.round(b * 200));
                    r1 = (int) (100 + Math.round(r * 100));
                    /*
                    if (maxA1<a1) {maxA1=a1; Log.d(TAG,"maxA1="+maxA1+",maxB1="+maxB1+",maxR1="+maxR1+",minA1="+minA1+",minB1="+minB1+",minR1="+minR1+", counter="+counter);}
                    if (maxB1<b1) {maxB1=b1; Log.d(TAG,"maxA1="+maxA1+",maxB1="+maxB1+",maxR1="+maxR1+",minA1="+minA1+",minB1="+minB1+",minR1="+minR1+", counter="+counter);}
                    if (maxR1<r1) {maxR1=r1; Log.d(TAG,"maxA1="+maxA1+",maxB1="+maxB1+",maxR1="+maxR1+",minA1="+minA1+",minB1="+minB1+",minR1="+minR1+", counter="+counter);}
                    if (minA1>a1) {minA1=a1; Log.d(TAG,"maxA1="+maxA1+",maxB1="+maxB1+",maxR1="+maxR1+",minA1="+minA1+",minB1="+minB1+",minR1="+minR1+", counter="+counter);}
                    if (minB1>b1) {minB1=b1; Log.d(TAG,"maxA1="+maxA1+",maxB1="+maxB1+",maxR1="+maxR1+",minA1="+minA1+",minB1="+minB1+",minR1="+minR1+", counter="+counter);}
                    if (minR1>r1) {minR1=r1; Log.d(TAG,"maxA1="+maxA1+",maxB1="+maxB1+",maxR1="+maxR1+",minA1="+minA1+",minB1="+minB1+",minR1="+minR1+", counter="+counter);}
                    // */
                    //if (a1>255 || a1<0 || b1<0 || b1>255 || r1<0 || r1>255)
                    //    Log.d(TAG,"borders! a1="+a1+", b1="+b1+", r1="+r1);
                    //Log.d(TAG,"хххххх: a1="+a1+", b1="+b1+", r1="+r1);
                    //a1 = (a1 & 0xFF);
                    //b1 = (b1 & 0xFF);
                    //r1 = (r1 & 0xFF);
                    if (a1!=b1 && Math.abs(a1+b1)<minA1plusB1) {
                        //автоподстройка
                        //Log.d(TAG,"Math.abs(a1+b1)="+Math.abs(a1+b1));
                        minA1plusB1 = Math.abs(a1+b1);
                    }

                    //Денойзер атонального шума. Log.d(TAG,"Критерии: ["+(Math.abs(a1+b1)>minA1plusB1)+"]-["+(a1<200)+"]-["+(b1>50)+"]");
                    //currentFPprobability = (int) Math.max(0,255-(Math.abs(a1+b1)>minA1plusB1?50:0)-(a1<200?100:0)-(b1>50?100:0)); //djd мои патентованные оценки
                    //currentFPprobability=(oldFPprobability+currentFPprobability)/2;
                    //рисуем
                    mycolor=255;//(oldFPprobability>100 && currentFPprobability>100?255:0);
                    //Log.d(TAG,"currentFPprobability="+currentFPprobability+",mycolor="+mycolor);
                    //poleP.setColor(Color.argb(mycolor,50,200,50));
                    //poleC.drawLine(numFrames - f, 100 + 2 * highFreq - 2 * lowFreq +oldA1,numFrames - f - 1, 100 + 2 * highFreq - 2 * lowFreq +a1, poleP);

                    //рисуем отмпечаток ---------------------
                    poleP.setColor(Color.argb(mycolor,100,255,100));
                    poleC.drawLine(numFrames - f    , highFreq-lowFreq+oldB1,
                                   numFrames - f - 1, highFreq-lowFreq+b1, poleP);

                    //poleC.drawLine(numFrames - f, 150 + 2 * highFreq - 2 * lowFreq +Math.abs(oldA1+oldB1),numFrames - f, 150 + 2 * highFreq - 2 * lowFreq +Math.abs(a1+b1), poleP);
                    //x=(int) 20+(a1/3);
                    //x1=(int) x+(b1/3);
                    /*
                    mycolor = (int) (60 + r1);
                    poleP.setColor(Color.argb(255,mycolor,mycolor,mycolor));
                    poleC.drawLine(numFrames - f - 1, 20 + 2 * highFreq - 2 * lowFreq, numFrames - f - 1, 20 + 2 * highFreq - 2 * lowFreq + (b1 / 2), poleP);
                    poleP.setColor(Color.argb(255,mycolor, (int) (mycolor*0.7),mycolor));
                    poleC.drawLine(numFrames - f - 1, 20 + 2 * highFreq - 2 * lowFreq + (b1 / 2), numFrames - f - 1, 20 + 2 * highFreq - 2 * lowFreq + (b1 / 2) + (a1 / 3), poleP);

                    int xx = (int) (Math.log10(1 + Math.abs(oldB1 - b1) * 10) * 20);
                    int xxx = (int) (Math.log10(1 + Math.abs(oldA1 - a1)) * 20);
                    mycolor = (int) (Math.abs(60 + oldR1 - r1) * 255);
                    poleP.setColor(Color.argb(255,mycolor,0,0));
                    poleC.drawLine(numFrames - f - 1, hh - 200-xx, numFrames - f - 1, hh - 200-xx - xxx, poleP);
                    */

                    //poleP.setColor(Color.argb(255,128,128,0));
                    //poleC.drawLine(numFrames - f - 1, hh - 100, numFrames - f - 1, hh - 100 - (oldR1+r1)/2, poleP);

                    oldA1 = a1;
                    oldB1 = b1;
                    oldR1 = r1;
                    oldFPprobability=currentFPprobability;

                    if (db!=null && db.isOpen()) {
                        if (recordMode>0) {
                            db.execSQL("insert into original"+recordMode+" (a,b,r) values (" + a1 + "," + b1 + "," + r1 + ")");
                            recordedFPs++;
                        } else {
                            //bd - тормоз!
                            //db.execSQL("insert into needle (a,b,r) values (" + a1 + "," + b1 + "," + r1 + ")");
                            if (needleFP==null) {
                                needleFP = new int[limitNeedle][3]; //a,b,r + сдвиг и его оценка
                                needleIndex=-1;
                            }
                            needleIndex=(needleIndex+1)%limitNeedle;
                            needleFP[needleIndex][0]=a1;
                            needleFP[needleIndex][1]=b1;
                            needleFP[needleIndex][2]=r1;
                            recognizeMe(true);
                        }
                    }


                    //Алгоритм №2
                    //что-то вроде побитового сравнения сдвиговых разниц. Имхо, плохой
                    /*
                    if (previousBarkArray==null){
                        previousBarkArray=new double[barkArray[f].length];
                    }
                    int rr=0;
                    for (int j = 2; j <  16; j++) {
                        rr= (int) (255*2*(barkArray[f][j][0]-barkArray[f][j+1][0]-previousBarkArray[j]+previousBarkArray[j+1]));
                        //Log.d(TAG,"rr["+j+"]="+rr);
                        mycolor= rr>0?255:0;//(int) Math.max(0,Math.min(255,rr));
                        poleP.setColor(Color.argb(255,mycolor,mycolor,0));
                        poleC.drawLine(numFrames - f-1, 3 * highFreq - 3 * lowFreq +j*3,numFrames - f - 1, 3 * highFreq - 3 * lowFreq +(j+1)*3, poleP);

                    }
                    for (int j = 2; j <= 16; j++) {
                        previousBarkArray[j] = barkArray[f][j][0];
                    }

                    //Алгоритм №3
                    // Split the FFT analysis into 4 filter bands. Take find the max value in each band and store the time
                    // The Hash value is computed from the max values.
                    int point1 = 0;
                    int point2 = 0;
                    int point3 = 0;
                    int point4 = 0;
                    int freq;
                    for (int j = lowFreq;j < highFreq; j++) {
                        freq = (int) Math.abs(absoluteSpectrogram[f][j]);
                        if (j < fftLength*0.25) {
                            if (Math.abs(absoluteSpectrogram[f][point1]) < freq || point1 == 0) {
                                point1 = j;
                            }
                        } else if (j < fftLength*0.5) {
                            if (Math.abs(absoluteSpectrogram[f][point2]) < freq || point2 == 0) {
                                point2 = j;
                            }
                        } else if (j < fftLength*0.75) {
                            if (Math.abs(absoluteSpectrogram[f][point3]) < freq || point3 == 0) {
                                point3 = j;
                            }
                        } else {
                            if (Math.abs(absoluteSpectrogram[f][point4]) < freq || point4 == 0) {
                                point4 = j;
                            }
                        }
                    }
                    int[] points = new int[4];
                    points[0] = point1;
                    points[1] = point2;
                    points[2] = point3;
                    points[3] = point4;
                    //fingerprint = new ArrayList<int[]>();
                    poleP.setColor(Color.argb(255,255,255,50));
                    for (int j = lowFreq;j < highFreq; j++) {
                        if (Math.abs(point1-j)<5 || Math.abs(point2-j)<5 || Math.abs(point3-j)<5 || Math.abs(point4-j)<5) {
                            poleC.drawPoint(numFrames - f - 1, 3.5f * highFreq - 4 * lowFreq + j / 5, poleP);
                        }

                    }
                    // */

                }
            }
            //pole.setImageBitmap(poleBM);
            pole.invalidate();
        }

    };

    void recognizeMe(boolean andDraw){
        int i;
        if (originalFP==null) {
            //creation
            originalFP = new int[3][][];
            poleK = new float[3];
            originalDuration = new int[3];
            equal=new float[3];
            diff=new float[3];
            maxEqual=new float[3];
            totalMaximum=new float[3];
            slidingMaximumIndex=new int[3];
            slidingMaximum=new float[3];
            maxIndex=new int[3];



            for (int t=1;t<=3;t++) {
                String q1 = "select * from original"+t+" order by id asc";
                Cursor c2 = db.rawQuery(q1, null);
                originalDuration[t-1] = c2.getCount() / (sampleRate / (fftSampleSize));
                if (c2.getCount() > 0 && c2 != null && c2.moveToFirst()) {
                    originalFP[t-1] = new int[c2.getCount()][6]; //a,b,r, сумма1,сумма2,average2
                    i = 0;
                    do {
                        originalFP[t-1][i][0] = c2.getInt(c2.getColumnIndex("a"));
                        originalFP[t-1][i][1] = c2.getInt(c2.getColumnIndex("b"));
                        originalFP[t-1][i][2] = c2.getInt(c2.getColumnIndex("r"));
                        originalFP[t-1][i][3] = 0;
                        originalFP[t-1][i][4] = 0;
                        originalFP[t-1][i][5] = 0;
                        i++;
                    } while (c2.moveToNext());
                    poleK[t-1] = (float) poleBMn[t-1].getWidth() / (float) (originalFP[t-1].length - limitNeedle);
                    Log.d(TAG, "Выбрано по образцу сэмплов #"+t+": " + c2.getCount() + " (" + originalFP[t-1].length + ") " + ", originalDuration1=" + originalDuration[t-1] + ", k=" + poleK[t-1]);
                }
                c2.close();
            }
        }
        ///////////////////////////////////////////

        /*
        if (needleFP==null)
            needleFP=new int[limitNeedle][3]; //a,b,r + сдвиг и его оценка
        String q1 = "select * from needle order by id desc limit "+limitNeedle;
        Cursor c2 = db.rawQuery(q1, null);
        if (c2.getCount() > 0 && c2 != null && c2.moveToLast()) {
            //needleFP=new int[c2.getCount()][3]; //a,b,r + сдвиг и его оценка
            i=0;
            do {
                needleFP[i][0] = c2.getInt(c2.getColumnIndex("a"));
                needleFP[i][1] = c2.getInt(c2.getColumnIndex("b"));
                needleFP[i][2] = c2.getInt(c2.getColumnIndex("r"));
                i++;
            } while (c2.moveToPrevious());
        }
        c2.close();
        */

        if (andDraw) {
            poleP.setColor(Color.argb(255, 255, 255, 255));
            for (int t=1;t<=3;t++) {
                poleCN[t - 1].drawRect(new Rect(0, 0, poleN[t - 1].getWidth(), 200), poleP);
            }
        }

        int jj,bonus;
        int wasRecognized=-1;
        for (int t=1;t<=3;t++) {
            if (originalFP[t-1] != null && needleFP != null) {
                //float oldMaxEqual=maxEqual[t-1];
                maxEqual[t-1] = 0;
                maxIndex[t-1] = -1;
                float maximum2=0.0f;
                int maximumI1=-1,maximumI2=-1;
                //poleP1.setColor(Color.argb(255,0,0,0));
                for (int shift = -limitNeedle/2; shift < originalFP[t-1].length + limitNeedle/2; shift++) {
                    equal[t - 1] = 0.0f;
                    diff[t - 1] = 0.0f;
                    for (int j = 0; j < needleFP.length; j++) {
                        i = j + shift;
                        jj = (needleIndex + 1 + j) % limitNeedle;
                        if (i >= 0 && i < originalFP[t - 1].length) {
                            /*
                            if (originalFP[t - 1][i][0] == needleFP[jj][0] && originalFP[t - 1][i][1] == needleFP[jj][1] && originalFP[t - 1][i][2] == needleFP[jj][2])
                                equal[t - 1]++;
                            else
                                diff[t - 1]++;// */
                            ///*
                            if (originalFP[t - 1][i][0] == needleFP[jj][0])
                                equal[t - 1]++;
                            else
                                diff[t - 1]++;
                            if (originalFP[t - 1][i][1] == needleFP[jj][1])
                                equal[t - 1]++;
                            else
                                diff[t - 1]++;
                            if (originalFP[t - 1][i][2] == needleFP[jj][2])
                                equal[t - 1]++;
                            else
                                diff[t - 1]++;
                            // */
                        }
                    }
                    //equal[t-1]=equal[t-1]/(diff[t-1]+1); //не надо , делает хуже
                    //if (diff[t-1]-equal[t-1]<10)...//не надо , делает хуже
                    if (equal[t - 1] > maxEqual[t - 1]) {
                        maxEqual[t - 1] = equal[t - 1];
                        maxIndex[t - 1] = shift;
                    }
                    if (andDraw && totalMaximum[t - 1] > 0) {
                        //Log.d(TAG,"diff[t-1]-equal[t-1]="+(diff[t-1]-equal[t-1]));
                        //poleCN[t-1].drawLine(poleK[t-1]*shift, 190, poleK[t-1]*shift, 190.0f-190.0f*((diff[t-1]-equal[t-1])/limitNeedle), poleP1);
                        //poleCN[t-1].drawLine(poleK[t-1]*shift, 190, poleK[t-1]*shift, 190.0f-(diff[t-1]-equal[t-1]), poleP1);
                        if (maxEqual[t - 1] - equal[t - 1] < 10)
                            poleP.setColor(Color.argb(255, 255, 0, 0));
                        else
                            poleP.setColor(Color.argb(255, 100, 100, 100));
                        /*poleCN[t-1].drawLine(poleK[t-1] * shift, 190, poleK[t-1] * shift
                                , 190 - 190.0f * (equal[t-1] / totalMaximum[t-1])
                                //, 190 - 190.0f * (equal[t-1] / totalMaximum[t-1]) * (20.0f / (20.0f + Math.abs(shift - slidingMaximumIndex[t-1])))
                                , poleP);// */
                    }
                    ///*
                    //вариант 2
                    if (shift >= 0 && shift < originalFP[t - 1].length && totalMaximum[t-1]>0) {
                        originalFP[t - 1][shift][4]+=100*equal[t-1]/totalMaximum[t-1];
                        //if (originalFP[t - 1][shift][4]>maximum) {
                        //    maximum = originalFP[t - 1][shift][4];
                        //    maximumI=shift;
                        //}
                    }// */
                    //////вариант 2
                }
                //вариант 2
                poleP.setColor(Color.argb(255, 255, 100, 100));
                ///*
                int area = limitNeedle/2;
                float average;
                int counter;
                for (int k = originalFP[t-1].length - 1; k > 0; k--) {
                    //originalFP[t-1][k][4] = (int) (originalFP[t-1][k-1 ][4]*0.98f);
                    //originalFP[t-1][k][4] = (int) (0.5f*(originalFP[t-1][k][4]+originalFP[t-1][k-1 ][4])/2);
                    originalFP[t - 1][k][4] = (int) (originalFP[t - 1][k - 1][4] * 0.995f);
                    if (originalFP[t - 1][k][4] < 0) {
                        originalFP[t - 1][k][4] = 0;
                    }
                    average = 0;
                    counter = 0;
                    for (int kk = k - area; kk < k + area; kk++)
                        if (kk >= 0 && kk < originalFP[t - 1].length) {
                            average += originalFP[t - 1][kk][4];
                            counter++;
                        }
                    if (counter > 0) {
                        average = average / (float)counter;
                        //averageTotal+=average;
                    }
                    originalFP[t - 1][k][5] = (int) (((float)originalFP[t - 1][k][4]-average)/(float)area);
                }
                //сбор максимумов, некоторые рядом
                int over100count=0,over150count=0;
                int over150maximumI=-2*originalFP[t-1].length;
                for (int k = 1; k<originalFP[t-1].length-1; k++) {
                    if (originalFP[t - 1][k][5] >= 100){
                        if (originalFP[t - 1][k][5] >= 150) {
                            over150count++;
                            over150maximumI=k;
                        } else if (Math.abs(k-over150maximumI)<3) {
                            over150count++;
                        } else
                            over100count++;
                    }
                }
                originalFP[t-1][0][4] = 0;
                for (int k = originalFP[t-1].length - 1; k > 0; k--) {
                    if (andDraw) {

                        if (over150count>0 && originalFP[t - 1][k][5]>=100) {
                            if (originalFP[t - 1][k][5] >= 150)
                                poleP.setColor(Color.argb(255, 255, 0, 0));
                            else
                                poleP.setColor(Color.argb(255, 100, 100, 100));
                            poleCN[t - 1].drawLine(poleK[t - 1] * k, 190, poleK[t - 1] * k, 190 - 0.5f * originalFP[t - 1][k][5], poleP);
                        }
                    }
                }

                if (over150count>0)
                    Log.d(TAG,"Распознан! #"+t+", over150count="+over150count+", over100count="+over100count);
                if (andDraw && over150count>0) {
                    if (over100count>1)
                        poleP.setColor(Color.argb(255,100,100,100));
                    else
                        poleP.setColor(Color.argb(255,0,0,255));
                    poleCN[t-1].drawRect(new Rect((int)(poleK[t-1] * over150maximumI-10)
                            , 0
                            , (int)(poleK[t-1] * over150maximumI+10)
                            , 50), poleP);
                }
                //////вариант 2

                if (maxEqual[t-1] > totalMaximum[t-1]) {
                    totalMaximum[t-1] = maxEqual[t-1];
                    //Log.d(TAG,"totalMaximum1="+totalMaximum1);
                }
                //int kk=originalFP1[0][3];
                final int oldSlidingMaximumIndex= slidingMaximumIndex[t-1];
                slidingMaximumIndex[t-1] = 0;
                slidingMaximum[t-1] = 0;
                for (int k = originalFP[t-1].length - 1; k > 0; k--) {
                    originalFP[t-1][k][3] =originalFP[t-1][k - 1][3] - 1;
                    if (originalFP[t-1][k][3] < 0)
                        originalFP[t-1][k][3] = 0;
                    else if (originalFP[t-1][k][3] > slidingMaximum[t-1]) {
                        slidingMaximum[t-1] = originalFP[t-1][k][3];
                        slidingMaximumIndex[t-1] = k;
                    }
                }
                originalFP[t-1][0][3] = 0;

                if (maxIndex[t-1] > 1 && maxIndex[t-1]<originalFP[t-1].length-2) {
                    //originalFP[t-1][maxIndex[t-1]-2][3] += 1;
                    //originalFP[t-1][maxIndex[t-1]-1][3] += 2;
                    originalFP[t-1][maxIndex[t-1]][3] += 3;
                    //originalFP[t-1][maxIndex[t-1]+1][3] += 2;
                    //originalFP[t-1][maxIndex[t-1]+2][3] += 1;
                    if (Math.abs(slidingMaximumIndex[t-1]-maxIndex[t-1])<1 && slidingMaximum[t-1]/limitNeedle>1) {
                        //всякие эмпирические нелинейности
                        //бонус за попадание в максимальную точку
                        bonus = 1;
                        //бонус за энергию
                        if (slidingMaximum[t-1]/limitNeedle>2)
                            bonus += 1;
                        //бонус за попадание в ту же старую точку
                        if (Math.abs(oldSlidingMaximumIndex - maxIndex[t - 1]) < 2)
                            bonus += 5;
                        originalFP[t-1][maxIndex[t-1]][3]+=bonus;
                        //Log.d(TAG, "bonus[" + (t - 1) + "]["+maxIndex[t-1]+"] = "+bonus);//slidingMaximumIndex[" + (t - 1) + "]="+slidingMaximumIndex[t - 1]+", =maxIndex="+maxIndex[t - 1] + ", oldSlidingMaximumIndex="+oldSlidingMaximumIndex+", slidingMaximum[" + (t - 1) + "]/limitNeedle=" + (slidingMaximum[t - 1] / limitNeedle));
                    }
                    //лимитер. Столько секунд будет откатываться назад
                    if (originalFP[t-1][maxIndex[t-1]][3]/limitNeedle>7) {
                        originalFP[t - 1][maxIndex[t - 1]][3] = 7 * limitNeedle;
                    }
                }
                if (slidingMaximum[t-1]/limitNeedle>3) {
                    wasRecognized=t-1;
                }
                if (false && andDraw) {
                    poleP.setColor(Color.argb(255, 255, 100, 0));
                    poleCN[t-1].drawLine(poleK[t-1] * maxIndex[t-1], 190, poleK[t-1] * maxIndex[t-1], 200, poleP);

                    poleP.setColor(Color.argb(255, 0, 157, 0));
                    for (int k = 0; k < originalFP[t-1].length; k++) {
                        if (originalFP[t-1][k][3] > 0) {
                            poleCN[t-1].drawLine(poleK[t-1] * k, 190, poleK[t-1] * k, 190 - 190 * originalFP[t-1][k][3] / limitNeedle, poleP);
                            poleCN[t-1].drawLine(poleK[t-1] * k, 190, poleK[t-1] * k, 200, poleP);
                            poleCN[t-1].drawLine(poleK[t-1] * k, 192, poleK[t-1] * k - 1, 200, poleP);
                            poleCN[t-1].drawLine(poleK[t-1] * k, 192, poleK[t-1] * k + 1, 200, poleP);
                            poleCN[t-1].drawLine(poleK[t-1] * k, 190, poleK[t-1] * k + 10, 200, poleP);
                            poleCN[t-1].drawLine(poleK[t-1] * k, 190, poleK[t-1] * k - 10, 200, poleP);
                        }
                    }
                    /*
                    poleP.setColor(Color.argb((int) Math.min(255,80* slidingMaximum[t-1] / limitNeedle), 0, 0, 255));
                    poleCN[t-1].drawRect(new Rect((int)(poleK[t-1] * slidingMaximumIndex[t-1]-10)
                            , (int)(200-70 * slidingMaximum[t-1] / limitNeedle)
                            , (int)(poleK[t-1] * slidingMaximumIndex[t-1]+10)
                            , 200), poleP);
                    // */
                }
            }
        }
        if (wasRecognized!=-1 && andDraw) {
            poleP.setTextSize(poleBMn[wasRecognized].getWidth()/10);
            poleP.setColor(Color.argb(255, 255, 255, 255));
            poleCN[wasRecognized].drawText("Распознан!",1,1+poleBMn[wasRecognized].getHeight()/2,poleP);
            poleP.setColor(Color.argb(255, 0, 0, 0));
            poleCN[wasRecognized].drawText("Распознан!",0,poleBMn[wasRecognized].getHeight()/2,poleP);
        }
    }

	public void start(Context context){
        //Log.d(TAG,"Recorder try starting... State="+state.toString());
        if (state == State.READY || state == State.INITIALIZING){

            //init Bark scale
            //соответствие 4000Гц к 300 к 17 полосам восприятия слуха
            //300 - asi is патентованная формула
            double z;
            Bark = new byte[300];
            for (int i=0;i<300;i++){
                z=i*15.625; //привеодим к 4000 гц
                z=Math.floor((26.81*z)/(1960+z))-0.53;
                if (z<2)
                    z=z+0.15*(2-z);
                else
                    if (z>20.1)
                        z=z+0.22*(z-20.1);
                Bark[i]=(byte) Math.floor(z);
                //Log.d(TAG,i+" - "+Bark[i]);
            }

            //and init graph
            pole = MyActivity.pole;
            ww=pole.getWidth(); if (ww==0) ww=200;
            hh=pole.getHeight(); if (hh==0) hh=256;
            poleBM = Bitmap.createBitmap(ww, hh, Bitmap.Config.RGB_565);
            poleBM1 = Bitmap.createBitmap(ww, hh, Bitmap.Config.RGB_565);
            poleC = new Canvas(poleBM);
            poleP = new Paint();
            poleP1 = new Paint();
            //poleC.drawBitmap(BitmapFactory.decodeResource(context.getResources(), R.drawable.logo), 0,0, null);
            pole.setImageBitmap(poleBM);

            poleN = new ImageView[3];
            poleBMn = new Bitmap[3];
            poleCN= new Canvas[3];

            poleN[0] = MyActivity.poleN1;
            poleBMn[0] = Bitmap.createBitmap(poleN[0].getWidth(), 200, Bitmap.Config.RGB_565);
            poleCN[0] = new Canvas(poleBMn[0]);
            poleN[0].setImageBitmap(poleBMn[0]);

            poleN[1] = MyActivity.poleN2;
            poleBMn[1] = Bitmap.createBitmap(poleN[1].getWidth(), 200, Bitmap.Config.RGB_565);
            poleCN[1] = new Canvas(poleBMn[1]);
            poleN[1].setImageBitmap(poleBMn[1]);

            poleN[2] = MyActivity.poleN3;
            poleBMn[2] = Bitmap.createBitmap(poleN[2].getWidth(), 200, Bitmap.Config.RGB_565);
            poleCN[2] = new Canvas(poleBMn[2]);
            poleN[2].setImageBitmap(poleBMn[2]);

            bufferFinregprintB=new byte[500];


            //open DB
            String dir = Environment.getExternalStorageDirectory().getPath();
            File dbfile = new File(dir+"/echolotFingerprints.db");
            db = SQLiteDatabase.openOrCreateDatabase(dbfile, null);
            if (!db.isOpen()){
                Log.d(TAG,"Its open? "  + db.isOpen());
                state = State.ERROR;
                return;
            }
            recordedFPs=0;
            if (recordMode>0) {
                try {
                    db.execSQL("drop table original"+recordMode);
                    Log.d(TAG, "удалена база данных original"+recordMode);
                } catch (SQLiteException exception) {
                    //Log.d(TAG,"ошибка удаления базы original");
                    ;//exception.printStackTrace();
                }
            } else {
                try {
                    db.execSQL("drop table needle");
                    //Log.d(TAG, "удалена база данных needle");
                } catch (SQLiteException exception) {
                    //Log.d(TAG,"ошибка удаления базы needle");
                    ;//exception.printStackTrace();
                }
            }
            if (recordMode>0)
                try {
                    db.execSQL("create table original"+recordMode+" ("
                            + "id integer primary key,"
                            + "a integer,"
                            + "b integer,"
                            + "r integer "
                            + ")");
                    //db.execSQL("CREATE INDEX original"+recordMode+"_a_idx ON original"+recordMode+" (a)");
                } catch (SQLiteException exception) {
                    if (recordMode==1) Log.e(TAG,"ошибка создания базы original"+recordMode+", "+exception.toString());
                    //exception.printStackTrace();
                }
            try {
                db.execSQL("create table needle ("
                        + "id integer primary key,"
                        + "a integer,"
                        + "b integer,"
                        + "r integer "
                        + ")");
                //db.execSQL("CREATE INDEX needle_a_idx ON needle (a)");
                //if (recordMode) Log.d(TAG, "создана новая база данных needle");
            } catch (SQLiteException exception) {
                if (recordMode==0) Log.e(TAG,"ошибка создания базы needle, "+exception.toString());
                //exception.printStackTrace();
            }


            //finally start rec
			audioRecorder.startRecording();
			audioRecorder.read(buffer, 0, buffer.length);
			state = State.RECORDING;

            //Log.d(TAG,"start() goes to record.");
		} else {
			Log.e(ExtAudioRecorder.class.getName(), "start() called on illegal state");
			state = State.ERROR;
		}
	}
	
	/**
	 * 
	 * 
	 *  Stops the recording, and sets the state to STOPPED.
	 * In case of further usage, a reset is needed.
	 * Also finalizes the wave file in case of uncompressed recording.
	 * 
	 */
	public void stop()
	{
		if (state == State.RECORDING){
			audioRecorder.stop();
			state = State.STOPPED;
            if (db!=null) {
                db.close();
                db=null;
            }
		} else {
			Log.e("djd", "stop() called on illegal state");
			state = State.ERROR;
		}
	}
	
	private short getShort(byte argB1, byte argB2) {
		return (short)(argB1 | (argB2 << 8));
	}
	
}
