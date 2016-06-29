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
    public boolean liteAlgo=true, andDraw=false;
    byte[] Bark;
    int ww,hh;
    int numFrames=1;
    double averAmp, lastAverageAmp =0;
    double maxIntens, minIntens, averIntens;
    SQLiteDatabase db;
    int mycolor;
    double perceptualDiff;
    double minValidAmp = 0.00000000001;
    double[][][] barkArray;
    int lowFreq,highFreq;
    float barkToFrequnits;

    float[] equal;
    float[] maxEqual;
    float[] totalMaximum;
    int[] slidingMaximumIndex;
    float[] slidingMaximum;
    int[] maxIndex;

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
    long oldA1,oldB1,oldR1;//,maxA1,maxB1,maxR1,minA1,minB1,minR1;
    long oldA12,oldB12,oldR12;
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
        }

        public void onPeriodicNotification(AudioRecord recorder) {
            //рисуем яркости частот, отпечаток, оцсилограмму
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

			audioRecorder.read(buffer, 0, buffer.length); // Fill buffer

            // get spectrogram's data
            final int numSamples = buffer.length;
            numFrames=numSamples/(fftSampleSize*2);  // WTF? numFrames:=1!! todo!
            if (numFrames>0) {
                //Сглаживаем в окне. Не влияет, но положено по теории.
                WindowFunction window = new WindowFunction();
                window.setWindowType("Hamming");//Rectangular");
                double[] win=window.generate(fftSampleSize);

                //забираем оцифровку из WAVE
                poleP.setColor(Color.GREEN);
                signals = new double[numFrames][];
                for (int f = 0; f < numFrames; f++) {
                    averAmp = 0.0d;
                    signals[f] = new double[fftSampleSize];
                    int startSample = f * fftSampleSize*2;
                    for (int n = 0; n < fftSampleSize; n++) {
                        signals[f][n] = getShort(buffer[startSample+n*2], buffer[startSample+n*2+1])*win[n];
                        averAmp += Math.abs(signals[f][n]);
                    }
                    //рисуем внизу амплитуду сигнала (скользящее)
                    averAmp = Math.min(200,0.2*averAmp/ fftSampleSize);
                    poleC.drawLine(numFrames - f - 1, (float) (hh-averAmp), numFrames - f , (float) (hh - lastAverageAmp), poleP);
                    lastAverageAmp= (lastAverageAmp+averAmp)/2;
                }

                //берем Фурье, переводим в реальные
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

                //нормализуем, отсекаем нерабочие частоты, компрессируем по яркости
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

                    for (int j = lowFreq; j < highFreq; j++) {
                        if (absoluteSpectrogram[f][j] < minValidAmp || maxIntens<=0)
                            spectrogram[f][j] = 0;
                        else {
                            //фильтр 1 - максимайзер
                            spectrogram[f][j] = (Math.log10(absoluteSpectrogram[f][j] / minIntens)) / perceptualDiff;    // 0..1
                            //варианты максимайзера, рабтают плохо.
                            //spectrogram[f][j] = (absoluteSpectrogram[f][j]-minIntens)/maxIntens;
                            //spectrogram[f][j] = absoluteSpectrogram[f][j]/maxIntens;

                            //фильтр 2 - гейт. ОЧЕНЬ НУЖЕН
                            if (spectrogram[f][j]<0.95)
                                spectrogram[f][j]=0;
                        }
                    }


                    ////Алгоритм №2, стоимость int(17)
                    //Отображаем в BARK. Патентованная шкала звукового восприятия. Разбиение в start()
                    for (int j = lowFreq; j < highFreq; j++) {
                        if (Bark[(int) (j*barkToFrequnits)] >= 0 && Bark[(int) (j*barkToFrequnits)] <= 16) {
                            barkArray[f][Bark[(int) (j*barkToFrequnits)]][0] += spectrogram[f][j];
                            barkArray[f][Bark[(int) (j*barkToFrequnits)]][1]++;
                        }
                    }
                    for (int j = 0; j < 17; j++) { // усредняем
                        if (barkArray[f][j][1] > 0)
                            barkArray[f][j][0] = barkArray[f][j][0] / barkArray[f][j][1];
                    }


                    //рисуем осцилляторы
                    poleP1.setColor(Color.argb(255, 255, 255, 255));
                    for (int j = lowFreq; j < highFreq; j++) {
                        if (maxIntens>0) {
                            //рисуем наверху интенсивности сигнала
                            mycolor=(int) (Math.min(100.0, 0.2 * 100.0 * absoluteSpectrogram[f][j] / averIntens));
                            poleP.setColor(Color.argb(255, (int) (Math.min(255.0, 2.0 * 255.0 * absoluteSpectrogram[f][j] / maxIntens)), mycolor, mycolor));
                            poleC.drawPoint(numFrames - f - 1, j - lowFreq, poleP);
                            //и яркие точки
                            if (spectrogram[f][j]>0)
                                poleC.drawPoint(numFrames - f - 1, j - lowFreq, poleP1);
                        }
                    }



                    //Алгоритм №1 стоимость byte(3)
                    //BARK+LS method - наклон, средний вектор между точек, ax+b+(r)
                    int a1;
                    int b1;
                    int r1;
                    {
                        double X_ = 0, Y_ = 0, YX_ = 0, X2_ = 0, Y2_ = 0, zz = 0, zz1 = 0, zr = 0, b, a, r;
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
                        //запакуем в 3 байта. Найдено экспериментально.
                        a1 = (int) (220 + Math.round(a * 2000));
                        b1 = (int) (50 + Math.round(b * 200));
                        r1 = (int) (100 + Math.round(r * 100));
                        //debug: проверим что помещаемся в байт
                        //if (a1>255 || a1<0 || b1<0 || b1>255 || r1<0 || r1>255)
                        //    Log.d(TAG,"Алгоритм №1 borders! a1="+a1+", b1="+b1+", r1="+r1);
                        a1 = (a1 & 0xFF);
                        b1 = (b1 & 0xFF);
                        r1 = (r1 & 0xFF);
                    }
                    //рисуем отпечаток ---------------------
                    poleP.setColor(Color.argb(255,100,255,100));
                    poleC.drawLine(numFrames - f    , 200+oldB1,
                                   numFrames - f - 1, 200+b1, poleP);
                    oldA1 = a1;
                    oldB1 = b1;
                    oldR1 = r1;


                    //пишем отпечаток в базу
                    if (recordMode>0) {
                        if (db!=null && db.isOpen()) {
                            //тут база тормозит, но увы
                            String s="insert into original" + recordMode
                                    + " (a,b,r"; //#1
                            for (int j = 0; j < 17; j++)
                                s+=", a"+j;
                            s+=") values (" + oldA1 + "," + oldB1 + "," + oldR1;
                            for (int j = 0; j < 17; j++) //#2
                                s+=", "+((int)(barkArray[f][j][0]*1000));
                            s+=")";
                            db.execSQL(s);
                            recordedFPs++;
                        }
                    } else {
                        if (needleFP==null) {
                            needleFP = new int[limitNeedle][3+17]; //a,b,r + a0..a16 (№2),
                            needleIndex=-1;
                        }
                        needleIndex=(needleIndex+1)%limitNeedle;
                        needleFP[needleIndex][0]= (int) oldA1;
                        needleFP[needleIndex][1]= (int) oldB1;
                        needleFP[needleIndex][2]= (int) oldR1;
                        for (int j = 0; j < 17; j++)
                            needleFP[needleIndex][3+j]= (int) (barkArray[f][j][0]*1000);
                        recognizeMe();
                    }


                    //Алгоритм №3
                    //все 16 запоминаем
                    //Алгоритм №4
                    //LS по всем частотам
                }
            }
            //pole.setImageBitmap(poleBM);
            pole.invalidate();
        }

    };

    void recognizeMe(){
        int i;

        poleP.setColor(Color.argb(255, 255, 255, 255));
        for (int t=1;t<=3;t++) {
            poleCN[t - 1].drawRect(new Rect(0, 0, poleN[t - 1].getWidth(), 200), poleP);
        }

        int jj;
        //по каждому сэмплу...
        for (int t=1;t<=3;t++) {
            if (originalFP[t-1] != null && needleFP != null) {
                maxEqual[t-1] = 0;
                maxIndex[t-1] = -1;
                //...сравниваем его с входящим звуком по всей длине
                for (int shift = -limitNeedle/2; shift < originalFP[t-1].length + limitNeedle/2; shift++) {
                    equal[t - 1] = 0.0f;
                    for (int j = 0; j < needleFP.length; j++) {
                        i = j + shift;
                        jj = (needleIndex + 1 + j) % limitNeedle;
                        if (i >= 0 && i < originalFP[t - 1].length) {
                            if (liteAlgo) {
                                //алгоритм №1 - векторизованный BARK, 3 величины
                                if (originalFP[t - 1][i][0] == needleFP[jj][0])
                                    equal[t - 1]++;
                                if (originalFP[t - 1][i][1] == needleFP[jj][1])
                                    equal[t - 1]++;
                                if (originalFP[t - 1][i][2] == needleFP[jj][2])
                                    equal[t - 1]++;
                            } else {
                                //алгоритм №2 - весь BARK, 17 величин
                                for (int j_ = 0; j_ < 17; j_++)
                                    if (needleFP[jj][3 + j_] > 0 && originalFP[t - 1][i][6 + j_] == needleFP[jj][3 + j_])
                                        equal[t - 1]++;
                            }
                        }
                    }
                    if (equal[t - 1] > maxEqual[t - 1]) {
                        maxEqual[t - 1] = equal[t - 1];
                        maxIndex[t - 1] = shift;
                    }
                    if (shift >= 0 && shift < originalFP[t - 1].length && totalMaximum[t-1]>0) {
                        originalFP[t - 1][shift][4]+=100*equal[t-1]/totalMaximum[t-1];
                    }
                }
                totalMaximum[t-1]=maxEqual[t - 1];

                poleP.setColor(Color.argb(255, 255, 100, 100));
                int area = limitNeedle/2; //усреднение на отрезках
                float average;
                float totalAverage=0;
                int counter;
                for (int k = originalFP[t-1].length - 1; k > 0; k--) {
                    originalFP[t - 1][k][4] = (int) (originalFP[t - 1][k - 1][4] * 0.995f); //самоочистка , строго 0.995. 0.9  не успевает собрать, нестабильное распознвание
                    if (originalFP[t - 1][k][4] < 0)
                        originalFP[t - 1][k][4] = 0;
                    average = 0;
                    counter = 0;
                    for (int kk = k - area; kk < k + area; kk++)
                        if (kk >= 0 && kk < originalFP[t - 1].length) {
                            average += originalFP[t - 1][kk][4];
                            counter++;
                        }
                    if (counter > 0) {
                        average = average / (float)counter;
                    }
                    originalFP[t - 1][k][5] = (int) (((float)originalFP[t - 1][k][4]-average)/(float)area);
                    if (originalFP[t - 1][k][5] >0)
                        totalAverage+=originalFP[t - 1][k][5];
                }
                originalFP[t-1][0][4] = 0;
                totalAverage=totalAverage/((float)originalFP[t-1].length);

                //сбор максимумов, некоторые рядом
                int over150count=0,over200count=0;
                int over200maximumI=-2*originalFP[t-1].length;
                int over200maximum=0;
                int over20timesCount=0;
                for (int k = 0; k<originalFP[t-1].length-1; k++) {
                    if (originalFP[t - 1][k][5] > 150){
                        if (originalFP[t - 1][k][5]/totalAverage>16)
                            over20timesCount++;
                        if (originalFP[t - 1][k][5] > 200) {
                            over200count++;
                            if (originalFP[t - 1][k][5]>over200maximum) {
                                over200maximum=originalFP[t - 1][k][5];
                                over200maximumI = k;
                            }
                        } else if (Math.abs(k-over200maximumI)<2) {
                            over200count++;
                        } else
                            over150count++;
                    }
                }


                poleP1.setColor(Color.argb(255, 0, 0, 100));
                poleCN[t - 1].drawLine(0,190 - 0.5f *17.0f*  totalAverage,poleBMn[t-1].getWidth(), 190 - 0.5f * 17.0f* totalAverage, poleP1);
                poleP.setColor(Color.argb(255, 255, 0, 0));
                if (andDraw){
                    for (int k=1;k<originalFP[t-1].length; k++) {
                        if (originalFP[t - 1][k-1][5]>0)
                            poleCN[t - 1].drawLine(poleK[t - 1] * (k-1),190 - 0.5f * originalFP[t - 1][k-1][5],poleK[t - 1] * k, 190 - 0.5f * originalFP[t - 1][k][5], poleP1);
                        if (originalFP[t - 1][k][5]>150)
                            poleCN[t - 1].drawLine(poleK[t - 1] * k, 190, poleK[t - 1] * k, 190 - 0.5f * originalFP[t - 1][k][5], poleP);
                    }
                }
                if (over20timesCount>0 && over200maximum>150 /*over200count>0*/) {
                    //Log.d(TAG, "Распознан!!! over20timesCount="+over20timesCount+", over200maximum="+over200maximum+", #" + t + ", over200maximum/totalAverage="+Math.round(over200maximum/totalAverage));//, over200count=" + over200count + ", over150count=" + over150count);
                    if (over150count > 5 || over200count>3)
                        poleP.setColor(Color.argb(255, 150, 100, 100)); //сомнительно
                    else
                        poleP.setColor(Color.argb(255, 255, 0, 0)); //уверенно
                    poleCN[t - 1].drawRect(new Rect((int) (poleK[t - 1] * over200maximumI - 10), 120, (int) (poleK[t - 1] * over200maximumI + 10), 200), poleP);

                    //рисуем тайминг
                    int k=over200maximumI*originalDuration[t-1]/originalFP[t-1].length;
                    poleCN[t - 1].drawText(""+(k/60)+":"+(k%60<10?"0":"")+(k%60),poleK[t - 1] * over200maximumI+5,50,poleP);
                }
            }
        }
    }

	public void start(Context context){
        //Log.d(TAG,"Recorder try starting... State="+state.toString());
        if (state == State.READY || state == State.INITIALIZING){

            //init Bark scale
            //соответствие 4000Гц к 300 к 17 полосам восприятия слуха
            //300 - as is патентованная формула
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

            poleP.setTextSize(poleBMn[0].getWidth()/10);



            //open DB
            String dir = Environment.getExternalStorageDirectory().getPath();
            File dbfile = new File(dir+"/echolotFingerprints.db");
            db = SQLiteDatabase.openOrCreateDatabase(dbfile, null);
            if (!db.isOpen()){
                //Log.d(TAG,"Its open? "  + db.isOpen());
                state = State.ERROR;
                return;
            }
            recordedFPs=0;
            if (recordMode>0) {
                //Режим записи.
                try {
                    db.execSQL("drop table original"+recordMode);
                    Log.d(TAG, "удалена база данных original"+recordMode);
                } catch (SQLiteException exception) {
                    //Log.d(TAG,"ошибка удаления базы original");
                    ;//exception.printStackTrace();
                }
                try {
                    String s="create table original"+recordMode+" ("
                            + "id integer primary key,"
                            + "a integer,"
                            + "b integer,"
                            + "r integer";
                    for (int j = 0; j < 17; j++)
                        s+=", a"+j+" integer";
                    s+=")";
                    db.execSQL(s);
                    //db.execSQL("CREATE INDEX original"+recordMode+"_a_idx ON original"+recordMode+" (a)");
                } catch (SQLiteException exception) {
                    if (recordMode==1) Log.e(TAG,"ошибка создания базы original"+recordMode+", "+exception.toString());
                    //exception.printStackTrace();
                }
                MyActivity.hint.setText("Идёт запись...");
            } else {
                //Режим распознавания
                originalFP = new int[3][][];
                poleK = new float[3];
                originalDuration = new int[3];
                equal=new float[3];
                maxEqual=new float[3];
                totalMaximum=new float[3];
                slidingMaximumIndex=new int[3];
                slidingMaximum=new float[3];
                maxIndex=new int[3];

                //выбираем образцы
                int i, amount=0;
                String hint="";
                for (int t=1;t<=3;t++) {
                    originalDuration[t-1] =0;
                    try {
                        String q1 = "select * from original"+t+" order by id asc";
                        Cursor c2 = db.rawQuery(q1, null);
                        originalDuration[t-1] = c2.getCount() / (sampleRate / (fftSampleSize));
                        if (originalDuration[t-1]<20)
                            hint+=" Сэмпл "+t+" слишком мал (<20s)";
                        if ((originalDuration[t-1]>400 && liteAlgo))
                            hint+=" Сэмпл "+t+" слишком велик (>400s)";
                        if ((originalDuration[t-1]>60 && !liteAlgo))
                            hint+=" Сэмпл "+t+" слишком велик (>60s)";
                        amount++;

                        if (c2.getCount() > 0 && c2 != null && c2.moveToFirst()) {
                            originalFP[t-1] = new int[c2.getCount()][6+17];
                            i = 0;
                            do {
                                originalFP[t-1][i][0] = c2.getInt(c2.getColumnIndex("a"));
                                originalFP[t-1][i][1] = c2.getInt(c2.getColumnIndex("b"));
                                originalFP[t-1][i][2] = c2.getInt(c2.getColumnIndex("r"));
                                originalFP[t-1][i][3] = 0;//не используется
                                originalFP[t-1][i][4] = 0; //относительный пик по мгновенному сравнению
                                originalFP[t-1][i][5] = 0; //абсолютный пик по сумме за несколько секунд
                                for (int j = 0; j < 17; j++)
                                    originalFP[t-1][i][6+j] = c2.getInt(c2.getColumnIndex("a"+j)); // алгоритм №2
                                i++;
                            } while (c2.moveToNext());
                            poleK[t-1] = (float) poleBMn[t-1].getWidth() / (float) (originalFP[t-1].length - limitNeedle);
                            Log.d(TAG, "Выбрано по образцу сэмплов #"+t+": " + c2.getCount() + " (" + originalFP[t-1].length + ") " + ", originalDuration1=" + originalDuration[t-1] + ", k=" + poleK[t-1]);
                        }
                        c2.close();
                    } catch (Exception e) {
                    }
                    if (originalDuration[t-1]==0) {
                        hint+=" Сэмпл "+t+" не записан.";
                        originalFP[t-1] = null;
                    }
                }
                if (amount==0) {
                    MyActivity.searchButton.setChecked(false);
                    MyActivity.hint.setText("Нет образцов! Запишите парочку");
                    return;
                } else {
                    MyActivity.hint.setText(hint+" Идентифицируем...");
                }

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
        MyActivity.hint.setText("Остановлен.");
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
