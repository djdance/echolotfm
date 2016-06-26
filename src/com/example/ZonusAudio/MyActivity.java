package com.example.ZonusAudio;

import android.app.Activity;
import android.os.Bundle;
import android.util.Log;
import android.view.*;
import android.widget.ImageView;
import android.widget.Toast;
import android.widget.ToggleButton;

public class MyActivity extends Activity implements View.OnClickListener {
    String TAG="djd";
    // detection parameters
    ExtAudioRecorder extAudioRecorder;
    public static ImageView pole;
    public static ImageView poleN1,poleN2,poleN3;
    ToggleButton searchButton,recButton1,recButton2,recButton3;

    @Override
    public void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.main);
        pole=(ImageView) findViewById(R.id.imageView);
        poleN1=(ImageView) findViewById(R.id.imageViewN1);
        poleN2=(ImageView) findViewById(R.id.imageViewN2);
        poleN3=(ImageView) findViewById(R.id.imageViewN3);
        searchButton=(ToggleButton) findViewById(R.id.searchToggleButton);
        searchButton.setOnClickListener(this);
        recButton1=(ToggleButton) findViewById(R.id.recToggleButton1);
        recButton1.setOnClickListener(this);
        recButton1.setTag(1);
        recButton2=(ToggleButton) findViewById(R.id.recToggleButton2);
        recButton2.setOnClickListener(this);
        recButton2.setTag(2);
        recButton3=(ToggleButton) findViewById(R.id.recToggleButton3);
        recButton3.setOnClickListener(this);
        recButton3.setTag(3);
        setTitle("Эхолот Распознавание");
    }

    void startRec(int recordMode){
        //Log.d(TAG, "starting...");
        if (extAudioRecorder!=null && extAudioRecorder.getState()== ExtAudioRecorder.State.RECORDING)
            stopRec();
        // Start recording
        extAudioRecorder = ExtAudioRecorder.getInstanse();
        extAudioRecorder.recordMode=recordMode;
        extAudioRecorder.start(getApplicationContext());
        //Log.d(TAG, "rec state=" + extAudioRecorder.getState().toString());
    }
    void stopRec(){
        //Log.d(TAG, "stopping...");
        if (extAudioRecorder!=null && extAudioRecorder.getState()== ExtAudioRecorder.State.RECORDING) {
            // Stop recording
            extAudioRecorder.stop();
            if (extAudioRecorder.recordMode>0)
                mytoast("Обработано сэмплов: "+extAudioRecorder.recordedFPs);
            extAudioRecorder=null;
            //Log.d(TAG, "ok");
        }
    }

    @Override
    protected void onPause() {
        super.onPause();
        stopRec();
    }

    @Override
    public void onClick(View v) {
        if (v==searchButton || v==recButton1 || v==recButton2 || v==recButton3){
            if (v==searchButton) {
                recButton1.setChecked(false);
                recButton2.setChecked(false);
                recButton3.setChecked(false);
            }else if (v==recButton1 || v==recButton2 || v==recButton3) {
                searchButton.setChecked(false);
                recButton1.setChecked(v==recButton1?((ToggleButton)v).isChecked():false);
                recButton2.setChecked(v==recButton2?((ToggleButton)v).isChecked():false);
                recButton3.setChecked(v==recButton3?((ToggleButton)v).isChecked():false);
            }

            if (((ToggleButton)v).isChecked())
                startRec(v.getTag()!=null?(Integer) (v.getTag()):0);
            else
                stopRec();
        }

    }

    void mytoast(final String s){
        runOnUiThread(new Runnable() {
            public void run() {
                Toast toast = Toast.makeText(getApplicationContext(), s, Toast.LENGTH_SHORT);
                toast.show();
            }
        });
    }

}
