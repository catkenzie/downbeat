# Downbeat
Real-time drum transcription: machine transforming a live drum-kit session into sheet music. The BPM is set by the drummer and the beat is recorded with a singular microphone which supports most basic setups that artists have at their disposal. 

## About the Project 
This collaborative project fullfilled capstone credit for the University of Michigan School of Engineering. 

Drum classification is implemented in C on the STM32 Nucleo F747ZI microprocessor and the note transcription is implemented in Python on the RPi4.


## Classification 

The drum classification algorithm uses non-negative matrix factorization (NNMF) to separate and determine which piece of the kit was hit. Each hit and its onset times are recorded and passed to the RPI for the transcription algorithm to determine the beat pattern. 

Isolated training samples of each drum piece are processed before runtime to minimize deviations in spectral content due to the acoustic environment, as well as any physical changes in the relative placement of the kit pieces and the microphone.


## Transcription

Pre-runtime the transcription algorithm takes the set BPM, determines all possible beat patterns, and weighs the likelihood of each pattern. During runtime everytime time a hit is detected, the hit spacing is calculated, rounding is accounted for and the beat pattern is created. 

The GUI is implemented in C using the Simple Directmedia Layer and Guido Music Notation library to render the musical notation as a PNG to display. 
